clc; clear; format short; tic;

%% === Load Input Data ===
mm = load('loaddata33bus.m');                % Bus data [Bus, Load_P, Load_Q, PV_Capacity] replace with 33 or 69 bus
ll = load('linedata33bus.m');                % Line data [From_Bus, To_Bus, R, X] replace with 33 or 69 bus
PVout  = load('PV_out_profile.m');           % Hourly PV output profile
L_prof = load('Residential_Load_Profile.m'); % Hourly load profile

%% === Define System Base Values ===
MVAb = 100;                                  % System base power (MVA)
KVb = 12.66;                                 % System base voltage (kV)
Zb = (KVb^2) / MVAb;                         % System base impedance (Ohm)

%% === Optimization Method Selection ===
opt = 'PSO_TS'; % Choose optimization algorithm: 'PSO', 'TS', or 'PSO_TS'

%% === Initialization Parameters ===
seed = round(sum(1000 * clock));             % Random seed based on current timestamp
rng(seed);

BESS_Number = 5;                             % Number of BESS units
Candidate_Buses = setdiff(mm(:,1)', 1);      % Candidate buses for BESS installation (excluding slack bus)
SoC_max = 0.9;                               % Maximum state of charge (per unit)
SoC_min = 0.2;                               % Minimum state of charge (per unit)
BESS_Demand = zeros(size(mm, 1), 1);         % Initial BESS demand (zero)
bus = size(mm, 1);                           % Total number of buses

stagnation_limit = 50;                       % Stagnation limit for optimization algorithms
cap = (sum(mm(:, 2)) / BESS_Number) * 1.25;  % Estimated capacity for BESS
upper_bound = round(cap / 50) * 50;          % Upper bound for BESS sizing (rounded to nearest 50)
lower_bound = -upper_bound;                  % Lower bound for BESS sizing
max_expected_Pwr_Losses = 0.05;              % Expected maximum power loss (p.u.)
max_expected_Vol_Dev = 0.05;                 % Expected maximum voltage deviation (p.u.)

%% === Initial Solution Settings ===

use_ranked_bus_guidance = true;               % Set to true to utilize a pre-generated initial solution
guided_fraction = 0.5;                        % 50% solution guided from ranked buses in initial iteration

if use_ranked_bus_guidance
    % Load bus ranking if file exists
    num_buses = size(mm, 1);
    ranking_filename = fullfile('results', sprintf('Bus_Ranked_%dbus.mat', num_buses));
    
    if exist(ranking_filename, 'file')
        ranking_data = load(ranking_filename);
        Bus_Ranked = ranking_data.Bus_Ranked';  % Ranked bus indices based on performance
        fprintf('>> Bus ranking loaded successfully from "%s"\n', ranking_filename);
    else
        warning('>> Ranked bus file "%s" not found. Proceeding without initial guidance.', ranking_filename);
        use_ranked_bus_guidance = false;
        Bus_Ranked = []; % No initial bus ranking used
    end
else
    Bus_Ranked = []; % No initial bus ranking used
end

%% === Preallocation for Result Matrices ===
Matrix_BESS_Demand    = zeros(size(mm, 1), 24);
Matrix_Voltages       = zeros(size(mm, 1), 24);
Matrix_Results        = zeros(24, 11);
obj_hourly            = zeros(24, 1);
fitness_iter          = zeros(24, 1000);
best_eval_hourly      = zeros(24, 1);
best_iter_hourly      = zeros(24, 1);

%% === Optimization Loop (24 Hours) ===
for hour = 1:24
    sel_pv = PVout(hour, 2);   % Selected PV scaling factor
    sel_lp = L_prof(hour, 2);  % Selected load scaling factor

    % Perform optimization according to selected method
    switch opt
        case 'PSO'
            [BESS_Location, BESS_Output, obj, fitness_history, fitness_eval, iter, eval_total, best_iter_idx, best_eval_idx] = ...
                Placement_Optimization_PSO(lower_bound, upper_bound, BESS_Number, Candidate_Buses, hour, mm, ll, sel_pv, sel_lp, ...
                MVAb, Zb, stagnation_limit, Bus_Ranked, use_ranked_bus_guidance, guided_fraction);
        case 'TS'
            [BESS_Location, BESS_Output, obj, fitness_history, fitness_eval, iter, eval_total, best_iter_idx, best_eval_idx] = ...
                Placement_Optimization_TS(lower_bound, upper_bound, BESS_Number, Candidate_Buses, hour, mm, ll, sel_pv, sel_lp, ...
                MVAb, Zb, stagnation_limit, Bus_Ranked, use_ranked_bus_guidance, guided_fraction);
        case 'PSO_TS'
            [BESS_Location, BESS_Output, obj, fitness_history, fitness_eval, iter, eval_total, best_iter_idx, best_eval_idx] = ...
                Placement_Optimization_PSO_TS(lower_bound, upper_bound, BESS_Number, Candidate_Buses, hour, mm, ll, sel_pv, sel_lp, ...
                MVAb, Zb, stagnation_limit, Bus_Ranked, use_ranked_bus_guidance, guided_fraction);
        otherwise
            error('Unsupported optimization method: %s', opt);
    end

    % Perform hourly load flow analysis
    BESS_Demand = zeros(size(mm, 1), 1);
    BESS_Demand(BESS_Location) = BESS_Output;
    [voltage, P_Loss_Kw, Q_Loss_KVAr, PL, QL, ld, node, branch, va, voltage_deviation, total_voltage_deviation] = ...
        HourlyLoadFlow(mm, ll, sel_pv, sel_lp, MVAb, Zb, BESS_Demand);


    % Store hourly results
    [minVoltage, minBus] = min(abs(voltage));
    [maxVoltage, maxBus] = max(abs(voltage));
    Matrix_Voltages(:, hour)     = voltage;
    Matrix_BESS_Demand(:, hour)  = BESS_Demand;
    Matrix_P_Nett(:, hour)       = ld(:, 2);
    Matrix_P_Loss(:, hour)       = P_Loss_Kw;
    Matrix_Q_Loss(:, hour)       = Q_Loss_KVAr;
    Matrix_Results(hour, :)      = [hour, sum(mm(:, 4) * sel_pv), 0, sum(mm(:, 2) * sel_lp), 0, sum(P_Loss_Kw), sum(Q_Loss_KVAr), maxVoltage, maxBus, minVoltage, minBus];
    obj_hourly(hour)             = obj;
    len = length(fitness_history);
    fitness_iter(hour, 1:len)    = fitness_history;
    best_eval_hourly(hour)       = best_eval_idx;
    best_iter_hourly(hour)       = best_iter_idx;
end

toc;
total_runtime = toc;

%% === Post-Processing ===
Matrix_Results(:, 3) = sum(Matrix_BESS_Demand, 1)';
Matrix_Results(:, 5) = sum(Matrix_P_Nett, 1)';
Matrix_Results(:, 6) = sum(Matrix_P_Loss, 1)';
Matrix_Results(:, 7) = sum(Matrix_Q_Loss, 1)';
Max_volt_dev = max(max(abs(Matrix_Voltages - 1.0)));
T_24h_obj             = sum(obj_hourly);
T_24h_effectives_eval = sum(best_eval_hourly);
T_24h_effectives_iter = sum(best_iter_hourly);
row_lengths = sum(fitness_iter ~= 0, 2);
max_length = max(row_lengths);
final_length = max(max_length, 300);
fitness_iter_hourly = fitness_iter(:, 1:final_length);

%% === Display and Save Results ===
BESS_Optimal_Location = placement_display_results(Matrix_BESS_Demand, fitness_iter_hourly, T_24h_obj, ...
    T_24h_effectives_eval, T_24h_effectives_iter, opt, total_runtime, Max_volt_dev, SoC_max, SoC_min, ...
    Matrix_Results, Matrix_Voltages, mm, BESS_Number, seed);

if ~exist('results', 'dir')
    mkdir('results');
end

filename = fullfile('results', ['BESS_Optimal_Location_' opt '_' num2str(bus) '.mat']);
save(filename, 'BESS_Optimal_Location', 'opt', 'bus');