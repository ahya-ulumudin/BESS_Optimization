clc; clear; format short; tic;

%%  LOAD NETWORK AND PROFILE DATA
mm     = load('loaddata33bus.m');              % Bus data (bus number, load data)chenge 33 or 69 bus
ll     = load('linedata33bus.m');              % Line data (connection, impedance)change 33 or 69 bus
PVout  = load('PV_out_profile.m');             % Hourly PV output profile
L_prof = load('Residential_Load_Profile.m');   % Hourly load profile

% Base system parameters
MVAb = 100;                                    % Base power (MVA)
KVb  = 12.66;                                  % Base voltage (kV)
Zb   = (KVb^2) / MVAb;                         % Base impedance (Ohm)

%% CONFIGURATION PARAMETERS
Bus_Placement = [14, 18, 24, 26, 31];           % Predefined BESS installation buses 33: (14, 18, 24, 26, 31) 69: (11, 26, 27, 50, 61) or manual provided
BESS_Number   = 5;                              % Number of BESS units
bus           = size(mm, 1);                    % Number of buses

% State of Charge Limits
SoC_max = 0.9;
SoC_min = 0.2;

% Random Seed Initialization
seed = round(sum(1000 * clock));                % Use current timestamp for reproducibility
rng(seed);

% Search Space Boundaries
cap = (sum(mm(:, 2)) / BESS_Number) * 1.25;
upper_bound = round(cap / 50) * 50;              % Rounded upper bound (multiple of 50)
lower_bound = -upper_bound;                      % Lower bound

% System Performance Limits
max_expected_Pwr_Losses = 0.05;                  % Acceptable maximum system losses
max_expected_Vol_Dev    = 0.05;                  % Acceptable maximum voltage deviation
BESS_RTE                = 0.9;                   % Round Trip Efficiency (90%)      
BESS_Eff                = sqrt(BESS_RTE);        % Charging/Discharging Efficiency

% Optimization Parameters
stagnation_limit = 100;                          % Stagnation limit for optimizer
opt = 'PSO_TS';                                      % Selected optimization method ('PSO', 'TS', or 'PSO_TS')

%% DATA STRUCTURES INITIALIZATION
BESS_Demand         = zeros(bus, 1);             % BESS injection per bus
Matrix_Results      = zeros(24, 11);             % Main results for each hour
Matrix_Voltages     = zeros(bus, 24);             % Bus voltages per hour
Matrix_BESS_Demand  = zeros(bus, 24);             % BESS dispatch per hour
Matrix_P_Nett       = zeros(bus, 24);             % Net active power per hour
Matrix_P_Loss       = zeros(size(ll,1), 24);      % Active power loss per hour
Matrix_Q_Loss       = zeros(size(ll,1), 24);      % Reactive power loss per hour
Matrix_Load         = zeros(bus, 24);             % Bus load per hour

%% LOAD PREVIOUS SOLUTION IF AVAILABLE
use_previous = false;  % Set 'true' to use previous result, 'false' for random initialization

filename = fullfile('results', sprintf('BESS_Demand_%dbus-%s.mat', bus, opt));

if use_previous && exist(filename, 'file')
    load(filename, 'BESS_Output', 'Bus_Placement');
    fprintf('>> Loaded previous BESS solution from %s\n', filename);
    initial_solution = reshape(BESS_Output, 1, []);  % Flatten to 1-D array
else
    if use_previous
        warning('>> Initial solution requested but file not found. Using random initialization.');
    else
        fprintf('>> Skipping initial solution. Using random initialization.\n');
    end
    initial_solution = [];
end

%% DEFINE OBJECTIVES FUCTION
objective_function = @(BESS_Output) Sizing_Objective( ...
    round(BESS_Output / 10) * 10, mm, ll, PVout, L_prof, MVAb, Zb, upper_bound, Bus_Placement, BESS_Eff);

%% CALL OPTIMIZER
switch opt
    case 'PSO'
        [BESS_Output, obj, fitness_history, iter, eval_total] = ...
            Sizing_Optimization_PSO(mm, lower_bound, upper_bound, BESS_Number, ...
            objective_function, stagnation_limit, initial_solution, BESS_Eff);

    case 'TS'
        [BESS_Output, obj, fitness_history, iter, eval_total] = ...
            Sizing_Optimization_TS(mm, lower_bound, upper_bound, BESS_Number, ...
            objective_function, stagnation_limit, initial_solution, BESS_Eff);

    case 'PSO_TS'
        [BESS_Output, obj, fitness_history, iter, eval_total] = ...
            Sizing_Optimization_PSO_TS(mm, lower_bound, upper_bound, BESS_Number, ...
            objective_function, stagnation_limit, initial_solution, BESS_Eff);

    otherwise
        error('Unsupported optimization method: %s', opt);
end

%% LOAD FLOW ANALYSIS OVER 24 HOURS
for hour = 1:24
    sel_pv = PVout(hour, 2);
    sel_lp = L_prof(hour, 2);

    % Update BESS Demand for Current Hour
    raw_out = round(BESS_Output(:, hour) / 10) * 10;   % kW per unit BESS
    BESS_Demand(:) = 0;
    
    pos_mask = raw_out >= 0;   % discharge
    neg_mask = raw_out <  0;   % charge
    if any(pos_mask)
        idx_pos = Bus_Placement(pos_mask);
        BESS_Demand(idx_pos) = raw_out(pos_mask) * BESS_Eff;
    end
    if any(neg_mask)
        idx_neg = Bus_Placement(neg_mask);
        BESS_Demand(idx_neg) = raw_out(neg_mask) / BESS_Eff;
    end

    % Perform Hourly Load Flow
    [voltage, P_Loss_Kw, Q_Loss_KVAr, PL, QL, ld, node, branch, va, voltage_deviation, total_voltage_deviation] = ...
        HourlyLoadFlow(mm, ll, sel_pv, sel_lp, MVAb, Zb, BESS_Demand);

    % Store Results
    [minVoltage, minBus] = min(abs(voltage));
    [maxVoltage, maxBus] = max(abs(voltage));

    Matrix_Voltages(:, hour)    = voltage;
    Matrix_BESS_Demand(:, hour) = BESS_Demand;
    Matrix_P_Nett(:, hour)      = ld(:, 2);
    Matrix_P_Loss(:, hour)      = P_Loss_Kw;
    Matrix_Q_Loss(:, hour)      = Q_Loss_KVAr;
    Matrix_Load(:, hour)        = mm(:, 2) * sel_lp;

    % Record to Main Results
    Matrix_Results(hour, :) = [ ...
        hour, ...
        sum(mm(:,4) * sel_pv), ... % Total generation
        sum(Matrix_BESS_Demand(:,hour)), ... % Total BESS demand
        sum(mm(:,2) * sel_lp), ... % Total load
        sum(Matrix_P_Nett(:,hour)), ... % Net system load
        sum(Matrix_P_Loss(:,hour)), ... % Active power losses
        sum(Matrix_Q_Loss(:,hour)), ... % Reactive power losses
        maxVoltage, ...
        maxBus, ...
        minVoltage, ...
        minBus];
end

toc;
total_runtime = toc;

%% === SAVE FINAL OUTPUT ===
results_dir = 'results';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

num_buses = size(mm, 1);
num_hours = size(BESS_Output, 2);
Max_Voltage_Deviation_Per_Bus = max(abs(Matrix_Voltages - 1.0), [], 2);

% Save Fitness History (.mat)
fitness_filename = sprintf('fitness_history_sizing_%d-bus-%s.mat', num_buses, opt);
save(fullfile(results_dir, fitness_filename), 'fitness_history');

% Save BESS Output (.mat)
bess_filename = sprintf('BESS_Demand_%dbus-%s.mat', num_buses, opt);
save(fullfile(results_dir, bess_filename), 'BESS_Output', 'Bus_Placement');

% Save Fitness History as .csv
fitness_csv = sprintf('fitness_history_sizing_%d-bus-%s.csv', num_buses, opt);
fitness_table = array2table(fitness_history(:), 'VariableNames', {'Fitness'});
writetable(fitness_table, fullfile(results_dir, fitness_csv));

% Save BESS Output as .csv
bess_csv = sprintf('BESS_Demand_%dbus-%s.csv', num_buses,  opt);
bess_output_table = array2table([Bus_Placement(:), BESS_Output], ...
    'VariableNames', ["Bus", strcat("Hour_", string(1:size(BESS_Output, 2)))]);
writetable(bess_output_table, fullfile(results_dir, bess_csv));

fprintf('\n Final results saved to: %s\n', results_dir);

%% DISPLAY FINAL RESULTS
display_results(Matrix_BESS_Demand, Bus_Placement, obj, max(Max_Voltage_Deviation_Per_Bus), SoC_max, SoC_min, ...
    Matrix_Results, Matrix_Voltages, mm, lower_bound, upper_bound, BESS_Output, opt, total_runtime, iter, eval_total, seed);
