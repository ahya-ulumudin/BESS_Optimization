%% =========================================================================
% 24-HOUR LOAD FLOW ANALYSIS WITHOUT BESS
% =========================================================================
clc; clear; format short; tic;

%% === Load System and Profile Data ===
mm = load('loaddata33bus.m');                 % Bus data [Bus, Load_P, Load_Q, PV_Capacity] choose 33 or 69 bus
ll = load('linedata33bus.m');                 % Line data [From_Bus, To_Bus, R, X] choose 33 or 69 bus
PVout  = load('PV_out_profile.m');            % PV output profile (hourly)
L_prof = load('Residential_Load_Profile.m');  % Load profile (hourly)

% Initialize fixed parameters
SoC_max = 0.9;                                % Maximum state of charge
SoC_min = 0.2;                                % Minimum state of charge
rng(42);                                      % Seed for reproducibility

% System base values
MVAb = 100;                                   % Base power (MVA)
KVb = 12.66;                                  % Base voltage (kV)
Zb = (KVb^2) / MVAb;                          % Base impedance (ohm)

% Preallocate matrices
Matrix_Results    = zeros(24, 11);             % Summary matrix: Hourly performance
Matrix_Voltages   = zeros(size(mm,1), 24);     % Bus voltages (per hour)
Matrix_BESS_Demand= zeros(size(mm,1), 24);     % BESS injection (per hour)
Matrix_P_Nett     = zeros(size(mm,1), 24);     % Net active power (per hour)
Matrix_P_Loss     = zeros(size(ll,1), 24);     % Active power loss per line
Matrix_Q_Loss     = zeros(size(ll,1), 24);     % Reactive power loss per line
Matrix_Load       = zeros(size(mm,1), 24);     % Load (per hour)
obj_hourly        = zeros(24, 1);              % Hourly objective values

% Initialize BESS variables (no BESS installed)
BESS_Demand    = zeros(size(mm,1), 1);
BESS_Locations = [];
BESS_Output    = [];

%% === 24-Hour Load Flow Simulation ===
for hour = 1:24
    % Select PV output and load scaling for the current hour
    sel_pv = PVout(hour, 2);                 
    sel_lp = L_prof(hour, 2);                 

    % Perform hourly load flow (without BESS contribution)
    [voltage, P_Loss_Kw, Q_Loss_KVAr, PL, QL, ld, node, branch, va, voltage_deviation, total_voltage_deviation] = ...
        HourlyLoadFlow(mm, ll, sel_pv, sel_lp, MVAb, Zb, BESS_Demand);

    % Evaluate objective performance
    obj = Placement_Objective(BESS_Locations, BESS_Output, mm, ll, sel_pv, sel_lp, MVAb, Zb);

    % Record hourly results
    [minVoltage, minBus] = min(abs(voltage));
    [maxVoltage, maxBus] = max(abs(voltage));

    Matrix_Voltages(:, hour)     = voltage;
    Matrix_BESS_Demand(:, hour)  = BESS_Demand;
    Matrix_P_Nett(:, hour)       = ld(:,2);
    Matrix_P_Loss(:, hour)       = P_Loss_Kw;
    Matrix_Q_Loss(:, hour)       = Q_Loss_KVAr;
    Matrix_Load  (:, hour)       = mm(:, 2) * sel_lp;

    Matrix_Results(hour, 1)  = hour;
    Matrix_Results(hour, 2)  = sum(mm(:,4)) * sel_pv;    % Total PV generation (kW)
    Matrix_Results(hour, 4)  = sum(mm(:,2)) * sel_lp;    % Total system load (kW)
    Matrix_Results(hour, 8)  = maxVoltage;               % Maximum voltage (p.u.)
    Matrix_Results(hour, 9)  = maxBus;                   % Bus index with maximum voltage
    Matrix_Results(hour,10)  = minVoltage;               % Minimum voltage (p.u.)
    Matrix_Results(hour,11)  = minBus;                   % Bus index with minimum voltage

    obj_hourly(hour) = obj;
end

%% === Post-Processing ===
Matrix_Results(:, 3) = sum(Matrix_BESS_Demand, 1)';    % Total hourly BESS injection (should be zero)
Matrix_Results(:, 5) = sum(Matrix_P_Nett, 1)';         % Total hourly net power
Matrix_Results(:, 6) = sum(Matrix_P_Loss, 1)';         % Total active power losses
Matrix_Results(:, 7) = sum(Matrix_Q_Loss, 1)';         % Total reactive power losses

Max_volt_dev = max(max(abs(Matrix_Voltages - 1.0)));   % Maximum voltage deviation (p.u.)
obj_24h = sum(obj_hourly);                             % Aggregate objective across 24 hours
total_BESS_Demand = sum(Matrix_BESS_Demand(:));        % Total BESS output (should be zero)

% ensure 'results' folder existing
if ~exist('results', 'dir')
    mkdir('results');
end

% Save .CSV file in results folder
num_buses = size(mm, 1);
output_filename = fullfile('results', sprintf('ranking_result_%dbus.csv', num_buses));
Bus_Ranked = bus_ranking(Matrix_Voltages, Matrix_P_Nett, mm, output_filename); % call bus ranking function

% Save .MAT file in results folder
matfile_name = fullfile('results', sprintf('Bus_Ranked_%dbus.mat', num_buses));
save(matfile_name, 'Bus_Ranked');

%% === Display Results and Generate Visualizations ===
display_results(obj_24h, Max_volt_dev, ...
                Matrix_Results, Matrix_Voltages, mm);

toc;


function display_results(obj_24h, Max_volt_dev, ...
                         Matrix_Results, Matrix_Voltages, mm)
% DISPLAY_RESULTS Displays and visualizes the 24-hour load flow analysis results.
%   Summarizes system performance metrics, displays key statistics, and generates
%   multiple plots for a comprehensive network evaluation.

    % === Display Hourly Load Flow Results ===
    disp('<strong>    DATA RESULTS:</strong>');
    column_names =  {'Hour', 'P_Gen_(kW)', 'BESS_Demand_(kW)', 'P_Load_(kW)', 'P_Nett_(kW)', ...
                     'P_Loss_(kW)', 'Q_Loss_(kVar)', 'V_Max_(p.u.)', 'V_Bus_Max', 'V_Min_(p.u.)', 'V_Bus_Min'};
    results_table = array2table(Matrix_Results, 'VariableNames', column_names);
    disp(results_table);
   
    num_buses = size(mm, 1);

    % === Calculate and Display Summary Statistics ===
    percentage_pv = (sum(mm(:, 4)) / sum(mm(:, 2))) * 100;
    percentage_losses = (sum(Matrix_Results(:, 6)) / sum(Matrix_Results(:, 4))) * 100;
    percentage_v_dev = Max_volt_dev / 1.0 * 100;
    BESS_Ratio = (sum(abs(Matrix_Results(:, 3)))) / (sum(Matrix_Results(:, 2)));
    pv_vs_load = (sum(Matrix_Results(:, 2)) / sum(Matrix_Results(:, 4))) * 100;
    bess_vs_load = (sum(abs(Matrix_Results(:, 3))) / sum(abs(Matrix_Results(:, 4)))) * 100;
    
    disp('<strong>24 HOUR LOAD FLOW SUMMARY:</strong>');
    fprintf('Total Load             : %.2f  kW\n', sum(mm(:, 2)));
    fprintf('Total Load Consumption : %.2f kWh\n', sum(Matrix_Results(:, 4)));
    fprintf('Total Installed PV     : %.2f  kW  (%.3f %%)\n', sum(mm(:, 4)), percentage_pv);
    fprintf('Total PV Supply        : %.2f kWh (%.3f %%)\n', sum(Matrix_Results(:, 2)), pv_vs_load);
    fprintf('Average P.Losses       : %.2f kW  (%.3f %%)\n', sum(Matrix_Results(:, 6)) / 24, percentage_losses);
    fprintf('Average V. Deviation   : %.3f p.u. (%.3f %%)\n', mean(mean(abs(Matrix_Voltages - 1))), (mean(mean(abs(Matrix_Voltages - 1))) / 1.0) * 100);
    fprintf('Max Voltage Deviation  : %.3f p.u. (%.3f %%)\n', Max_volt_dev, percentage_v_dev);
    fprintf('Max Bus Voltage        : %.2f p.u.\n', max(Matrix_Voltages(:)));
    fprintf('Min Bus Voltage        : %.2f p.u.\n', min(Matrix_Voltages(:)));
    fprintf('\n24 Hour Total Objectives Value: %.2f\n', obj_24h);

    fig_position = [100, 100, 1920, 1080];
    % === Plot Hourly Trends for Power and Voltage ===
    hours = Matrix_Results(:, 1)';
    figure('Position', fig_position);
    
    % Plot power-related trends (left y-axis)
    yyaxis left
    p1 = plot(hours, Matrix_Results(:, 4), '--', 'Color', [0.0, 0.0, 1.0], 'LineWidth', 3, 'DisplayName', 'Total P Load (kW)');
    hold on;
    p2 = plot(hours, Matrix_Results(:, 2), '--', 'Color', [1.0, 0.5, 0.0], 'LineWidth', 3, 'DisplayName', 'Total PV Gen (kW)');
    p3 = plot(hours, Matrix_Results(:, 5), '-o', 'Color', [0.6, 0.0, 0.0], 'LineWidth', 3, 'DisplayName', 'Total P Net (kW)');
    p4 = plot(hours, Matrix_Results(:, 6), '-v', 'Color', [1.0, 0.0, 0.0], 'LineWidth', 3, 'DisplayName', 'P Loss (kW)');
    ylabel('Power (kW)', 'FontSize', 18, 'FontWeight', 'bold');
    xlabel('Hour', 'FontSize', 18, 'FontWeight', 'bold');
    ylim([-3000 7500]);
    grid on;
    
    % Plot voltage trends (right y-axis)
    yyaxis right
    p5 = plot(hours, Matrix_Results(:, 8), '-+', 'Color', [0.0, 0.9, 0.7], 'LineWidth', 3, 'DisplayName', 'V Max (p.u.)');
    p6 = plot(hours, Matrix_Results(:, 10), '-x', 'Color', [0.7, 0.9, 0.0], 'LineWidth', 3, 'DisplayName', 'V Min (p.u.)');
    ylabel('Voltage (p.u.)', 'FontSize', 18, 'FontWeight', 'bold');
    ylim([0.35 1.15]);
    grid on;
    
    xticks(hours);
    set(gca, 'FontSize', 15, 'FontWeight', 'normal');
    legend([p1, p2, p3, p4, p5, p6], 'Location', 'southoutside', 'FontSize', 18, 'NumColumns', 2);
    title(['24-Hour Load Flow Overview Trend - IEEE ', num2str(num_buses), '-Bus (Baseline)'], 'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'none');
    hold off;

    % === Plot Voltage Distribution per Bus ===
    figure('Position', fig_position);
    hold on;
    
    v_min = min(Matrix_Voltages, [], 2);
    v_max = max(Matrix_Voltages, [], 2);
    v_avg = mean(Matrix_Voltages, 2);
    v_range = v_max - v_min;
    bus_numbers = 1:length(v_avg);
    
    bar_stack = bar(bus_numbers, [v_min, v_range], 'stacked');
    bar_stack(1).FaceColor = 'none';
    bar_stack(1).EdgeColor = 'none';
    bar_stack(2).FaceColor = [0.6 0.8 1];
    bar_stack(2).EdgeColor = 'none';
    
    errorbar(bus_numbers, v_avg, v_avg - v_min, v_max - v_avg, ...
        'k', 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 6, 'LineWidth', 1.5, 'CapSize', 5);
    
    yline(1.00, '--', 'Nominal', 'Color', [0.1 0.1 0.1], 'FontSize', 18, 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    yline(1.05, '--', 'Upper Limit', 'Color', [1.0 0.3 0.3], 'FontSize', 18, 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    yline(0.95, '--', 'Lower Limit', 'Color', [1.0 0.3 0.3], 'FontSize', 18, 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    
    violated = find(any(Matrix_Voltages < 0.95 | Matrix_Voltages > 1.05, 2));
    if ~isempty(violated)
        scatter(violated, v_avg(violated), 60, 'r', 'filled', 'DisplayName', 'Voltage Violation');
    end
    
    xlabel('Bus Number', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Voltage (p.u.)', 'FontSize', 18, 'FontWeight', 'bold');
    xticks(bus_numbers);
    ylim([0.90 1.10]);
    grid on;
    %set(gca, 'FontSize', 15, 'FontWeight', 'normal', 'XGrid', 'on', 'YGrid', 'on', 'Layer', 'top');
    set(gca, 'FontSize', 15);
    title(['Bus Voltage Swings Range — IEEE ', num2str(num_buses), '-Bus System'], 'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'none');
    xtickangle(90);
    hold off;

    % === Plot Power Loss, Load, and PV Output ===
    figure('Position', fig_position);
    plot(hours, Matrix_Results(:,6), '-r', 'LineWidth', 2.5); hold on;
    plot(hours, Matrix_Results(:,4), '--b', 'LineWidth', 1.8);
    plot(hours, Matrix_Results(:,2), '--', 'Color', [1 0.5 0], 'LineWidth', 1.8);
    legend('Power Loss (kW)', 'Load (kW)', 'PV Output (kW)', 'Location', 'northwest');
    xlabel('Hour', 'FontSize', 18);
    ylabel('Power (kW)', 'FontSize', 18);
    title(sprintf('Power Loss vs Load and PV – IEEE %d-Bus Distribution System', num_buses), 'FontSize', 20);
    xticks(hours);
    ylim([0 4000]);
    set(gca, 'FontSize', 15);
    grid on;

    % === Plot Voltage Heatmap ===
    figure('Position', fig_position);
    imagesc(Matrix_Voltages');
    colormap(jet);
    colorbar;
    caxis([0.9 1.1]);
    xlabel('Bus Number', 'FontSize', 18);
    ylabel('Hour', 'FontSize', 18);
    title(sprintf('Hourly Voltage Profile Heatmap (p.u.) – IEEE %d-Bus Distribution System', num_buses), 'FontSize', 20);
    set(gca, 'XTick', 1:size(Matrix_Voltages,1));
    set(gca, 'YTick', 1:24);
    xtickangle(90);
    set(gca, 'FontSize', 15);
    grid on;

    % === Plot Reverse Power Flow and Net Power ===
    figure('Position', fig_position);
    t = Matrix_Results(:, 1);
    netP = Matrix_Results(:, 5);
    pv = Matrix_Results(:, 2);
    loadP = Matrix_Results(:, 4);
    
    plot(t, netP, '-', 'Color', [0.6 0 0], 'LineWidth', 2.5, 'DisplayName', 'Net Power (kW)'); hold on;
    reverse_idx = netP < 0;
    scatter(t(reverse_idx), netP(reverse_idx), 80, 'r', 'filled', 'DisplayName', 'Reverse Power Flow');
    plot(t, pv, '--', 'Color', [1 0.5 0], 'LineWidth', 1.8, 'DisplayName', 'PV Output (kW)');
    plot(t, loadP, '--b', 'LineWidth', 1.8, 'DisplayName', 'Load (kW)');
    xlabel('Hour', 'FontSize', 18);
    ylabel('Power (kW)', 'FontSize', 18);
    title(sprintf('Net Power and Reverse Power Flow – IEEE %d-Bus Distribution System', num_buses), 'FontSize', 20);
    xticks(t);
    legend('Location', 'northwest', 'FontSize', 18);
    set(gca, 'FontSize', 15);
    grid on;
    hold off;

    figure('Position', fig_position); hold on;
    v_min = min(Matrix_Voltages, [], 1);
    v_max = max(Matrix_Voltages, [], 1);
    v_avg = mean(Matrix_Voltages, 1);
    
    plot(hours, v_min, '-xr', 'LineWidth', 2);
    plot(hours, v_max, '-ob', 'LineWidth', 2);
    plot(hours, v_avg, '-g', 'LineWidth', 2.2);
    
    yline(1.05, 'r--', 'LineWidth', 1.2);
    yline(0.95, 'r--', 'LineWidth', 1.2);
    yline(1.00, 'k-.', 'LineWidth', 1.2);
    
    xlabel('Hour', 'FontSize', 14);
    ylabel('Voltage (p.u.)', 'FontSize', 14);
    title(sprintf('Voltage Profile Trend — IEEE %d-Bus Distribution System', num_buses), 'FontSize', 20);
    ylim([0.9 1.1]);
    xticks(1:24);
    legend({'Min Voltage', 'Max Voltage', 'Avg Voltage', 'Upper Limit', 'Lower Limit', 'Nominal'}, ...
           'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 15);
    set(gca, 'FontSize', 15);
    grid on;
    hold off;

    %% === Save All Figures Automatically ===
    fig_dir = fullfile('figures_summary', sprintf('baseline_%dbus', num_buses));
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
    
    figHandles = findall(groot, 'Type', 'figure');
    figHandles = flipud(figHandles);  % Ensure save order matches creation
    
    for i = 1:length(figHandles)
        fig = figHandles(i);
        if isvalid(fig)
            fig_name = sprintf('figure_%02d.png', i);
            exportgraphics(fig, fullfile(fig_dir, fig_name), 'Resolution', 300);
            close(fig);  % Optional: auto-close after saving
        end
    end
    fprintf('\nFigures saved to folder:\n- %s\n', fig_dir);

    %% === Save Matrix_Results to CSV ===
    results_matrix_table = array2table(Matrix_Results, ...
        'VariableNames', {'Hour', 'P_Gen_(kW)', 'BESS_Demand_(kW)', 'P_Load_(kW)', ...
                          'P_Nett_(kW)', 'P_Loss_(kW)', 'Q_Loss_(kVar)', ...
                          'V_Max_(p.u.)', 'V_Bus_Max', 'V_Min_(p.u.)', 'V_Bus_Min'});
    
    % Ensure 'results_summary' folder existing
    if ~exist('results_summary', 'dir')
        mkdir('results_summary');
    end
    
    % Save .CSV File load flow results summary to 'results_summary' folder
    results_matrix_filename = fullfile('results_summary', ...
        sprintf('Load_Flow_Summary_Baseline_%dbus.csv', num_buses));
    
    writetable(results_matrix_table, results_matrix_filename);
    
    fprintf('\n Matrix_Results saved to:\n  - %s\n', results_matrix_filename);

end