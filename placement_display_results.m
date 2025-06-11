function BESS_Optimal_Location = placement_display_results(Matrix_BESS_Demand, fitness_iter_hourly, T_24h_obj,...
         T_24h_effectives_eval, T_24h_effectives_iter, opt, total_runtime, Max_volt_dev, SoC_max, SoC_min, Matrix_Results, Matrix_Voltages, mm, BESS_Number, seed)
    %% Display Optimal BESS Location and Size
    % disp('<strong>    OPTIMAL BESS LOCATION AND SIZE:</strong>');
    % Selected_Bus = Matrix_BESS_Demand(Bus_Placement, :);
    % BESS_max_out = max(abs(Selected_Bus), [], 2);
    % BESS_max_cap = (max(cumsum(Selected_Bus, 2), [], 2) - min(cumsum(Selected_Bus, 2), [], 2)) / (SoC_max - SoC_min);
    % %BESS_max_cap = abs(sum(Selected_Bus, 2)) / (SoC_max - SoC_min);
    % pz_size = table(Bus_Placement', BESS_max_out, BESS_max_cap, ...
    %     'VariableNames', {'Bus_Number', 'Max_Output_(kW)', 'Max_Capacity_(kWh)'});
    % disp(pz_size);
    num_buses = size(mm, 1);
    % Display Data Results
    disp('<strong>    DATA RESULTS:</strong>');
    column_names =  {'Hour', 'P_Gen_(kW)', 'BESS_Demand_(kW)', 'P_Load_(kW)', 'P_Nett_(kW)', ...
                    'P_Loss_(kW)', 'Q_Loss_(kVar)', 'V_Max_(p.u.)', 'V_Bus_Max', 'V_Min_(p.u.)', 'V_Bus_Min'};
    results_table = array2table(Matrix_Results, 'VariableNames', column_names);
    disp(results_table);

    %% Calculate and Display Statistics
    percentage_pv = (sum(mm(:, 4)) / sum(mm(:, 2))) * 100;
    %percentage_losses = (sum(Matrix_Results(:, 6)) / 24) / sum(mm(:, 4)) * 100;
    percentage_losses = (sum(Matrix_Results(:, 6)) / sum(Matrix_Results(:, 4))) * 100;
    percentage_v_dev = Max_volt_dev / 1.0 * 100;
    BESS_Ratio = (sum(abs(Matrix_Results(:, 3)))) / (sum(Matrix_Results(:, 2)));
    pv_vs_load = ((sum(Matrix_Results(:, 2))) / (sum(Matrix_Results(:, 4)))) * 100;
    bess_vs_load = ((sum(abs(Matrix_Results(:, 3)))) / (sum(abs(Matrix_Results(:, 4))))) * 100;
    disp('<strong>24 HOUR LOAD FLOW SUMMARY:</strong>');
    fprintf('Total Load               : %.2f  kW\n', sum(mm(:, 2)));
    fprintf('Total Load Consumption   : %.2f kWh\n', (sum(Matrix_Results(:, 4))));
    fprintf('Total Installed PV       : %.2f  kW  (%.2f %%)\n', sum(mm(:, 4)), percentage_pv);
    fprintf('Total PV Supply          : %.2f kWh (%.2f %%)\n', (sum(Matrix_Results(:, 2))), pv_vs_load);
    fprintf('BESS Total Abs Supply    : %.2f kWh (%.2f %%)\n', (sum(abs(Matrix_Results(:, 3)))), bess_vs_load);
    fprintf('PV:BESS Supply Ratio = 1 : %.2f \n', BESS_Ratio);
    fprintf('Average P Losses         : %.2f kW   (%.2f %%)\n', (sum(Matrix_Results(:, 6)) / 24), percentage_losses);
    fprintf('Max Bus Voltage Deviation: %.3f p.u.  (%.2f %%)\n', Max_volt_dev, percentage_v_dev);
    fprintf('Average V. Deviation     : %.3f p.u.  (%.2f %%)\n', mean(mean(abs(1 - Matrix_Voltages))), (mean(mean(abs(1 - Matrix_Voltages)))/1.0)*100);
    fprintf('Max Bus Voltage          : %.3f p.u.\n', max(Matrix_Voltages(:)));
    fprintf('Min Bus Voltage          : %.3f p.u.\n', min(Matrix_Voltages(:)));
    fprintf('\n24 Hour Total Objectives Value     : %.2f', T_24h_obj); % Print the final objective value
    fprintf('\n24 Hour Total Effectives Iteration : %d', T_24h_effectives_iter); 
    fprintf('\n24 Hour Total Effectives Eval Count: %d \n', T_24h_effectives_eval); 

    % Sum absolute values of BESS demand for each bus over 24 hours
    BESS_Absolute_Sum = sum(abs(Matrix_BESS_Demand), 2); % Sum across 24 hours
    
    % Find the top sort buses with the highest absolute demand
    [Top_Values, Top_Indices] = maxk(BESS_Absolute_Sum, BESS_Number);
    
    % Sort buses in ascending order
    [Sorted_Bus, Sort_Index] = sort(Top_Indices); % Sort bus numbers
    Sorted_Values = Top_Values(Sort_Index);       % Rearrange values accordingly
    
    % Store result as a matrix
    BESS_Optimal_Location = Sorted_Bus';
    
    % Display the sorted results
    disp('Buses with Highest Absolute BESS Demand (Sorted by Bus Number):');
    disp(BESS_Optimal_Location);

    %% PLOT TREND
    % === PLOT OVERVIEW TREND ===
    hours = Matrix_Results(:, 1)';
    fig = figure;
    set(fig, 'Units', 'pixels', 'Position', [100, 100, 1920, 1080]);  % Full HD
    
    % Left y-axis: Power
    yyaxis left
    p1 = plot(hours, Matrix_Results(:, 4), '--', 'Color', [0.0, 0.0, 1.0], 'LineWidth', 3, 'DisplayName', 'Total P Load (kW)');
    hold on;
    p2 = plot(hours, Matrix_Results(:, 2), '--', 'Color', [1.0, 0.5, 0.0], 'LineWidth', 3, 'DisplayName', 'Total PV Gen (kW)');
    p3 = plot(hours, Matrix_Results(:, 6), '-v', 'Color', [1.0, 0.0, 0.0], 'LineWidth', 3, 'DisplayName', 'Total Losses (kW)');
    p4 = plot(hours, Matrix_Results(:, 5), '-o', 'Color', [0.6, 0.0, 0.0], 'LineWidth', 3, 'DisplayName', 'Total P Net (kW)');
    p5 = plot(hours, Matrix_Results(:, 3), '-*', 'Color', [0.2, 0.7, 0.3], 'LineWidth', 3, 'DisplayName', 'Total BESS Output (kW)');
    ylabel('Power (kW)', 'FontSize', 18, 'FontWeight', 'bold');
    xlabel('Hour', 'FontSize', 18, 'FontWeight', 'bold');
    ylim([-3000 7500]);
    grid on;
    
    % Right y-axis: Voltage
    yyaxis right
    p6 = plot(hours, Matrix_Results(:, 8), '-+', 'Color', [0.0, 0.9, 0.7], 'LineWidth', 3, 'DisplayName', 'V Max (p.u.)'); 
    p7 = plot(hours, Matrix_Results(:, 10), '-x', 'Color', [0.7, 0.9, 0.0], 'LineWidth', 3, 'DisplayName', 'V Min (p.u.)'); 
    ylabel('Voltage (p.u.)', 'FontSize', 18, 'FontWeight', 'bold');
    ylim([0.35 1.15]);
    grid on;
    
    xticks(hours);
    set(gca, 'FontSize', 15, 'FontWeight', 'bold');  % Axis tick font
    legend([p1, p2, p3, p4, p5, p6, p7], 'Location', 'southoutside', ...
           'FontSize', 16, 'NumColumns', 2);
    
    % Title with opt method
    title(['24-Hour Load Flow Overview Trend - IEEE ', num2str(num_buses), ...
           '-Bus | Optimizer: ', opt], 'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'none');
    
    hold off;
    
    % === BUS VOLTAGE DISTRIBUTION  ===
    fig = figure;
    set(fig, 'Units', 'pixels', 'Position', [100, 100, 1920, 1080], 'Color', 'w');
    hold on;
    
    % Calculate voltage stats for each bus
    v_min = min(Matrix_Voltages, [], 2);
    v_max = max(Matrix_Voltages, [], 2);
    v_avg = mean(Matrix_Voltages, 2);
    v_range = v_max - v_min;
    bus_numbers = 1:length(v_avg);
    
    % Plot voltage range using stacked bars (transparent base + shaded range)
    bar_stack = bar(bus_numbers, [v_min, v_range], 'stacked');
    bar_stack(1).FaceColor = 'none';  % Transparent base
    bar_stack(1).EdgeColor = 'none';
    bar_stack(2).FaceColor = [0.6 0.8 1];  % Soft blue for voltage range
    bar_stack(2).EdgeColor = 'none';
    
    % Plot average voltage with error bars
    errorbar(bus_numbers, v_avg, v_avg - v_min, v_max - v_avg, ...
        'k', 'LineStyle', 'none', 'Marker', 'o', ...
        'MarkerSize', 6, 'LineWidth', 1.5, 'CapSize', 5);
    
    % Add voltage reference lines (nominal ±5%)
    yline(1.00, 'Color', [0.1 0.1 0.1], 'LineWidth', 1.5, 'LineStyle', '--', 'Label', 'Nominal', 'LabelHorizontalAlignment', 'left','FontSize', 14);
    yline(1.05, 'Color', [1.0 0.3 0.3], 'LineWidth', 1.5, 'LineStyle', '--', 'Label', 'Upper Limit', 'LabelHorizontalAlignment', 'left','FontSize', 14);
    yline(0.95, 'Color', [1.0 0.3 0.3], 'LineWidth', 1.5, 'LineStyle', '--', 'Label', 'Lower Limit', 'LabelHorizontalAlignment', 'left','FontSize', 14);
    
    % Highlight buses with violations (below 0.95 or above 1.05 at any hour)
    violated = find(any(Matrix_Voltages < 0.95 | Matrix_Voltages > 1.05, 2));
    if ~isempty(violated)
        scatter(violated, v_avg(violated), 60, 'r', 'filled', 'DisplayName', 'Voltage Violation');
    end
    
    % Axis settings
    xlabel('Bus Number', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Voltage (p.u.)', 'FontSize', 18, 'FontWeight', 'bold');
    xticks(bus_numbers);
    ylim([0.90 1.10]);
    grid on;
    set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'XGrid', 'on', 'YGrid', 'on', 'Layer', 'top');
    
    % Add title with optimization method
    title(['Bus Voltage Profile — IEEE ', num2str(num_buses), ...
           '-Bus System | Optimizer: ', opt], ...
           'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'none');
    
    hold off;
    
    % === BESS DEMAND TREND: Active Buses Only Over 24 Hours ===
    fig = figure;
    set(fig, 'Units', 'pixels', 'Position', [100, 100, 1920, 1080], 'Color', 'w');
    hold on;
    
    [num_buses, ~] = size(Matrix_BESS_Demand);
    legend_entries = {};  % To store bus labels for active buses
    plotted_lines = [];   % To store plot handles for legend
    
    for bus = 1:num_buses
        demand = Matrix_BESS_Demand(bus, :);
        if any(demand ~= 0)  % Only plot if the bus has non-zero demand
            p = plot(1:24, demand, ...
                'LineWidth', 1.8, ...
                'DisplayName', ['Bus ', num2str(bus)]);
            plotted_lines(end+1) = p;
            legend_entries{end+1} = ['Bus ', num2str(bus)];
        end
    end
    
    % Axis labeling
    xlabel('Hour', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('BESS Demand (kW)', 'FontSize', 14, 'FontWeight', 'bold');
    title(['BESS Demand Trend Over 24 Hours | Optimizer: ', opt], ...
          'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');
    
    % Axis settings
    xlim([1, 24]);
    ylim([-950, 950]);
    xticks(1:24);
    grid on;
    
    % Show legend only for active buses
    legend(plotted_lines, legend_entries, ...
           'Location', 'eastoutside', 'FontSize', 10);
    
    set(gca, 'FontSize', 12, 'FontWeight', 'normal', ...
             'XGrid', 'on', 'YGrid', 'on');
    hold off;

    % === PLOT: Total Absolute BESS Demand per Bus) ===
    fig = figure;
    set(fig, 'Units', 'pixels', 'Position', [100, 100, 1920, 1080], 'Color', 'w');
    
    % Define bar colors: blue for default, red for selected buses
    bar_colors = repmat([0.2, 0.6, 1], length(bus_numbers), 1); % Blue
    selected_buses = BESS_Optimal_Location;
    selected_indices = ismember(bus_numbers, selected_buses);
    bar_colors(selected_indices, :) = repmat([1, 0.2, 0.2], sum(selected_indices), 1); % Red
    
    % Create the bar chart
    b = bar(bus_numbers, BESS_Absolute_Sum, 'FaceColor', 'flat', 'BarWidth', 0.75);
    b.CData = bar_colors;
    
    % Axis labeling
    xlabel('Bus Number', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Total Absolute BESS Demand (kW)', 'FontSize', 14, 'FontWeight', 'bold');
    title(['Total Absolute BESS Demand per Bus (24 Hours) | Optimizer: ', opt], ...
          'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'none');
    
    % Axis settings
    xlim([min(bus_numbers) - 1, max(bus_numbers) + 1]);
    ylim([0, max(BESS_Absolute_Sum) * 1.1]);
    xticks(bus_numbers);
    set(gca, 'FontSize', 12, 'XGrid', 'on', 'YGrid', 'on', 'Layer', 'top');
    
    % Add value labels on top of each bar
    for i = 1:length(bus_numbers)
        val = BESS_Absolute_Sum(i);
        if val > 0
            text(bus_numbers(i), val + 0.02 * max(BESS_Absolute_Sum), ...
                num2str(round(val, 1)), ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 10, 'Color', 'k');
        end
    end
    
    % Add legend only if selected buses exist
    if any(selected_indices)
        hold on;
        dummy = bar(nan, nan, 'FaceColor', [1, 0.2, 0.2]);
        legend(dummy, 'Projected BESS Location', 'Location', 'northeast', ...
               'FontSize', 14);
        hold off;
    end


    % === PLOT: Fitness Heatmap (Generation vs Hour) ===
    fig = figure;
    set(fig, 'Units', 'pixels', 'Position', [100, 100, 1920, 1080], 'Color', 'w');
    
    imagesc(fitness_iter_hourly);
    caxis([0 40]); % max scale at 50
    colormap(jet(32));
    colorbar;
    
    % Axis labeling
    xlabel('Generation', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Hour', 'FontSize', 14, 'FontWeight', 'bold');
    title(['Fitness Heatmap (Iterations vs Hours) | Optimizer: ', opt], ...
          'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');
    
    % Tick settings
    [num_hours, num_generations] = size(fitness_iter_hourly);
    set(gca, 'YTick', 1:num_hours);  % One tick per hour
    
    % Adjust x-ticks
    tick_interval = max(1, round(num_generations / 10));
    x_tick_vals = tick_interval:tick_interval:num_generations;
    set(gca, 'XTick', x_tick_vals);
    
    % Improve axis look
    set(gca, 'FontSize', 12, 'FontWeight', 'normal');
    axis tight;

    % Define figure resolution
    fig_position = [100, 100, 1920, 1080];
    
    % === POWER LOSS vs LOAD and PV ===
    figure('Position', fig_position);
    plot(hours, Matrix_Results(:,6), '-r', 'LineWidth', 2.5); hold on; % power loss
    plot(hours, Matrix_Results(:,4), '--b', 'LineWidth', 1.8); % load profile
    plot(hours, Matrix_Results(:,2), '--', 'Color', [1 0.5 0], 'LineWidth', 1.8); % pv out
    legend('Power Loss (kW)', 'Load (kW)', 'PV Output (kW)', 'Location', 'northwest','FontSize', 14);
    xlabel('Hour', 'FontSize', 14);
    ylabel('Power (kW)', 'FontSize', 14);
    title(sprintf('Power Loss vs Load and PV – IEEE %d-Bus Distribution System', num_buses), 'FontSize', 16);
    xticks(hours);
    ylim([0 4000]);
    grid on;
    
    % === VOLTAGE HEATMAP ===
    figure('Position', fig_position);
    imagesc(Matrix_Voltages');
    colormap(jet);
    colorbar;
    caxis([0.9 1.1]);
    xlabel('Bus Number', 'FontSize', 14);
    ylabel('Hour', 'FontSize', 14);
    title(sprintf('Hourly Voltage Profile Heatmap (p.u.) – IEEE %d-Bus Distribution System', num_buses), 'FontSize', 16);
    set(gca, 'XTick', 1:size(Matrix_Voltages,1));
    set(gca, 'YTick', 1:24);
    grid on;
    
    % === REVERSE POWER FLOW OVERLAY on NET POWER ===
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
    xlabel('Hour', 'FontSize', 14);
    ylabel('Power (kW)', 'FontSize', 14);
    title(sprintf('Net Power and Reverse Power Flow – IEEE %d-Bus Distribution System', num_buses), 'FontSize', 16);
    xticks(t);
    legend('Location', 'northwest', 'FontSize', 14);
    grid on;
    hold off;
    
    % === VOLTAGE PROFILE TREND with LIMIT SHADING ===
    figure('Position', fig_position); hold on;
    v_min = min(Matrix_Voltages, [], 1);
    v_max = max(Matrix_Voltages, [], 1);
    v_avg = mean(Matrix_Voltages, 1);
    
    plot(hours, v_min, '-xr', 'LineWidth', 2, 'DisplayName', 'Min Voltage (p.u.)');
    plot(hours, v_max, '-ob', 'LineWidth', 2, 'DisplayName', 'Max Voltage (p.u.)');
    plot(hours, v_avg, '-g', 'LineWidth', 2.2, 'DisplayName', 'Avg Voltage (p.u.)');
    
    yline(1.05, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Upper Limit (1.05 p.u.)');
    yline(0.95, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Lower Limit (0.95 p.u.)');
    yline(1.00, 'k-.', 'LineWidth', 1.2, 'DisplayName', 'Nominal Voltage (1.00 p.u.)');
    
    xlabel('Hour', 'FontSize', 14);
    ylabel('Voltage (p.u.)', 'FontSize', 14);
    title(sprintf('Voltage Profile Trend – IEEE %d-Bus Distribution System', num_buses), 'FontSize', 16);
    ylim([0.9 1.1]);
    xticks(1:24);
    grid on;
    legend('Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 12);
    hold off;

    %% === Save Summary Data ===
    % Create folder if it doesn't exist
    results_dir = 'results_summary';
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end
    
    % Create result structure
    Summary = struct();
    Summary.Objective_Value           = T_24h_obj;
    Summary.Total_Load                = sum(mm(:, 2));
    Summary.Total_Load_Consumption    = sum(Matrix_Results(:, 4));
    Summary.Total_PV_Installed        = sum(mm(:, 4));
    Summary.Total_PV_Supply           = sum(Matrix_Results(:, 2));
    Summary.Total_BESS_Demand         = sum(abs(Matrix_Results(:, 3)));
    Summary.BESS_to_Load_Percent      = bess_vs_load;
    Summary.PV_to_Load_Percent        = pv_vs_load;
    Summary.P_Loss_Avg                = sum(Matrix_Results(:, 6)) / 24;
    Summary.P_Loss_Percent            = percentage_losses;
    Summary.Avg_Volt_Dev              = mean(mean(abs(1 - Matrix_Voltages)));
    Summary.Max_Voltage_Dev           = Max_volt_dev;
    Summary.Voltage_Dev_Percent       = percentage_v_dev;
    Summary.Max_Voltage               = max(Matrix_Voltages(:));
    Summary.Min_Voltage               = min(Matrix_Voltages(:));
    Summary.Total_Effective_Iter      = T_24h_effectives_iter;
    Summary.Total_Effective_Eval      = T_24h_effectives_eval;
    Summary.BESS_Location             = join(string(BESS_Optimal_Location), ' | ');
    Summary.Total_Selected_BESS_Demand_per_Bus = join(string(Sorted_Values), ' | ');
    Summary.Total_Runtime_Seconds     = total_runtime;
    Summary.seed                      = seed;
    Summary.Timestamp                 = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    
    % Save as table
    summary_table = struct2table(Summary);
    csv_name = sprintf('summary_optimal_placement__%s_%dbus_%s.csv', opt, num_buses, datestr(now, 'yyyymmdd_HHMMSS'));
    csv_path = fullfile(results_dir, csv_name);
    writetable(summary_table, csv_path);

    fprintf('\nSummary data saved to:\n  -> %s\n', csv_path);

    %% === Save All Figures Automatically ===
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    fig_dir = fullfile('figures_summary', sprintf('Placement_%s_%dbus_%s', opt, num_buses, timestamp));
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    figHandles = findall(groot, 'Type', 'figure');
    figHandles = flipud(figHandles);  % Ensure correct order

    for i = 1:length(figHandles)
        fig = figHandles(i);
        if isgraphics(fig, 'figure')
            fig_name = fullfile(fig_dir, sprintf('figure_%02d.png', i));
            try
                saveas(fig, fig_name);
                close(fig);  % opsional: close figure after saved
            catch ME
                fprintf('Warning: Failed to save figure %d (%s)\n', i, ME.message);
            end
        end
    end

    fprintf('All figures saved to:\n  -> %s\n', fig_dir);
end