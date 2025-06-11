function display_results(Matrix_BESS_Demand, Bus_Placement, obj, Max_volt_dev, SoC_max, SoC_min, ...
                         Matrix_Results, Matrix_Voltages, mm, lower_bound, upper_bound, BESS_Output, opt, total_runtime, iter, eval_total, seed)
    %% Display Optimal BESS Location and Size
    disp('<strong>    OPTIMAL BESS LOCATION AND SIZE:</strong>');
    Selected_Bus = Matrix_BESS_Demand(Bus_Placement, :);
    BESS_max_out = max(abs(Selected_Bus), [], 2);
    BESS_max_cap = (max(cumsum(Selected_Bus, 2), [], 2) - min(cumsum(Selected_Bus, 2), [], 2)) / (SoC_max - SoC_min);
    BESS_Cap_ratio = ((sum(BESS_max_cap))/(sum(Matrix_Results(:, 4))))*100;
    pz_size = table(Bus_Placement', BESS_max_out, BESS_max_cap, ...
        'VariableNames', {'Bus_Number', 'Max_Output_(kW)', 'Optimal_Capacity_(kWh)'});
    disp(pz_size);
    percent_demand = (cumsum(-Matrix_BESS_Demand(Bus_Placement, :), 2) ./ BESS_max_cap)*100;
    max_values = abs((max(percent_demand, [], 2))-(SoC_max*100));
    balance_OP = percent_demand + max_values;
    
    % Average unbalance
    charging  = sum(-min(BESS_Output, 0), 2);
    discharge = sum(max(BESS_Output, 0), 2);
    unbalance = abs(charging - discharge) / (charging + discharge + eps) * 100;
    ebp = mean (unbalance);
    Unbalance_Percent = max(ebp);
    
    % Rata-rata unbalance antar unit
    Unbalance_kW = mean(charging - discharge);

    % Display Data Results
    disp('<strong>    DATA RESULTS:</strong>');
    column_names =  {'Hour', 'P_Gen_(kW)', 'BESS_Demand_(kW)', 'P_Load_(kW)', 'P_Nett_(kW)', ...
                    'P_Loss_(kW)', 'Q_Loss_(kVar)', 'V_Max_(p.u.)', 'V_Bus_Max', 'V_Min_(p.u.)', 'V_Bus_Min'};
    results_table = array2table(Matrix_Results, 'VariableNames', column_names);
    disp(results_table);
    num_buses = size(mm, 1);

    %% Calculate and Display Statistics
    percentage_pv = (sum(mm(:, 4)) / sum(mm(:, 2))) * 100;
    %percentage_losses = (sum(Matrix_Results(:, 6)) / 24) / sum(mm(:, 4)) * 100;
    percentage_losses = (sum(Matrix_Results(:, 6)) / sum(Matrix_Results(:, 4))) * 100;
    percentage_v_dev = Max_volt_dev / 1.0 * 100;
    BESS_Sup_Ratio = (sum(abs(Matrix_Results(:, 3)))) / (sum(Matrix_Results(:, 2)));
    pv_vs_load = ((sum(Matrix_Results(:, 2))) / (sum(Matrix_Results(:, 4)))) * 100;
    bess_vs_load = ((sum(abs(Matrix_Results(:, 3)))) / (sum(abs(Matrix_Results(:, 4))))) * 100;
    
    % Change rate penalty per BESS unit (cyclic 24 hour)
    BESS_shifted = [BESS_Output(:,24), BESS_Output(:,1:23)];
    delta_output = abs(diff(BESS_shifted, 1, 2));
    avg_Demand_rate = mean(mean(delta_output, 2));
    max_Demand_rate = max(max(delta_output, [], 2));
    
    % Tampilkan hasil
    
    %fprintf('by BESS on bus #%d between hour %d and %d\n', unit_idx, hour_idx, hour_idx+1);
    disp('<strong>24 HOUR LOAD FLOW SUMMARY:</strong>');
    fprintf('Total Load Capacity     : %.2f  kW\n', sum(mm(:, 2)));
    fprintf('Total Load Consumption  : %.2f kWh\n', (sum(Matrix_Results(:, 4))));
    fprintf('Total Total PV Capacity : %.2f  kW  (%.2f %%)\n', sum(mm(:, 4)), percentage_pv);
    fprintf('Total Total PV Supply   : %.2f kWh (%.2f %%)\n', (sum(Matrix_Results(:, 2))), pv_vs_load);
    fprintf('BESS Total Capacity     : %.2f kWh (%.2f %%)\n', sum(BESS_max_cap), BESS_Cap_ratio);
    fprintf('BESS Total Abs Supply   : %.2f kWh (%.2f %%)\n', (sum(abs(Matrix_Results(:, 3)))), bess_vs_load);
    fprintf('BESS Demand Unbalance   : %.2f kWh    (%.2f %%)\n', Unbalance_kW, Unbalance_Percent);
    fprintf('Max Ramp Rate           : %.2f kW/min \n', (max_Demand_rate/60));
    fprintf('Average Ramp Rate       : %.2f kW/min \n', (avg_Demand_rate/60));
    fprintf('PV:BESS Ratio =       1 : %.2f kWh\n', BESS_Sup_Ratio);
    fprintf('Average P. Losses       : %.2f kW   (%.2f %%)\n', (sum(Matrix_Results(:, 6)) / 24), percentage_losses);
    fprintf('Average V. Deviation    : %.3f p.u.  (%.2f %%)\n', mean(mean(abs(1 - Matrix_Voltages))), (mean(mean(abs(1 - Matrix_Voltages)))/1.0)*100);
    fprintf('Max Voltage Deviation   : %.3f p.u.  (%.2f %%)\n', Max_volt_dev, percentage_v_dev);
    fprintf('Max Bus Voltage         : %.3f p.u.\n', max(Matrix_Voltages(:)));
    fprintf('Min Bus Voltage         : %.3f p.u.\n', min(Matrix_Voltages(:)));

    fprintf('\n24 Hour Total Effectives Iterations : %.2f \n', iter); % Print Total Effectives Evaluation
    fprintf('24 Hour Total Effectives Evaluations: %.2f \n', eval_total); % Print Total Effectives Evaluation
    fprintf('24 Hour Total Objectives Value: %.2f \n', obj); % Print Total Objectives Value

    %% PLOT TREND
    % === COMBINED PLOT: BESS Demand (Left) & SoC (Right with Clean Rotated Labels) ===
    fig = figure;
    set(fig, 'Units', 'pixels', 'Position', [100, 100, 1920, 1080], 'Color', 'w');
    
    num_BESS = length(Bus_Placement);
    
    for i = 1:num_BESS
        % --- DEMAND subplot (left column) ---
        subplot(num_BESS, 2, 2*i - 1);
        bar(1:24, Matrix_BESS_Demand(Bus_Placement(i), :), ...
            'FaceColor', [0.2 0.6 0.9], 'EdgeColor', 'k');
        title(sprintf('BESS #%d (Bus %d) — Demand', i, Bus_Placement(i)), ...
              'FontSize', 14, 'FontWeight', 'bold');
        ylabel('Demand (kW)', 'FontSize', 12, 'FontWeight', 'bold');
        ylim([-1000, 1000]);
        yticks(-1000:250:1000);
        xticks(1:24);
        yline(0, 'k--');
        grid on;
        if i == num_BESS
            xlabel('Hour', 'FontSize', 12, 'FontWeight', 'bold');
        else
            set(gca, 'XTickLabel', []);
        end
    
        % --- SoC subplot (right column) ---
        subplot(num_BESS, 2, 2*i);
        soc_data = balance_OP(i, :);
        b = bar(1:24, soc_data, 'FaceColor', [0.2 0.8 0.3], 'EdgeColor', 'k');
        title(sprintf('BESS #%d (Bus %d) — SoC', i, Bus_Placement(i)), ...
              'FontSize', 14, 'FontWeight', 'bold');
        ylabel('SoC (%)', 'FontSize', 12, 'FontWeight', 'bold');
        ylim([0 100]);
        yticks(0:20:100);
        xticks(1:24);
        grid on;
    
        % === Show SoC values inside the bars with vertical orientation ===
        xtips = b.XData;
        ytips = b.YData;
        labels = string(round(ytips, 1));
        for j = 1:length(xtips)
            % Only display if height is enough
            if ytips(j) > 10
                text(xtips(j), ytips(j) - 12, labels(j), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontSize', 8, ...       % Smaller font
                    'Rotation', 90, ...      % Vertical label
                    'Color', 'k');           % Inside the bar
            end
        end
    
        if i == num_BESS
            xlabel('Hour', 'FontSize', 12, 'FontWeight', 'bold');
        else
            set(gca, 'XTickLabel', []);
        end
    end
    
    % === Main title ===
    sgtitle(['BESS Demand and SoC per Unit (24 Hours) for IEEE ', num2str(num_buses), ...
            '-Bus | Optimizer: ', opt],'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'none');

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
    title(['24-Hour Load Flow Overview Trend IEEE ', num2str(num_buses), ...
           '-Bus | Optimizer: ', opt], 'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'none');
    
    hold off;
     
    % === PLOT DEMAND TREND ===
    fig = figure;
    set(fig, 'Units', 'pixels', 'Position', [100, 100, 1920, 1080]);  % Set size to 1920x1080
    hold on;
    
    % Plot each BESS demand
    for i = 1:length(Bus_Placement)
        bus = Bus_Placement(i);
        plot(1:24, Matrix_BESS_Demand(bus, :), 'LineWidth', 3, ...
             'DisplayName', ['BESS #', num2str(i), ' (Bus ', num2str(bus), ')']);
    end
    
    % Labels and axis
    xlabel('Hour', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('BESS Demand (kW)', 'FontSize', 18, 'FontWeight', 'bold');
    ylim([-1000, 1000]);
    xlim([1 24]);
    xticks(1:24);
    grid on;
    
    num_buses = size(Matrix_BESS_Demand, 1);

    % Title and legend
    title(sprintf('BESS Demand Trend for IEEE %d Bus | Optimizer: %s', num_buses, opt), ...
          'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'none');
    
    legend('show', 'Location', 'southoutside', 'FontSize', 16);
    set(gca, 'FontSize', 15, 'FontWeight', 'bold');  % Tick font
    
    hold off;

    % === GROUPED BESS DEMAND SCHEDULE ===
    fig = figure;
    set(fig, 'Units', 'pixels', 'Position', [100, 100, 1920, 1080], 'Color', 'w');
    
    % Extract BESS demand data for selected buses
    selected_BESS_Demand = Matrix_BESS_Demand(Bus_Placement, :);
    num_BESS = size(selected_BESS_Demand, 1);
    
    % Calculate symmetric Y-axis limit (rounded up to nearest multiple of 200)
    y_max = ceil(max(abs(selected_BESS_Demand(:))) / 200) * 200;
    ytick_vals = -y_max:200:y_max;
    
    % Create grouped bar chart (transposed for hours as categories)
    b = bar(1:24, selected_BESS_Demand', 'grouped');
    colormap(parula(num_BESS));  % Apply a soft color map for better readability
    
    % Axis labels
    xlabel('Hour', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('BESS Demand (kW)', 'FontSize', 18, 'FontWeight', 'bold');
    
    % Set axis limits and ticks
    ylim([-y_max, y_max]);
    yticks(ytick_vals);
    xticks(1:24);
    xlim([0.4 24.6]);  % Add horizontal margin to prevent edge clipping
    grid on;
    
    % Style the grid
    set(gca, 'GridAlpha', 0.3, 'GridLineStyle', '--');  % Light dashed grid
    yline(0, 'k--', 'LineWidth', 1.5);  % Zero reference line
    
    % Add vertical guide lines between hours
    hold on;
    for h = 1.5:1:23.5
        xline(h, ':k', 'LineWidth', 0.8);  % Thin dotted separator lines
    end
    hold off;
    
    % Create legend showing each BESS with its bus number
    legend(arrayfun(@(i) sprintf('BESS #%d (Bus %d)', i, Bus_Placement(i)), ...
           1:num_BESS, 'UniformOutput', false), ...
           'Location', 'southoutside', 'FontSize', 14, ...
           'NumColumns', min(4, num_BESS));  % Use up to 4 columns for space efficiency
    
    % Add descriptive figure title with optimization method
    title(['Grouped BESS Demand Schedule (24 Hours) | Optimizer: ', opt], ...
          'FontSize', 20, 'FontWeight', 'bold', 'Interpreter', 'none');
    
    % Set font size for axes ticks
    set(gca, 'FontSize', 15, 'FontWeight', 'bold');
    
    % === BUS VOLTAGE DISTRIBUTION ===
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

    % Define figure resolution
    fig_position = [100, 100, 1920, 1080];
    
    % === POWER LOSS vs LOAD and PV ===
    figure('Position', fig_position);
    plot(hours, Matrix_Results(:,6), '-r', 'LineWidth', 2.5); hold on;
    plot(hours, Matrix_Results(:,4), '--b', 'LineWidth', 1.8);
    plot(hours, Matrix_Results(:,2), '--', 'Color', [1 0.5 0], 'LineWidth', 1.8);
    legend('Power Loss (kW)', 'Load (kW)', 'PV Output (kW)', 'Location', 'northwest','FontSize', 14);
    xlabel('Hour', 'FontSize', 14);
    ylabel('Power (kW)', 'FontSize', 14);
    title(['Power Loss vs Load and PV — IEEE ', num2str(num_buses), '-Bus | Optimizer: ', opt], ...
    'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');
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
    title(['Voltage Heatmap — IEEE ', num2str(num_buses), '-Bus | Optimizer: ', opt], ...
    'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');
    set(gca, 'XTick', 1:size(Matrix_Voltages,1));
    set(gca, 'YTick', 1:24);
    grid on;
    
    % === REVERSE POWER FLOW OVERLAY ===
    figure('Position', fig_position);
    t = Matrix_Results(:, 1);
    netP = Matrix_Results(:, 5);
    pv = Matrix_Results(:, 2);
    loadP = Matrix_Results(:, 4);
    
    plot(t, netP, '-', 'Color', [0.6 0 0], 'LineWidth', 2.5); hold on;
    reverse_idx = netP < 0;
    scatter(t(reverse_idx), netP(reverse_idx), 80, 'r', 'filled');
    plot(t, pv, '--', 'Color', [1 0.5 0], 'LineWidth', 1.8);
    plot(t, loadP, '--b', 'LineWidth', 1.8);
    xlabel('Hour', 'FontSize', 14);
    ylabel('Power (kW)', 'FontSize', 14);
    title(['Net Power and Reverse Flow — IEEE ', num2str(num_buses), '-Bus | Optimizer: ', opt], ...
    'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');
    xticks(t);
    legend('Net Power', 'Reverse Flow', 'PV Output', 'Load', 'Location', 'northwest', 'FontSize', 14);
    grid on;
    hold off;
    
    % === VOLTAGE PROFILE TREND ===
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
    title(['Voltage Profile Trend — IEEE ', num2str(num_buses), '-Bus | Optimizer: ', opt], ...
    'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');
    ylim([0.9 1.1]);
    xticks(1:24);
    legend({'Min Voltage', 'Max Voltage', 'Avg Voltage', 'Upper Limit', 'Lower Limit', 'Nominal'}, ...
           'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 12);
    grid on;
    hold off;

     %% === Save Summary Data ===
    ResultsSummary = struct();
    ResultsSummary.Objective_Value     = obj;
    ResultsSummary.Total_Load          = sum(mm(:,2));
    ResultsSummary.Total_PV            = sum(Matrix_Results(:,2));
    ResultsSummary.Total_BESS          = sum(abs(Matrix_Results(:,3)));
    ResultsSummary.PV_to_Load_Percent  = pv_vs_load;
    ResultsSummary.BESS_Sup_Percent    = bess_vs_load;
    ResultsSummary.Unbalance_Percent   = Unbalance_Percent;
    ResultsSummary.Unbalance_kW        = Unbalance_kW;
    ResultsSummary.BESS_Cap_Percent    = BESS_Cap_ratio;
    ResultsSummary.Ploss_Avg           = sum(Matrix_Results(:,6)) / 24;
    ResultsSummary.Ploss_Percent       = percentage_losses;
    ResultsSummary.Max_Volt_Dev        = Max_volt_dev;
    ResultsSummary.Avg_Volt_Dev        = mean(mean(abs(1 - Matrix_Voltages)));
    ResultsSummary.Ramp_Rate_Max       = max_Demand_rate / 60;
    ResultsSummary.Ramp_Rate_Avg       = avg_Demand_rate / 60;
    ResultsSummary.Bus_Placement       = Bus_Placement(:)';
    ResultsSummary.BESS_Max_Capacity   = BESS_max_cap(:)';
    ResultsSummary.BESS_Max_Output     = BESS_max_out(:)';
    ResultsSummary.Total_Runtime_Seconds= total_runtime;
    ResultsSummary.Iter_Total          = iter;
    ResultsSummary.Eval_Total          = eval_total;
    ResultsSummary.seed                = seed;
    ResultsSummary.Time_stamp          = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    % Save folder
    results_dir = 'results_summary';
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end

    % Use optimization method and timestamp as filename
    num_buses = size(mm, 1);
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    csv_filename = fullfile(results_dir, sprintf('summary_optimal_sizing_%s_%dbus_%s.csv', opt, num_buses, timestamp));

    % Convert to table and save as .csv
    SummaryTable = struct2table(ResultsSummary);
    writetable(SummaryTable, csv_filename);
    
    fprintf('\n Results data saved to:\n  - %s\n', csv_filename);

    %% === Save Matrix_Results to CSV ===
    results_matrix_table = array2table(Matrix_Results, ...
        'VariableNames', {'Hour', 'P_Gen_(kW)', 'BESS_Demand_(kW)', 'P_Load_(kW)', ...
                          'P_Nett_(kW)', 'P_Loss_(kW)', 'Q_Loss_(kVar)', ...
                          'V_Max_(p.u.)', 'V_Bus_Max', 'V_Min_(p.u.)', 'V_Bus_Min'});
    
    results_matrix_filename = fullfile(results_dir, ...
        sprintf('summary_matrix_results_%s_%dbus_%s.csv', opt, num_buses, timestamp));
    writetable(results_matrix_table, results_matrix_filename);
    
    fprintf('\n Matrix Results saved to:\n  - %s\n', results_matrix_filename);

        %% === Save All Figures Automatically ===
    fig_dir = fullfile('figures_summary', sprintf('sizing_%s_%dbus_%s', opt, num_buses, timestamp));
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
    
    fprintf('\n All Figures saved to:\n- %s\n', fig_dir);
end