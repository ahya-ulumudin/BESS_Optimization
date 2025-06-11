function obj = Sizing_Objective(BESS_Output, mm, ll, PVout, L_prof, MVAb, Zb, upper_bound, Bus_Placement)
% SIZING_OBJECTIVE Evaluate the composite objective function for BESS sizing and 24-hour scheduling optimization.

    %% DECLARE GLOBAL VARIABLE
    global F_Voltage_Correction_24h F_Losses_Correction_24h F_Reverse_Power_24h
    global F_BESS_Throughput_24h F_Energy_Unbalance_24h F_Ramp_Rate_24h    

    %% PARAMETER INITIALIZATION
    num_buses = size(mm, 1);
    max_expected_Pwr_Losses = 0.05;    % Acceptable system losses threshold
    max_expected_Vol_Dev    = 0.05;    % Acceptable voltage deviation threshold

    % Weighting factors for each objective component
    w_voltage_24h    = 2.5;
    w_losses_24h     = 2.5;
    w_reverse_24h    = 2.5;
    w_throughput_24h = 0.5;
    w_ramp_24h       = 1.5;
    w_unbalance_24h  = 0.5;

    % Initialize storage matrices
    Matrix_Voltages    = zeros(num_buses, 24);
    Matrix_P_Loss      = zeros(size(ll,1), 24);
    Matrix_P_Nett      = zeros(num_buses, 24);
    Matrix_Load        = zeros(num_buses, 24);
    Matrix_BESS_Demand = zeros(num_buses, 24);

    %% LOAD FLOW CALCULATION  FOR 24 HOURS
    for hour = 1:24
        sel_pv = PVout(hour, 2);
        sel_lp = L_prof(hour, 2);

        % Generate BESS injection profile for the hour
        BESS_Demand = zeros(num_buses, 1);
        BESS_Demand(Bus_Placement) = BESS_Output(:, hour);

        % Perform hourly load flow analysis
        [voltage, P_Loss_Kw, ~, ~, ~, ld, ~, ~, ~, ~, ~] = ...
            HourlyLoadFlow(mm, ll, sel_pv, sel_lp, MVAb, Zb, BESS_Demand);

        % Record hourly results
        Matrix_Voltages(:, hour)    = voltage;
        Matrix_P_Loss(:, hour)      = P_Loss_Kw;
        Matrix_P_Nett(:, hour)      = ld(:,2);
        Matrix_Load(:, hour)        = mm(:,2) * sel_lp;
        Matrix_BESS_Demand(:, hour) = BESS_Demand;
    end

    %% PENALTY CALCULATION

    % Voltage Stability Penalty (mean and maximum deviation)
    bus_voltage_range            = max(Matrix_Voltages, [], 2) - min(Matrix_Voltages, [], 2);
    bus_max_voltage_deviation    = max(abs(Matrix_Voltages - 1.0), [], 2);
    F_Voltage_Correction_24h     = (0.2 * mean(bus_voltage_range) + 0.8 * max(bus_voltage_range)) / max_expected_Vol_Dev * 100 ...
                                 + (0.2 * mean(bus_max_voltage_deviation) + 0.8 * max(bus_max_voltage_deviation)) / max_expected_Vol_Dev * 100;

    % Power Loss Penalty (mean and peak loss percentage)
    P_Loss_Percentage       = sum(Matrix_P_Loss, 1) ./ sum(Matrix_Load, 1) * 100;
    F_Losses_Correction_24h = ((0.2 * mean(P_Loss_Percentage) + 0.8 * max(P_Loss_Percentage)) / max_expected_Pwr_Losses);

    % Reverse Power Flow Penalty (percentage of hours with reverse flow)
    Max_reverse_power_perbus = ((max(sum(Matrix_P_Nett<0, 2)))/24) * 100;
    Reverse_Flow_Ratio       = sum(abs(Matrix_P_Nett(Matrix_P_Nett < 0)))/sum(sum(abs(Matrix_P_Nett))) * 100;
    Reverse_Flow_Penalty     = any(sum(sum(Matrix_P_Nett, 1) < 0)) * 100;
    F_Reverse_Power_24h      = Reverse_Flow_Ratio + Reverse_Flow_Penalty;

    % BESS Throughput Efficiency Penalty
    F_BESS_Throughput_24h = (sum(sum(abs(Matrix_BESS_Demand))) / sum(sum(abs(Matrix_Load)))) * 100;

    % Energy Unbalance Penalty (charge-discharge mismatch)
    total_charging         = sum(-min(BESS_Output, 0), 2);
    total_discharging      = sum(max(BESS_Output, 0), 2);
    unbalance_ratio        = abs(total_charging - total_discharging) ./ (total_charging + total_discharging + eps) * 100;
    F_Energy_Unbalance_24h = (1 + max(unbalance_ratio))^2;

    % Ramp Rate Penalty (limit BESS output variations)
    max_ramp        = upper_bound * 0.2;
    BESS_shifted    = [BESS_Output(:,24), BESS_Output(:,1:23)];
    ramp_change     = abs(diff(BESS_shifted, 1, 2));
    ramp_violation  = max(0, ramp_change - max_ramp);
    ramp_penalty    = (0.5 * mean(ramp_violation, 2) + 0.5 * max(ramp_violation, [], 2));
    F_Ramp_Rate_24h = mean(ramp_penalty) / max_ramp * 100;

    %% FINAL COMPOSITE OBJECTIVE CALCULATION
    total_weight_24h = w_voltage_24h + w_losses_24h + w_reverse_24h + w_throughput_24h + w_ramp_24h + w_unbalance_24h;
    
    obj = (w_voltage_24h    * F_Voltage_Correction_24h + ...
           w_losses_24h     * F_Losses_Correction_24h + ...
           w_reverse_24h    * F_Reverse_Power_24h + ...
           w_throughput_24h * F_BESS_Throughput_24h + ...
           w_ramp_24h       * F_Ramp_Rate_24h + ...
           w_unbalance_24h  * F_Energy_Unbalance_24h) / total_weight_24h;

end