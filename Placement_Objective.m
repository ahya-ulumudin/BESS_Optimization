function obj = Placement_Objective(BESS_Locations, BESS_Output, mm, ll, sel_pv, sel_lp, MVAb, Zb, BESS_Eff)
% PLACEMENT_OBJECTIVE Evaluate the objective function for initial BESS placement optimization (snapshot 1 hour).

    %% Initialize
    num_buses = size(mm, 1);
    max_expected_Pwr_Losses = 0.05;
    max_expected_Vol_Dev = 0.05;

    % Weights
    w_voltage_1h    = 2.5;
    w_losses_1h     = 2.5;
    w_reverse_1h    = 2.5;
    w_throughput_1h = 0.5;

    %% Generate BESS Active Power Injection
    BESS_Demand = zeros(num_buses, 1);
    BESS_Demand(BESS_Locations) = BESS_Output;
    pos_mask = BESS_Output >= 0;   % unit yang discharge
    neg_mask = BESS_Output <  0;   % unit yang charge
    if any(pos_mask)
        idx_pos = BESS_Locations(pos_mask);
        BESS_Demand(idx_pos) = BESS_Output(pos_mask) * BESS_Eff;
    end
    if any(neg_mask)
        idx_neg = BESS_Locations(neg_mask);
        BESS_Demand(idx_neg) = BESS_Output(neg_mask) / BESS_Eff;
    end

    %% Perform Load Flow
    [voltage, P_Loss_Kw, ~, ~, ~, ld, ~, ~, ~, ~, ~] = ...
        HourlyLoadFlow(mm, ll, sel_pv, sel_lp, MVAb, Zb, BESS_Demand);

    %% Penalty Calculations
    total_load = sum(mm(:,2)) * sel_lp;
    total_loss = sum(P_Loss_Kw);
    net_power_flow = ld(:,2);

    % --- Voltage Stability Penalty
    voltage_deviation = abs(voltage - 1.0);
    F_Voltage_Correction_1h = (0.2 * mean(voltage_deviation) + 0.8 * max(voltage_deviation)) / max_expected_Vol_Dev * 100;

    % --- Power Loss Penalty
    Percent_loss = (total_loss / (total_load + eps)) * 100;
    F_Losses_Correction_1h = (0.2 * Percent_loss + 0.8 * Percent_loss) / max_expected_Pwr_Losses * 100;

    % --- Reverse Power Flow Penalty
    reverse_flow_bus = net_power_flow < 0;           % bus experiencing reverse
    Reverse_Bus_Penalty = (sum(abs(net_power_flow(reverse_flow_bus))) / (sum(abs(net_power_flow)) + eps)) * 100;
    System_Reverse_Penalty = (sum(net_power_flow) < 0) * 100;
    F_Reverse_Power_1h = Reverse_Bus_Penalty + System_Reverse_Penalty;

    % --- BESS Throughput Penalty
    F_BESS_Throughput_1h = (sum(abs(BESS_Demand)) / (total_load + eps)) * 100;

    %% Final Composite Objective
    total_weight_1h = w_voltage_1h + w_losses_1h + w_reverse_1h + w_throughput_1h;
    
    obj = ((w_voltage_1h   * F_Voltage_Correction_1h + ...
           w_losses_1h     * F_Losses_Correction_1h + ...
           w_reverse_1h    * F_Reverse_Power_1h + ...
           w_throughput_1h * F_BESS_Throughput_1h) / total_weight_1h)/100;

end
