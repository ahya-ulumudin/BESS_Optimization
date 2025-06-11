function Bus_Ranked = bus_ranking(Matrix_Voltages, Matrix_P_Nett, mm, output_filename)
%BUS_RANKING Ranks buses based on voltage deviation, net power fluctuation, and PV capacity.
%
% This function evaluates each bus based on three criteria:
% (1) Voltage Deviation, (2) Net Power Fluctuation, and (3) Distributed Generation Capacity,
% normalizes the values (based on the total sum per metric), computes a combined total score,
% ranks the buses accordingly, and exports the ranking table to a CSV file.

    %% === Data Dimension ===
    num_buses = size(Matrix_Voltages, 1);

    %% === Calculate Voltage Deviation for Each Bus ===
    v_min = min(Matrix_Voltages, [], 2);
    v_max = max(Matrix_Voltages, [], 2);
    vol_deviation = v_max - v_min; % Voltage range per bus over 24 hours

    % Normalize Voltage Deviation by total deviation sum
    sum_vol_deviation = sum(vol_deviation);
    vol_dev_score = vol_deviation / max(sum_vol_deviation, eps); % Avoid division by zero

    %% === Calculate Net Power Fluctuation for Each Bus ===
    p_min = min(Matrix_P_Nett, [], 2);
    p_max = max(Matrix_P_Nett, [], 2);
    net_pwr_fluctuation = p_max - p_min; % Net power fluctuation per bus

    % Normalize Net Power Fluctuation by total fluctuation sum
    sum_net_pwr_fluctuation = sum(net_pwr_fluctuation);
    net_pwr_fluct_score = net_pwr_fluctuation / max(sum_net_pwr_fluctuation, eps); % Avoid division by zero

    %% === Extract Distributed Generation (DG) Capacity ===
    dg_capacity = mm(:, 4); % Installed PV capacity per bus

    % Normalize DG Capacity by total installed capacity
    total_dg_capacity = sum(dg_capacity);
    dg_cap_score = dg_capacity / max(total_dg_capacity, eps); % Avoid division by zero

    %% === Calculate Total Score per Bus ===
    total_score = vol_dev_score + net_pwr_fluct_score + dg_cap_score;

    %% === Create the Ranking Table ===
    RankData = table;
    RankData.Bus = (1:num_buses)';                      % Bus number
    RankData.Vol_Dev = vol_deviation;                   % Voltage deviation magnitude
    RankData.Vol_Dev_Score = vol_dev_score;              % Normalized voltage deviation score
    RankData.Net_Pwr_Fluct = net_pwr_fluctuation;        % Net power fluctuation magnitude
    RankData.Net_Pwr_Fluct_Score = net_pwr_fluct_score;  % Normalized net power fluctuation score
    RankData.DG_Cap = dg_capacity;                      % Distributed Generation installed capacity
    RankData.DG_Cap_Score = dg_cap_score;                % Normalized DG capacity score
    RankData.Total_Score = total_score;                  % Combined total score

    % Sort buses by Total_Score in descending order
    RankData = sortrows(RankData, 'Total_Score', 'descend');

    % Assign ranking based on sorted order
    RankData.Rank = (1:height(RankData))';
    RankData = movevars(RankData, 'Rank', 'Before', 'Bus');

    %% === Export the Result to CSV File ===
    writetable(RankData, output_filename);

    %% === Prepare Output Bus Ranking ===
    Bus_Ranked = RankData.Bus;  % Output only the sorted bus numbers

    %% === Display Summary Information ===
    fprintf('>> Bus ranking successfully saved to "%s"\n', output_filename);

end