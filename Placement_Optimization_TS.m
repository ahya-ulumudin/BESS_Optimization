function [BESS_Location, BESS_Output, obj, fitness_history, fitness_eval, iter, eval_total, best_iter_idx, best_eval_idx] = ...
    Placement_Optimization_TS(lower_bound, upper_bound, BESS_Number, Candidate_Buses, ...
    hour, mm, ll, sel_pv, sel_lp, MVAb, Zb, stagnation_limit, Bus_Ranked, use_ranked_bus_guidance, guided_fraction, BESS_Eff)
% PLACEMENT_OPTIMIZATION_TS_R1
% Tabu Search for BESS placement with aspiration criterion and soft restart on stagnation.

    %% === Parameters ===
    max_iter = 300;
    num_neighbors = 100;
    num_variables = BESS_Number * 2;
    num_candidates = length(Candidate_Buses);
    tabu_tenure = min(40, max(20, round(0.85 * max_iter)));
    gaussian_std_start = 50;
    gaussian_std_end = 5;
    tolerance = 1e-6;
    soft_restart_threshold = round(0.5 * stagnation_limit);

    %% === Initialization ===
    tabu_list = zeros(tabu_tenure, num_variables);
    eval_total = 0;
    eval_best = 0;
    best_iter_idx = 0;
    best_eval_idx = 0;
    fitness_eval = [];

    selected_buses = sort(randsample(Candidate_Buses, BESS_Number));
    selected_outputs = round((lower_bound + (upper_bound - lower_bound) * rand(1, BESS_Number)) / 10) * 10;
    current_solution = [selected_buses, selected_outputs];
    best_solution = current_solution;

    obj = Placement_Objective(selected_buses, selected_outputs, mm, ll, sel_pv, sel_lp, MVAb, Zb, BESS_Eff);
    best_obj = obj;
    prev_best_obj = best_obj;
    eval_total = eval_total + 1;
    fitness_eval(end+1) = best_obj;
    fitness_history = zeros(1, max_iter);
    stagnation_count = 0;
    termination_logged = false;

    %% === Main Loop ===
    for k = 1:max_iter
        current_gaussian_std = gaussian_std_start - (gaussian_std_start - gaussian_std_end) * (k / max_iter);
        best_neighbor_obj = inf;
        best_neighbor = [];
        is_aspired = false;

        for n = 1:num_neighbors
            neighbor = current_solution;

            if use_ranked_bus_guidance && n <= round(guided_fraction * num_neighbors)
                selected_buses = Bus_Ranked(1:BESS_Number);
                selected_outputs = round((lower_bound + (upper_bound - lower_bound) * rand(1, BESS_Number)) / 10) * 10;
                selected_outputs = max(lower_bound, min(upper_bound, selected_outputs));
                neighbor = [selected_buses(:)', selected_outputs];
            else
                for i = 1:BESS_Number
                    if rand < 0.5
                        neighbor(i) = Candidate_Buses(randi(num_candidates));
                    end
                end
                loc_part = unique(neighbor(1:BESS_Number));
                while length(loc_part) < BESS_Number
                    new_bus = Candidate_Buses(randi(num_candidates));
                    if ~ismember(new_bus, loc_part)
                        loc_part(end+1) = new_bus;
                    end
                end
                neighbor(1:BESS_Number) = sort(loc_part);

                for i = BESS_Number+1:num_variables
                    if rand < 0.5
                        neighbor(i) = neighbor(i) + round(normrnd(0, current_gaussian_std));
                        neighbor(i) = round(min(max(neighbor(i), lower_bound), upper_bound) / 10) * 10;
                    end
                end
            end

            loc = neighbor(1:BESS_Number);
            out = neighbor(BESS_Number+1:end);
            fval = Placement_Objective(loc, out, mm, ll, sel_pv, sel_lp, MVAb, Zb, BESS_Eff);
            eval_total = eval_total + 1;

            is_in_tabu = ismember(neighbor, tabu_list, 'rows');
            if is_in_tabu && fval < best_obj
                is_in_tabu = false;
                is_aspired = true;
            end

            if ~is_in_tabu && fval < best_neighbor_obj
                best_neighbor_obj = fval;
                best_neighbor = neighbor;
            end

            fitness_eval(end+1) = min(fval, fitness_eval(end));
        end

        if isempty(best_neighbor)
            fprintf('[Hour %02d] [TS] Iteration %3d: No valid neighbor found. Skipping.\n', hour, k);
            continue;
        end

        current_solution = best_neighbor;
        tabu_list = [tabu_list(2:end, :); best_neighbor];

        if best_neighbor_obj < best_obj
            best_solution = best_neighbor;
            best_obj = best_neighbor_obj;
            eval_best = eval_total;
            best_eval_idx = eval_total;
            best_iter_idx = k;
        end

        delta = abs(best_obj - prev_best_obj);
        if delta < tolerance
            stagnation_count = stagnation_count + 1;
        else
            stagnation_count = 0;
        end
        prev_best_obj = best_obj;

        fitness_history(k) = best_obj;

        if mod(k, 25) == 0 || k == 1
            % fprintf('[Hour %02d] [TS] Iteration %3d | Best Fitness: %.6f | Evaluations: %d\n', hour, k, best_obj, eval_total);
        end

        %% === Soft Restart ===
        if stagnation_count == soft_restart_threshold
            % fprintf('[Hour %02d] [TS] Soft restart triggered at iteration %3d due to partial stagnation.\n', hour, k);
            selected_buses = sort(randsample(Candidate_Buses, BESS_Number));
            selected_outputs = round((lower_bound + (upper_bound - lower_bound) * rand(1, BESS_Number)) / 10) * 10;
            current_solution = [selected_buses, selected_outputs];
            tabu_list = zeros(tabu_tenure, num_variables);
            % stagnation_count = 0;

            obj = Placement_Objective(selected_buses, selected_outputs, mm, ll, sel_pv, sel_lp, MVAb, Zb, BESS_Eff);
            fitness_eval(end+1) = obj;
            eval_total = eval_total + 1;
            continue;
        end

        %% === Stop on Full Stagnation ===
        if stagnation_count >= stagnation_limit
            fprintf('[Hour %02d] [TS] Terminated due to stagnation at iteration %3d | Final Fitness: %.6f\n', hour, k, best_obj);
            termination_logged = true;
            break;
        end
    end

    if ~termination_logged
        fprintf('[Hour %02d] [TS] Completed all iterations. Final Fitness: %.6f\n', hour, best_obj);
    end

    %% === Output ===
    BESS_Location = best_solution(1:BESS_Number);
    BESS_Output   = best_solution(BESS_Number+1:end);
    obj           = best_obj;
    iter          = best_iter_idx;
end
