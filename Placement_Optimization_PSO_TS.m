function [BESS_Location, BESS_Output, obj, fitness_history, fitness_eval, iter, eval_total, best_iter_idx, best_eval_idx] = ...
    Placement_Optimization_PSO_TS(lower_bound, upper_bound, BESS_Number, Candidate_Buses, ...
    hour, mm, ll, sel_pv, sel_lp, MVAb, Zb, stagnation_limit, Bus_Ranked, use_ranked_bus_guidance, guided_fraction)
% PLACEMENT_OPTIMIZATION_PSO_TS_R1
% Hybrid PSO-TS optimization for BESS placement and sizing

    %% === PSO-TS Parameters ===
    pop_size = 100;
    max_gen = 150;
    ts_max_iter = 30;
    ts_interval = 30;
    w_max = 0.7; w_min = 0.4;
    num_neighbors = 100;
    num_candidates = length(Candidate_Buses);
    num_variables = BESS_Number * 2;
    tolerance = 1e-6;

    %% === Initialization ===
    guided_pop = round(guided_fraction * pop_size);

    position = zeros(pop_size, num_variables);
    velocity = zeros(pop_size, num_variables);
    fitness = zeros(1, pop_size);
    fitness_eval = [];
    fitness_history = [];
    eval_total = 0;
    stagnation_count = 0;
    termination_logged = false;
    best_iter_idx = 1;
    best_eval_idx = 1;

    %% === Initial Population Generation ===
    for i = 1:pop_size
        if use_ranked_bus_guidance && i <= guided_pop
            selected_buses = Bus_Ranked(1:BESS_Number);
            selected_outputs = lower_bound + (upper_bound - lower_bound) * rand(1, BESS_Number);
            selected_outputs = round(selected_outputs / 10) * 10;
            selected_outputs = max(lower_bound, min(upper_bound, selected_outputs));
            position(i,:) = [selected_buses, selected_outputs];
        else
            loc_idx = randi(num_candidates, 1, BESS_Number);
            selected_buses = Candidate_Buses(loc_idx);
            selected_outputs = lower_bound + (upper_bound - lower_bound) * rand(1, BESS_Number);
            selected_outputs = round(selected_outputs / 10) * 10;
            selected_outputs = max(lower_bound, min(upper_bound, selected_outputs));
            position(i,:) = [selected_buses, selected_outputs];
        end

        [loc, out] = decode_particle(position(i,:), BESS_Number, Candidate_Buses, lower_bound, upper_bound);
        fitness(i) = Placement_Objective(loc, out, mm, ll, sel_pv, sel_lp, MVAb, Zb);
        eval_total = eval_total + 1;
        fitness_eval(end+1) = fitness(i);
    end

    %% === Initialize pbest and gbest ===
    pbest = position;
    pbest_fitness = fitness;
    [obj, gbest_idx] = min(fitness);
    gbest = position(gbest_idx,:);
    fitness_history(end+1) = obj;

    %% === Main PSO-TS Loop ===
    for gen = 1:max_gen
        w = w_max - ((w_max - w_min) * gen / max_gen);
        c1 = 2.0 - 1.0 * (gen / max_gen);
        c2 = 1.0 + 1.5 * (gen / max_gen);

        for i = 1:pop_size
            r1 = rand(1, num_variables);
            r2 = rand(1, num_variables);

            velocity(i,:) = w * velocity(i,:) + c1 * r1 .* (pbest(i,:) - position(i,:)) + c2 * r2 .* (gbest - position(i,:));
            position(i,:) = max(lower_bound, min(upper_bound, position(i,:) + velocity(i,:)));
            position(i,1:BESS_Number) = round(position(i,1:BESS_Number));

            [loc, out] = decode_particle(position(i,:), BESS_Number, Candidate_Buses, lower_bound, upper_bound);
            fit = Placement_Objective(loc, out, mm, ll, sel_pv, sel_lp, MVAb, Zb);
            eval_total = eval_total + 1;
            fitness_eval(end+1) = fit;

            if fit < pbest_fitness(i)
                pbest(i,:) = position(i,:);
                pbest_fitness(i) = fit;
                if fit < obj
                    gbest = position(i,:);
                    obj = fit;
                    best_iter_idx = gen;
                    best_eval_idx = eval_total;
                end
            end
        end

        fitness_history(end+1) = obj;

        if length(fitness_history) > 1
            delta = abs(fitness_history(end) - fitness_history(end-1));
            stagnation_count = (delta < tolerance) * (stagnation_count + 1);
        end

        if stagnation_count >= stagnation_limit
            termination_logged = true;
            break;
        end

        if mod(gen, ts_interval) == 0
            [refined_gbest, refined_obj, ts_eval, ts_hist, ts_terminated] = TabuSearch_LocalRefine( ...
                gbest, obj, ts_max_iter, lower_bound, upper_bound, BESS_Number, Candidate_Buses, ...
                mm, ll, sel_pv, sel_lp, MVAb, Zb, num_neighbors, stagnation_limit, hour, tolerance);

            eval_total = eval_total + ts_eval;

            for k = 1:length(ts_hist)
                fitness_history(end+1) = min(ts_hist(k), fitness_history(end));
                if length(fitness_history) > 1
                    delta = abs(fitness_history(end) - fitness_history(end-1));
                    stagnation_count = (delta < tolerance) * (stagnation_count + 1);
                end
                if stagnation_count >= stagnation_limit
                    termination_logged = true;
                    break;
                end
            end

            if refined_obj < obj
                gbest = refined_gbest;
                obj = refined_obj;
            end

            if termination_logged
                break;
            end
        end
    end

    if termination_logged
        fprintf('[Hour %02d] [PSO-TS] Terminated at iteration %d | Final Fitness: %.6f\n', hour, gen, obj);
    else
        fprintf('[Hour %02d] [PSO-TS] Completed | Final Fitness: %.6f\n', hour, obj);
    end

    [BESS_Location, BESS_Output] = decode_particle(gbest, BESS_Number, Candidate_Buses, lower_bound, upper_bound);
    iter = length(fitness_history);
end

function [best_solution, best_obj, eval_count, fitness_history, terminated] = TabuSearch_LocalRefine( ...
    start_solution, start_obj, max_iter, ...
    lower_bound, upper_bound, BESS_Number, Candidate_Buses, ...
    mm, ll, sel_pv, sel_lp, MVAb, Zb, num_neighbors, stagnation_limit, hour, tolerance)
% TABUSEARCH_LOCALREFINE
% Local refinement using Tabu Search algorithm.

    num_variables = length(start_solution);
    current_solution = start_solution;
    best_solution = start_solution;
    best_obj = start_obj;
    prev_best_obj = best_obj;
    eval_count = 0;
    stagnation_count = 0;
    terminated = false;
    fitness_history = zeros(1, max_iter);

    tabu_tenure = min(40, max(20, round(0.85 * max_iter)));
    tabu_list = zeros(tabu_tenure, num_variables);

    soft_restart_threshold = round(0.5 * stagnation_limit);

    for iter = 1:max_iter
        current_gaussian_std = max(5, 30 - 25 * (iter / max_iter));
        best_neighbor_obj = inf;
        best_neighbor = [];

        for n = 1:num_neighbors
            neighbor = current_solution;

            for i = 1:BESS_Number
                if rand < 0.5
                    neighbor(i) = Candidate_Buses(randi(length(Candidate_Buses)));
                end
            end

            loc_part = unique(neighbor(1:BESS_Number));
            while length(loc_part) < BESS_Number
                new_bus = Candidate_Buses(randi(length(Candidate_Buses)));
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

            % === Evaluate fitness ===
            [loc, out] = decode_particle(neighbor, BESS_Number, Candidate_Buses, lower_bound, upper_bound);
            neighbor_obj = Placement_Objective(loc, out, mm, ll, sel_pv, sel_lp, MVAb, Zb);
            eval_count = eval_count + 1;

            % === Tabu and Aspiration Check ===
            is_in_tabu = ismember(neighbor, tabu_list, 'rows');
            if is_in_tabu && neighbor_obj < best_obj
                is_in_tabu = false;  % Aspiration criterion triggered
            end

            if ~is_in_tabu && neighbor_obj < best_neighbor_obj
                best_neighbor_obj = neighbor_obj;
                best_neighbor = neighbor;
            end
        end

        if isempty(best_neighbor)
            continue;
        end

        current_solution = best_neighbor;
        tabu_list = [tabu_list(2:end, :); best_neighbor];

        if best_neighbor_obj < best_obj
            best_solution = best_neighbor;
            prev_best_obj = best_obj;
            best_obj = best_neighbor_obj;

            delta = abs(best_obj - prev_best_obj);
            if delta < tolerance
                stagnation_count = stagnation_count + 1;
            else
                stagnation_count = 0;
            end
        else
            stagnation_count = stagnation_count + 1;
        end

        prev_best_obj = best_obj;
        fitness_history(iter) = best_obj;

        %% === Soft Restart if Partial Stagnation ===
        if stagnation_count == soft_restart_threshold
            % fprintf('[Hour %02d] [TS-Refine] Soft restart triggered at iteration %d\n', hour, iter);
            new_buses = sort(randsample(Candidate_Buses, BESS_Number));
            new_outputs = round((lower_bound + (upper_bound - lower_bound) * rand(1, BESS_Number)) / 10) * 10;
            new_outputs = max(lower_bound, min(upper_bound, new_outputs));
            current_solution = [new_buses, new_outputs];
            tabu_list = zeros(tabu_tenure, num_variables);
            % stagnation_count = 0;

            [loc, out] = decode_particle(current_solution, BESS_Number, Candidate_Buses, lower_bound, upper_bound);
            best_obj = Placement_Objective(loc, out, mm, ll, sel_pv, sel_lp, MVAb, Zb);
            eval_count = eval_count + 1;
            fitness_history(iter) = best_obj;

            continue;
        end

        %% === Early Termination ===
        if stagnation_count >= stagnation_limit
            terminated = true;
            break;
        end
    end

    fitness_history = fitness_history(1:iter);
end

function [locations, outputs] = decode_particle(particle, BESS_Number, Candidate_Buses, lower_bound, upper_bound)
% DECODE_PARTICLE

    locations = round(particle(1:BESS_Number));

    for i = 1:BESS_Number
        if isnan(locations(i)) || ~isfinite(locations(i)) || locations(i) <= 0
            locations(i) = Candidate_Buses(randi(length(Candidate_Buses)));
        else
            [~, idx] = min(abs(Candidate_Buses - locations(i)));
            if isempty(idx)
                locations(i) = Candidate_Buses(randi(length(Candidate_Buses)));
            else
                locations(i) = Candidate_Buses(idx);
            end
        end
    end

    locations = unique(locations);
    while length(locations) < BESS_Number
        new_bus = Candidate_Buses(randi(length(Candidate_Buses)));
        if ~ismember(new_bus, locations)
            locations(end+1) = new_bus;
        end
    end
    locations = sort(locations);

    outputs = particle(BESS_Number+1:end);
    outputs = max(lower_bound, min(upper_bound, round(outputs / 10) * 10));
end