function [BESS_Output, obj, fitness_history, iter, eval_total] = ...
    Sizing_Optimization_PSO_TS(mm, lower_bound, upper_bound, BESS_Number, objective_function, stagnation_limit, initial_solution, BESS_Eff)
% SIZING_OPTIMIZATION_PSO_TS_TRY
% Hybrid Particle Swarm Optimization (PSO) and Tabu Search (TS)
% for optimal 24-hour BESS demand scheduling in distribution systems.

    %% === Algorithm Parameters ===
    max_gen  = 500;
    pop_size = 100;
    w_max = 0.7; w_min = 0.4;
    tolerance = 1e-6;
    num_variables = BESS_Number * 24;

    %% === Initialization ===
    position = round((lower_bound + (upper_bound - lower_bound) .* rand(pop_size, num_variables)) / 10) * 10;
    if exist('initial_solution', 'var') && ~isempty(initial_solution)
        if numel(initial_solution) == num_variables
            position(1,:) = initial_solution;
            fprintf('Initial solution injected into the first particle (PSO-TS)\n');
        else
            warning('Initial solution size mismatch. Skipping initial solution injection.');
        end
    end

    velocity = zeros(pop_size, num_variables);
    fitness = arrayfun(@(i) objective_function(reshape(position(i,:), BESS_Number, 24)), 1:pop_size);
    eval_total = pop_size;

    pbest = position;
    pbest_fitness = fitness;
    [obj, gbest_idx] = min(fitness);
    gbest = position(gbest_idx, :);

    fitness_history = [];
    stagnation_counter = 0;
    termination_logged = false;
    total_ts_iterations = 0;

    fprintf('Initial Global Best Fitness: %.6f\n', obj);

    %% === Main PSO Optimization Loop ===
    for gen = 1:max_gen
        w = w_max - ((w_max - w_min) * gen / max_gen);
        c1 = 2.0 - 1.0 * gen / max_gen;
        c2 = 1.0 + 1.5 * gen / max_gen;

        for i = 1:pop_size
            r1 = rand(1, num_variables);
            r2 = rand(1, num_variables);

            velocity(i,:) = w * velocity(i,:) + ...
                            c1 * r1 .* (pbest(i,:) - position(i,:)) + ...
                            c2 * r2 .* (gbest - position(i,:));
            position(i,:) = round((position(i,:) + velocity(i,:)) / 10) * 10;
            position(i,:) = max(lower_bound, min(upper_bound, position(i,:)));

            fit = objective_function(reshape(position(i,:), BESS_Number, 24));
            eval_total = eval_total + 1;

            if fit < pbest_fitness(i)
                pbest(i,:) = position(i,:);
                pbest_fitness(i) = fit;
                if fit < obj
                    gbest = position(i,:);
                    obj = fit;
                end
            end
        end

        fitness_history(end+1) = obj;

        %% === Stagnation Check ===
        if length(fitness_history) > 1
            delta = abs(fitness_history(end) - fitness_history(end-1));
            if delta < tolerance
                stagnation_counter = stagnation_counter + 1;
            else
                stagnation_counter = 0;
            end
        end

        if stagnation_counter >= stagnation_limit
            fprintf('[PSO-TS] Terminated early due to stagnation at PSO generation %3d | Fitness: %.6f\n', gen, obj);
            termination_logged = true;
            break;
        end

        %% === Tabu Search Refinement every 50 generations ===
        if mod(gen, 50) == 0
            fprintf('PSO Generation %3d | Best Fitness: %.6f | Evaluations: %d\n', gen, obj, eval_total);

            [refined_gbest, refined_obj, ts_eval, ts_hist, ts_iter] = ...
                TabuSearch_LocalRefine(gbest, obj, 50, lower_bound, upper_bound, BESS_Number, objective_function, tolerance, BESS_Eff);

            eval_total = eval_total + ts_eval;
            total_ts_iterations = total_ts_iterations + ts_iter;

            %% Reset stagnation if TS improves the solution
            if refined_obj < obj
                gbest = refined_gbest;
                obj = refined_obj;
                stagnation_counter = 0;
            end

            %% === Append TS history and re-check stagnation
            for k = 2:length(ts_hist)
                fitness_history(end+1) = min(ts_hist(k), fitness_history(end));
                if length(fitness_history) > 1
                    delta = abs(fitness_history(end) - fitness_history(end-1));
                    if delta < tolerance
                        stagnation_counter = stagnation_counter + 1;
                    else
                        stagnation_counter = 0;
                    end
                end
                if stagnation_counter >= stagnation_limit
                    fprintf('[PSO-TS] Terminated early after TS refinement at generation %3d | Fitness: %.6f\n', gen, refined_obj);
                    termination_logged = true;
                    break;
                end
            end

            if termination_logged
                break;
            end
        end
    end

    if ~termination_logged
        fprintf('[PSO-TS] Final Fitness after %d generations: %.6f\n', gen, obj);
    end

    %% === Final Output ===
    BESS_Output = reshape(round(gbest / 10) * 10, BESS_Number, 24);
    iter = length(fitness_history);

    %% === Convergence Plot ===
    figure;
    plot(1:iter, fitness_history, 'b-', 'LineWidth', 2);
    xlabel('Iteration', 'FontSize', 14);
    ylabel('Objective Value', 'FontSize', 14);
    title('Convergence Trend of Hybrid PSO–TS Sizing Optimization', 'FontSize', 16, 'FontWeight', 'bold');
    grid on;

    % === Highlight last fitness value ===
    x_last = iter;
    y_last = fitness_history(end);

    hold on;
    plot(x_last, y_last, 'ro', 'MarkerFaceColor', 'r');  % Red dot marker
    text(x_last, y_last, sprintf('  %.2f', y_last), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 12, 'FontWeight', 'bold', ...
        'Color', 'red');
    hold off;

end

function [refined_gbest, refined_obj, ts_eval, ts_hist, ts_iter] = ...
    TabuSearch_LocalRefine(gbest, obj, max_iter, lower_bound, upper_bound, BESS_Number, objective_function, tolerance, BESS_Eff)
% TABUSEARCH_LOCALREFINE
% Tabu Search local refinement with aspiration criterion and soft restart.

    %% === Initialization ===
    num_variables = BESS_Number * 24;
    current_solution = gbest;
    current_obj = obj;

    refined_gbest = current_solution;
    refined_obj = current_obj;
    ts_hist = refined_obj;
    ts_eval = 0;
    ts_iter = 0;
    stagnation_counter = 0;
    termination_logged = false;

    stagnation_limit = 30;
    soft_restart_threshold = round(0.5 * stagnation_limit);

    tabu_tenure = min(40, max(20, round(0.85 * max_iter))); % Adaptive tabu list size
    tabu_list = zeros(tabu_tenure, num_variables);
    tabu_ptr = 1;

    num_neighbors = 100;

    %% === Main Tabu Search Loop ===
    for iter = 1:max_iter
        ts_iter = ts_iter + 1;
        current_gaussian_std = max(5, 30 - 25 * (iter / max_iter));

        best_neighbor_obj = inf;
        best_neighbor = [];

        % === Generate Neighbors ===
        for n = 1:num_neighbors
            neighbor = current_solution;

            % Mutate random subset
            idx = randperm(num_variables, randi([1, round(0.1 * num_variables)]));
            perturb = round(normrnd(0, current_gaussian_std, size(idx)));
            neighbor(idx) = neighbor(idx) + perturb;

            % Quantize and limit
            neighbor = round(neighbor / 10) * 10;
            neighbor = max(lower_bound, min(upper_bound, neighbor));

            % Evaluate
            neighbor_obj = objective_function(reshape(neighbor, BESS_Number, 24));
            ts_eval = ts_eval + 1;

            % Tabu Check
            is_tabu = any(all(abs(tabu_list - neighbor) < 1e-6, 2));

            % === Aspiration Criterion ===
            if is_tabu && neighbor_obj < refined_obj
                is_tabu = false;
            end

            % Update best neighbor
            if ~is_tabu && neighbor_obj < best_neighbor_obj
                best_neighbor = neighbor;
                best_neighbor_obj = neighbor_obj;
            end
        end

        % === No Valid Neighbor ===
        if isempty(best_neighbor)
            continue;
        end

        % === Move to Best Neighbor ===
        current_solution = best_neighbor;
        current_obj = best_neighbor_obj;

        tabu_list(tabu_ptr, :) = best_neighbor;
        tabu_ptr = mod(tabu_ptr, size(tabu_list, 1)) + 1;

        % === Update Global Best if Improved ===
        if current_obj < refined_obj
            refined_gbest = current_solution;
            refined_obj = current_obj;
        end

        ts_hist(end+1) = refined_obj;

        %% === Stagnation Tracking ===
        if length(ts_hist) > 1
            delta = abs(ts_hist(end) - ts_hist(end-1));
            if delta < tolerance
                stagnation_counter = stagnation_counter + 1;
            else
                stagnation_counter = 0;
            end
        end

        %% === Soft Restart if Partial Stagnation ===
        if stagnation_counter == soft_restart_threshold
            fprintf('TS Soft restart triggered at iteration %d due to partial stagnation.\n', iter);
            current_solution = round((lower_bound + (upper_bound - lower_bound) .* rand(1, num_variables)) / 10) * 10;
            tabu_list = zeros(tabu_tenure, num_variables);
            stagnation_counter = 0;

            current_obj = objective_function(reshape(current_solution, BESS_Number, 24));
            ts_eval = ts_eval + 1;
            ts_hist(end+1) = current_obj;

            continue;
        end

        %% === Early Termination ===
        if stagnation_counter >= stagnation_limit
            fprintf('TS Terminated early due to stagnation at TS iteration %3d | Best Fitness: %.6f\n', iter, refined_obj);
            termination_logged = true;
            break;
        end
    end

    %% === Final Logging ===
    if ~termination_logged
        fprintf('TS Refinement %d  | Best Fitnes : %.6f | Evaluations: %d\n', ts_iter, refined_obj, ts_eval);
    end
end
