function [best_solution, obj, fitness_history, iter, eval_total] = ...
    Sizing_Optimization_TS(mm, lower_bound, upper_bound, BESS_Number, ...
    objective_function, stagnation_limit, initial_solution)
% SIZING_OPTIMIZATION_TS
% Tabu Search optimization for 24-hour BESS sizing and scheduling with
% aspiration criterion, soft restart, and stable stagnation logic.

    %% === Parameters ===
    max_iter = 1000;
    num_neighbors = 100;
    tabu_tenure = min(40, max(20, round(0.85 * max_iter)));
    gaussian_std_start = 50;
    gaussian_std_end = 5;
    num_variables = BESS_Number * 24;
    tolerance = 1e-6;
    soft_restart_threshold = round(0.5 * stagnation_limit);

    %% === Initialization ===
    if exist('initial_solution', 'var') && ~isempty(initial_solution) && numel(initial_solution) == num_variables
        current_solution = initial_solution;
        fprintf('Initial solution injected into Tabu Search.\n');
    else
        current_solution = round((lower_bound + (upper_bound - lower_bound) .* rand(1, num_variables)) / 10) * 10;
    end

    best_solution = current_solution;
    obj = objective_function(reshape(current_solution, BESS_Number, 24));
    eval_total = 1;
    iter = 0;

    tabu_list = zeros(tabu_tenure, num_variables);
    fitness_history = zeros(1, max_iter);
    stagnation_count = 0;

    %% === Main TS Loop ===
    for k = 1:max_iter
        iter = iter + 1;
        prev_obj = obj;
        current_gaussian_std = gaussian_std_start - ...
            (gaussian_std_start - gaussian_std_end) * (k / max_iter);

        best_neighbor_obj = inf;
        best_neighbor = [];

        %% === Generate Neighbors ===
        for n = 1:num_neighbors
            neighbor = current_solution + normrnd(0, current_gaussian_std, 1, num_variables);
            neighbor = round(neighbor / 10) * 10;
            neighbor = max(lower_bound, min(upper_bound, neighbor));

            fval = objective_function(reshape(neighbor, BESS_Number, 24));
            eval_total = eval_total + 1;

            is_tabu = any(all(abs(tabu_list - neighbor) < 1e-6, 2));

            % Aspiration Criterion
            if is_tabu && fval < best_solution
                is_tabu = false;
            end

            if ~is_tabu && fval < best_neighbor_obj
                best_neighbor = neighbor;
                best_neighbor_obj = fval;
            end
        end

        %% === No Valid Neighbor ===
        if isempty(best_neighbor)
            fprintf('Iteration %3d: No valid neighbor found. Skipping.\n', k);
            fitness_history(k) = obj;
            continue;
        end

        %% === Move to Best Neighbor ===
        current_solution = best_neighbor;
        tabu_list = [tabu_list(2:end, :); best_neighbor];

        if best_neighbor_obj < obj
        best_solution = best_neighbor;
        obj = best_neighbor_obj;
        stagnation_count = 0;
    else
        delta = abs(obj - prev_obj);
        if delta < tolerance
            stagnation_count = stagnation_count + 1;
        else
            stagnation_count = 0;
        end
    end
    
    prev_obj = obj;

        fitness_history(k) = obj;

        if mod(k, 50) == 0 || k == 1
            fprintf('TS Iteration %3d | Best Fitness: %.6f | Evaluations: %d\n', k, obj, eval_total);
        end

        %% === Soft Restart ===
        if stagnation_count == soft_restart_threshold
            fprintf('Soft restart triggered at iteration %3d.\n', k);
            current_solution = round((lower_bound + (upper_bound - lower_bound) .* rand(1, num_variables)) / 10) * 10;
            tabu_list = zeros(tabu_tenure, num_variables);
            obj_restart = objective_function(reshape(current_solution, BESS_Number, 24));
            eval_total = eval_total + 1;
            continue;
        end

        %% === Early Termination
        if stagnation_count >= stagnation_limit
            fprintf('Stagnation detected at iteration %3d.\n', k);
            break;
        end
    end

    %% === Final Output ===
    best_solution = reshape(best_solution, BESS_Number, 24);

    % Plot convergence
    figure;
    plot(1:iter, fitness_history(1:iter), 'b-', 'LineWidth', 2);
    xlabel('Iteration', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Best Objective Value', 'FontSize', 12, 'FontWeight', 'bold');
    title('Convergence Trend of TS Sizing Optimization', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;

    % Highlight final fitness
    x_last = iter;
    y_last = fitness_history(iter);
    text(x_last, y_last, sprintf('  %.2f', y_last), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 12, 'FontWeight', 'bold', ...
        'Color', 'red');
    hold on;
    plot(x_last, y_last, 'ro', 'MarkerFaceColor', 'r');
    hold off;
end