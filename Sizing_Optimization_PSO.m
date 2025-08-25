function [BESS_Output, obj, fitness_history, iter, eval_total] = Sizing_Optimization_PSO(mm, lower_bound, upper_bound, BESS_Number, objective_function, stagnation_limit, initial_solution, BESS_Eff)
% SIZING_OPTIMIZATION_PSO
% Particle Swarm Optimization (PSO) for 24-hour BESS demand scheduling optimization.

    %% === PSO Parameters ===
    pop_size = 100;                   % Population size
    max_gen = 1000;                    % Maximum number of generations
    w_max = 0.7; w_min = 0.4;          % Inertia weight range
    tolerance = 1e-6;                 % Stagnation detection threshold
    num_variables = BESS_Number * 24; % Number of decision variables
    v_max = (upper_bound - lower_bound) / 2; % Maximum velocity

    %% === Particle Initialization ===
    position = lower_bound + (upper_bound - lower_bound) .* rand(pop_size, num_variables);
    velocity = zeros(pop_size, num_variables);
    fitness = zeros(1, pop_size);
    eval_total = 0;

    % Inject initial solution if provided
    if exist('initial_solution', 'var') && ~isempty(initial_solution)
        if numel(initial_solution) == num_variables
            position(1, :) = initial_solution;
            fprintf('Initial solution injected into the first particle.\n');
        else
            warning('Initial solution size mismatch. Injection skipped.');
        end
    end

    % Initial fitness evaluation
    for i = 1:pop_size
        rounded_pos = round(position(i, :) / 10) * 10;
        fitness(i) = objective_function(reshape(rounded_pos, BESS_Number, 24));
        eval_total = eval_total + 1;
    end

    % Initialize personal bests and global best
    pbest = position;
    pbest_fitness = fitness;
    [obj, gbest_idx] = min(fitness);
    gbest = position(gbest_idx, :);

    fitness_history = obj;
    stagnation_counter = 0;
    iter = 0;
    prev_best_obj = obj;

    %% === Main PSO Iterative Loop ===
    for gen = 1:max_gen
        iter = gen;

        % Adaptive inertia weight and acceleration coefficients
        w = w_max - ((w_max - w_min) * (gen / max_gen)^2); % Quadratic decreasing
        c1 = 2.0 - 1.0 * (gen / max_gen);                  % Cognitive component
        c2 = 1.0 + 1.5 * (gen / max_gen);                  % Social component

        for i = 1:pop_size
            r1 = rand(1, num_variables);
            r2 = rand(1, num_variables);

            % Update velocity and position
            velocity(i, :) = w * velocity(i, :) ...
                           + c1 * r1 .* (pbest(i, :) - position(i, :)) ...
                           + c2 * r2 .* (gbest - position(i, :));
            velocity(i, :) = max(-v_max, min(v_max, velocity(i, :)));

            position(i, :) = position(i, :) + velocity(i, :);
            position(i, :) = max(lower_bound, min(upper_bound, position(i, :)));

            % Fitness evaluation
            rounded_pos = round(position(i, :) / 10) * 10;
            fit = objective_function(reshape(rounded_pos, BESS_Number, 24));
            eval_total = eval_total + 1;
            fitness(i) = fit;

            % Update personal and global bests
            if fit < pbest_fitness(i)
                pbest(i, :) = position(i, :);
                pbest_fitness(i) = fit;
                if fit < obj
                    gbest = position(i, :);
                    obj = fit;
                end
            end
        end

        fitness_history(end+1) = obj;

        % Progress Logging
        if mod(gen, 50) == 0 || gen == 1
            fprintf('PSO Generation %3d | Best Fitness: %.6f | Evaluations: %d\n', gen, obj, eval_total);
        end

        % Stagnation detection
        if length(fitness_history) > 1
            delta = abs(fitness_history(end) - fitness_history(end-1));
            if delta < tolerance
                stagnation_counter = stagnation_counter + 1;
            else
                stagnation_counter = 0;
            end
        end

        % Early termination on stagnation
        if stagnation_counter >= stagnation_limit
            fprintf('Early stopping: stagnation detected at generation %d (Δ < %.1e for %d iterations).\n', gen, tolerance, stagnation_limit);
            break;
        end
    end

    %% === Final Output ===
    BESS_Output = reshape(round(gbest / 10) * 10, BESS_Number, 24);

    %% === Convergence Plot ===
    figure;
    plot(0:iter, fitness_history, 'b-', 'LineWidth', 2);
    xlabel('Generation', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Objective Value', 'FontSize', 14, 'FontWeight', 'bold');
    title('Convergence Trend of PSO Sizing Optimization', 'FontSize', 16, 'FontWeight', 'bold');
    grid on;

    % === Highlight last fitness value ===
    x_last = iter;
    y_last = fitness_history(end);  % fitness_history has length (iter + 1)

    hold on;
    plot(x_last, y_last, 'ro', 'MarkerFaceColor', 'r');  % Red circle marker
    text(x_last, y_last, sprintf('  %.2f', y_last), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 12, 'FontWeight', 'bold', ...
        'Color', 'red');
    hold off;

end
