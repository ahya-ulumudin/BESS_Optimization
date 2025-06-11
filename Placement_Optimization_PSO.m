function [BESS_Location, BESS_Output, obj, fitness_history, fitness_eval, iter, eval_total, best_iter_idx, best_eval_idx] = ...
    Placement_Optimization_PSO(lower_bound, upper_bound, BESS_Number, Candidate_Buses, ...
    hour, mm, ll, sel_pv, sel_lp, MVAb, Zb, stagnation_limit, Bus_Ranked, use_ranked_bus_guidance, guided_fraction)
% PLACEMENT_OPTIMIZATION_PSO_R1
% Performs Particle Swarm Optimization for BESS placement and sizing

    %% === PSO Parameters ===
    max_gen  = 300;                          % Maximum number of generations
    pop_size = 100;                          % Population size (number of particles)
    w_max = 0.7;                             % Maximum inertia weight
    w_min = 0.4;                             % Minimum inertia weight
    num_candidates = length(Candidate_Buses);% Number of candidate buses
    num_variables = BESS_Number * 2;         % Total decision variables (locations + outputs)
    v_max = (upper_bound - lower_bound) / 2; % Maximum velocity constraint
    tolerance = 1e-6;                        % Tolerance for stagnation detection

    %% === Global Evaluation Tracking Initialization ===
    eval_total = 0;
    eval_best = 0;
    best_iter_idx = 0;
    best_eval_idx = 0;
    fitness_eval = [];

    %% === Particle Position and Velocity Initialization ===
    position = zeros(pop_size, num_variables);
    velocity = zeros(pop_size, num_variables);
    fitness = zeros(1, pop_size);
    guided_pop = round(guided_fraction * pop_size);

    for i = 1:pop_size
        if use_ranked_bus_guidance && i <= guided_pop
            % === Guided Initialization using Bus_Ranked ===
            Ranked_Buses = Bus_Ranked;  % initial_solution is Bus_Ranked
    
            % Pick top BESS_Number buses and transpose to row
            selected_buses = Ranked_Buses(1:BESS_Number);
    
            % Generate random demand for selected buses
            selected_outputs = lower_bound + (upper_bound - lower_bound) * rand(1, BESS_Number);
            selected_outputs = round(selected_outputs / 10) * 10;
            selected_outputs = max(lower_bound, min(upper_bound, selected_outputs));
    
            % Assign to position
            position(i, :) = [selected_buses, selected_outputs];
    
        else
            % === Random Initialization ===
            loc_idx = randi(num_candidates, 1, BESS_Number);
            locations = Candidate_Buses(loc_idx);
    
            outputs = lower_bound + (upper_bound - lower_bound) * rand(1, BESS_Number);
            outputs = round(outputs / 10) * 10;
            outputs = max(lower_bound, min(upper_bound, outputs));
    
            position(i, :) = [locations, outputs];
        end
    end

    %% === Initial Fitness Evaluation ===
    for i = 1:pop_size
        loc = round(position(i, 1:BESS_Number));
        loc = Candidate_Buses(min(max(loc, 1), num_candidates));
        loc = unique(loc);
        while length(loc) < BESS_Number
            new_bus = Candidate_Buses(randi(num_candidates));
            if ~ismember(new_bus, loc)
                loc(end+1) = new_bus;
            end
        end
        loc = sort(loc);

        out = round(position(i, BESS_Number+1:end) / 10) * 10;
        out = max(lower_bound, min(upper_bound, out));

        fitness(i) = Placement_Objective(loc, out, mm, ll, sel_pv, sel_lp, MVAb, Zb);
        eval_total = eval_total + 1;

        [temp_best, ~] = min(fitness);
        fitness_eval(end+1) = temp_best;
    end

    %% === Initialize Personal Bests (pbest) and Global Best (gbest) ===
    pbest = position;
    pbest_fitness = fitness;
    [gbest_fitness, gbest_idx] = min(fitness);
    gbest = position(gbest_idx, :);
    prev_gbest_fitness = gbest_fitness;
    eval_best = eval_total;
    best_eval_idx = eval_total;
    best_iter_idx = 0;

    fitness_history = zeros(1, max_gen);
    stagnation_count = 0;
    termination_logged = false;

    %% === Main PSO Iterative Loop ===
    for gen = 1:max_gen
        % --- Adaptive inertia weight and acceleration coefficients ---
        w = w_max - ((w_max - w_min) * (gen / max_gen)^2); % Quadratic decreasing inertia
        c1 = 2.0 - 1.0 * (gen / max_gen);                  % Cognitive component decreases
        c2 = 1.0 + 1.5 * (gen / max_gen);                  % Social component increases

        % --- Particle Update Loop ---
        for i = 1:pop_size
            r1 = rand(1, num_variables);
            r2 = rand(1, num_variables);

            % Update velocity and clamp
            velocity(i,:) = w * velocity(i,:) ...
                          + c1 .* r1 .* (pbest(i,:) - position(i,:)) ...
                          + c2 .* r2 .* (gbest - position(i,:));
            velocity(i,:) = max(-v_max, min(v_max, velocity(i,:)));

            % Update position and clamp
            position(i,:) = position(i,:) + velocity(i,:);
            position(i,:) = max(lower_bound, min(upper_bound, position(i,:)));

            % Decode the particle
            loc = round(position(i, 1:BESS_Number));
            loc = Candidate_Buses(min(max(loc, 1), num_candidates));
            loc = unique(loc);
            while length(loc) < BESS_Number
                new_bus = Candidate_Buses(randi(num_candidates));
                if ~ismember(new_bus, loc)
                    loc(end+1) = new_bus;
                end
            end
            loc = sort(loc);

            out = round(position(i, BESS_Number+1:end) / 10) * 10;
            out = max(lower_bound, min(upper_bound, out));

            % Fitness evaluation
            fit = Placement_Objective(loc, out, mm, ll, sel_pv, sel_lp, MVAb, Zb);
            eval_total = eval_total + 1;

            % Update pbest and gbest
            if fit < pbest_fitness(i)
                pbest(i,:) = position(i,:);
                pbest_fitness(i) = fit;
                if fit < gbest_fitness
                    gbest = position(i,:);
                    gbest_fitness = fit;
                    eval_best = eval_total;
                    best_eval_idx = eval_total;
                    best_iter_idx = gen;
                end
            end
            fitness_eval(end+1) = gbest_fitness;
        end

        % Record best fitness at generation
        fitness_history(gen) = gbest_fitness;

        % --- Stagnation Detection ---
        delta = abs(gbest_fitness - prev_gbest_fitness);
        if delta < tolerance
            stagnation_count = stagnation_count + 1;
        else
            stagnation_count = 0;
        end
        prev_gbest_fitness = gbest_fitness;

        % --- Early Termination ---
        if stagnation_count >= stagnation_limit
            fprintf('[Hour %02d] [PSO] Terminated due to stagnation at generation %3d | Final fitness: %.6f\n', hour, gen, gbest_fitness);
            termination_logged = true;
            break;
        end
    end

    if ~termination_logged
        fprintf('[Hour %02d] [PSO] Final Fitness: %.6f\n', hour, gbest_fitness);
    end

    %% === Decode Final Best Solution ===
    BESS_Location = round(gbest(1:BESS_Number));
    BESS_Location = Candidate_Buses(min(max(BESS_Location, 1), num_candidates));
    BESS_Location = unique(BESS_Location);
    while length(BESS_Location) < BESS_Number
        new_bus = Candidate_Buses(randi(num_candidates));
        if ~ismember(new_bus, BESS_Location)
            BESS_Location(end+1) = new_bus;
        end
    end
    BESS_Location = sort(BESS_Location);

    BESS_Output = round(gbest(BESS_Number+1:end) / 10) * 10;
    BESS_Output = max(lower_bound, min(upper_bound, BESS_Output));

    obj = gbest_fitness;
    iter = gen;
end