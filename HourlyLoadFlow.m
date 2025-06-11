function [voltage, P_Loss_Kw, Q_Loss_KVAr, PL, QL, ld, node, branch, va, voltage_deviation, total_voltage_deviation] = ...
    HourlyLoadFlow(mm, ll, sel_pv, sel_lp, MVAb, Zb, BESS_Demand)            
        % Load data adjustment for the current hour
        ld = zeros(size(mm, 1), 3); % Initialize load data
        ld(:, 1) = mm(:, 1);
        ld(:, 2) = (mm(:, 2) * sel_lp) - ((mm(:, 4) * sel_pv) + BESS_Demand); % P net = Pload - Pgen
        ld(:, 3) = (mm(:, 3) * sel_lp) - mm(:, 5); % Q net = Qload - Qgen
        
        branch = length(ll); % Number of branches (lines)
        node = size(ld,1); % Number of nodes (buses)
        f = 0; % Initialize variable for tracking
        d = 0; % Initialize variable for tracking
        
        % Per unit Values Calculation
        for i = 1:branch
            R(i, 1) = (ll(i, 4)) / Zb; % Calculate per unit resistance
            X(i, 1) = (ll(i, 5)) / Zb; % Calculate per unit reactance
        end
        
        for i = 1:node
            P(i, 1) = ((ld(i, 2)) / (1000 * MVAb)); % Calculate per unit active power
            Q(i, 1) = ((ld(i, 3)) / (1000 * MVAb)); % Calculate per unit reactive power
        end
        
        % Create Connectivity Matrix
        C = zeros(branch, node); % Initialize connectivity matrix with zeros
        for i = 1:branch
            a = ll(i, 2); % Starting bus of the line
            b = ll(i, 3); % Ending bus of the line
            for j = 1:node
                if a == j
                    C(i, j) = -1; % Indicate starting bus in the matrix
                end
                if b == j
                    C(i, j) = 1; % Indicate ending bus in the matrix
                end
            end
        end
    
        % Find End Nodes
        e = 1; % Initialize counter for end nodes
        for i = 1:node
            d = 0; % Reset end node flag
            for j = 1:branch
                if C(j, i) == -1
                    d = 1; % Mark as not an end node if it has a starting connection
                end
            end
            if d == 0
                endnode(e, 1) = i; % Store end node
                e = e + 1; % Increment end node counter
            end
        end
    
        h = length(endnode); % Number of end nodes
        for j = 1:h
            e = 2; % Initialize index for path storage
            f = endnode(j, 1); % Current end node
            for s = 1:node
                if (f ~= 1) % If not the root bus
                    k = 1;  
                    for i = 1:branch
                        if (C(i, f) == 1 && k == 1)
                            f = i; % Move to parent node
                            k = 2; % Mark as visited
                        end
                    end
                    k = 1;
                    for i = 1:node
                        if (C(f, i) == -1 && k == 1)
                            f = i; % Find child node
                            g(j, e) = i; % Store child node
                            e = e + 1; % Increment index
                            k = 3; % Mark as visited
                        end            
                    end
                end
            end
        end
    
        for i = 1:h
            g(i, 1) = endnode(i, 1); % Store end node in the path matrix
        end
    
        w = length(g(1, :)); % Width of the path matrix
        for i = 1:h
            j = 1;
            for k = 1:node 
                for t = 1:w
                    if g(i, t) == k
                        g(i, t) = g(i, j); % Shift elements in the path
                        g(i, j) = k; % Set current node
                        j = j + 1; % Increment index
                    end
                end
            end
        end
    
        % Build the adjacency matrix
        for k = 1:branch
            e = 1;
            for i = 1:h
                for j = 1:w-1
                    if (g(i, j) == k) 
                        if g(i, j + 1) ~= 0
                            adjb(k, e) = g(i, j + 1); % Store adjacent node
                            e = e + 1;
                        else
                            adjb(k, 1) = 0; % No adjacent nodes
                        end
                    end
                end
            end
        end
    
        % Clean up the adjacency matrix
        for i = 1:branch-1
            for j = h:-1:1
                for k = j:-1:2
                    if adjb(i, j) == adjb(i, k-1)
                       adjb(i, j)  = 0; % Remove duplicate edges
                    end
                end
            end
        end
    
        x  = length(adjb(:, 1)); % Number of branches in the adjacency matrix
        ab = length(adjb(1, :)); % Number of nodes in the adjacency matrix
    
        % Shift and renumber the adjacency matrix
        for i = 1:x
            for j = 1:ab
                if adjb(i, j) == 0 && j ~= ab
                    if adjb(i, j + 1) ~= 0
                        adjb(i, j) = adjb(i, j + 1);
                        adjb(i, j + 1) = 0; % Shift elements left
                    end
                end
                if adjb(i, j) ~= 0
                   adjb(i, j)  = adjb(i, j) - 1; % Renumber nodes to zero-based indexing
                end
            end
        end
    
        for i = 1:x-1
            for j = 1:ab
                adjcb(i, j) = adjb(i + 1, j); % Create a compact adjacency matrix
            end
        end
        b = length(adjcb); % Length of the compact adjacency matrix
    
        % Voltage Calculation Program
        for i = 1:node
            vol_bus(i, 1) = 1.0; % Initialize bus voltages to 1.0 p.u.
        end
    
        % Iteratively calculate node voltages
        for s = 1:10
            for i = 1:node
                I_node(i, 1) = conj(complex(P(i, 1), Q(i, 1))) / (vol_bus(i, 1)); % Calculate load current
            end
            
            for i = 1:branch
                I_branch(i, 1) = I_node(i + 1, 1); % Branch current
            end
            
            xy = length(adjcb(1, :)); % Number of connections in the compact adjacency matrix
            for i = branch-1:-1:1
                for k = 1:xy
                    if adjcb(i, k) ~= 0
                        u = adjcb(i, k); % Get connected node
                        I_branch(i, 1) = I_branch(i, 1) + I_branch(u, 1); % Sum currents from connected nodes
                    end
                end      
            end
    
            for i = 2:node
                g = 0; % Initialize flag for voltage update
                for a = 1:b 
                    if xy > 1
                        if size(adjcb,2) >= 2 && adjcb(a, 2) == i - 1 
                            u = adjcb(a, 1);
                            vol_bus(i, 1) = (vol_bus(u, 1)) - (I_branch(i - 1, 1) * complex(R(i - 1, 1), X(i - 1, 1)));
                            g = 1; % Voltage updated from connected node
                        end
                        if size(adjcb,2) >= 3 && adjcb(a, 3) == i - 1 
                            u = adjcb(a, 1);
                            vol_bus(i, 1) = (vol_bus(u, 1)) - (I_branch(i - 1, 1) * complex(R(i - 1, 1), X(i - 1, 1)));
                            g = 1; % Voltage updated from connected node
                        end
                    end
                end
                if g == 0
                    if i > 1
                    vol_bus(i, 1) = (vol_bus(i - 1, 1)) - (I_branch(i - 1, 1) * complex(R(i - 1, 1), X(i - 1, 1))); % Update voltage from previous node
                end
                end
            end
        end
        Sr = (1:node)';
        endnode = zeros(node,1);
        g = zeros(node,node);
        adjb = zeros(branch,node);
        adjcb = zeros(branch-1,node); % Create a sequence for bus numbers for table xls

        % Power loss calculation
        PL = 0; % Active power loss
        QL = 0; % Reactive power loss
        
        for f = 1:branch
            Pl(f, 1) = (abs(I_branch(f, 1))^2) * R(f, 1); % Active power loss for branch f
            Ql(f, 1) = (abs(I_branch(f, 1))^2) * X(f, 1); % Reactive power loss for branch f
            PL = PL + Pl(f, 1); % Accumulate total active power loss
            QL = QL + Ql(f, 1); % Accumulate total reactive power loss
        end
    
        % Convert losses from per unit to kilowatts (kW) and kilovolt-amperes reactive (kVAR)
        P_Loss_Kw   = (Pl) * 100000; % Active power loss in kW
        Q_Loss_KVAr = (Ql) * 100000; % Reactive power loss in kVAR
        PL = (PL) * 100000; % Total active power loss in kW
        QL = (QL) * 100000; % Total reactive power loss in kVAR

        % Store calculated voltage and angle in a matrix
        vol_bus_p = [abs(vol_bus) angle(vol_bus) * 180 / pi]; % Voltage magnitude and angle
        % Populate 'va' array with voltage data
        for i=1:node
            va = zeros(node, 3);
        va(i,2:3)=vol_bus_p(i,1:2); % Store voltage magnitude and angle
        end
        for i=1:node
            va(i,1)=i; % Assign bus number
        end
        va; %

        % Extract voltage data
        voltage = vol_bus_p(:, 1);
        V_angle = vol_bus_p(:, 2) * (pi / 180); % Voltage angle in radians
        
        % Calculate the total active and reactive power losses
        sum(P_Loss_Kw); % Total active power loss
        sum(Q_Loss_KVAr); % Total reactive power loss
        
        % Store the total losses in the last index of the losses arrays
        P_Loss_Kw(branch, 1) = PL; % Store total active power loss
        Q_Loss_KVAr(branch, 1) = QL; % Store total reactive power loss

        % Calculate total absolute voltage deviation
        voltage_deviation = (voltage - 1.0); % Absolute deviation of voltage from nominal
        total_voltage_deviation = sum(abs(voltage_deviation)); % Sum of all deviations
end
