function output = parameter_initialization(cell, SOC)
    num_of_breakpoints = length(SOC);

    switch cell
        case "EVE280"
            %% OCV
            OCV_init = 3.2 * ones(num_of_breakpoints, 1); 
            OCV_lb = 2.4 * ones(num_of_breakpoints, 1); 
            OCV_ub = 3.75 * ones(num_of_breakpoints, 1); 
            
            %% Resistances
            R = 2e-4 * ones(num_of_breakpoints, 1); 
            R_lb = 1e-5 * ones(num_of_breakpoints, 1); 
            R_ub = 1e-3 * ones(num_of_breakpoints, 1); 

            [R0_init, R1_init, R2_init, R3_init] = deal(R);
            [R0_lb, R1_lb, R2_lb, R3_lb] = deal(R_lb);
            [R0_ub, R1_ub, R2_ub, R3_ub] = deal(R_ub);

            %% Time constants
            %v4
            tau1_init = 2 * ones(num_of_breakpoints, 1);
            tau1_lb = 1 * ones(num_of_breakpoints, 1); 
            tau1_ub = 10 * ones(num_of_breakpoints, 1);
            
            tau2_init = 50 * ones(num_of_breakpoints, 1); 
            tau2_lb = 20 * ones(num_of_breakpoints, 1); 
            tau2_ub = 80 * ones(num_of_breakpoints, 1); 
            
            tau3_init = 1000 * ones(num_of_breakpoints, 1);  %1000
            tau3_lb = 500 * ones(num_of_breakpoints, 1);  %800
            tau3_ub = 1800 * ones(num_of_breakpoints, 1);  %1200

        case"SamsungE35"

            R = 3.5e-5 * ones(num_of_breakpoints, 1); 
            R_lb = 1e-6 * ones(num_of_breakpoints, 1); 
            R_ub = 1e-3 * ones(num_of_breakpoints, 1); 

            [R0_init, R1_init, R2_init, R3_init] = deal(R);
            [R0_lb, R1_lb, R2_lb, R3_lb] = deal(R_lb);
            [R0_ub, R1_ub, R2_ub, R3_ub] = deal(R_ub);

            OCV_init = 3.75 * ones(num_of_breakpoints, 1);
            OCV_lb = 2.5 * ones(num_of_breakpoints, 1);
            OCV_ub = 4.2 * ones(num_of_breakpoints, 1);
            
            tau1_init = 1 * ones(num_of_breakpoints, 1);
            tau1_lb = 0.1 * ones(num_of_breakpoints, 1);
            tau1_ub = 10 * ones(num_of_breakpoints, 1);
            
            tau2_init = 50 * ones(num_of_breakpoints, 1);
            tau2_lb = 100 * ones(num_of_breakpoints, 1);
            tau2_ub = 200 * ones(num_of_breakpoints, 1);
            
            tau3_init = 1000 * ones(num_of_breakpoints, 1);
            tau3_lb = 500 * ones(num_of_breakpoints, 1);
            tau3_ub = 2000 * ones(num_of_breakpoints, 1);
    end

    param_initialization = table(SOC, R0_init, R0_lb, R0_ub, R1_init, R1_lb, R1_ub, R2_init, R2_lb, R2_ub, R3_init, R3_lb, R3_ub, tau1_init, tau1_lb, tau1_ub, tau2_init, tau2_lb, tau2_ub, tau3_init, tau3_lb, tau3_ub, OCV_init, OCV_lb, OCV_ub);

    output = param_initialization; 
end

