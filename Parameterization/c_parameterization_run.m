% Notes:
% Script assumes data current is such that +ve = charge & -ve = discharge. Corrected under "Make sure data is in the correct format"

clear

% Folders
folder_current = cd; 
folder_functions = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\Scripts\Parameterization\functions"; 
folder_data1 = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\7 - Pulse Test";
folder_data2 = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\11 - Drive Cycle (0p8C)\"; % where experimental data is saved
folder_data3 = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\6 - Drive cycle\";
folder_data4 = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\10 - Characterization\";
folder_result = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\Results\"; % where results from this script is saved
addpath(folder_current, folder_functions, folder_data1, folder_data2, folder_data3, folder_data4, genpath(folder_result));

% Test settings
cell = "EVE280";
deltaT = 0.1; 
switch cell
    case "EVE280"
        Vmax = 3.65; 
        Vmin = 2.5;
end

% Set data to optimize on and polynomial degrees to try 
dataSet = ["PUL_25degC_0p8_15min";];
polySet = ["12"; "13"; "14";"15";"16";"17";"18";"19";"20";"21";]; % ["7";"8";"9";"10";"11";"15";"18";"19";"21";];

% Other optimization settings
optimMethod = "PSO"; % Can use "LS", "PSO" or "GA"
constrainOCV = 0;   % Whether OCV is strictly decreasing

%%
for dataNum = 1:height(dataSet)

    % Load initial guess
    load("param_init.mat", "param"); % SOC from 1 to 0
    param = parameter_initialization(cell, param.SOC);

    for polyNum = 1:height(polySet)
        %% Load test profile & OCV polynomial (if applicable)
        
        clear segments
        
        % Load data set to optimize on
        dataSample = dataSet(dataNum); % Name of the data optimized on
        load(dataSample + ".mat", "meas")
        dataType = "driveCycle";

        % Select data type & delete unnecessary data
        if dataType == "pulse" 
            %% Pulse is pre divided into segments
            segments(:, 3) = [pulse.segment2]; 
            SOC = param.SOC;
        else
            %% Divide into segments by SOC

            % Delete unnecessary fields
            try
                fields = {'TimeStamp','StepTime','Procedure','Wh','Power','Battery_Temp_degC'};
                meas = rmfield(meas, fields);
            catch
                % Do nothing
            end

            % Entire profile to be optimized
            try
                meas_optim = struct2table(meas);
            catch
                meas_optim = meas;
            end

            % Delete inf & NaN
            meas_optim = meas_optim( ~any( isnan( meas_optim.Time ) | isinf( meas_optim.Time ), 2 ),: );
            meas_optim = meas_optim( ~any( isnan( meas_optim.Current ) | isinf( meas_optim.Current ), 2 ),: );
            meas_optim = meas_optim( ~any( isnan( meas_optim.Voltage ) | isinf( meas_optim.Voltage ), 2 ),: );

            % Delete anything after Vmin is reached
            ind = find(meas_optim.Voltage <= Vmin, 1);
            if ~isempty(ind)
                meas_optim = meas_optim(1:ind, :);
            end

            % Make sure data is in the correct format
            Q = abs(meas_optim.Ah(end)-meas_optim.Ah(1));
            meas_optim.Current = -meas_optim.Current; % Assumes input data current has opposite sign convension
            meas_optim.Ah = meas_optim.Ah - meas_optim.Ah(1);
            meas_optim.SOC = 1 - (-meas_optim.Ah)/Q;
            meas_optim.Time = meas_optim.Time - meas_optim.Time(1);

            % Divide data into segments by SOC
            SOC = param.SOC;
            validateattributes(SOC, "double", "decreasing")
            for i = 1:height(SOC)-1
                i1 = find(meas_optim.SOC <= SOC(i), 1);
                i2 = find(meas_optim.SOC <= SOC(i+1), 1);
                iStart = min(i1, i2); 
                iEnd = max(i1, i2); 
                if isempty(iEnd)
                    iEnd = height(meas_optim);
                    segments{i, 1} = [SOC(i), SOC(i+1)];
                    segments{i, 2} = [iStart, iEnd];
                    segments{i, 3} = meas_optim(iStart:iEnd, :);
                    break
                end
                segments{i, 1} = [SOC(i), SOC(i+1)];
                segments{i, 2} = [iStart, iEnd];
                segments{i, 3} = meas_optim(iStart:iEnd, :);
            end

            % Plot segments for validation
            figure("WindowStyle", "docked"); hold on
            for i = 1:height(segments)
                plot(segments{i, 3}.Time, segments{i, 3}.Voltage)
            end

            % Report error is segments doesn't exist
            if isempty(segments{1, 3})
                error("CHECK SEGMENTS")
            end

        end

        % Select OCV polynomial
        polyDeg = polySet(polyNum);
        switch polyDeg
            case "0"
                useOCVpolynomial = 0;
                coeff = 0;
            otherwise 
                useOCVpolynomial = 1;
                load("OCVpoly_" + polyDeg + "_25degC.mat", "coeff") % OCV coefficients (7, 8, 9, 10, 11, 15, 18, 19, 21)
        end
        
        clear meas meas_optim

        %% Optimize each segment

        for i = 1:height(segments)
            %% Load ECM parameters for segment i
        
            % Parameters for segment i
            if i > 1

                % Load parameters from prev segment
                filename = string(folder_result) + "Parameters_Optimized_" + dataSample + "_" + optimMethod + "_Poly" + polyDeg + "_" + string(i-1) + ".mat";
                load(filename)
                initStates = endStates; 

                % Set parameters for current segment
                param_SOC = [param.SOC(i-1); param.SOC(i)];
                param_prev = [param.OCV_opt(i-1); param.R0_opt(i-1); param.R1_opt(i-1); param.R2_opt(i-1); param.R3_opt(i-1); param.tau1_opt(i-1); param.tau2_opt(i-1); param.tau3_opt(i-1)]; 
                param_opt = [param.OCV_init(i); param.R0_init(i); param.R1_init(i); param.R2_init(i); param.R3_init(i); param.tau1_init(i); param.tau2_init(i); param.tau3_init(i)]; 
                lb = [param.OCV_lb(i); param.R0_lb(i); param.R1_lb(i); param.R2_lb(i); param.R3_lb(i); param.tau1_lb(i); param.tau2_lb(i); param.tau3_lb(i)]; 
                ub = [param.OCV_ub(i); param.R0_ub(i); param.R1_ub(i); param.R2_ub(i); param.R3_ub(i); param.tau1_ub(i); param.tau2_ub(i); param.tau3_ub(i)]; 
                
                if constrainOCV
                    % Whether OCV is constrained to be strictly decreasing
                    ub(1) = min(param.OCV_ub(i), param.OCV_opt(i-1)); 
                end

                if ub(1) <= lb(1)
                    error("Check here"); 
                end

            elseif i == 1
                param_SOC = [param.SOC(1); param.SOC(2)];
                param_opt = [param.OCV_init(1); param.R0_init(1); param.R1_init(1); param.R2_init(1); param.R3_init(1); param.tau1_init(1); param.tau2_init(1); param.tau3_init(1);...
                                param.OCV_init(2); param.R0_init(2); param.R1_init(2); param.R2_init(2); param.R3_init(2); param.tau1_init(2); param.tau2_init(2); param.tau3_init(2)]; 
                lb = [param.OCV_lb(1); param.R0_lb(1); param.R1_lb(1); param.R2_lb(1); param.R3_lb(1); param.tau1_lb(1); param.tau2_lb(1); param.tau3_lb(1);...
                        param.OCV_lb(2); param.R0_lb(2); param.R1_lb(2); param.R2_lb(2); param.R3_lb(2); param.tau1_lb(2); param.tau2_lb(2); param.tau3_lb(2)]; 
                ub = [param.OCV_ub(1); param.R0_ub(1); param.R1_ub(1); param.R2_ub(1); param.R3_ub(1); param.tau1_ub(1); param.tau2_ub(1); param.tau3_ub(1);...
                        param.OCV_ub(2); param.R0_ub(2); param.R1_ub(2); param.R2_ub(2); param.R3_ub(2); param.tau1_ub(2); param.tau2_ub(2); param.tau3_ub(2);]; 

                % Initial states
                initStates.I_RC = [0, 0, 0];
                initStates.Vt = segments{i, 3}.Voltage(1);
                initStates.SOC = segments{1, 3}.SOC(1);
            end

            %% Optimize

            seg_optim = segments{i, 3};

            BatteryModel = @Model_3RC; 
            
            if i > 1
                fun = @(param_opt)parameterization_Objective(param_opt, param_prev, param_SOC, BatteryModel, seg_optim, initStates, Q, deltaT, optimMethod, coeff);
            elseif i == 1
                fun = @(param_opt)parameterization_Objective(param_opt, [], param_SOC, BatteryModel, seg_optim, initStates, Q, deltaT, optimMethod, coeff);
            end

            x = parameterization_Optimize(fun, param_opt, lb, ub, optimMethod); 

            %% Plot Vt before & after optimization

            if i > 1
                parameters_before.SOC = param_SOC; 
                parameters_before.OCV = [param_prev(1); param_opt(1)];
                parameters_before.R0 = [param_prev(2); param_opt(2)];
                parameters_before.R1 = [param_prev(3); param_opt(3)];
                parameters_before.R2 = [param_prev(4); param_opt(4)];
                parameters_before.R3 = [param_prev(5); param_opt(5)];
                parameters_before.T1 = [param_prev(6); param_opt(6)];
                parameters_before.T2 = [param_prev(7); param_opt(7)];
                parameters_before.T3 = [param_prev(8); param_opt(8)];

                parameters_after.SOC = param_SOC; 
                parameters_after.OCV = [param_prev(1); x(1)];
                parameters_after.R0 = [param_prev(2); x(2)];
                parameters_after.R1 = [param_prev(3); x(3)];
                parameters_after.R2 = [param_prev(4); x(4)];
                parameters_after.R3 = [param_prev(5); x(5)];
                parameters_after.T1 = [param_prev(6); x(6)];
                parameters_after.T2 = [param_prev(7); x(7)];
                parameters_after.T3 = [param_prev(8); x(8)];
            elseif i == 1
                parameters_before.SOC = param_SOC;
                parameters_before.OCV = [param_opt(1); param_opt(9)];
                parameters_before.R0 = [param_opt(2); param_opt(10)];
                parameters_before.R1 = [param_opt(3); param_opt(11)];
                parameters_before.R2 = [param_opt(4); param_opt(12)];
                parameters_before.R3 = [param_opt(5); param_opt(13)];
                parameters_before.T1 = [param_opt(6); param_opt(14)];
                parameters_before.T2 = [param_opt(7); param_opt(15)];
                parameters_before.T3 = [param_opt(8); param_opt(16)];

                parameters_after.SOC = param_SOC; 
                parameters_after.OCV = [x(1); x(9)];
                parameters_after.R0 = [x(2); x(10)];
                parameters_after.R1 = [x(3); x(11)];
                parameters_after.R2 = [x(4); x(12)];
                parameters_after.R3 = [x(5); x(13)];
                parameters_after.T1 = [x(6); x(14)];
                parameters_after.T2 = [x(7); x(15)];
                parameters_after.T3 = [x(8); x(16)];
            end

            [Vt_before, ~, ~] = Model_3RC(seg_optim.Current, initStates, parameters_before, Q, deltaT, coeff); 
            [Vt_after, ~, endStates] = Model_3RC(seg_optim.Current, initStates, parameters_after, Q, deltaT, coeff); 

            % Plot the measured and simulated data.
            figName = 'Pulse ' + string(i);
            figure('Name', figName, 'WindowStyle', 'docked'); hold on
            plot(seg_optim.Time, seg_optim.Voltage, 'color', '#0072BD')
            plot(seg_optim.Time, Vt_before, 'color', '#D95319');
            plot(seg_optim.Time, Vt_after, 'color', '#77AC30');
            title(append('Simulated and Measured Responses Before Estimation (Pulse ', string(i), ')'))
            legend('Measured Vt', 'Before Optim', 'After Optim');

            %% Update parameters

            % Add optimized parameters to table "param"
            if i > 1
                param.OCV_opt(i+1) = x(1);
                param.R0_opt(i+1) = x(2);
                param.R1_opt(i+1) = x(3);
                param.R2_opt(i+1) = x(4);
                param.R3_opt(i+1) = x(5);
                param.tau1_opt(i+1) = x(6);
                param.tau2_opt(i+1) = x(7);
                param.tau3_opt(i+1) = x(8);
            elseif i == 1
                param.OCV_opt(i) = x(1);
                param.R0_opt(i) = x(2);
                param.R1_opt(i) = x(3);
                param.R2_opt(i) = x(4);
                param.R3_opt(i) = x(5);
                param.tau1_opt(i) = x(6);
                param.tau2_opt(i) = x(7);
                param.tau3_opt(i) = x(8);

                param.OCV_opt(i+1) = x(9);
                param.R0_opt(i+1) = x(10);
                param.R1_opt(i+1) = x(11);
                param.R2_opt(i+1) = x(12);
                param.R3_opt(i+1) = x(13);
                param.tau1_opt(i+1) = x(14);
                param.tau2_opt(i+1) = x(15);
                param.tau3_opt(i+1) = x(16);
            end

            % Save optimized parameters
            filename = string(folder_result) + "Parameters_Optimized_" + dataSample + "_" + optimMethod + "_Poly" + polyDeg + "_" + string(i) + ".mat";
            save(filename, 'param', 'initStates', 'endStates')

        end

    end
end
