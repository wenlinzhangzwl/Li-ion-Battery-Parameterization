clear

% Folders
folder_current = cd; 
folder_project = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design';
folder_functions = append(folder_project, '\Scripts\Parameterization\functions\'); 
folder_data = append(folder_project, '\Test Results'); % where experimental data is saved
folder_result = append(folder_project, '\Test Results\7 - Pulse Test\Results\'); % where results from this script is saved
addpath(folder_current);
addpath(folder_functions);
addpath(genpath(folder_data));
addpath(folder_result);

% Load initial guess from pulse discharge
load('param_Pulse0p1.mat')

% Set upper & lower bounds
param.OCV_ub =  param.OCV_init*1.005;
param.OCV_lb =  param.OCV_init*0.998;
param.R0_ub =  param.R0_init*1.5;
param.R0_lb =  param.R0_init*0.5;
param.R1_ub =  param.R1_init*1.5;
param.R1_lb =  param.R1_init*0.5;
param.R2_ub =  param.R2_init*1.5;
param.R2_lb =  param.R2_init*0.5;
param.R3_ub =  param.R3_init*1.5;
param.R3_lb =  param.R3_init*0.5;
param.tau1_init(:) =  2;
param.tau1_ub(:) =  3;
param.tau1_lb(:) =  1;
param.tau2_init(:) =  50;
param.tau2_ub(:) =  60;
param.tau2_lb(:) =  40;
param.tau3_ub =  param.tau3_init*1.5;
param.tau3_lb =  param.tau3_init*0.5;

% Validate OCV ub is greater than OCV lb 
OCVdiff = param.OCV_ub - param.OCV_lb; 
errorInd = find(OCVdiff<=0); 
if exist(errorInd)
    error('Check ub & lb of OCV (see errorInd)');
end

% Set data to optimize on and polynomial degrees to try 
% dataSet = ["UDDS"; "US06"; "HWFET"; "MIXED"; "NEDC"; "Pulse"];
% polySet = ["0";];

% dataSet = ["US06"];
% polySet = ["7";"8";"9";"10";"11";"15";"18";"19";"21";];

dataSet = ["Pulse0p1"];
polySet = ["0"];
% reStart = 'yes'; 

for dataNum = 1:height(dataSet)
    
    clearvars -except dataNum dataSet folder_current folder_data_DriveCycle folder_data_OCVtest folder_data_pulseTest...
                      folder_functions folder_result param polySet reStart
    
    for polyNum = 1:height(polySet)
        %% Set up input data & parameters
        
        clear segments
        
        % Select data set to optimize on
        dataSample = dataSet(dataNum);
        switch dataSample
            case "UDDS"
                load('UDDS25.mat');
                dataType = "driveCycle";
            case "US06"
                load('US0625.mat');
                dataType = "driveCycle";
            case "HWFET"
                load('HWFET25.mat');
                dataType = "driveCycle";
            case "MIXED"
                load('MIX25.mat');
                dataType = "driveCycle";
            case "NEDC"
                load('NEDC23.mat');
                dataType = "driveCycle";
            case "Pulse0p8"
                load('PROCESSED_PUL25dch.mat');
                dataType = 'pulse';
            case "Pulse0p1"
                load('PROCESSED_PUL_0p1_discharge.mat');
                dataType = 'pulse';
        end

        % Select data type
        switch dataType
            case "driveCycle"
                % Delete unnecessary fields
                fields = {'TimeStamp','StepTime','Procedure','Wh','Power','Battery_Temp_degC'};
                meas = rmfield(meas, fields); % Delete unnecessary fields
                
                % Entire profile to be optimized
                meas_optim = struct2table(meas);

                % Delete last entry
                meas_optim = meas_optim(1:end-1, :);

                % Delete inf & NaN
                meas_optim = meas_optim( ~any( isnan( meas_optim.Time ) | isinf( meas_optim.Time ), 2 ),: );
                meas_optim = meas_optim( ~any( isnan( meas_optim.Current ) | isinf( meas_optim.Current ), 2 ),: );
                meas_optim = meas_optim( ~any( isnan( meas_optim.Voltage ) | isinf( meas_optim.Voltage ), 2 ),: );

                % Calculate parameters
                ind = find(meas_optim.Voltage <=2.5, 1);
                if ~isempty(ind)
                    meas_optim = meas_optim(1:ind, :);
                end
                Q = abs(meas_optim.Ah(end)-meas_optim.Ah(1));
                meas_optim.Current = -meas_optim.Current;
                meas_optim.Ah = meas_optim.Ah - meas_optim.Ah(1);
                meas_optim.SOC = 1 - (-meas_optim.Ah)/Q;
                meas_optim.Time = meas_optim.Time - meas_optim.Time(1);

                % Divide data into segments by SOC
                SOC = param.SOC;
                SOC_flip = flip(SOC);
                for i = 1:height(SOC_flip)-1
                    i1 = find(meas_optim.SOC <= SOC_flip(i), 1);
                    i2 = find(meas_optim.SOC <= SOC_flip(i+1), 1);
                    iStart = min(i1, i2); 
                    iEnd = max(i1, i2); 
                    if isempty(iEnd)
                        iEnd = height(meas_optim);
                        segments{i, 1} = [SOC_flip(i), SOC_flip(i+1)];
                        segments{i, 2} = [iStart, iEnd];
                        segments{i, 3} = meas_optim(iStart:iEnd, :);
                        break
                    end
                    segments{i, 1} = [SOC_flip(i), SOC_flip(i+1)];
                    segments{i, 2} = [iStart, iEnd];
                    segments{i, 3} = meas_optim(iStart:iEnd, :);
                end
            case 'pulse'
                segments(:, 3) = [pulse.segment2]; 
                clear meas_resampled meas_t pulse
                SOC = param.SOC;
        end
        
        if isempty(segments{1, 3})
            error("CHECK SEGMENTS")
        end
        
        % Select OCV polynomial
        polyDeg = polySet(polyNum);
        switch polyDeg
            case "0"
                useOCVpolynomial = 0;
                coeff = 0;
            otherwise 
                useOCVpolynomial = 1;
                load("C20CCCV_OCVSOC" + polyDeg + ".mat") % OCV coefficients (7, 8, 9, 10, 11, 15, 18, 19, 21)
                clear coeff_chg coeff_dch OCV_SOC
        end
        
        clear meas meas_optim fields
        
        %% Optimize with drive cycle
        
        if exist("reStart", "var")
            clear param SOC
            load('param_opt_Pulse_Poly0_26.mat') %**********************************************CHANGE
            param = flip(param);
            startInd = 27; %**********************************************CHANGE
        else
            startInd = 1; 
        end

        % Set optimization model
        model = 'Parameterization_DriveCycle';
        open(model)
        
        param = flip(param);
        for i = startInd:height(segments)
            param_prev = flip(param(1:i, :));
            param_current = param(i+1, :);
        
            %% Set Model Parameters
            % Parameter look up tables with initial guesses
            if i ~= 1
                SOC_prev    = param_prev.SOC;
                OCV_prev    = param_prev.OCV_opt;
                R0_prev     = param_prev.R0_opt;
                R1_prev     = param_prev.R1_opt;
                R2_prev     = param_prev.R2_opt;
                R3_prev     = param_prev.R3_opt;
                tau1_prev   = param_prev.tau1_opt;
                tau2_prev   = param_prev.tau2_opt;
                tau3_prev   = param_prev.tau3_opt;

                SOC_current  = param_current.SOC;
                OCV_current  = param_current.OCV_init;
                R0_current   = param_current.R0_init;
                R1_current   = param_current.R1_init;
                R2_current   = param_current.R2_init;
                R3_current   = param_current.R3_init;
                tau1_current = param_current.tau1_init;
                tau2_current = param_current.tau2_init;
                tau3_current = param_current.tau3_init;
            elseif i == 1
                SOC_prev    = param_prev.SOC;
                OCV_prev    = param_prev.OCV_init;
                R0_prev     = param_prev.R0_init;
                R1_prev     = param_prev.R1_init;
                R2_prev     = param_prev.R2_init;
                R3_prev     = param_prev.R3_init;
                tau1_prev   = param_prev.tau1_init;
                tau2_prev   = param_prev.tau2_init;
                tau3_prev   = param_prev.tau3_init;

                SOC_current  = param_current.SOC;
                OCV_current  = param_current.OCV_init;
                R0_current   = param_current.R0_init;
                R1_current   = param_current.R1_init;
                R2_current   = param_current.R2_init;
                R3_current   = param_current.R3_init;
                tau1_current = param_current.tau1_init;
                tau2_current = param_current.tau2_init;
                tau3_current = param_current.tau3_init;

                % Initial states
                init_I1 = 0; 
                init_I2 = 0;
                init_I3 = 0; 
                init_Q = Q;
                init_Vt = [0, 0, 0];
                init_SOC = segments{1, 3}.SOC(1);
            end

            % Data
            seg_optim = segments{i, 3};
            seg_optim.Time = seg_optim.Time - seg_optim.Time(1);
            Time = seg_optim.Time;

            % Model inputs
            t = [Time Time];
            meas_I = [Time seg_optim.Current];
            meas_Vt = [Time seg_optim.Voltage];
            clear seg_optim

            % Simulation settings
            stepSize = 0.1;
            runTime = Time(end) - Time(1);

            %% Define the Estimation Experiment

            % Create an experiment object
            Exp = sdo.Experiment(model);

            % Define the signal to estimate
            sim_Vt           = Simulink.SimulationData.Signal;
            sim_Vt.Name      = 'Vt_sim';
            sim_Vt.BlockPath = append(model, '/StateUpdate');
            sim_Vt.PortType  = 'outport';
            sim_Vt.PortIndex = 2;
            sim_Vt.Values    = timeseries(meas_Vt(:,2),meas_Vt(:,1));
            Exp.OutputData   = sim_Vt;
            clear sim_Vt

            % Add the initial state from the model & set them to non-tunable
            Exp.InitialStates = sdo.getStateFromModel(model);
            for j = 1:5
                Exp.InitialStates(j).Free = 0;
            end

            %% Compare the Measured Output and the Initial Simulated Output

            % Create a simulation using the experiment
            Simulator = createSimulator(Exp);
            Simulator = sim(Simulator);

            % Search for the Vt_sim signal in the logged simulation data.
            SimLog = find(Simulator.LoggedData,get_param(model,'SignalLoggingName'));
            signal_Vt = find(SimLog,'Vt_sim');

            % Plot the measured and simulated data.
            figName = 'Pulse ' + string(i) + '(before)';
            figure('Name', figName, 'WindowStyle', 'docked'); hold on
            plot(meas_Vt(:, 1), meas_Vt(:, 2), '.-')
            plot(signal_Vt.Values.Time,signal_Vt.Values.Data, '.--');
            title(append('Simulated and Measured Responses Before Estimation (Pulse ', string(i), ')'))
            legend('Measured Vt', 'Simulated Vt');
            clear signal_Vt

            %% Specify the Parameters to Estimate & Their Upper/Lower Bounds

            % Get parameters to estimate ('p') & set upper/lower bounds
            if i ~= 1

                p = sdo.getParameterFromModel(model,{'R0_current','R1_current','R2_current','R3_current','tau1_current','tau2_current','tau3_current',...
                                                         'OCV_current'});
                % R0_current
                p(1).Minimum = param_current.R0_lb;        
                p(1).Maximum = param_current.R0_ub;
                % R1_current
                p(2).Minimum = param_current.R1_lb;       
                p(2).Maximum = param_current.R1_ub;
                % R2_current
                p(3).Minimum = param_current.R2_lb;       
                p(3).Maximum = param_current.R2_ub; 
                % R3_current
                p(4).Minimum = param_current.R3_lb;   
                p(4).Maximum = param_current.R3_ub; 
                % tau1_current
                p(5).Minimum = param_current.tau1_lb;   
                p(5).Maximum = param_current.tau1_ub; 
                % tau2_current
                p(6).Minimum = param_current.tau2_lb;  
                p(6).Maximum = param_current.tau2_ub; 
                % tau3_current
                p(7).Minimum = param_current.tau3_lb;  
                p(7).Maximum = param_current.tau3_ub; 
                % OCV_current
                p(8).Minimum = param_current.OCV_lb;
                p(8).Maximum = min(param_prev.OCV_opt(1), param_current.OCV_ub);
                if p(8).Maximum <= p(8).Minimum
        %             p(8).Maximum = p(8).Minimum; 
                    pause
                end

            elseif i == 1

                p = sdo.getParameterFromModel(model,{'R0_prev','R1_prev','R2_prev','R3_prev','tau1_prev','tau2_prev','tau3_prev',...
                                                         'R0_current','R1_current','R2_current','R3_current','tau1_current','tau2_current','tau3_current',...
                                                         'OCV_prev','OCV_current',});      
                % R0_prev
                p(1).Minimum = param_prev.R0_lb;        
                p(1).Maximum = param_prev.R0_ub;
                % R1_prev
                p(2).Minimum = param_prev.R1_lb;      
                p(2).Maximum = param_prev.R1_ub;
                % R2_prev
                p(3).Minimum = param_prev.R2_lb;      
                p(3).Maximum = param_prev.R2_ub;
                % R3_prev
                p(4).Minimum = param_prev.R3_lb;  
                p(4).Maximum = param_prev.R3_ub;
                % tau1_prev
                p(5).Minimum = param_prev.tau1_lb;   
                p(5).Maximum = param_prev.tau1_ub;
                % tau2_prev
                p(6).Minimum = param_prev.tau2_lb;  
                p(6).Maximum = param_prev.tau2_ub;
                % tau3_prev
                p(7).Minimum = param_prev.tau3_lb;  
                p(7).Maximum = param_prev.tau3_ub;

                % R0_current
                p(8).Minimum = param_current.R0_lb;        
                p(8).Maximum = param_current.R0_ub;
                % R1_current
                p(9).Minimum = param_current.R1_lb;       
                p(9).Maximum = param_current.R1_ub;
                % R2_current
                p(10).Minimum = param_current.R2_lb;       
                p(10).Maximum = param_current.R2_ub; 
                % R3_current
                p(11).Minimum = param_current.R3_lb;   
                p(11).Maximum = param_current.R3_ub; 
                % tau1_current
                p(12).Minimum = param_current.tau1_lb;   
                p(12).Maximum = param_current.tau1_ub; 
                % tau2_current
                p(13).Minimum = param_current.tau2_lb;  
                p(13).Maximum = param_current.tau2_ub; 
                % tau3_current
                p(14).Minimum = param_current.tau3_lb;  
                p(14).Maximum = param_current.tau3_ub; 

                % OCV_prev
                p(15).Minimum = param_prev.OCV_lb;
                p(15).Maximum = param_prev.OCV_ub;
                % OCV_current
                p(16).Minimum = param_current.OCV_lb;
                p(16).Maximum = param_current.OCV_ub;

            end

            %% Define the Estimation Objective Function & Estimate the Parameters

            estFcn = @(v) ElectricalModel_estimation_Objective(v, Simulator, Exp, model);

            optOptions = sdo.OptimizeOptions(...
                'OptimizedModel',Simulator,...
                'Method','lsqnonlin',...
                'UseParallel',true);
            optOptions.MethodOptions.FunctionTolerance = 0.01;
            optOptions.MethodOptions.OptimalityTolerance = 0.01;
            optOptions.MethodOptions.StepTolerance = 0.1;
            optOptions.MethodOptions.MaxIterations = 10;

            % Estimate the parameters.
            vOpt = sdo.optimize(estFcn, p, optOptions);

            %% Compare the Measured Output and the Final Simulated Output

            % Update the experiments with the estimated parameter values.
            Exp = setEstimatedValues(Exp,vOpt);

            % Simulate the model using the updated experiment
            Simulator   = createSimulator(Exp,Simulator);
            Simulator   = sim(Simulator);
            SimLog      = find(Simulator.LoggedData,get_param(model,'SignalLoggingName'));
            signal_Vt   = find(SimLog,'Vt_sim');
            signal_SOC  = find(SimLog,'SOC_sim');
            signal_I1   = find(SimLog,'I1_sim');
            signal_I2   = find(SimLog,'I2_sim');
            signal_I3   = find(SimLog,'I3_sim');

            % Plot the measured and simulated data.
            figName = 'Pulse ' + string(i) + '(after)';
            figure('Name', figName, 'WindowStyle', 'docked');
            plot(meas_Vt(:, 1), meas_Vt(:, 2), '.-', ...
                signal_Vt.Values.Time,signal_Vt.Values.Data,'.--');
            title(append('Simulated and Measured Responses After Estimation (Pulse ', string(i), ')'))
            legend('Measured Vt', 'Simulated Vt');

            %% Update the Model Parameter Values 

            % Update the model with the estimated values.
            sdo.setValueInModel(model,vOpt(1:end));

            % Add optimized parameters to table |param|
            if i ~= 1 && i ~= height(param)
                param.R0_opt(i+1) = vOpt(1).Value(1);
                param.R1_opt(i+1) = vOpt(2).Value(1);
                param.R2_opt(i+1) = vOpt(3).Value(1);
                param.R3_opt(i+1) = vOpt(4).Value(1);
                param.tau1_opt(i+1) = vOpt(5).Value(1);
                param.tau2_opt(i+1) = vOpt(6).Value(1);
                param.tau3_opt(i+1) = vOpt(7).Value(1);
                param.OCV_opt(i+1) = vOpt(8).Value(1);

            elseif i == 1 % first pulse
                param.R0_opt(i) = vOpt(1).Value(1);
                param.R1_opt(i) = vOpt(2).Value(1);
                param.R2_opt(i) = vOpt(3).Value(1);
                param.R3_opt(i) = vOpt(4).Value(1);
                param.tau1_opt(i) = vOpt(5).Value(1);
                param.tau2_opt(i) = vOpt(6).Value(1);
                param.tau3_opt(i) = vOpt(7).Value(1);

                param.R0_opt(i+1) = vOpt(8).Value(1);
                param.R1_opt(i+1) = vOpt(9).Value(1);
                param.R2_opt(i+1) = vOpt(10).Value(1);
                param.R3_opt(i+1) = vOpt(11).Value(1);
                param.tau1_opt(i+1) = vOpt(12).Value(1);
                param.tau2_opt(i+1) = vOpt(13).Value(1);
                param.tau3_opt(i+1) = vOpt(14).Value(1);

                param.OCV_opt(i) = vOpt(15).Value(1);
                param.OCV_opt(i+1) = vOpt(16).Value(1);

            elseif i == height(param) % last pulse
                param.R0_opt(i+1) = vOpt(1).Value(1);
                param.R1_opt(i+1) = vOpt(2).Value(1);
                param.R2_opt(i+1) = vOpt(3).Value(1);
                param.R3_opt(i+1) = vOpt(4).Value(1);
                param.tau1_opt(i+1) = vOpt(5).Value(1);
                param.tau2_opt(i+1) = vOpt(6).Value(1);
                param.tau3_opt(i+1) = vOpt(7).Value(1);

                param.R0_opt(i+2) = vOpt(8).Value(1);
                param.R1_opt(i+2) = vOpt(9).Value(1);
                param.R2_opt(i+2) = vOpt(10).Value(1);
                param.R3_opt(i+2) = vOpt(11).Value(1);
                param.tau1_opt(i+2) = vOpt(12).Value(1);
                param.tau2_opt(i+2) = vOpt(13).Value(1);
                param.tau3_opt(i+2) = vOpt(14).Value(1);

                param.OCV_opt(i+1) = vOpt(15).Value(1);
                param.OCV_opt(i+2) = vOpt(16).Value(1);
            end

            % Update initial states for next pulse
            if i <= height(param) -1
                sim_time = signal_I1.Values.Time;
                sim_SOC = signal_SOC.Values.Data;
                sim_I1 = signal_I1.Values.Data;
                sim_I2 = signal_I2.Values.Data;
                sim_I3 = signal_I3.Values.Data;

                init_I1 = sim_I1(end);
                init_I2 = sim_I2(end);
                init_I3 = sim_I3(end);
                init_SOC = sim_SOC(end);
            end

            % Save optimized parameters
            filename = string(folder_result) + "param_opt_" + dataSample + "_Poly" + polyDeg + "_" + string(i);
            save(filename, 'param', 'SOC', 'init_I1', 'init_I2', 'init_I3', 'init_SOC')

            % Clear unnecessary variables
            clear seg_optim SDOSimTest_Log signal_Vt 

        end
        
        close all
    end
end
