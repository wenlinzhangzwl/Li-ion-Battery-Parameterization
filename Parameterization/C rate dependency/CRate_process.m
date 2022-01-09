clear
% clearvars -except folder_current folder_data folder_result folder_fcn param data Q

%% Get data & test conditions
folder_current = cd; 
folder_data = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\5 - C Rate Dependancy Tests\';            % where experimental data is saved
folder_result = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\5 - C Rate Dependancy Tests\Results\';  % where results from this script is saved
folder_fcn = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Model\Electric Model\Parameterization\functions\';      % where functions used in this script is saved
addpath(folder_current, folder_data, folder_result, folder_fcn);

% Load data
cd(folder_data);
highC = load('HighC25.mat');
midC = load('MidC25.mat');
lowC = load('LowC25.mat');
cd(folder_current);

% Load model parameters
load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Discharge\Results\param_opt')

%% Process data 
data = {highC, midC, lowC};

for i = 1:length(data)
    
    % Calculate SOC
    data{1, i} = struct2table(data{1, i}.meas);
    data{1, i}.Time = data{1, i}.Time - data{1, i}.Time(1);
    data{1, i}.SOC = 1 - (-data{1, i}.Ah)/Q;
    
    % Delete duplicated data points
    iDelete = [];
    for j = 1:length(data{1, i}.Time)-1
        if data{1, i}.Time(j) >= data{1, i}.Time(j+1)
            iDelete = [iDelete; j];
        end
    end
    data{1, i}(iDelete,:) = [];
    validateattributes(data{1, i}.Time, {'double'}, {'increasing'});
    
    % Divide into segments by SOC
    segNum = find(param.SOC <= min(data{1, i}.SOC(end)), 1);
    for j = 1:segNum-1
        iStart = find(data{1, i}.SOC <= param.SOC(j), 1);
        iEnd = find(data{1, i}.SOC <= param.SOC(j+1), 1);
        
        if isempty(iEnd)
            iEnd = height(data{1, i});
        end
        
        segment{j, 1} = [param.SOC(j), param.SOC(j+1)];
        segment{j, 2} = [iStart, iEnd];
        segment{j, 3} = data{1, i}(iStart:iEnd, :);
    end
    data{2, i} = segment; 
end

%% Optimize parameters for SOC 90-85
optimizedParam = [];

for profileNum = 1:length(data)
    
    % segment to be optimziaed & prior segments
    segment_opt = data{2, profileNum}{end, 3};
    segment_pri = [];
    for i = 1:height(data{2, profileNum}) - 1
        segment_pri = [segment_pri; data{2, profileNum}{i, 3}];
    end
    
    %% Set model parameters
    
    ind = find(param.SOC <= segment_opt.SOC(1), 1); 
    param_prev = flip(param(1:ind-1, :), 1);
    param_current = param(ind, :);

    % Parameter look up tables
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
    OCV_current  = param_current.OCV_opt;
    R0_current   = param_current.R0_opt;
    R1_current   = param_current.R1_opt;
    R2_current   = param_current.R2_opt;
    R3_current   = param_current.R3_opt;
    tau1_current = param_current.tau1_opt;
    tau2_current = param_current.tau2_opt;
    tau3_current = param_current.tau3_opt;
    
    % Initial states
    init_I1 = 0; 
    init_I2 = 0;
    init_I3 = 0; 
    init_SOC = segment_pri.SOC(1);

    %% Get initial conditions for optimization
    
    % Input data
    Time = segment_pri.Time;
    t = [Time Time];                        %#ok<NASGU>
    meas_I = [Time segment_pri.Current];    %#ok<NASGU>
    meas_Vt = [Time segment_pri.Voltage];   %#ok<NASGU>

    % Simulation settings
    stepSize = 0.01;
    runTime = Time(end) - Time(1);

    % Simulate the model
    open('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Model\Electric Model\C Rate Dependency\ElectricalModel_estimation_CRate')
    model = 'ElectricalModel_estimation_CRate';
    simout_init = sim(model);

    % Log the end of each parameter as initial conditions for optimization
    init_I1 = simout_init.I1(end);
    init_I2 = simout_init.I2(end);
    init_I3 = simout_init.I3(end);
    init_SOC = simout_init.SOC(end);
    
    % Update input data
    Time = segment_opt.Time - segment_opt.Time(1);
    t = [Time Time];
    meas_I = [Time segment_opt.Current];
    meas_Vt = [Time segment_opt.Voltage];
    
    %% Define the estimation exp & simulate with initial guess

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

    % Add the initial state from the model & set them to non-tunable
    Exp.InitialStates = sdo.getStateFromModel(model);
    for j = 1:5
        Exp.InitialStates(j).Free = 0;
    end

    % Create a simulation using the experiment
    Simulator = createSimulator(Exp);
    Simulator = sim(Simulator);

    % Search for the Vt_sim signal in the logged simulation data.
    SimLog = find(Simulator.LoggedData,get_param(model,'SignalLoggingName'));
    signal_Vt_sim = find(SimLog,'Vt_sim');
    signal_Vt_exp = find(SimLog,'Vt_exp');

    rmse_Vt = sqrt(mean((signal_Vt_sim.Values.Data - signal_Vt_exp.Values.Data).^2));

    % Plot the measured and simulated data.
    figName = 'Profile' + string(profileNum) + '(before)';
    figure('Name', figName, 'WindowStyle', 'docked'); hold on; grid on
    plot(meas_Vt(:, 1), meas_Vt(:, 2), '.-')
    plot(signal_Vt_sim.Values.Time,signal_Vt_sim.Values.Data, '.--');
    title(append('Simulated and Measured Responses Before Estimation (Pulse ', string(i), ')'))
    legend('Measured Vt', 'Simulated Vt');
    annotation('textbox', [0.14, 0.82, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt))

    %% Specify parameters to estimate & their bounds

    % Get parameters to estimate ('p') & set upper/lower bounds
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
    p(8).Maximum = param_prev.OCV_opt(1);

    % Group the model parameters and initial states to be estimated together.
    % Only has parameters in this case, no initial states to be estimated
    v = p;

    %% Define the Estimation Objective Function & Estimate the Parameters

    estFcn = @(v) ElectricalModel_estimation_Objective(v,Simulator,Exp, model);

    optOptions = sdo.OptimizeOptions(...
        'OptimizedModel',Simulator,...
        'Method','lsqnonlin',...
        'UseParallel',true);
    optOptions.MethodOptions.FunctionTolerance = 0.05;
    optOptions.MethodOptions.OptimalityTolerance = 0.01;
    optOptions.MethodOptions.StepTolerance = 0.1;
    optOptions.MethodOptions.MaxIterations = 10;

    % Estimate the parameters.
    vOpt = sdo.optimize(estFcn, v, optOptions);

    %% Compare the Measured Output and the Final Simulated Output

    % Update the experiments with the estimated parameter values.
    Exp = setEstimatedValues(Exp,vOpt);

    % Simulate the model using the updated experiment
    Simulator   = createSimulator(Exp,Simulator);
    Simulator   = sim(Simulator);
    SimLog      = find(Simulator.LoggedData,get_param(model,'SignalLoggingName'));
    signal_Vt_sim = find(SimLog,'Vt_sim');
    signal_Vt_exp = find(SimLog,'Vt_exp');
    signal_SOC  = find(SimLog,'SOC_sim');

    % Error
    rmse_Vt = sqrt(mean((signal_Vt_sim.Values.Data - signal_Vt_exp.Values.Data).^2));

    % Plot the measured and simulated data.
    figName = 'Pulse ' + string(i) + '(after)';
    figure('Name', figName, 'WindowStyle', 'docked');
    plot(meas_Vt(:, 1), meas_Vt(:, 2), '.-', ...
        signal_Vt_sim.Values.Time,signal_Vt_sim.Values.Data,'.--');
    title(append('Simulated and Measured Responses After Estimation (Pulse ', string(i), ')'))
    legend('Measured Vt', 'Simulated Vt');
    annotation('textbox', [0.14, 0.82, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt))

    %% Update the Model Parameter Values 

    % Update the model with the estimated values.
    sdo.setValueInModel(model,vOpt(1:end));

    % Add optimized parameters to table |param|
    optimizedParam.R0_opt(profileNum, 1) = vOpt(1).Value(1);
    optimizedParam.R1_opt(profileNum, 1) = vOpt(2).Value(1);
    optimizedParam.R2_opt(profileNum, 1) = vOpt(3).Value(1);
    optimizedParam.R3_opt(profileNum, 1) = vOpt(4).Value(1);
    optimizedParam.tau1_opt(profileNum, 1) = vOpt(5).Value(1);
    optimizedParam.tau2_opt(profileNum, 1) = vOpt(6).Value(1);
    optimizedParam.tau3_opt(profileNum, 1) = vOpt(7).Value(1);
    optimizedParam.OCV_opt(profileNum, 1) = vOpt(8).Value(1);
  
end

optimizedParam = struct2table(optimizedParam);

% Visualize optimized parameters
figure;
subplot(1, 8, 1); bar(optimizedParam.R0_opt); title('R0'); ylabel('ohm')
subplot(1, 8, 2); bar(optimizedParam.R1_opt); title('R1'); ylabel('ohm')
subplot(1, 8, 3); bar(optimizedParam.R2_opt); title('R2'); ylabel('ohm')
subplot(1, 8, 4); bar(optimizedParam.R3_opt); title('R3'); ylabel('ohm')
subplot(1, 8, 5); bar(optimizedParam.tau1_opt); title('tau1'); ylabel('s')
subplot(1, 8, 6); bar(optimizedParam.tau2_opt); title('tau2'); ylabel('s')
subplot(1, 8, 7); bar(optimizedParam.tau3_opt); title('tau3'); ylabel('s')
subplot(1, 8, 8); bar(optimizedParam.OCV_opt); title('OCV'); ylabel('V')