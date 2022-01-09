clear
% % clearvars -except pulse
% 
% %% Load data from the experiment
% 
% % Load files
% folder_current = cd; 
% folder_data   = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Test\';                % where experimental data is saved
% folder_result = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\';        % where results from this script is saved
% folder_fcn    = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Model\Electric Model\Parameterization\functions\'; % where functions used in this script is saved
% addpath(folder_current);
% addpath(folder_data);
% addpath(folder_result);
% addpath(folder_fcn);
% 
% cd(folder_data);
% % load('PUL25dch_meas_t.mat'); % Trimmed data
% load('PUL25dch_pulse.mat');  % Info of each pulse
% cd(folder_current);
% 
% % Measured data (previously downsampled)
% meas_resampled = [];
% for i = 1:height(pulse)
%     meas_resampled = [meas_resampled; pulse.segment2{i}];
% end
% meas_resampled = unique(meas_resampled, 'rows');
% 
% meas_Vt = meas_resampled.Voltage; 
% meas_I = -meas_resampled.Current; 
% meas_Time = meas_resampled.Time;
% validateattributes(meas_Time, {'double'}, {'increasing'});
% 
% % Create a Battery.PulseSequency object with all measurement data
% psObj = Battery.PulseSequence;
% psObj.ModelName = 'BatteryEstim3RC_PTBS';
% addData(psObj, meas_Time, meas_Vt, meas_I);
% psObj.plot();
% 
% % Identify the pulses within the data set
% psObj.createPulses(...
%     'CurrentOnThreshold',0.025,... %minimum current magnitude to identify pulse events
%     'NumRCBranches',3,... %how many RC pairs in the model
%     'RCBranchesUse2TimeConstants',false,... %do RC pairs have different time constant for discharge and rest?
%     'PreBufferSamples',10,... %how many samples to include before the current pulse starts
%     'PostBufferSamples',50); %how many samples to include after the next pulse starts
% psObj.plotIdentifiedPulses();
% 
% %%  Estimate Parameters
% Params = psObj.Parameters;
% 
% % Set R0 constraints and initial guesses
% Params.R0(:) = 2.5e-4;
% Params.R0Min(:) = 1e-4;
% Params.R0Max(:) = 1e-3;
% 
% % Estimate Em
% Params.Em(1, :) = 3.2;
% Params.Em(1, 1) = 2.5;
% Params.Em(1, end) = 3.65;
% 
% Params.EmMin(1, :) = 2.4;
% Params.EmMax(1, :) = 3.6;
% 
% % Set Tx constraints and initial guesses
% % [Params.Tx(1, 1:26),    Params.Tx(1, 27:37)]    = deal(2,   2); 
% % [Params.TxMin(1, 1:26), Params.TxMin(1, 27:37)] = deal(0.1, 0.1);
% % [Params.TxMax(1, 1:26), Params.TxMax(1, 27:37)] = deal(10,   10);
% Params.Tx(1,:)      = 2;
% Params.TxMin(1,:)   = 0.1;
% Params.TxMax(1,:)   = 10;
% 
% Params.Tx(2,:)      = 50;
% Params.TxMin(2,:)   = 30;
% Params.TxMax(2,:)   = 100;
% 
% Params.Tx(3, :)     = 1500; 
% Params.TxMin(3, :)  = 500;
% Params.TxMax(3, :)  = 3600*2; %don't set this bigger than the relaxation time available
% 
% % Set Rx constraints and initial guesses
% Params.Rx(1,:)      = 5e-4;
% Params.RxMin(1,:)   = 5e-5;
% Params.RxMax(1,:)   = 1e-3;
% 
% Params.Rx(2,:)      = 5e-4;
% Params.RxMin(2,:)   = 5e-5;
% Params.RxMax(2,:)   = 1e-3;
% 
% Params.Rx(3, :)     = 5e-4; 
% Params.RxMin(3, :)  = 5e-5;
% Params.RxMax(3, :)  = 1e-3;
% 
% % [Params.Rx(1, 1:36),    Params.Rx(1, 37:37)] = deal(1e-4, 1.6e-4);
% % [Params.RxMin(1, 1:36), Params.RxMin(1, 37:37)] = deal(1e-4, 1e-4);
% % [Params.RxMax(1, 1:36), Params.RxMax(1, 37:37)] = deal(1.4e-4, 1e-3);
% % 
% % [Params.Rx(2, 1:11),    Params.Rx(2, 12:32),     Params.Rx(2, 33:37)   ] = deal(2e-4, 2e-4, 2e-4);
% % [Params.RxMin(2, 1:11), Params.RxMin(2, 12:32),  Params.RxMin(2, 33:37)] = deal(1e-4, 1e-4, 1e-4);
% % [Params.RxMax(2, 1:11), Params.RxMax(2, 12:32),  Params.RxMax(2, 33:37)] = deal(1e-3, 3e-4, 1e-3);
% % 
% % [Params.Rx(3, 1:11),    Params.Rx(3, 12:37)    ] = deal(5e-4, 5e-4);
% % [Params.RxMin(3, 1:11), Params.RxMin(3, 12:37) ] = deal(1e-4, 1e-4);
% % [Params.RxMax(3, 1:11), Params.RxMax(3, 12:37) ] = deal(0.01, 1e-3);
% 
% % Update parameters
% psObj.Parameters = Params;
% psObj.plotLatestParameters();
% 
% % Estimate initial R0 values
% psObj.estimateInitialEmR0(...
%     'SetEmConstraints',true,... %Update EmMin or EmMax values based on what we learn here
%     'EstimateEm',true,... %Keep this on to perform Em estimates
%     'EstimateR0',true); %Keep this on to perform R0 estimates
% 
% % Plot results
% psObj.plotLatestParameters();
% 
% % Get initial Tx (Tau) values
% psObj.estimateInitialTau(...
%     'UpdateEndingEm',false,... %Keep this on to update Em estimates at the end of relaxations, based on the curve fit
%     'ShowPlots',true,... %Set this true if you want to see plots while this runs
%     'ReusePlotFigure',true,... %Set this true to overwrite the plots in the same figure
%     'UseLoadData',false,... %Set this true if you want to estimate Time constants from the load part of the pulse, instead of relaxation
%     'PlotDelay',0.5); %Set this to add delay so you can see the plots 
% 
% % % Plot results
% psObj.plotLatestParameters(); %See what the parameters look like so far
% % psObj.plotSimulationResults(); %See what the result looks like so far
% 
% % Get initial Em and Rx values using a linear system approach - pulse by pulse
% psObj.estimateInitialEmRx(...
%     'IgnoreRelaxation',false,... %Set this true if you want to ignore the relaxation periods during this step
%     'ShowPlots',true,...  %Set this true if you want to see plots while this runs
%     'ShowBeforePlots',true,... %Set this true if you want to see the 'before' value on the plots
%     'PlotDelay',0.5,... %Set this to add delay so you can see the plots 
%     'EstimateEm',true,... %Set this true to allow the optimizer to change Em further in this step
%     'RetainEm',true,... %Set this true keep any changes made to Em in this step
%     'EstimateR0',true,... %Set this true to allow the optimizer to change R0 further in this step
%     'RetainR0',true); %Set this true keep any changes made to R0 in this step
% 
% Simulink.SimulationData.ModelLoggingInfo
% 
% % Plot results
% psObj.plotLatestParameters(); %See what the parameters look like so far
% psObj.plotSimulationResults(); %See what the result looks like so far
% 
% % Update OCV at SOC = 0
% psObj.Parameters.Em(1) = 2.5;
% psObj.Parameters.Em(end) = 3.65;
% 
% % Update R0 at SOC = 0 & SOC = 1
% psObj.Parameters.R0(1) = 4.72e-4;
% psObj.Parameters.R0(end) = 4.72e-4;
% 
% psParam = psObj.Parameters;
% 
% % Export initial guess
% param.SOC = psParam.SOC';
% [param.R0_init, param.R1_init,      param.R2_init,      param.R3_init,      param.tau1_init,    param.tau2_init,    param.tau3_init,    param.OCV_init] = deal(...
%  psParam.R0',   psParam.Rx(1, :)',  psParam.Rx(2, :)',  psParam.Rx(3, :)',  psParam.Tx(1, :)',  psParam.Tx(2, :)',  psParam.Tx(3, :)',  psParam.Em');
% 
% [param.R0_lb,       param.R1_lb,            param.R2_lb,            param.R3_lb, param.tau1_lb, param.tau2_lb, param.tau3_lb, param.OCV_lb] = deal(...
%  psParam.R0Min',        psParam.RxMin(1, :)',   psParam.RxMin(2, :)',   psParam.RxMin(3, :)', ...
%  psParam.TxMin(1, :)',  psParam.TxMin(2, :)',   psParam.TxMin(3, :)',   psParam.Em'*0.995);
% 
% [param.R0_ub, param.R1_ub, param.R2_ub, param.R3_ub, param.tau1_ub, param.tau2_ub, param.tau3_ub, param.OCV_ub] = deal(...
%  psParam.R0Max',        psParam.RxMax(1, :)',   psParam.RxMax(2, :)',   psParam.RxMax(3, :)', ...
%  psParam.TxMax(1, :)',  psParam.TxMax(2, :)',   psParam.TxMax(3, :)',   psParam.Em'*1.005);
% 
% param = struct2table(param);





clear
load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_init.mat')

% Set upper & lower bounds
bound = 'constrained';
switch bound
    case 'unconstrained'
        param.OCV_ub =  param.OCV_init*1.005;
        param.OCV_lb =  param.OCV_init*0.998;

        param.R0_ub =  param.R0_init*2;
        param.R0_lb =  param.R0_init*0.5;
        param.R1_ub =  param.R1_init*2;
        param.R1_lb =  param.R1_init*0.5;
        param.R2_ub =  param.R2_init*2;
        param.R2_lb =  param.R2_init*0.5;
        param.R3_ub =  param.R3_init*2;
        param.R3_lb =  param.R3_init*0.5;

        param.tau1_ub =  param.tau1_init*2;
        param.tau1_lb =  param.tau1_init*0.5;
        param.tau2_ub =  param.tau2_init*2;
        param.tau2_lb =  param.tau2_init*0.5;
        param.tau3_ub =  param.tau3_init*2;
        param.tau3_lb =  param.tau3_init*0.5;

        param.R0_ub(36) =  param.R0_init(36)*1.1;
        param.R0_lb(36) =  param.R0_init(36)*0.9;
        param.R1_ub(36) =  param.R1_init(36)*1.1;
        param.R1_lb(36) =  param.R1_init(36)*0.9;
        param.R2_ub(36) =  param.R2_init(36)*1.1;
        param.R2_lb(36) =  param.R2_init(36)*0.9;
        param.R3_ub(36) =  param.R3_init(36)*1.1;
        param.R3_lb(36) =  param.R3_init(36)*0.9;

        param.tau1_ub(36) =  param.tau1_init(36)*1.1;
        param.tau1_lb(36) =  param.tau1_init(36)*0.9;
        param.tau2_ub(36) =  param.tau2_init(36)*1.1;
        param.tau2_lb(36) =  param.tau2_init(36)*0.9;
        param.tau3_ub(36) =  param.tau3_init(36)*1.1;
        param.tau3_lb(36) =  param.tau3_init(36)*0.9;
    case 'constrained'
        param.OCV_ub =  param.OCV_init*1.005;
        param.OCV_lb =  param.OCV_init*0.998;

        param.R0_ub =  param.R0_init*2;
        param.R0_lb =  param.R0_init*0.5;
        param.R1_ub =  param.R1_init*2;
        param.R1_lb =  param.R1_init*0.5;
        param.R2_ub =  param.R2_init*2;
        param.R2_lb =  param.R2_init*0.5;
        param.R3_ub =  param.R3_init*2;
        param.R3_lb =  param.R3_init*0.5;

        param.tau1_init(:) =  2;
        param.tau1_ub(:) =  3;
        param.tau1_lb(:) =  1;
        
        param.tau2_init(:) =  50;
        param.tau2_ub(:) =  60;
        param.tau2_lb(:) =  40;
        
        param.tau3_ub =  param.tau3_init*2;
        param.tau3_lb =  param.tau3_init*0.5;
end

load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Test\PUL25chg_R0.mat')

% Add R0_chg to param & flip parameters
param.R0_init_chg = abs(param_chg(:, 2));
param.R0_ub_chg = param.R0_init_chg*2;
param.R0_lb_chg = param.R0_init_chg*0.5;
param.R0_ub_chg(36) = param.R0_init_chg(36)*1.1;
param.R0_lb_chg(36) = param.R0_init_chg(36)*0.9;

folder_current = cd; 
folder_data   = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Test\';                % where experimental data is saved
folder_result = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\';        % where results from this script is saved
folder_fcn    = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Model\Electric Model\Parameterization\functions\'; % where functions used in this script is saved
folder_OCV    = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\'; % where OCV polynomial coefficients are saved
addpath(folder_current);
addpath(folder_data);
addpath(folder_result);
addpath(folder_fcn);
addpath(folder_OCV);

load('C20CCCV_OCVSOC10_new.mat');  % OCV coefficients (8, 9, 10, 11, 15, 18, 19, 21)
clear coeff_chg coeff_dch OCV_SOC

%% Set up input data & parameters

% Data sample
meas_US06 = load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\8 - Capacity test & US06\US0625.mat');
% meas_UDDS = load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\6 - Drive cycle\UDDS25.mat');
% meas_MIX = load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\6 - Drive cycle\MIX25.mat');
% meas_NEDC = load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\6 - Drive cycle\NEDC23.mat');
% meas_pulse = load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Test\PUL25dch_pulse.mat');

dataType = 'driveCycle';

switch dataType
    case 'driveCycle'
        meas_optim = struct2table(meas_US06.meas);

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
            Q = abs(meas_optim.Ah(end)-meas_optim.Ah(1));
        end
        meas_optim.Current = -meas_optim.Current;
        meas_optim.Ah = meas_optim.Ah - meas_optim.Ah(1);
        meas_optim.SOC = 1 - (-meas_optim.Ah)/Q;
        meas_optim.Time = meas_optim.Time - meas_optim.Time(1);

        % Divide data into segments by SOC
        SOC = param.SOC;
        SOC_flip = flip(SOC);
        for i = 1:height(SOC_flip)-1
            iStart = find(meas_optim.SOC <= SOC_flip(i), 1);
            iEnd = find(meas_optim.SOC <= SOC_flip(i+1), 1);
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
        segments(:, 3) = [meas_pulse.pulse.segment2]; 
        for i = 1:height(segments)
            segments{i, 3}.Current = -segments{i, 3}.Current;
        end
        SOC = param.SOC;
end

clear meas_optim meas_UDDS meas_US06 meas_MIX param_chg meas_pulse meas_NEDC

%% Optimize with drive cycle

% clear param SOC
% load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_6')
% param = flip(param);

% Set optimization model
model = 'ElectricalModel_estimation_DriveCycle';
useOCVpolynomial = 0; 
open(model)

param = flip(param);
for i = 1:height(segments)
    % Relevant parameters
    param_prev = flip(param(1:i, :));
    param_current = param(i+1, :);
    
    %% Set Model Parameters
    % Parameter look up tables with initial guesses
    if i ~= 1
        SOC_prev    = param_prev.SOC;
        OCV_prev    = param_prev.OCV_opt;
        R0_prev     = [param_prev.R0_opt param_prev.R0_opt_chg];
        R1_prev     = param_prev.R1_opt;
        R2_prev     = param_prev.R2_opt;
        R3_prev     = param_prev.R3_opt;
        tau1_prev   = param_prev.tau1_opt;
        tau2_prev   = param_prev.tau2_opt;
        tau3_prev   = param_prev.tau3_opt;

        SOC_current  = param_current.SOC;
        OCV_current  = param_current.OCV_init;
        R0_current   = [param_current.R0_init param_current.R0_init_chg];
        R1_current   = param_current.R1_init;
        R2_current   = param_current.R2_init;
        R3_current   = param_current.R3_init;
        tau1_current = param_current.tau1_init;
        tau2_current = param_current.tau2_init;
        tau3_current = param_current.tau3_init;
    else
        SOC_prev    = param_prev.SOC;
        OCV_prev    = param_prev.OCV_init;
        R0_prev     = [param_prev.R0_init param_prev.R0_init_chg];
        R1_prev     = param_prev.R1_init;
        R2_prev     = param_prev.R2_init;
        R3_prev     = param_prev.R3_init;
        tau1_prev   = param_prev.tau1_init;
        tau2_prev   = param_prev.tau2_init;
        tau3_prev   = param_prev.tau3_init;

        SOC_current  = param_current.SOC;
        OCV_current  = param_current.OCV_init;
        R0_current   = [param_current.R0_init param_current.R0_init_chg];
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
    sim_Vt.BlockPath = 'ElectricalModel_estimation_DriveCycle/StateUpdate';
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
        p(1).Minimum = [param_current.R0_lb, param_current.R0_lb_chg];        
        p(1).Maximum = [param_current.R0_ub, param_current.R0_ub_chg];
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
            break
        end
        
    elseif i == 1
        
        p = sdo.getParameterFromModel(model,{'R0_prev','R1_prev','R2_prev','R3_prev','tau1_prev','tau2_prev','tau3_prev',...
                                                 'R0_current','R1_current','R2_current','R3_current','tau1_current','tau2_current','tau3_current',...
                                                 'OCV_prev','OCV_current',});      
        % R0_prev
        p(1).Minimum = [param_prev.R0_lb, param_prev.R0_lb_chg];        
        p(1).Maximum = [param_prev.R0_ub, param_prev.R0_ub_chg];
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
        p(8).Minimum = [param_current.R0_lb, param_current.R0_lb_chg];        
        p(8).Maximum = [param_current.R0_ub, param_current.R0_ub_chg];
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
        param.R0_opt_chg(i+1) = vOpt(1).Value(2);
        param.R1_opt(i+1) = vOpt(2).Value(1);
        param.R2_opt(i+1) = vOpt(3).Value(1);
        param.R3_opt(i+1) = vOpt(4).Value(1);
        param.tau1_opt(i+1) = vOpt(5).Value(1);
        param.tau2_opt(i+1) = vOpt(6).Value(1);
        param.tau3_opt(i+1) = vOpt(7).Value(1);
        param.OCV_opt(i+1) = vOpt(8).Value(1);
        
    elseif i == 1 % first pulse
        param.R0_opt(i) = vOpt(1).Value(1);
        param.R0_opt_chg(i) = vOpt(1).Value(2);
        param.R1_opt(i) = vOpt(2).Value(1);
        param.R2_opt(i) = vOpt(3).Value(1);
        param.R3_opt(i) = vOpt(4).Value(1);
        param.tau1_opt(i) = vOpt(5).Value(1);
        param.tau2_opt(i) = vOpt(6).Value(1);
        param.tau3_opt(i) = vOpt(7).Value(1);
                
        param.R0_opt(i+1) = vOpt(8).Value(1);
        param.R0_opt_chg(i+1) = vOpt(8).Value(2);
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
        param.R0_opt_chg(i+1) = vOpt(1).Value(2);
        param.R1_opt(i+1) = vOpt(2).Value(1);
        param.R2_opt(i+1) = vOpt(3).Value(1);
        param.R3_opt(i+1) = vOpt(4).Value(1);
        param.tau1_opt(i+1) = vOpt(5).Value(1);
        param.tau2_opt(i+1) = vOpt(6).Value(1);
        param.tau3_opt(i+1) = vOpt(7).Value(1);

        param.R0_opt(i+2) = vOpt(8).Value(1);
        param.R0_opt_chg(i+2) = vOpt(8).Value(2);
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
    save(append(folder_result, 'param_opt_', string(i)), 'param', 'SOC', 'init_I1', 'init_I2', 'init_I3', 'init_SOC')
    
    % Clear unnecessary variables
    clear seg_optim SDOSimTest_Log signal_Vt 
    
end
