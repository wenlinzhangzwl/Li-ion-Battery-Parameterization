clear

%% Load data from the experiment

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

% Load data
load('PROCESSED_PUL_0p1_discharge.mat'); 
validateattributes(meas_t.Time, {'double'}, {'increasing'});
validateattributes(meas_resampled.Time, {'double'}, {'increasing'});

% Create a Battery.PulseSequency object with all measurement data
psObj = Battery.PulseSequence;
psObj.ModelName = 'BatteryEstim3RC_PTBS';
addData(psObj, meas_resampled.Time, meas_resampled.Voltage, meas_resampled.Current); % The MATLAB functions take negative current as discharge
psObj.plot();

% Identify the pulses within the data set
psObj.createPulses(...
    'CurrentOnThreshold',0.025,... %minimum current magnitude to identify pulse events
    'NumRCBranches',3,... %how many RC pairs in the model
    'RCBranchesUse2TimeConstants',false,... %do RC pairs have different time constant for discharge and rest?
    'PreBufferSamples',10,... %how many samples to include before the current pulse starts
    'PostBufferSamples',50); %how many samples to include after the next pulse starts
psObj.plotIdentifiedPulses();

%%  Estimate Parameters
Params = psObj.Parameters;

% Set R0 constraints and initial guesses
Params.R0(:) = 2.5e-4;
Params.R0Min(:) = 1e-4;
Params.R0Max(:) = 1e-3;

% Estimate Em
Params.Em(1, :) = 3.2;
Params.Em(1, 1) = 2.5;
Params.Em(1, end) = 3.65;

Params.EmMin(1, :) = 2.4;
Params.EmMax(1, :) = 3.6;

% Set Tx constraints and initial guesses
Params.Tx(1,:)      = 2;
Params.TxMin(1,:)   = 1;
Params.TxMax(1,:)   = 3;

Params.Tx(2,:)      = 50;
Params.TxMin(2,:)   = 40;
Params.TxMax(2,:)   = 60;

Params.Tx(3, :)     = 1000; 
Params.TxMin(3, :)  = 100;
Params.TxMax(3, :)  = 3600*2; %don't set this bigger than the relaxation time available

% Set Rx constraints and initial guesses
Params.Rx(1,:)      = 2.5e-4;
Params.RxMin(1,:)   = 1e-4;
Params.RxMax(1,:)   = 5e-4;

Params.Rx(2,:)      = 3e-4;
Params.RxMin(2,:)   = 1e-5;
Params.RxMax(2,:)   = 1e-3;

Params.Rx(3, :)     = 5e-4; 
Params.RxMin(3, :)  = 1e-5;
Params.RxMax(3, :)  = 8e-4;

% [Params.Rx(1, 1:36),    Params.Rx(1, 37:37)] = deal(1e-4, 1.6e-4);
% [Params.RxMin(1, 1:36), Params.RxMin(1, 37:37)] = deal(1e-4, 1e-4);
% [Params.RxMax(1, 1:36), Params.RxMax(1, 37:37)] = deal(1.4e-4, 1e-3);
% 
% [Params.Rx(2, 1:11),    Params.Rx(2, 12:32),     Params.Rx(2, 33:37)   ] = deal(2e-4, 2e-4, 2e-4);
% [Params.RxMin(2, 1:11), Params.RxMin(2, 12:32),  Params.RxMin(2, 33:37)] = deal(1e-4, 1e-4, 1e-4);
% [Params.RxMax(2, 1:11), Params.RxMax(2, 12:32),  Params.RxMax(2, 33:37)] = deal(1e-3, 3e-4, 1e-3);
% 
% [Params.Rx(3, 1:11),    Params.Rx(3, 12:37)    ] = deal(5e-4, 5e-4);
% [Params.RxMin(3, 1:11), Params.RxMin(3, 12:37) ] = deal(1e-4, 1e-4);
% [Params.RxMax(3, 1:11), Params.RxMax(3, 12:37) ] = deal(0.01, 1e-3);

% Update parameters
psObj.Parameters = Params;
psObj.plotLatestParameters();

% Estimate initial R0 values
psObj.estimateInitialEmR0(...
    'SetEmConstraints',true,... %Update EmMin or EmMax values based on what we learn here
    'EstimateEm',true,... %Keep this on to perform Em estimates
    'EstimateR0',true); %Keep this on to perform R0 estimates

% Plot results
psObj.plotLatestParameters();

% Get initial Tx (Tau) values
psObj.estimateInitialTau(...
    'UpdateEndingEm',false,... %Keep this on to update Em estimates at the end of relaxations, based on the curve fit
    'ShowPlots',true,... %Set this true if you want to see plots while this runs
    'ReusePlotFigure',true,... %Set this true to overwrite the plots in the same figure
    'UseLoadData',false,... %Set this true if you want to estimate Time constants from the load part of the pulse, instead of relaxation
    'PlotDelay',0.5); %Set this to add delay so you can see the plots 

% % Plot results
psObj.plotLatestParameters(); %See what the parameters look like so far
% psObj.plotSimulationResults(); %See what the result looks like so far

% Get initial Em and Rx values using a linear system approach - pulse by pulse
psObj.estimateInitialEmRx(...
    'IgnoreRelaxation',false,... %Set this true if you want to ignore the relaxation periods during this step
    'ShowPlots',true,...  %Set this true if you want to see plots while this runs
    'ShowBeforePlots',true,... %Set this true if you want to see the 'before' value on the plots
    'PlotDelay',0.5,... %Set this to add delay so you can see the plots 
    'EstimateEm',true,... %Set this true to allow the optimizer to change Em further in this step
    'RetainEm',true,... %Set this true keep any changes made to Em in this step
    'EstimateR0',true,... %Set this true to allow the optimizer to change R0 further in this step
    'RetainR0',true); %Set this true keep any changes made to R0 in this step

Simulink.SimulationData.ModelLoggingInfo

% Plot results
psObj.plotLatestParameters(); %See what the parameters look like so far
psObj.plotSimulationResults(); %See what the result looks like so far

% Update OCV at SOC = 0
psObj.Parameters.Em(1) = 2.5;
psObj.Parameters.Em(end) = 3.65;

% Update R0 at SOC = 0 & SOC = 1
psObj.Parameters.R0(1) = 4.72e-4;
psObj.Parameters.R0(end) = 4.72e-4;

psParam = psObj.Parameters;

% Export initial guess
param.SOC = psParam.SOC';
[param.R0_init, param.R1_init,      param.R2_init,      param.R3_init,      param.tau1_init,    param.tau2_init,    param.tau3_init,    param.OCV_init] = deal(...
 psParam.R0',   psParam.Rx(1, :)',  psParam.Rx(2, :)',  psParam.Rx(3, :)',  psParam.Tx(1, :)',  psParam.Tx(2, :)',  psParam.Tx(3, :)',  psParam.Em');

[param.R0_lb,       param.R1_lb,            param.R2_lb,            param.R3_lb, param.tau1_lb, param.tau2_lb, param.tau3_lb, param.OCV_lb] = deal(...
 psParam.R0Min',        psParam.RxMin(1, :)',   psParam.RxMin(2, :)',   psParam.RxMin(3, :)', ...
 psParam.TxMin(1, :)',  psParam.TxMin(2, :)',   psParam.TxMin(3, :)',   psParam.Em'*0.995);

[param.R0_ub, param.R1_ub, param.R2_ub, param.R3_ub, param.tau1_ub, param.tau2_ub, param.tau3_ub, param.OCV_ub] = deal(...
 psParam.R0Max',        psParam.RxMax(1, :)',   psParam.RxMax(2, :)',   psParam.RxMax(3, :)', ...
 psParam.TxMax(1, :)',  psParam.TxMax(2, :)',   psParam.TxMax(3, :)',   psParam.Em'*1.005);

param = struct2table(param);

try
    validateattributes(param.OCV_init, {'double'}, {'increasing'})
catch
    errorInd = [];
    for i = 2:height(param)
        if param.OCV_init(i) <= param.OCV_init(i-1)
            errorInd = [errorInd; i];
        end
    end
    error("Check if param.OCV_init is strictly increasing (see errorInd)")
end