clear

%% Load data from the experiment

% Folders
folder_current = cd; 
folder_project = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design";
folder_functions = folder_project + "\Scripts\Parameterization\functions\"; 
folder_data = folder_project + "\Test Results"; % where experimental data is saved
folder_result = folder_project + "\Test Results\7 - Pulse Test\Results\"; % where results from this script is saved to
addpath(folder_current);
addpath(folder_functions);
addpath(genpath(folder_data));
addpath(folder_result);

% Load data
load('PROCESSED_PUL25dch.mat'); 
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
    'PostBufferSamples',10); %how many samples to include after the next pulse starts
psObj.plotIdentifiedPulses();

%%  Estimate Parameters
Params = psObj.Parameters;

R_init = 2e-4; 
R_lb = 1e-5; 
R_ub = 1e-3; 

% Set R0 constraints and initial guesses
Params.R0(:) = R_init;
Params.R0Min(:) = R_lb;
Params.R0Max(:) = R_ub;

% Estimate Em
Params.Em(1, :) = 3.2;
% Params.Em(1, 1) = 2.5;
% Params.Em(1, end) = 3.65;
Params.EmMin(1, :) = 2.4;
Params.EmMax(1, :) = 3.75;

% Set Tx constraints and initial guesses
Params.Tx(1,:)      = 1;
Params.TxMin(1,:)   = 0.1;
Params.TxMax(1,:)   = 4;

Params.Tx(2,:)      = 50;
Params.TxMin(2,:)   = 20;
Params.TxMax(2,:)   = 80;

Params.Tx(3, :)     = 1500; 
Params.TxMin(3, :)  = 1000;
Params.TxMax(3, :)  = 2000;

% Params.Tx(1,:)      = 2;
% Params.TxMin(1,:)   = 0.1;
% Params.TxMax(1,:)   = 10;
% 
% Params.Tx(2,:)      = 50;
% Params.TxMin(2,:)   = 10;
% Params.TxMax(2,:)   = 500;
% 
% Params.Tx(3, :)     = 1000; 
% Params.TxMin(3, :)  = 500;
% Params.TxMax(3, :)  = 2000;

% Set Rx constraints and initial guesses
Params.Rx(1,:)      = R_init;
Params.RxMin(1,:)   = R_lb;
Params.RxMax(1,:)   = R_ub;

Params.Rx(2,:)      = R_init;
Params.RxMin(2,:)   = R_lb;
Params.RxMax(2,:)   = R_ub;

Params.Rx(3, :)     = R_init; 
Params.RxMin(3, :)  = R_lb;
Params.RxMax(3, :)  = R_ub;

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



%% Export initial guess
psParam = psObj.Parameters;

% % Update OCV at SOC = 0
% psParam.Em(1) = 2.5;
% psParam.Em(end) = 3.65;
% 
% % Update R0 at SOC = 0 & SOC = 1
% psParam.R0(1) = 4.72e-4;
% psParam.R0(end) = 4.72e-4;

% Round OCV to 0.01 to avoid make sure it's non-decreasing
psParam.EmMin = round(psObj.Parameters.EmMin, 2);
psParam.EmMax = round(psObj.Parameters.EmMax, 2);
psParam.Em = round(psObj.Parameters.Em, 2);

param.SOC = psParam.SOC';

[param.R0_init, param.R1_init,      param.R2_init,      param.R3_init,      param.tau1_init,    param.tau2_init,    param.tau3_init,    param.OCV_init] = deal(...
 psParam.R0',   psParam.Rx(1, :)',  psParam.Rx(2, :)',  psParam.Rx(3, :)',  psParam.Tx(1, :)',  psParam.Tx(2, :)',  psParam.Tx(3, :)',  psParam.Em');

[param.R0_lb,       param.R1_lb,            param.R2_lb,            param.R3_lb, param.tau1_lb, param.tau2_lb, param.tau3_lb, param.OCV_lb] = deal(...
 psParam.R0Min',        psParam.RxMin(1, :)',   psParam.RxMin(2, :)',   psParam.RxMin(3, :)', ...
 psParam.TxMin(1, :)',  psParam.TxMin(2, :)',   psParam.TxMin(3, :)',   psParam.EmMin');

[param.R0_ub, param.R1_ub, param.R2_ub, param.R3_ub, param.tau1_ub, param.tau2_ub, param.tau3_ub, param.OCV_ub] = deal(...
 psParam.R0Max',        psParam.RxMax(1, :)',   psParam.RxMax(2, :)',   psParam.RxMax(3, :)', ...
 psParam.TxMax(1, :)',  psParam.TxMax(2, :)',   psParam.TxMax(3, :)',   psParam.EmMax');

param = struct2table(param);

%% Plot the parameters & their bounds
figure('WindowStyle', 'docked');

ax1 = subplot(2, 4, 1); hold on; 
plot(param.SOC, param.R0_init, '.-'); plot(param.SOC, param.R0_lb, '.-'); plot(param.SOC, param.R0_ub, '.-'); 
title('R0'); grid on; xlabel('SOC'); ylabel('R [Ohm]');

ax2 = subplot(2, 4, 2); hold on; 
plot(param.SOC, param.R1_init, '.-'); plot(param.SOC, param.R1_lb, '.-'); plot(param.SOC, param.R1_ub, '.-'); 
title('R1'); grid on; xlabel('SOC'); ylabel('R [Ohm]');

ax3 = subplot(2, 4, 3); hold on; 
plot(param.SOC, param.R2_init, '.-'); plot(param.SOC, param.R2_lb, '.-'); plot(param.SOC, param.R2_ub, '.-'); 
title('R2'); grid on; xlabel('SOC'); ylabel('R [Ohm]');

ax4 = subplot(2, 4, 4); hold on; 
plot(param.SOC, param.R3_init, '.-'); plot(param.SOC, param.R3_lb, '.-'); plot(param.SOC, param.R3_ub, '.-'); 
title('R3'); grid on; xlabel('SOC'); ylabel('R [Ohm]');

ax5 = subplot(2, 4, 5); hold on; 
plot(param.SOC, param.tau1_init, '.-'); plot(param.SOC, param.tau1_lb, '.-'); plot(param.SOC, param.tau1_ub, '.-'); 
title('tau1'); grid on; xlabel('SOC'); ylabel('tau [s]');

ax6 = subplot(2, 4, 6); hold on; 
plot(param.SOC, param.tau2_init, '.-'); plot(param.SOC, param.tau2_lb, '.-'); plot(param.SOC, param.tau2_ub, '.-'); 
title('tau2'); grid on; xlabel('SOC'); ylabel('tau [s]');

ax7 = subplot(2, 4, 7); hold on; 
plot(param.SOC, param.tau3_init, '.-'); plot(param.SOC, param.tau3_lb, '.-'); plot(param.SOC, param.tau3_ub, '.-');  
title('tau3'); grid on; xlabel('SOC'); ylabel('tau [s]');

ax8 = subplot(2, 4, 8); hold on; 
plot(param.SOC, param.OCV_init, '.-'); plot(param.SOC, param.OCV_lb, '.-'); plot(param.SOC, param.OCV_ub, '.-'); 
title('OCV'); grid on; xlabel('SOC'); ylabel('OCV [V]');

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8], 'x')

%% Make sure OCV is strictly increasing
try
    validateattributes(param.OCV_init, {'double'}, {'nondecreasing'})
catch
    errorInd1 = [];
    for i = 2:height(param)
        if param.OCV_init(i) <= param.OCV_init(i-1)
            errorInd1 = [errorInd1; i];
        end
    end
    error("Check if param.OCV_init is strictly increasing (see errorInd1)")
end

% Validate OCV ub is greater than OCV lb 
OCVdiff = param.OCV_ub - param.OCV_lb; 
errorInd2 = find(OCVdiff<=0, 1); 
if ~isempty(errorInd2)
    error('Check ub & lb of OCV (see errorInd2)');
end