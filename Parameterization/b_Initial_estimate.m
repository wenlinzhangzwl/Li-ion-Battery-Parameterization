clear

%% Load data from the experiment

% Folders
folder_current = cd; 
folder_functions = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\Scripts\Parameterization\functions"; 
folder_data = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\10 - Characterization\"; % where experimental data is saved
folder_result = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\Results\"; % where results from this script is saved to
addpath(folder_current, folder_functions, folder_data, folder_result);

% Test settings
cell = "EVE280";
switch cell
    case "EVE280"
        Vmax = 3.65; 
        Vmin = 2.5;
end

% Load data
file = "PUL_25degC_0p8_1min";
filename = file +".mat"; 
load(filename, 'meas'); 
try
    meas = struct2table(meas);
end
meas.Current = -meas.Current; % Convert to convension: -ve current = charge, +ve current = discharge

% Add in an artifical data point for the function to recognize the first pulse
if meas.Current(1) ~= 0
    fake_data = meas(1, :);
    fake_data.Current = 0; 
    fake_data.Step = 99; 

    meas.Time = meas.Time + meas.Time(2);
    meas = [fake_data; meas];
end

% Delete data where time is not strictly increasing
time_diff = meas.Time(2:end) - meas.Time(1:end-1);
ind = find(time_diff <= 0) + 1; 
meas(ind, :) = []; 
validateattributes(meas.Time, {'double'}, {'increasing'});

% Create a Battery.PulseSequency object with all measurement data
psObj = Battery.PulseSequence;
psObj.ModelName = 'BatteryEstim3RC_PTBS';
addData(psObj, meas.Time, meas.Voltage, meas.Current); % The MATLAB functions take +ve current as discharge
psObj.plot();

% Identify the pulses within the data set
psObj.createPulses(...
    'CurrentOnThreshold',0.1,... %minimum current magnitude to identify pulse events
    'NumRCBranches',3,... %how many RC pairs in the model
    'RCBranchesUse2TimeConstants',false,... %do RC pairs have different time constant for discharge and rest?
    'PreBufferSamples',10,... %how many samples to include before the current pulse starts
    'PostBufferSamples',10); %how many samples to include after the next pulse starts
psObj.plotIdentifiedPulses();

%%  Estimate Parameters

Params = psObj.Parameters;

% Initialize parameters & ub/lb
param_initialization = parameter_initialization(cell, Params.SOC');
[Params.Em, Params.EmMin, Params.EmMax] = deal(param_initialization.OCV_init', param_initialization.OCV_lb', param_initialization.OCV_ub');
[Params.R0, Params.Rx(1, :), Params.Rx(2, :), Params.Rx(3, :)] = deal(param_initialization.R0_init', param_initialization.R1_init', ...
                                                                      param_initialization.R2_init', param_initialization.R3_init');
[Params.R0Min, Params.RxMin(1, :), Params.RxMin(2, :), Params.RxMin(3, :)] = deal(param_initialization.R0_lb');
[Params.R0Max, Params.RxMax(1, :), Params.RxMax(2, :), Params.RxMax(3, :)] = deal(param_initialization.R0_ub');
[Params.Tx(1, :), Params.Tx(2, :), Params.Tx(3, :)] = deal(param_initialization.tau1_init', param_initialization.tau2_init', param_initialization.tau3_init');
[Params.TxMin(1, :), Params.TxMin(2, :), Params.TxMin(3, :)] = deal(param_initialization.tau1_lb', param_initialization.tau2_lb', param_initialization.tau3_lb');
[Params.TxMax(1, :), Params.TxMax(2, :), Params.TxMax(3, :)] = deal(param_initialization.tau1_ub', param_initialization.tau2_ub', param_initialization.tau3_ub');

% Update parameters
psObj.Parameters = Params;
psObj.plotLatestParameters();

% Estimate initial R0 values
psObj.estimateInitialEmR0(...
    'SetEmConstraints',false,... %Update EmMin or EmMax values based on what we learn here
    'EstimateEm',true,... %Keep this on to perform Em estimates
    'EstimateR0',true); %Keep this on to perform R0 estimates

% Plot results
psObj.plotLatestParameters();

% Get initial Tx (Tau) values
psObj.estimateInitialTau(...
    'UpdateEndingEm',true,... %Keep this on to update Em estimates at the end of relaxations, based on the curve fit
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

% Round OCV to 0.01 to make sure it's non-decreasing
psParam.EmMin = round(psObj.Parameters.EmMin, 2);
psParam.EmMax = round(psObj.Parameters.EmMax, 2);
psParam.Em = round(psObj.Parameters.Em, 2);

param.SOC = psParam.SOC';

[param.R0_init, param.R1_init,      param.R2_init,      param.R3_init,      param.tau1_init,    param.tau2_init,    param.tau3_init,    param.OCV_init] = deal(...
 psParam.R0',   psParam.Rx(1, :)',  psParam.Rx(2, :)',  psParam.Rx(3, :)',  psParam.Tx(1, :)',  psParam.Tx(2, :)',  psParam.Tx(3, :)',  psParam.Em');

[param.R0_lb,           param.R1_lb,            param.R2_lb,            param.R3_lb, param.tau1_lb, param.tau2_lb, param.tau3_lb, param.OCV_lb] = deal(...
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


% Save parameters
filename = folder_result + "Parameters_Layered_" + file + ".mat"; 
save(filename, "param", "-mat")