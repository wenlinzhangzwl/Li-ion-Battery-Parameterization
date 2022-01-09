clear

%% Get data & test conditions

% Folders
folder_current = cd; 
folder_project = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design';
folder_functions = append(folder_project, '\Scripts\Parameterization\functions\'); 
folder_testResults = append(folder_project, '\Test Results');
folder_result = append(folder_project, '\Test Results\7 - Pulse Test\Results\');  % where results from this script is saved
addpath(folder_current);
addpath(folder_functions);
addpath(genpath(folder_testResults));
addpath(folder_result);

% Load data
capacityTest = load('CAP25.mat');
Q = -capacityTest.meas.Ah(end);
load('PUL_0p1_discharge.mat'); % Pulse discharge test at 0.8C

%% Divide data into pulses

% 0.8C Pulse Test
% % Step numbers *******************************************CHANGE
% iRelax0 =   16;     % Rest before the start of the test
% iLoad1 =    18;     % 10 charge/discharge pulses from 100% to 90% SOC
% iRelax1 =   19;     % rest in between
% iLoad2 =    22;     % 16 charge/discharge pulses from 90% to 10% SOC
% iRelax2 =   23;     % rest in between
% iLoad3 =    26;     % 10 charge/discharge pulses from 10% to 0% SOC
% iRelax3 =   27;     % rest in between
% iRelax4 =   29;     % rest after test reached Vmin
% steps = [iRelax0 iLoad1 iRelax1 iLoad2 iRelax2 iLoad3 iRelax3 iRelax4];

% 0.1C Pulse Test
% Step numbers *******************************************CHANGE
iRelax0 =   14;     % Rest before the start of the test
iLoad1 =    16;     % 10 charge/discharge pulses from 100% to 90% SOC
iRelax1 =   17;     % rest in between
iLoad2 =    20;     % 16 charge/discharge pulses from 90% to 10% SOC
iRelax2 =   21;     % rest in between
iLoad3 =    24;     % 10 charge/discharge pulses from 10% to 0% SOC
iRelax3 =   25;     % rest in between
iRelax4 =   27;     % rest after test reached Vmin
steps = [iRelax0 iLoad1 iRelax1 iLoad2 iRelax2 iLoad3 iRelax3 iRelax4];

% Correct current sign, negative = charge, positive = charge
meas.Current = -meas.Current; 

% Delete unnecessary fields
meas_t = struct2table(meas);
meas_t.TimeStamp = [];
meas_t.StepTime = [];
meas_t.Procedure = [];
meas_t.Wh = [];
meas_t.Power = [];
meas_t.Battery_Temp_degC = [];

% Trim unnecessary data at beginning & end
iBegin = find(meas.Step == min(steps), 1, 'first'); 
iEnd = find(meas.Step == max(steps), 1, 'last'); 
meas_t = meas_t(iBegin:iEnd, :);
meas_t = meas_t(meas_t.Step <= max(steps), :);
meas_t.Time = meas_t.Time - meas_t.Time(1);
meas_t.Ah = meas_t.Ah - meas_t.Ah(1);

% Remove data where time is not strictly increasing
[~,ia, ~] = unique(meas_t.Time);
meas_t = meas_t(ia,:);
validateattributes(meas_t.Time, {'double'}, {'increasing'})

% Calculate SOC
dch = 1; % Whether it is charging or discharging *******************************************CHANGE
if dch == 1
    meas_t.SOC = 1 - (-meas_t.Ah)/Q;
else
    meas_t.SOC = 0 + (meas_t.Ah)/Q;
end

% Divide data into pulses
ind = getInd(meas_t, steps);
pulse = getPulse(meas_t, ind);

% Plot segment1
figure
hold on; 
plot(meas_t.Time/3600, meas_t.Voltage, '.');
for i = 1:height(pulse)
    t = pulse.segment1{i,1}.Time/3600;
    Voltage = pulse.segment1{i,1}.Voltage;
    plot(t, Voltage, '*');
end
title('segment1 (for time constant curve fitting)'); xlabel('Time [h]'); ylabel('Voltage'); grid on
legend('data', 'sample')
hold off

% Plot segment2
figure
hold on; 
plot(meas_t.Time, meas_t.Voltage, '.');
for i = 1:height(pulse)
    t = pulse.segment2{i,1}.Time;
    Voltage = pulse.segment2{i,1}.Voltage;
    plot(t, Voltage, 'o');
end
title('segment2 (for optimization)'); xlabel('Time [s]'); ylabel('Voltage'); grid on
legend('data', 'sample')
hold off

% Combine downsampled data for optimization
meas_resampled = [];
for i = 1:height(pulse)
    meas_resampled = [meas_resampled; pulse.segment2{i}];
end
meas_resampled = unique(meas_resampled, 'rows');

%% Find initial conditions for each pulse
% 
% % Find initial conditions
% cd(folder_fcn);
% param = getInitialcondition(meas_t, pulse);
% cd(folder_current);
% % for i = 1:height(pulse)
% %     pulse.param(i) = {param(i, :)};
% % end
% 
% % Plot parameters for validation
% SOC = [1; pulse.SOC; 0];
% 
% figure;
% ax1 = subplot(4, 1, 1); plot(SOC, param.R0_init); grid on; title('R0')
% ax2 = subplot(4, 1, 2); plot(SOC, param.R1_init); grid on; title('R1')
% ax3 = subplot(4, 1, 3); plot(SOC, param.R2_init); grid on; title('R2')
% ax4 = subplot(4, 1, 4); plot(SOC, param.R3_init); grid on; title('R3')
% linkaxes([ax1,ax2,ax3,ax4], 'x')
% 
% % figure; hold on
% % plot(SOC, param.OCV_init)
% % plot(SOC, param.OCV_ub)
% % grid on; legend('init (lb)', 'ub'); title('OCV')
% 
% % figure; hold on
% % plot(SOC, param.tau1_init);
% % plot(SOC, param.tau2_init);
% % plot(SOC, param.tau3_init);
% % hold off; grid on; legend('tau1', 'tau2', 'tau3')

%% Save results

% Remove raw data from 'pulse' & only leave downsampled ones for optimization
pulse = removevars(pulse, 'segment1_t');
pulse = removevars(pulse, 'segment1');
pulse = removevars(pulse, 'segment2_t');
pulse = removevars(pulse, 'param');
pulse = removevars(pulse, 'Optim');

msgbox('Save the following parameters: meas_t, meas_resampled, pulse, Q')