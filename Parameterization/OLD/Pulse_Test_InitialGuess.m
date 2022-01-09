clear
% clearvars -except meas meas_t steps ind pulse

%% Get data & test conditions
folder_current = cd; 
addpath(folder_current);
folder_data = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Discharge\';               % where experimental data is saved
folder_result = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Discharge\Results\';     % where results from this script is saved
folder_fcn = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Model\Electric Model\Parameterization\functions\'; % where functions used in this script is saved

% Load data
cd(folder_data);
load('PUL25dch.mat');
cd(folder_current);
Q = -meas.Ah(end);

% Test conditions
temp = 25; 
dch = 1; % Whether it is charging or discharging

%% Divide data into pulses
% Input step numbers
iRelax0 =   16;     % Rest before the start of the test
iLoad1 =    18;     % 10 charge/discharge pulses from 100% to 90% SOC
iRelax1 =   19;     % rest in between
iLoad2 =    22;     % 16 charge/discharge pulses from 90% to 10% SOC
iRelax2 =   23;     % rest in between
iLoad3 =    26;     % 10 charge/discharge pulses from 10% to 0% SOC
iRelax3 =   27;     % rest in between
iRelax4 =   29;     % rest after test reached Vmin
steps = [iRelax0 iLoad1 iRelax1 iLoad2 iRelax2 iLoad3 iRelax3 iRelax4];

% Delete data before the first pulse
meas_t = struct2table(meas);
iBegin = find(meas.Step == min(steps), 1, 'first'); 
iEnd = find(meas.Step == max(steps), 1, 'last'); 
meas_t = meas_t(iBegin:iEnd, :);
meas_t = meas_t(meas_t.Step <= max(steps), :);
meas_t.Time = meas_t.Time - meas_t.Time(1);
meas_t.Ah = meas_t.Ah - meas_t.Ah(1);

% Delete data where the time is not strictly increasing
iDelete = [];
for i = 1:length(meas_t.Time)-1
    if meas_t.Time(i) >= meas_t.Time(i+1)
        iDelete = [iDelete; i];
    end
end
meas_t(iDelete,:) = [];

% Calculate SOC
if dch == 1
    meas_t.SOC = 1 - (-meas_t.Ah)/Q;
else
    meas_t.SOC = 0 + (meas_t.Ah)/Q;
end

% Divide data into pulses ****************************
cd(folder_fcn);
ind = getInd(meas_t, steps);
pulse = getPulse(meas_t, ind);
cd(folder_current);

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
pulse = removevars(pulse, 'segment2_t');
pulse = removevars(pulse, 'Optim');

Time = meas_t.Time; 
validateattributes(Time, {'double'}, {'increasing'})

msgbox('Save the following parameters: pulse, Q, temp')