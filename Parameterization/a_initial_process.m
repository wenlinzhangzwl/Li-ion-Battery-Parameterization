clear

%% Get data & test conditions

% Folders
folder_current = cd; 
folder_functions = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\Scripts\Parameterization\functions";
folder_data = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\10 - Characterization\";
% folder_result = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\7 - Pulse Test\";  % where results from this script is saved
addpath(folder_current, folder_functions, folder_data);

% Load data
% load("C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\8 - Capacity test & US06\Results\Capacity_25degC_Sept23_2021.mat");
load("C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\10 - Characterization\Results\Capacity_25degC_March22_2022.mat");
load('PUL_25degC_0p8_30min.mat');
meas.Current = -meas.Current; % Assume +ve current = discharge

% Plot data to identify steps
figure; 
ax(1) = subplot(5, 1, 1); plot(meas.Time, meas.Voltage, '.-'); title("Voltage"); xlabel("Time [s]"); ylabel("Voltage [V]"); grid on
ax(2) = subplot(5, 1, 2); plot(meas.Time, meas.Current, '.-'); title("Current"); xlabel("Time [s]"); ylabel("Current [A]"); grid on
ax(3) = subplot(5, 1, 3); plot(meas.Time, meas.Ah, '.-'); title("Ah"); xlabel("Time [s]"); ylabel("Ah"); grid on
ax(4) = subplot(5, 1, 4); plot(meas.Time, meas.Battery_Temp_degC, '.-'); title("Temperature"); xlabel("Time [s]"); ylabel("Ah"); grid on
ax(5) = subplot(5, 1, 5); plot(meas.Time, meas.Step, '.-'); title("Step"); xlabel("Time [s]"); ylabel("Ah"); grid on; ylim([10,50]);
linkaxes(ax, 'x')

%% Divide data into pulses

% Step numbers *******************************************CHANGE
iRelax0 =   99;   % Rest before the start of the test. Needed for MATLAB toolbox to recognize the first pulse. Script will ddd in an artifical data point if set to 99
iLoad1 =    15;     % 10 charge/discharge pulses from 100% to 90% SOC
iRelax1 =   16;     % rest in between
iLoad2 =    19;     % 16 charge/discharge pulses from 90% to 10% SOC
iRelax2 =   20;     % rest in between
iLoad3 =    23;     % 10 charge/discharge pulses from 10% to 0% SOC
iRelax3 =   24;     % rest in between
iRelax4 =   26;     % rest after test reached Vmin
steps = [iRelax0 iLoad1 iRelax1 iLoad2 iRelax2 iLoad3 iRelax3 iRelax4];

meas_t = struct2table(meas);

% Add in an artifical data point for MATLAB script to recognize the first pulse
if steps(1) == 99
    fake_data = meas_t(1, :);
    fake_data.Current = 0; 
    fake_data.Step = 99; 

    meas_t.Time = meas_t.Time + meas_t.Time(2);
    meas_t = [fake_data; meas_t];
    
end

% Correct current sign, negative = charge, positive = charge
meas.Current = -meas.Current; 


% Trim unnecessary data at beginning & end
% iBegin = find(meas.Step == min(steps), 1, 'first');
% if isempty(iBegin)
%     iBegin = 1; 
% end
% iEnd = find(meas.Step == max(steps), 1, 'last'); 
% meas_t = meas_t(iBegin:iEnd, :);
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

meas = meas_t;

% msgbox('Save the following parameters: meas_t, meas_resampled, pulse, Q')
msgbox('Save the following parameters: meas_t, meas')