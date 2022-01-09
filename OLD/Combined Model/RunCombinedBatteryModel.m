clear

%% Load input current data from the test
folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Prismatic Cell\Test Results\1 - Characterization test 1\';
file = append(folder, '12-11-20_08.01 1336_Charge3_HPPC.mat');
load(file);
samplingTime = 10;
eta25 = 0.99;

% Calculate other model parameters
Q = abs(meas.Ah(end) - meas.Ah(1));                                         % Total capacity
meas.SOC = 1 - (-meas.Ah)/Q;                                                % SOC

global Batt;
Batt = struct('RecordingTime', meas.Time(4000:36000), ...
              'I', -meas.Current(4000:36000),...
              'SOC_Actual', meas.SOC(4000:36000),...
              'V_Actual', meas.Voltage(4000:36000),...
              'samplingTime', samplingTime,...
              'Q', Q,...
              'eta', eta25);

% Runs battery model
[Batt.SOC_Actual, Batt.V_Actual] = CombinedBatteryModel(Batt);


%% Ploting
% figure
% % subplot(3,1,1)
% plot(Batt.RecordingTime, Batt.SOC_Actual * 100, 'b');
% xlabel('time[s]')
% ylabel('SOC [%]')
% title('SOC')
% grid minor
% hold on
% plot(meas.Time(4000:36000), meas.SOC(4000:36000)*100, 'r');
% legend('sim', 'exp')


figure % subplot(3,1,2)
plot(Batt.RecordingTime, Batt.V_Actual, 'b')
xlabel('time[s]')
ylabel('TerminalVoltage [V]')
title('Voltage')
grid minor
hold on
plot(meas.Time(4000:36000), meas.Voltage(4000:36000), 'r');
legend('sim', 'exp')


% figure % subplot(3,1,3)
% plot (Batt.RecordingTime, Batt.I, 'b')
% xlabel('Time[s]')
% ylabel('Current [I]')
% title('Current')
% grid minor
% hold on
% plot(meas.Time, meas.Current, 'r');
% legend('sim', 'exp')