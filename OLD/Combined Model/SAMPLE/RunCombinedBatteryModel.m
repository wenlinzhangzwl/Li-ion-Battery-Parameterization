clear

%% Load input current data from the test
% folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Prismatic Cell\Test Results\1 - Characterization test 1\';
% file = append(folder, '12-11-20_08.01 1336_Charge3_HPPC.mat');
% load(file);
% global Batt; 
% Batt.RecordingTime = meas.Time; 
% Batt.I = -meas.Current; 

file = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Prismatic Cell\Model\Electric Model\Combined Model\SAMPLE\UDDS_50.csv';
Data = xlsread('UDDS_50.csv', 'A69:AD13771');
global Batt;
Batt.RecordingTime          = Data(:,1);
Batt.I                      = -Data(:,2);

% Runs battery model
[Batt.SOC_Actual, Batt.V_Actual] = CombinedBatteryModel(Batt.I, Batt.RecordingTime);


%% Ploting
figure
subplot(3,1,1)
plot(Batt.RecordingTime/60, Batt.SOC_Actual * 100);
xlabel('time')
ylabel('SOC [%]')
grid minor

subplot(3,1,2)
plot(Batt.RecordingTime/60, Batt.V_Actual)
ylabel('TerminalVoltage [V]')
grid minor

subplot(3,1,3)
plot (Batt.RecordingTime/60, Batt.I)
xlabel('Time')
ylabel('Current [I]')
grid minor