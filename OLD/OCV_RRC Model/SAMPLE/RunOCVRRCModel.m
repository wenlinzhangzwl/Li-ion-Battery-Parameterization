clear

%% Load and process data from the experiment

% Load files
folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Prismatic Cell\Test Results\1 - Characterization test 1\';
file = append(folder, '12-11-20_08.01 1336_Charge3_HPPC.mat');
load(file);
file = append(folder,'OCV_SOC.mat');
load(file)

% Model parameters
samplingTime = 10;
eta = 0.99;
Q = abs(meas.Ah(end) - meas.Ah(1));     % Total capacity
meas.SOC = 1 - (-meas.Ah)/Q;            % SOC

% Run battery model
global Batt;
Batt = struct('t', meas.Time(4000:36000), ...
              'Current', -meas.Current(4000:36000),...
              'SOC_exp', meas.SOC(4000:36000),...
              'V_exp', meas.Voltage(4000:36000),...
              'samplingTime', samplingTime,...
              'Q', abs(meas.Ah(36000) - meas.Ah(4000)),...
              'eta', eta,...
              'OCVSOC', OCV_SOC);
[Batt.SOC_sim, Batt.V_sim] = OCVRRCModel(Batt);

%% Save simulation data
clearvars -except Batt
meas.Time = Batt.t; 
meas.Voltage = Batt.V_sim;
meas.Current = -Batt.Current; 
meas.SOC = Batt.SOC_sim;
meas.Q = Batt.Q;
save('meas.mat', 'meas')

%% Ploting
% figure
% plot(Batt.t, Batt.SOC_sim * 100, 'b');
% xlabel('time[s]')
% ylabel('SOC [%]')
% title('SOC')
% grid minor
% hold on
% plot(Batt.t, Batt.SOC_exp * 100, 'r');
% legend('sim', 'exp')
% 
% figure
% plot(Batt.t, Batt.V_sim, 'b')
% xlabel('time[s]')
% ylabel('TerminalVoltage [V]')
% title('Voltage')
% grid minor
% hold on
% plot(Batt.t, Batt.V_exp, 'r');
% legend('sim', 'exp')

% figure
% plot (Batt.t, Batt.Current, 'b')
% xlabel('Time[s]')
% ylabel('Current [I]')
% title('Current')
% grid minor