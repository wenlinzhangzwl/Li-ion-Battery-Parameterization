clearvars -except Batt meas

%% Load and process data from the experiment

% Load files
currentFolder = cd; 
addpath(currentFolder);
folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\4 - Characterization test 3\';
cd(folder);
% load('12-11-20_08.01 1336_Charge3_HPPC.mat'); 
load('OCV_SOC.mat');
load('PUL25dchOptimResult');


% Run battery model
global Batt;
% Batt = struct('Time', meas.Time, ...
%               'Current', -meas.Current,...
%               'SOC_exp', meas.SOC,...
%               'Vt_exp', meas.Voltage,...
%               'Q', Q, ...
%               'OCVSOC', OCV_SOC);
[Batt.V_sim, Batt.SOC_sim] = OCVRRCModel(Batt);
clearvars -except Batt meas OptimResult

%% Save simulation data

% meas.Time = Batt.Time; 
% meas.Voltage = Batt.V_sim;
% meas.Current = -Batt.Current; 
% meas.SOC = Batt.SOC_sim;
% meas.Q = Batt.Q;
% save('meas.mat', 'meas')

%% Ploting
% figure
% plot(Batt.Time, Batt.SOC_sim * 100, 'b');
% xlabel('time[s]')
% ylabel('SOC [%]')
% title('SOC')
% grid minor
% hold on
% plot(Batt.Time, Batt.SOC_exp * 100, 'r');
% legend('sim', 'exp')

figure
plot(Batt.Time, Batt.V_sim, 'b')
xlabel('time[s]')
ylabel('TerminalVoltage [V]')
title('Voltage')
grid minor
hold on
plot(Batt.Time, Batt.Voltage_exp, 'r');
hold off;
legend('sim', 'exp')
% xlim([8950,11000]);

% hold on
% plot(Batt.Time, Batt.V_exp);
% legend('sim1', 'exp', 'sim2')

% figure
% plot (Batt.Time, Batt.Current, 'b')
% xlabel('Time[s]')
% ylabel('Current [I]')
% title('Current')
% grid minor