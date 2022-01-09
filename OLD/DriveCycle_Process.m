clear
% clearvars -except meas meas_t

%% Get data & test parameters
current_folder = cd;
data_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\6 - Drive cycle\';
cd(data_folder);
addpath(current_folder);
filename = 'UDDS25';
file = append(data_folder, filename, '.mat');
load(file)

%% Load info to run simulation
sim_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\';

sim_deg = '19'; %********************************************CHANGE

% parameters
% filename_param = append(sim_folder, 'characterization_', string(sim_deg));
% load(filename_param);
Q = 279;
results_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\Aug 18\';
num = 33;
load(append(results_folder, 'param_opt19_', string(num), '.mat'))

sim_case = 'psParam_v0_optim';
switch sim_case
    case 19
        num = 33;
        load(append('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\Aug 18\param_opt19_', string(num), '.mat'))

        R0 = param.R0_init;
        R1 = param.R1_init;
        R2 = param.R2_init;
        R3 = param.R3_init;
        tau1 = param.tau1_init;
        tau2 = param.tau2_init;
        tau3 = param.tau3_init;
        
%         R0(1:num+1) = param.R0_opt(1:num+1);
%         R1(1:num+1) = param.R1_opt(1:num+1);
%         R2(1:num+1) = param.R2_opt(1:num+1);
%         R3(1:num+1) = param.R3_opt(1:num+1);
%         tau1(1:num+1) = param.tau1_opt(1:num+1);
%         tau2(1:num+1) = param.tau2_opt(1:num+1);
%         tau3(1:num+1) = param.tau3_opt(1:num+1);

    case 'psParam_v0'
        load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\psParam_v0.mat')
        SOC = psParam.SOC;
        R0 = psParam.R0;
        R1 = psParam.Rx(1, :);
        R2 = psParam.Rx(2, :);
        R3 = psParam.Rx(3, :);
        tau1 = psParam.Tx(1, :);
        tau2 = psParam.Tx(2, :);
        tau3 = psParam.Tx(3, :);
    case 'psParam_v0_optim'
        load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\psParam_v0_optim.mat')
        SOC = psParam.SOC;
        R0 = psParam.R0;
        R1 = psParam.Rx(1, :);
        R2 = psParam.Rx(2, :);
        R3 = psParam.Rx(3, :);
        tau1 = psParam.Tx(1, :);
        tau2 = psParam.Tx(2, :);
        tau3 = psParam.Tx(3, :);
end

% OCV curve
filename_OCV = append(sim_folder, 'C20CCCV_OCVSOC', string(sim_deg));
load(filename_OCV);

% Trim & process data
meas_t = struct2table(meas);
ind = find(meas.Voltage <=2.5, 1);
meas_t = meas_t(1:ind, :);
meas_t.Ah = meas_t.Ah - meas_t.Ah(1);
meas_t.Time = meas_t.Time - meas_t.Time(1);
meas_t.SOC = 1 - (-meas_t.Ah)/Q;
% figure; plot(meas_t.Time, meas_t.Voltage); grid on

% Inputs
stepSize = 0.1;
time = meas_t.Time; 
t = [time time];
I = [time meas_t.Current]; 
runTime = meas_t.Time(end);
SOC_init = meas_t.SOC(1);
I1_init = 0;
I2_init = 0;
I3_init = 0;

% Experimental results
Vt_exp = [time meas_t.Voltage];
SOC_exp = [time meas_t.SOC];


%% Run simulation
model = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Model\Electric Model\Parameterization\ElectricalModel.slx';
sim_result = sim(model);

%% Simulation results

% Results
t_sim = sim_result.Vt.time;
Vt_sim = sim_result.Vt.signals(1).values;
Vt_exp2 = sim_result.Vt.signals(2).values;
SOC_exp2 = sim_result.SOC.signals(1).values;
SOC_sim = sim_result.SOC.signals(2).values;

% Errors
rmse_Vt = sqrt(mean((Vt_exp2 - Vt_sim).^2));

% % Save results
% sim_resultsFolder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\5 - Drive cycle\Results\';
% sim_resultsFilename = append(sim_resultsFolder, 'DriveCycle', string(sim_deg));
% save(sim_resultsFilename, 't_sim', 'Vt_sim', 'Vt_exp2', 'SOC_exp2', 'SOC_sim');

%% Plot results

% Vt drivecycle vs Time
figure; hold on
plot(t_sim, Vt_exp2, '-');
plot(t_sim, Vt_sim, '--');
title(append('Vt drivecycle ', string(sim_deg))); xlabel('Time (s)'); ylabel('Voltage (V)'); legend('exp', 'sim'); grid on;
annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt))

% % Vt_err drivecycle vs Time
% figure; hold on
% plot(t_sim, Vt_exp2-Vt_sim);
% title(append('Vt_err drivecycle ', string(sim_deg)), 'Interpreter', 'none'); xlabel('Time (s)'); ylabel('Voltage (V)'); grid on;
% annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt))

% Voltage vs SOC
figure; hold on
plot(SOC_sim, Vt_exp2, '-');
plot(SOC_sim, Vt_sim, '--');
title(append('Vt drivecycle ', string(sim_deg))); xlabel('SOC'); ylabel('Voltage (V)'); legend('exp', 'sim'); grid on;
annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt))

% % Voltage_err vs SOC
% figure; hold on
% plot(SOC_sim, Vt_exp2-Vt_sim);
% title(append('Vt_err drivecycle ', string(sim_deg)), 'Interpreter', 'none'); xlabel('SOC'); ylabel('Voltage (V)'); grid on;
% annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt))

% % SOC_exp vs SOC_sim
% figure; hold on
% plot(t_sim, SOC_exp2);
% plot(t_sim, SOC_sim, '--');
% title(append('SOC drivecycle ', string(sim_deg))); xlabel('SOC'); ylabel('Voltage (V)'); legend('exp', 'sim'); grid on;

% % Plot resistances
% figure; 
% subplot(2, 2, 1); plot(SOC, R0*1000); title('R0'); grid on; xlabel('SOC'); ylabel('R [mOhm]')
% subplot(2, 2, 2); plot(SOC, R1*1000); title('R1'); grid on; xlabel('SOC'); ylabel('R [mOhm]')
% subplot(2, 2, 3); plot(SOC, R2*1000); title('R2'); grid on; xlabel('SOC'); ylabel('R [mOhm]')
% subplot(2, 2, 4); plot(SOC, R3*1000); title('R3'); grid on; xlabel('SOC'); ylabel('R [mOhm]')
% 
% % Plot time constants
% figure; 
% subplot(3, 1, 1); plot(SOC, tau1); title('tau1'); grid on; xlabel('SOC'); ylabel('tau [s]')
% subplot(3, 1, 2); plot(SOC, tau2); title('tau2'); grid on; xlabel('SOC'); ylabel('tau [s]')
% subplot(3, 1, 3); plot(SOC, tau3); title('tau3'); grid on; xlabel('SOC'); ylabel('tau [s]')
 

