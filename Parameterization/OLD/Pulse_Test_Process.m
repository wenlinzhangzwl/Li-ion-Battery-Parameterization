% clear
clearvars -except meas meas_t steps ind pulse

%% Get data & test parameters
currentFolder = cd; 
addpath(currentFolder);
data_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Test\';
results_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\';
function_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Model\Electric Model\Parameterization\functions';
addpath(function_folder)

cd(data_folder);
% load('PUL25chg.mat');
cd(currentFolder);

% Known parameters
temp = 25; 
Q = abs(meas.Ah(end));
dch = 0; % Whether it is charging or discharging

% %% Divide data into pulses
% % % Input step numbers (Discharging)
% % iRelax0 =   16;     % Rest before the start of the test
% % iLoad1 =    18;     % 10 charge/discharge pulses from 100% to 90% SOC
% % iRelax1 =   19;     % rest in between
% % iLoad2 =    22;     % 16 charge/discharge pulses from 90% to 10% SOC
% % iRelax2 =   23;     % rest in between
% % iLoad3 =    26;     % 10 charge/discharge pulses from 10% to 0% SOC
% % iRelax3 =   27;     % rest in between
% % iRelax4 =   29;     % rest after test reached Vmin
% 
% % Input step numbers (Charging)
% iRelax0 =   13;     % Rest before the start of the test
% iLoad1 =    15;     % 10 charge/discharge pulses from 100% to 90% SOC
% iRelax1 =   16;     % rest in between
% iLoad2 =    19;     % 16 charge/discharge pulses from 90% to 10% SOC
% iRelax2 =   20;     % rest in between
% iLoad3 =    23;     % 10 charge/discharge pulses from 10% to 0% SOC
% iRelax3 =   24;     % rest in between
% iRelax4 =   26;     % rest after test reached Vmin
% 
% steps = [iRelax0 iLoad1 iRelax1 iLoad2 iRelax2 iLoad3 iRelax3 iRelax4];
% 
% % Delete data before the first pulse
% meas_t = struct2table(meas);
% iBegin = find(meas.Step == min(steps), 1, 'first'); 
% iEnd = find(meas.Step == max(steps), 1, 'last'); 
% meas_t = meas_t(iBegin:iEnd, :);
% meas_t = meas_t(meas_t.Step <= max(steps), :);
% meas_t.Time = meas_t.Time - meas_t.Time(1);
% meas_t.Ah = meas_t.Ah - meas_t.Ah(1);
% meas_t.Current = -meas_t.Current;
% 
% % Delete data where the time is not strictly increasing
% iDelete = [];
% for i = 1:length(meas_t.Time)-1
%     if meas_t.Time(i) >= meas_t.Time(i+1)
%         iDelete = [iDelete; i];
%     end
% end
% meas_t(iDelete,:) = [];
% 
% % Calculate SOC
% if dch == 1
%     meas_t.SOC = 1 - (-meas_t.Ah)/Q;
% else
%     meas_t.SOC = 0 + (meas_t.Ah)/Q;
% end
% 
% figure; 
% ax1 = subplot(3, 1, 1); plot(meas_t.Time, meas_t.Current)
% ax2 = subplot(3, 1, 2); plot(meas_t.Time, meas_t.Voltage)
% ax3 = subplot(3, 1, 3); plot(meas_t.Time, meas_t.SOC)
% 
% % Divide data into pulses ****************************
% ind = getInd(meas_t, steps);
% pulse = getPulse(meas_t, ind);
% 
% % Plot segment1
% figure
% hold on; 
% plot(meas_t.Time/3600, meas_t.Voltage, '.');
% for i = 1:height(pulse)
%     t = pulse.segment1{i,1}.Time/3600;
%     Voltage = pulse.segment1{i,1}.Voltage;
%     plot(t, Voltage, '*');
% end
% title('segment1 (for time constant curve fitting)'); xlabel('Time [h]'); ylabel('Voltage'); grid on
% legend('data', 'sample')
% hold off
% 
% % Plot segment2
% figure
% hold on; 
% plot(meas_t.Time, meas_t.Voltage, '.');
% for i = 1:height(pulse)
%     t = pulse.segment2{i,1}.Time;
%     Voltage = pulse.segment2{i,1}.Voltage;
%     plot(t, Voltage, 'o');
% end
% title('segment2 (for optimization)'); xlabel('Time [s]'); ylabel('Voltage'); grid on
% legend('data', 'sample')
% hold off

%% Find initial conditions for each pulse

% Find initial conditions
param = getInitialcondition(meas_t, pulse, dch);
% for i = 1:height(pulse)
%     pulse.param(i) = {param(i, :)};
% end

% Plot parameters for validation
SOC = param.SOC;

figure; 
plot(SOC, param.R0_init, '.-'); grid on; title('R0')

% figure;
% ax1 = subplot(4, 1, 1); plot(SOC, param.R0_init); grid on; title('R0')
% ax2 = subplot(4, 1, 2); plot(SOC, param.R1_init); grid on; title('R1')
% ax3 = subplot(4, 1, 3); plot(SOC, param.R2_init); grid on; title('R2')
% ax4 = subplot(4, 1, 4); plot(SOC, param.R3_init); grid on; title('R3')
% linkaxes([ax1,ax2,ax3,ax4], 'x')

figure; hold on
plot(SOC, param.OCV_init)
% plot(SOC, param.OCV_ub)
% grid on; legend('init (lb)', 'ub'); title('OCV')

% figure; hold on
% plot(SOC, param.tau1_init);
% plot(SOC, param.tau2_init);
% plot(SOC, param.tau3_init);
% hold off; grid on; legend('tau1', 'tau2', 'tau3')

% % Remove raw data from 'pulse' & only leave downsampled ones for optimization
% pulse = removevars(pulse, 'segment1_t');
% pulse = removevars(pulse, 'segment2_t');
% pulse = removevars(pulse, 'Optim');
% 
% Time = meas_t.Time; 
% validateattributes(Time, {'double'}, {'increasing'})

param_chg = [param.SOC param.R0_init];

%% Run ga()
% global Batt
% Batt.Q = Q;
% % Batt.coeff = coeff;
% % Batt.OCV_SOC = OCV_SOC;
% 
% % Run GA optimizer for each segment
% for i = 1:height(pulse)
%     %% Define GA optimization options
% 
%     % Initial guess & upper/lower bounds
%     init_optim = table2array(param(i, 1:8));
%     lb_optim = table2array(param(i, 9:16));
%     ub_optim = table2array(param(i, 17:24));
% 
%     % Define options
%     PopulationSize_Data = 100; 
%     Generations_Data = 100; 
%     options = optimoptions('ga');  
%     options = optimoptions(options, 'PopInitRange', [lb_optim; ub_optim],...
%                                     'PopulationSize', PopulationSize_Data,...
%                                     'Generations', Generations_Data,...
%                                     'InitialPopulation', init_optim,...
%                                     'Display', 'off',...
%                                     'PlotFcns', {@gaplotbestf @gaplotbestindiv},...
%                                     'Vectorized', 'off',...
%                                     'UseParallel', false,...
%                                     'MaxStallGenerations', 10);
% 
%     nvars = 8; 
%     
%     %% Run ga()
%     Batt.Time = pulse.segment1{i,1}.Time;
%     Batt.Voltage_exp = pulse.segment1{i,1}.Voltage;
%     Batt.Current = pulse.segment1{i,1}.Current;
%     Batt.SOC_exp = pulse.segment1{i,1}.SOC;
% 
%     [results, fval, exitflag, output, population, score] = ga(@RunOCVRRC3Model, nvars, [], [], [], [], lb_optim, ub_optim, [], [], options);
%     
%     %% Record results
%     Optim = struct('results', results,... 
%                    'fval', fval,...
%                    'exitflag', exitflag,...
%                    'output', output,...
%                    'population', population,...
%                    'score', score,...
%                    'SOC', pulse{i, 'SOC'},...
%                    'RMSE_Vt', Batt.RMSE_Vt,...
%                    'RMSE_SOC', Batt.RMSE_SOC);
%     pulse{i,'Optim'} = {Optim};
%     
% end
% 
% % Display & save optimization results
% fieldNum = 8;
% OptimResult = [];
% for i = 1:height(pulse)
%     oOptimResult = struct2table(pulse.Optim{i,1}, 'AsArray',true);
%     OptimResult = [OptimResult; oOptimResult];
% end
% 
% % Export 'pulse' & 'OptimResult' tables for later use
% ofilename = append(folder, filename, '_pulse.mat');
% save(ofilename, 'pulse');
% ofilename = append(folder, filename, '_OptimResult.mat');
% save(ofilename, 'OptimResult');
% 
% %% Plot optimization results
% 
% % Get parameters
% R0 = OptimResult.results(:, 1);
% R1 = OptimResult.results(:, 2);
% R2 = OptimResult.results(:, 3);
% R3 = OptimResult.results(:, 4);
% tau1 = OptimResult.results(:, 5);
% tau2 = OptimResult.results(:, 6);
% tau3 = OptimResult.results(:, 7);
% OCV = OptimResult.results(:, 8);
% 
% % Plot resistances
% figure; 
% hold on
% plot(OptimResult.SOC, R1*1000);
% plot(OptimResult.SOC, R2*1000);
% plot(OptimResult.SOC, R3*1000);
% plot(OptimResult.SOC, R0*1000);
% hold off
% legend('R1', 'R2', 'R3', 'R0');
% title('Resistances'); grid minor; xlabel('SOC [%]'); ylabel('R [mOhm]')
% 
% % Plot time constants
% figure; 
% hold on
% plot(OptimResult.SOC, tau1);
% plot(OptimResult.SOC, tau2);
% plot(OptimResult.SOC, tau3);
% hold off
% legend('tau1', 'tau2', 'tau3');
% title('tau'); grid minor; xlabel('SOC [%]'); ylabel('tau [s]')

%% Functions

function [OBJ] = RunOCVRRC3Model(Opt_Param)
% Objective function for GA optimization
    
    %% Run Model with Optimization Parameters
    global Batt;   
    Time = Batt.Time;
    Vt_exp = Batt.Voltage_exp;
    SOC_exp = Batt.SOC_exp; 
    
    [Vt_sim, SOC_sim] = OCVRRC3Model(Opt_Param);

    %% Calculate the Terminal Voltage errors
    error_Vt = Vt_exp - Vt_sim;
    error_SOC = SOC_exp - SOC_sim;
    Batt.RMSE_Vt = sqrt(mean(error_Vt.^2))*1000;
    Batt.RMSE_SOC = sqrt(mean(error_SOC.^2))*1000;

    OBJ = Batt.RMSE_Vt;
    
%     % DEBUGGING
%     figure; hold on
%     plot(Batt.Time, Batt.Voltage_exp);
%     plot(Batt.Time, Vt_sim);
%     legend('exp', 'sim')
end

function [Vt, SOC] = OCVRRC3Model(param)
% Simulates battery with parameters from optimization

    % Optimization Parameters
    R0 = param(1);
    R1 = param(2);
    R2 = param(3); 
    R3 = param(4);
    tau1 = param(5);
    tau2 = param(6);
    tau3 = param(7);
    OCV = param(8);

    % Battery Fixed/Known Parameters 
    global Batt
    Q = Batt.Q*3600;           % Capacity [Amp*s]
    Current = Batt.Current;    % Current {Amp]
    
    %% Run Actual Battery Model
    SOC = [Batt.SOC_exp(1)];
    Vt  = [Batt.Voltage_exp(1)];
    Irc1 = [0];
    Irc2 = [0];
    Irc3 = [0];

for k = 2:length(Current)

   oSOC = SOC(k-1);
   deltaT = Batt.Time(k) - Batt.Time(k-1);
   
   % Current through RC branchs
    oIrc1 =  ( 1 - exp(-deltaT/tau1) ) * Current(k) + exp(-deltaT/tau1) * Irc1(k-1);
    oIrc2 =  ( 1 - exp(-deltaT/tau2) ) * Current(k) + exp(-deltaT/tau2) * Irc2(k-1);
    oIrc3 =  ( 1 - exp(-deltaT/tau3) ) * Current(k) + exp(-deltaT/tau3) * Irc3(k-1);
    
    % Calculate current Vt & SOC
    oVt   = OCV + (R0  * Current(k)) + R1 * oIrc1 + R2 * oIrc2 + R3 * oIrc3;
    oSOC  = oSOC + (deltaT / Q)* Current(k);
    
    % Update parameters
    SOC = [SOC; oSOC];
    Vt  = [Vt; oVt];
    Irc1 = [Irc1; oIrc1];
    Irc2 = [Irc2; oIrc2];
    Irc3 = [Irc2; oIrc3];
    
end

end

