% % Calculates electrical parameters at different SOC
% clearvars -except meas
% 
% %% Get data & test parameters
% current_folder = cd;
% folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\4 - Characterization test 3\';
% cd(folder);
% addpath(current_folder);
% filename = 'PUL25chg';
% file = append(folder, filename, '.mat');
% % load(file)
% load('C20CCCV_OCVSOC.mat');            % SOC OCV relationship
% 
% % Known parameters
% temp = 25; 
% Q = 279.29;
% 
% %% Divide data into different steps
% % Input step numbers
% iDch1 = 12;     % 10 discharge pulses from 100% to 90% SOC
% iPau1 = 13;     % 1 hr rest in between
% iDch2 = 16;     % 16 discharge pulses from 90% to 10% SOC
% iPau2 = 17;     % 4 hr rest in between
% iDch3 = 20;     % 10 discharge pulses from 10% to 0% SOC
% iPau3 = 21;     % 1 hr rest in between
% steps = [iDch1 iPau1 iDch2 iPau2 iDch3 iPau3];
% 
% % Delete data before the first discharge pulse
% meas_t = struct2table(meas);
% meas_t = meas_t(meas.Step >= iDch1, :);
% meas_t.Ah = meas_t.Ah - meas_t.Ah(1);
% 
% % Calculate SOC
% meas_t.SOC = meas_t.Ah / Q;
% 
% % Find indices of each step in meas
%     % Returns matrix 'ind': 
%     % column1 = step number; 
%     % column2 = all indices of this step number
%     % column3 = indices of segments of this step number 
% pulse = table('Size', [length(steps), 3],...
%                'VariableTypes', ["double","double", "double"],...
%                'VariableNames',["stepNum", "indices_all", "indices_seg"]);    
% for i = 1:length(steps)
%     ind{i,1} = steps(i);
%     ind{i,2} = find(meas_t.Step == steps(i));
%     indices = ind{i,2};             % Indices of data points with this step number
%     segments = [indices(1) 0];      % Start & end points of each segment with this start number
%     p = 1;                          % Counter
%     
%     for j = 1:length(indices)
%         if j == length(indices)
%             segments(p,2) = indices(j);
%         elseif indices(j+1) ~= indices(j)+1
%             segments(p,2) = indices(j);
%             segments(p+1,1) = indices(j+1);
%             p = p + 1;
%         end
%     end
%     
%     ind{i,3} = segments; 
% end
% 
% %% Find each pulse-rest segment
% fields = string(fieldnames(meas));
% 
% % table containing info of each segment
% % column1 = step number, 
% % column2 = number of pulse in that step, 
% % column3 = beginning & end indices of pulse
% % column4 = beginning & end indices of corresponding rest    
% pulse = table('Size', [length(steps)/2, 6],...
%                'VariableTypes', ["cell","double", "cell", "cell", "double", "cell"],...
%                'VariableNames',["stepNum", "pulseNum", "segment_t", "segment", "SOC", "Optim"]);    
% 
% m = 1;
% n = 1;
% for i = 1:2:length(steps)
%     ind1 = ind{i,3};
%     ind2 = ind{i+1,3};
%     
%     % Check if pulse is cut off because Vmin was reached
%     if length(ind1) == length(ind2)
%         len = length(ind1);
%     else
%         len = length(ind2);
%     end
%     
%     % Form 'pulses' table
%     for j = 1:len-1
%         % Segment information
%         pulse.stepNum(n) = {[ind{i,1}; ind{i+1,1}]};
%         pulse.pulseNum(n) = j;
% 
%         % Divide into segments
%         iBeg = ind1(j, 1);
%         iEnd = ind1(j+1, 2);
%         segment_t = meas_t(iBeg:iEnd, :);
%         pulse.segment_t(n) = {segment_t};
%         
%         % Resample points from segment
%         div1 = 100;
%         div2 = 100;
%         div3 = 200;
%         inc1 = 1;
%         inc2 = 30;
%         inc3 = fix((height(segment_t) - 2*div1*inc1 - 2*div2*inc2) / (div3));
%         div = (div1 + div2)*2 + div3 + 1; % Total number of points sampled
%         p = 1;
%         for k = 1:div
%             if or( k<=div1, k>div1+div2 && k<=2*div1+div2)
%                 segment(k,:) = segment_t(p,:);
%                 p = p + inc1;
%             elseif or( k>div1 && k<=div1+div2, k>2*div1+div2 && k<=2*(div1+div2) ) 
%                 segment(k,:) = segment_t(p,:);
%                 p = p + inc2;
%             elseif k == div
%                 segment(k,:) = segment_t(end,:); 
%             else
%                 segment(k,:) = segment_t(p,:);
%                 p = p + inc3; 
%             end
%         end
%         pulse.segment(n) = {segment};
%                 
%         % SOC
%         pulse.SOC(n) = segment_t.SOC(2);
%         
%         n = n + 1;
%     end
%     
%     m = m + 1; 
% end
% 
% % Plot original & resampled points
% figure
% hold on; 
% plot(meas.Time/3600, meas.Voltage, '.');
% for i = 1:height(pulse)
%     t = pulse.segment{i,1}.Time/3600;
%     Voltage = pulse.segment{i,1}.Voltage;
%     plot(t, Voltage, '*');
% end
% title('Resampling'); xlabel('Time [h]'); ylabel('Voltage'); grid on
% legend('data', 'sample')
% hold off

for i = 1:height(pulse)
     seg = pulse.segment_t{i};
     OCV(i, 1) = seg.Voltage(end);
end

%% Define GA optimization options

% Initialization
R0_init = 0.25e-3;
R1_init = 0.5e-3;
R2_init = 0.5e-3; 
R3_init = 0.5e-3; 
C1_init = 5e3; 
C2_init = 5e3; 
C3_init = 5e3; 
InitialPopulation_Data = [R0_init, R1_init, R2_init, R3_init, C1_init, C2_init, C3_init];

% Lower and upper bounds of varialbes
R0_lb = 0.1e-3; R0_ub = 2.5e-3; 
R1_lb = 0.1e-3; R1_ub = 5e-3; 
R2_lb = 0.1e-3; R2_ub = 5e-3; 
R3_lb = 0; R3_ub = 2.5e-3;
C1_lb = 1e3; C1_ub = 1000e3;
C2_lb = 1e3; C2_ub = 1000e3;
C3_lb = 0; C3_ub = 1000e3;
lb = [R0_lb, R1_lb, R2_lb, R3_lb, C1_lb, C2_lb, C3_lb];
ub = [R0_ub, R1_ub, R2_ub, R3_ub, C1_ub, C2_ub, C3_lb];

% Define options
PopulationSize_Data = 100; 
Generations_Data = 100; 
options = optimoptions('ga');  
options = optimoptions(options, 'PopInitRange', [lb; ub],...
                                'PopulationSize', PopulationSize_Data,...
                                'Generations', Generations_Data,...
                                'InitialPopulation', InitialPopulation_Data,...
                                'Display', 'off',...
                                'PlotFcns', {@gaplotbestf @gaplotbestindiv},...
                                'Vectorized', 'off',...
                                'UseParallel', false,...
                                'MaxStallGenerations', 10);

nvars = 7; 

%% Run ga()
global Batt
Batt.Q = Q;
Batt.coeff = coeff;
Batt.OCV_SOC = OCV_SOC;

% Run GA optimizer for each segment
for i = 1:height(pulse)
    Batt.Time = pulse.segment{i,1}.Time;
    Batt.Voltage_exp = pulse.segment{i,1}.Voltage;
    Batt.Current = pulse.segment{i,1}.Current;
    Batt.SOC_exp = pulse.segment{i,1}.SOC;

    [results, fval, exitflag, output, population, score] = ga(@RunOCVRRC3Model, nvars, [], [], [], [], lb, ub, [], [], options);
    
    Optim = struct('results', results,... 
                   'fval', fval,...
                   'exitflag', exitflag,...
                   'output', output,...
                   'population', population,...
                   'score', score,...
                   'SOC', pulse{i, 'SOC'},...
                   'RMSE_Vt', Batt.RMSE_Vt,...
                   'RMSE_SOC', Batt.RMSE_SOC);
    pulse{i,'Optim'} = {Optim};
end

% Display & save optimization results
fieldNum = 8;
OptimResult = [];
for i = 1:height(pulse)
    oOptimResult = struct2table(pulse.Optim{i,1}, 'AsArray',true);
    OptimResult = [OptimResult; oOptimResult];
end

% Export 'pulse' & 'OptimResult' tables for later use
ofilename = append(folder, filename, '_pulse.mat');
save(ofilename, 'pulse');
ofilename = append(folder, filename, '_OptimResult.mat');
save(ofilename, 'OptimResult');

%% Plot optimization results

% Get parameters
R0 = OptimResult.results(:, 1);
R1 = OptimResult.results(:, 2);
R2 = OptimResult.results(:, 3);
R3 = OptimResult.results(:, 4);
C1 = OptimResult.results(:, 5);
C2 = OptimResult.results(:, 6);
C3 = OptimResult.results(:, 7);

% Plot resistances
figure; 
hold on
plot(OptimResult.SOC*100, R1*1000);
plot(OptimResult.SOC*100, R2*1000);
plot(OptimResult.SOC*100, R3*1000);
plot(OptimResult.SOC*100, R0*1000);
hold off
legend('R1', 'R2', 'R3', 'R0');
title('Resistances'); grid minor; xlabel('SOC [%]'); ylabel('R [mOhm]')

% Plot capacitances
figure; 
hold on
plot(OptimResult.SOC*100, C1/1000);
plot(OptimResult.SOC*100, C2/1000);
plot(OptimResult.SOC*100, C3/1000);
hold off
legend('C1', 'C2', 'C3');
title('Capacitances'); grid minor; xlabel('SOC [%]'); ylabel('C [kF]')

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
    C1 = param(5);
    C2 = param(6);
    C3 = param(7);
    tau1 = R1 * C1;
    tau2 = R2 * C2;
    tau3 = R3 * C3;

    % Battery Fixed/Known Parameters 
    global Batt
    Q = Batt.Q*3600;           % Capacity [Amp*s]
    Current = Batt.Current;    % Current {Amp]
    
    %% Run Actual Battery Model
    SOC = [Batt.SOC_exp(1)];
%     OCV = [polyval(Batt.coeff, SOC(1))];
% A = Batt.OCV_SOC(:, 1);
% B = Batt.OCV_SOC(:, 2);
% oSOC = Batt.SOC_exp(1);
% OCV = [interp1(A, B, oSOC)];
OCV = [Batt.Voltage_exp(end)];
    Vt  = [Batt.Voltage_exp(1)];
    Irc1 = [0];
    Irc2 = [0];
    Irc3 = [0];

for k = 2:length(Current)

   oSOC = SOC(k-1);
%    oOCV = polyval(Batt.coeff, oSOC);
% oOCV = interp1(A, B, oSOC);
oOCV = Batt.Voltage_exp(end);
   deltaT = Batt.Time(k) - Batt.Time(k-1);
   
   % Current through RC branchs
    oIrc1 =  ( 1 - exp(-deltaT/tau1) ) * Current(k) + exp(-deltaT/tau1) * Irc1(k-1);
    oIrc2 =  ( 1 - exp(-deltaT/tau2) ) * Current(k) + exp(-deltaT/tau2) * Irc2(k-1);
    oIrc3 =  ( 1 - exp(-deltaT/tau3) ) * Current(k) + exp(-deltaT/tau3) * Irc3(k-1);
    
    % Calculate current Vt & SOC
    oVt   = oOCV + (R0  * Current(k)) + R1 * oIrc1 + R2 * oIrc2 + R3 * oIrc3;
    oSOC  = oSOC + (deltaT / Q)* Current(k);
    
    % Update parameters
    SOC = [SOC; oSOC];
    Vt  = [Vt; oVt];
    OCV = [OCV; oOCV];
    Irc1 = [Irc1; oIrc1];
    Irc2 = [Irc2; oIrc2];
    Irc3 = [Irc2; oIrc3];
    
end

%     % DEBUGGING
%     global Batt
%     Batt.SOC_sim = OCV; 
%     Batt.Voltage_sim = OCV; 
%     Batt.OCV_sim = OCV; 
%     Batt.Irc1 = Irc1;
%     Batt.Irc2 = Irc2;
%     Batt.Irc3 = Irc3;
end