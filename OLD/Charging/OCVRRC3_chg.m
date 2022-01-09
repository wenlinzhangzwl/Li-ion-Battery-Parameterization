clearvars -except meas

%% Load and process data from the experiment

% Load files
currentFolder = cd; 
addpath(currentFolder);
folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\4 - Characterization test 3\';
cd(folder);
% load('PUL25chg.mat');           % Raw data
load('C20CCCV_OCVSOC.mat');     % SOC OCV relationship
% load('PUL25dch_OptimResult');   % GA optimization results

% Trim data so it starts at desired point
stepStart = 12;
stepEnd = 21;
meas_t = struct2table(meas);
meas_t = meas_t(and(meas.Step >= stepStart, meas.Step <= stepEnd), :);

% Cell parameters
Q = 279.29;
temp = 25; 
meas_t.Ah = meas_t.Ah - meas_t.Ah(1);
meas_t.SOC = meas_t.Ah/Q;

% Resample data from 10Hz to 1Hz
measResampled = meas_t(1:10:end, :);

%% Run battery model
global Batt;
Batt = struct('Time', measResampled.Time, ...
              'Current', measResampled.Current,...
              'SOC_exp', measResampled.SOC,...
              'Voltage_exp', measResampled.Voltage,...
              'Q', Q, ...
              'coeff', coeff,...
              'OCV_SOC', OCV_SOC,...
              'OptimResult', OptimResult);
[Batt.Voltage_sim, Batt.SOC_sim] = OCVRRCModel(Batt);

%% Post processing - SOC
% Error
Batt.errSOC = Batt.SOC_exp - Batt.SOC_sim;
Batt.RMSE_SOC = sqrt(mean(Batt.errSOC.^2));

% Plot SOC_exp & SOC_sim
figure
plot(Batt.Time/3600, Batt.SOC_sim * 100);
hold on; 
plot(Batt.Time/3600, Batt.SOC_exp * 100);
xlabel('time[h]'); ylabel('SOC [%]'); title('SOC - chg'); grid minor
legend('sim', 'exp')

%% Post processing - Voltage
% Error
Batt.errV = Batt.Voltage_exp - Batt.Voltage_sim;
errV = Batt.errV( ~any( isnan( Batt.errV ) | isinf( Batt.errV ), 2 ),: );
Batt.RMSE_V = sqrt(mean(rmmissing(Batt.errV.^2)));

% Plot Voltage_sim & Voltage_exp 
figure
hold on
plot(Batt.Time/3600, Batt.Voltage_exp, 'r');
plot(Batt.Time/3600, Batt.Voltage_sim, 'b');
xlabel('Time [h]'); ylabel('Terminal Voltage [V]'); title('Voltage - chg'); grid minor
legend('exp', 'sim')

% Plot Voltage error
figure
plot(Batt.Time/3600, Batt.errV*1000);
xlabel('Time [h]'); ylabel('Error [mV]'); title('Voltage Error - chg'); grid minor

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
title('Resistances - chg'); grid minor; xlabel('SOC [%]'); ylabel('R [mOhm]')

% Plot capacitances
figure; 
hold on
plot(OptimResult.SOC*100, C1/1000);
plot(OptimResult.SOC*100, C2/1000);
plot(OptimResult.SOC*100, C3/1000);
hold off
legend('C1', 'C2', 'C3');
title('Capacitances - chg'); grid minor; xlabel('SOC [%]'); ylabel('C [kF]')

%% Save simulation results
% filename = append(folder, 'OCVRRC3_PUL25chg.mat');
% save(filename, 'Batt', 'OptimResult');

%% Debugging
% figure
% hold on
% plot(Batt.Time/3600, Batt.Voltage_exp);
% plot(Batt.Time/3600, Batt.OCV_sim);
% hold off
% legend('exp', 'sim')

%% Function
function [Vt, SOC] = OCVRRCModel(Batt)

    % Battery Fixed/Known Parameters 
    Q = Batt.Q*3600;           % Capacity [Amp*s]
    Current = Batt.Current;    % Current {A]

    SOC = [Batt.SOC_exp(1)];
    Vt  = [polyval(Batt.coeff, SOC(1))];
%     OCV = [polyval(Batt.coeff, SOC(1))];
A = Batt.OCV_SOC(:, 1);
B = Batt.OCV_SOC(:, 2);
oSOC = Batt.SOC_exp(1);
OCV = [pchip(A, B, oSOC)];
    Irc1 = [0];
    Irc2 = [0];
    Irc3 = [0];

    for k = 2:length(Current)
        % Calculate current step parameters
        oSOC = SOC(k-1);
%         oOCV = polyval(Batt.coeff, oSOC);
oOCV = pchip(A, B, oSOC);
        deltaT = Batt.Time(k) - Batt.Time(k-1);

        % Avoids extrapolation
        param_t = table2array(Batt.OptimResult(:, 'results'));
        SOC_t = table2array(Batt.OptimResult(:, 'SOC'));
        if oSOC > max(SOC_t)
            oSOC2 = max(SOC_t);
        elseif oSOC < min(SOC_t)
            oSOC2 = min(SOC_t);
        else
            oSOC2 = oSOC;
        end
        
        % Interpolate for model parameters based on SOC
        R0 = interp1(SOC_t, param_t(:, 1), oSOC2, 'linear');
        R1 = interp1(SOC_t, param_t(:, 2), oSOC2, 'linear');
        R2 = interp1(SOC_t, param_t(:, 3), oSOC2, 'linear');
        R3 = interp1(SOC_t, param_t(:, 4), oSOC2, 'linear');
        C1 = interp1(SOC_t, param_t(:, 5), oSOC2, 'linear');
        C2 = interp1(SOC_t, param_t(:, 6), oSOC2, 'linear');
        C3 = interp1(SOC_t, param_t(:, 1), oSOC2, 'linear');
%         results = [8e-4, 9e-5, 9e-5, 9e-5, 8e4, 8e4, 8e4];
%         R0 = results(1);
%         R1 = results(2);
%         R2 = results(3);
%         R3 = results(4);
%         C1 = results(5);
%         C2 = results(6);
%         C3 = results(7);

        tau1 = R1 * C1;
        tau2 = R2 * C2;
        tau3 = R3 * C3;

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
    
    % Debugging
    global Batt
    Batt.OCV_sim = OCV; 
    Batt.Irc1 = Irc1;
    Batt.Irc2 = Irc2;
    Batt.Irc3 = Irc3;
    beep
end
