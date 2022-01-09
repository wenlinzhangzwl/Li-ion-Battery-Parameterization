% Gets SOC-OCV Relationship and Cell Capacity (Q) from C/20 CC-CV test
clear

%% Get data & test parameters
current_folder = cd;
% folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\Samsung INR18650-35E Cell\Testing\Cell2 3.5 Ah\';
folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\'; 
addpath(folder); 
addpath(current_folder);

cd(folder);
filename = 'C20CCCV';
file = append(folder, filename, '.mat');
load(file)

% % Samsung INR18650-35E
% iChg = 14; 
% iDch = 16; 
% temp = 23; 
% Q = 3.1043;
% Vmax = 4.2;
% Vmin = 2.65; 

% EVE 280
iChg = 15; 
iDch = 13; 
temp = 25;
Q = 277.2393; 
Vmax = 3.65; 
Vmin = 2.5; 

%% Divide data into charge & discharge

ind_chg = find(meas.Step == iChg);
ind_dch = find(meas.Step == iDch);

meas_t = struct2table(meas);
chg_t = meas_t(ind_chg, :);
dch_t = meas_t(ind_dch, :);

% % Remove charging data after Vmax ia reached
% ind = find(chg_t.Voltage >= Vmax, 2); 
% chg_t = chg_t(1:ind(2));

% Force Ah to be 0 at the beginning
chg_t.Ah = chg_t.Ah - chg_t.Ah(1); 
dch_t.Ah = dch_t.Ah - dch_t.Ah(1); 

% Calculate SOC during discharging & charging
% Q = max(abs(dch_t.Ah)); 
dch_t.SOC = 1 - (-dch_t.Ah)/Q;
chg_t.SOC = 0 + chg_t.Ah/Q;

% Remove repeating SOC values
[~,ia, ~] = unique(chg_t.SOC);
chg_t = chg_t(ia,:);
[~,ia, ~] = unique(dch_t.SOC);
dch_t = dch_t(ia,:);

% % Plot for validation
% figure; 
% subplot(3, 1, 1); hold on; grid on
% plot(chg_t.Time, chg_t.Voltage, '.'); 
% plot(dch_t.Time, dch_t.Voltage, '.'); 
% subplot(3, 1, 2); hold on; grid on
% plot(chg_t.Time, chg_t.Current, '.'); 
% plot(dch_t.Time, dch_t.Current, '.'); 
% subplot(3, 1, 3); hold on; grid on
% plot(chg_t.Time, chg_t.SOC, '.'); 
% plot(dch_t.Time, dch_t.SOC, '.'); 

%% Calculating resistance from instaneous voltage drop
% meas_t = struct2table(meas);
% 
% % Find cell resistance 
% % (dividing by 5 because sampling time is too long, doesn't reflect instaneous voltage)
% R0_dch1 = (pau2.Voltage(1)-dch.Voltage(end)) / 5 / (-dch.Current(end));  % 0% SOC discharge
% R0_dch2 = (pau1.Voltage(end)-dch.Voltage(1)) / 5 / (-dch.Current(1));    % 100% SOC discharge
% R0_chg1 = (chg.Voltage(1)-pau2.Voltage(end)) / 5 / chg.Current(1);     % 0% SOC charge
% R0_chg2 = (chg.Voltage(end)-pau3.Voltage(1)) / 5 / chg.Current(end);   % 100% SOC charge
% 
% % Find ressitance at each timestamp assuming it changes linearly from 0% to 100% SOC
% dch.R = (R0_dch2 - R0_dch1) * dch.SOC + R0_dch1;
% chg.R = (R0_chg2 - R0_chg1) * chg.SOC + R0_chg1;
% 
% % Adjust voltages by removing R0*i
% dch.Vadj = dch.Voltage - dch.R .* dch.Current; 
% chg.Vadj = chg.Voltage - chg.R .* chg.Current;

%% Calculate average OCV

% Approximate average OCV
SOC = (0:0.0001:1)';
OCVdch = interp1([dch_t.SOC], [dch_t.Voltage], SOC, 'linear');
OCVchg = interp1([chg_t.SOC], [chg_t.Voltage], SOC, 'linear');
OCVavg = (OCVdch+OCVchg)/2;

% Resample n data points and the beginning & end of OCVavg (n+2 points in total)
OCV_SOC = rmmissing([SOC, OCVavg]);
SOCfit = OCV_SOC(:, 1);
OCVfit = OCV_SOC(:, 2);

OCV_SOC_dch = rmmissing([SOC, OCVdch]);
SOCfit_dch = OCV_SOC_dch(:, 1);
OCVfit_dch = OCV_SOC_dch(:, 2);

OCV_SOC_chg = rmmissing([SOC, OCVchg]);
SOCfit_chg = OCV_SOC_chg(:, 1);
OCVfit_chg = OCV_SOC_chg(:, 2);

%% Determine optimal polynomial degree
% Polynomial coefficients returned as variable 'fit'
deg = 7;
for i = 1:10
    fit(i).deg = deg;
    [fit(i).p, fit(i).s] = polyfit(SOCfit, OCVfit, deg); 
    [fit(i).OCVpoly] = polyval(fit(i).p, SOC, fit(i).s);
    fit(i).err = OCVavg - fit(i).OCVpoly;
    err = rmmissing(fit(i).err);
    fit(i).RMSE = sqrt(mean(err.^2));
    deg = deg + 1; 
end

figure; hold on
polydeg = [];
deg = deg - length(fit);
for i = 1:length(fit)
    plot(deg, fit(i).RMSE*1000, '*');
    polydeg = [polydeg fit(i).deg];
    deg = deg + 1; 
end
polydeg = string(polydeg);
polydeg = 'deg ' + polydeg;
% legend(polydeg(1), polydeg(2), polydeg(3), polydeg(4), polydeg(5), polydeg(6), polydeg(7), polydeg(8), polydeg(9), polydeg(10));
grid on
title('RMSE Error'); xlabel('Polynomial Degree'); ylabel('Polynomial Fit RMSE Error [mV]');

%% Fit polynomials & plot error
deg = 21;

% Fit polynomial to OCV curves
[coeff, s] = polyfit(SOCfit, OCVfit, deg); 
[coeff_dch, s_dch] = polyfit(SOCfit_dch, OCVfit_dch, deg); 
[coeff_chg, s_chg] = polyfit(SOCfit_chg, OCVfit_chg, deg);

[OCVpoly] = polyval(coeff, SOC, s);
[OCVpoly_dch] = polyval(coeff_dch, SOC, s_dch);
[OCVpoly_chg] = polyval(coeff_chg, SOC, s_chg);

err = OCVavg - OCVpoly;
err2 = rmmissing(err);
RMSE = sqrt(mean(err2.^2));

% Plot polynomial fit
figure; hold on
plot(SOC, OCVavg,...    % original
    SOC, OCVpoly);      % fitted polynomial        
% plot(SOCfit, OCVfit, '*');          % data points
legend('exp', 'fit (' + string(deg) + 'th deg)'); grid on
xlabel('SOC'); ylabel('Voltage [V]'); title('OCV-SOC Curve')

% Plot polynomial error
figure;
plot(SOC, err*1000);
grid minor; title('Error - (' + string(deg) + 'th deg)'); xlabel('SOC'); ylabel('Error [mV]');

%% Plot and export OCV-SOC curve

% % Plot SOC-OCV curves
% figure; hold on;
% plot(dch_t.SOC, dch_t.Voltage);
% plot(chg_t.SOC, chg_t.Voltage); 
% plot(SOC,OCVavg);
% % plot(SOC, OCVpoly);
% grid on; title('OCV-SOC'); legend('OCVdch', 'OCVchg', 'OCVavg'); 
% % grid on; title('OCV-SOC'); legend('OCVdch', 'OCVchg', 'OCVavg', 'OCVpolyfit'); 
% xlabel('SOC'); ylabel('OCV (V)'), %xlim([0,1]); ylim([2.4,4]);

% Export OCV-SOC relationship
filename = append(folder, 'Results\', filename, '_OCVSOC', string(deg), '.mat');
save(filename, 'OCV_SOC', 'coeff');
