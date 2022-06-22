% Gets SOC-OCV Relationship and Cell Capacity (Q) from C/20 CC-CV test

clear
close all

saveCoeff = 0;  %!!!!!!!!!!!!!!!!!!!!CHANGE TO SAVE/NOT SAVE 

folder_current = cd;
folder_data = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\4 - Characterization test 3\"; 
folder_results = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\Results\OCV Polynomials\";
addpath(folder_current, folder_data, folder_results); 

% Test info
temp = 25; 
cycler = "Digatron";
cell = "EVE280";
iDch = 13;  % Set iChg = 1 for D&V data
iChg = 15;  % Set iDch = 0 for D&V data
if cell == "Samsung_INR18650_35E"
    % Samsung INR18650-35E
    Vmax = 4.2;
    Vmin = 2.65;
    Q_nominal = 3.4; % Alghough minimum Q specificed by spec sheet is 3.35Ah
elseif cell == "EVE280"
    % EVE 280
    Vmax = 3.65; 
    Vmin = 2.5; 
    Q_nominal = 280; 
end

% Divide data into discharge & charge
if cycler == "Digatron"
% For data from Digatron

    % Load data
    filename = "OCV_25degC_0p05C.mat";
    load(filename, "meas");
    meas.Current = -meas.Current;  % Convert to the convension: -ve = charge, +ve = discharge
    meas.Ah = -meas.Ah; 

    % Divide data into chg_t & dch_t
    ind_chg = find(meas.Step == iChg);
    ind_dch = find(meas.Step == iDch);
    
    try
        meas = struct2table(meas); % Convert to table if it's not
    catch
        % Do nothing
    end
    meas_chg = meas(ind_chg, :);
    meas_dch = meas(ind_dch, :);

elseif cycler == "DnV"
    % For data from D&V:
    deltaT = 1; 
    
    filename = "OCVdch25.mat";
    load(folder_data + filename);
    meas_dch = struct2table(meas);
    meas_dch.Step = zeros(height(meas_dch), 1); 
    
    filename = "OCVchg25.mat";
    load(folder_data + filename);
    meas_chg = struct2table(meas);
    meas_chg.Time = meas_chg.Time + meas_dch.Time(end) + deltaT;
    meas_chg.Step = ones(height(meas_chg), 1); 

end

% Force Ah to be 0 at the beginning
meas_chg.Ah = meas_chg.Ah - meas_chg.Ah(1); 
meas_dch.Ah = meas_dch.Ah - meas_dch.Ah(1); 

% Calculate SOC during discharging & charging
Q_dch = abs(max(meas_dch.Ah)-min(meas_dch.Ah)); 
meas_dch.SOC = 1 - meas_dch.Ah/Q_dch;
meas_chg.SOC = 0 - meas_chg.Ah/Q_dch;

% Delete meas_chg data after SOC = 1 is reached
ind = find(meas_chg.SOC>1);
meas_chg(ind, :) = []; 

% Remove repeating SOC values
[~,ia, ~] = unique(meas_chg.SOC);
meas_chg = meas_chg(ia,:);
[~,ia, ~] = unique(meas_dch.SOC);
meas_dch = meas_dch(ia,:);

% Plot for validation
figure("WindowStyle", "docked"); 
subplot(3, 1, 1); hold on; grid on
plot(meas_chg.Time, meas_chg.Voltage, '.'); 
plot(meas_dch.Time, meas_dch.Voltage, '.'); 
title("Voltage"); xlabel("Time [s]"); ylabel("Voltage [V]")
subplot(3, 1, 2); hold on; grid on
plot(meas_chg.Time, meas_chg.Current, '.'); 
plot(meas_dch.Time, meas_dch.Current, '.'); 
title("Current"); xlabel("Time [s]"); ylabel("Current [A]")
subplot(3, 1, 3); hold on; grid on
plot(meas_chg.Time, meas_chg.SOC, '.'); 
plot(meas_dch.Time, meas_dch.SOC, '.'); 
title("SOC"); xlabel("Time [s]"); ylabel("SOC")

%% Calculate average OCV

% SOC where the OCVs are interpolated
SOC_min = max(min(meas_chg.SOC), min(meas_dch.SOC));
SOC_max = min(max(meas_chg.SOC), max(meas_dch.SOC));
SOC = (SOC_min:0.001:SOC_max)';

% Approximate average OCV
OCVdch = interp1([meas_dch.SOC], [meas_dch.Voltage], SOC, 'linear');
OCVchg = interp1([meas_chg.SOC], [meas_chg.Voltage], SOC, 'linear');
OCVavg = (OCVdch+OCVchg)/2;

% Plot for validation
figure("WindowStyle", "docked"); hold on
plot(SOC, OCVdch, '.-')
plot(SOC, OCVchg, '.-')
plot(SOC, OCVavg, '.-')
legend("dch", "chg", "avg"); grid on; xlabel("SOC"); ylabel("Voltage")

%% Determine optimal polynomial degree
% Polynomial coefficients returned as variable 'fit'
deg = 5;
for i = 1:17
    fit(i).deg = deg;
    [fit(i).p, fit(i).s] = polyfit(SOC, OCVavg, deg); 
    [fit(i).OCVpoly] = polyval(fit(i).p, SOC, fit(i).s);
    fit(i).err = OCVavg - fit(i).OCVpoly;
    err = rmmissing(fit(i).err);
    fit(i).RMSE = sqrt(mean(err.^2));
    deg = deg + 1; 
end
fit = struct2table(fit); 

figure("WindowStyle", "docked"); hold on
plot(fit.deg, fit.RMSE*1000, '*-')
grid on
title('RMSE of Fitted Polynomial'); xlabel('Polynomial Degree'); ylabel('RMSE [mV]');

%% Fit polynomials & plot error
deg = 5;

for i = 1:17
    % Fit polynomial to OCV curves
    [coeff, s] = polyfit(SOC, OCVavg, deg); 
    [coeff_dch, s_dch] = polyfit(SOC, OCVdch, deg); 
    [coeff_chg, s_chg] = polyfit(SOC, OCVchg, deg);
    
    [OCVpoly] = polyval(coeff, SOC, s);
    [OCVpoly_dch] = polyval(coeff_dch, SOC, s_dch);
    [OCVpoly_chg] = polyval(coeff_chg, SOC, s_chg);
    
    err = OCVavg - OCVpoly;
    err2 = rmmissing(err);
    RMSE = sqrt(mean(err2.^2));
    
    % Plot polynomial fit
    figure("WindowStyle", "docked");
    ax1 = subplot(2, 1, 1); hold on
    plot(SOC, OCVavg);     % original
    plot(SOC, OCVpoly);    % fitted polynomial        
    legend('exp', 'fit (' + string(deg) + 'th deg)'); grid on
    xlabel('SOC'); ylabel('Voltage [V]'); title('OCV-SOC Curve')
    ax2 = subplot(2, 1, 2);
    plot(SOC, err*1000);
    grid minor; title('Error - (' + string(deg) + 'th deg)'); xlabel('SOC'); ylabel('Error [mV]');
    
    % Export OCV-SOC relationship
    if saveCoeff == 1
        filename = folder_results + "OCVpoly_" + string(deg) + "_" + string(temp) + "degC.mat";
        save(filename, "coeff");
    end

    deg = deg + 1; 
end