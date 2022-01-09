% Gets SOC-OCV Relationship and Cell Capacity (Q) from C/20 CC-CV test

clear

%% Get data & test parameters
current_folder = cd;
folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\3 - Characterization test 2\';
cd(folder);
addpath(current_folder);
file = append(folder, '01-05-21_10.29 1344_TS001243_C20DisCh.mat');
load(file)

%% Divide data into charge & discharge

iPau = find(meas.Current == 0); % Indices of paused steps (one after discharging & one after charging)

% Find first paused step (after the discharge step)
i = 1;
while iPau(i+1) - iPau(i) == 1
    i = i+1;
end

% Find indices of each step 
% Step 1: Discharge (iDch1, iDch2)
% Step 2: Pause (iPau1, iPau2)
% Step 3: Charge (iChg1, iChg2)
% Step 4: Pause (iPau3, iPau4)
iPau1 = iPau(1); 
iPau2 = iPau(i); 
iPau3 = iPau(i+1); 
iPau4 = iPau(end); 
iDch1 = 1; 
iDch2 = iPau1-1;
iChg1 = iPau2+1; 
iChg2 = iPau3-1;

% Separate data
fields = string(fieldnames(meas));
meas_t = struct2table(meas);

pau1_t = meas_t(iPau1:iPau2, :);
pau1 = table2struct(pau1_t,'ToScalar',true);

dch_t = meas_t(iDch1:iDch2, :);
dch = table2struct(dch_t,'ToScalar',true);

pau2_t = meas_t(iPau3:iPau4, :);
pau2 = table2struct(pau2_t,'ToScalar',true);

chg_t = meas_t(iChg1:iChg2, :);
chg = table2struct(chg_t,'ToScalar',true);

 clearvars -except folder meas pau1 pau2 dch chg

%% Find static parameters
% Force Ah to start at 0 Ah for each step
dch.Ah = dch.Ah - dch.Ah(1);
pau1.Ah = pau1.Ah - pau1.Ah(1);
chg.Ah = chg.Ah - chg.Ah(1);
pau2.Ah = pau2.Ah - pau2.Ah(1);

% Coulombic efficiency at 25C
eta = abs(dch.Ah(end)/chg.Ah(end));   % Coulombic efficiency at 25C
chg.Ah = chg.Ah * eta;                % Adjust charge Ah per eta (not all the charge goes to the battery due to inefficiencies

% Capacity at 25C
Q = 280; 

%% Find OCV-SOC relationship
% SOC during discharging & charging
dch.SOC = 1 - abs(dch.Ah)/Q;
chg.SOC = 0 + chg.Ah/Q;

% Create look up table
SOC = (0:0.005:1)';
OCVdch = interp1([dch.SOC], [dch.Voltage], SOC, 'linear');%,'extrap');
OCVchg = interp1([chg.SOC], [chg.Voltage], SOC, 'linear');%,'extrap');
OCV = (OCVdch+OCVchg)/2;
OCV_SOC = struct('SOC',SOC, 'OCV',OCV);

% Plot SOC-OCV curves
figure; plot(dch.SOC, dch.Voltage);
hold on; plot(chg.SOC, chg.Voltage); 
hold on; plot(OCV_SOC.SOC,OCV_SOC.OCV);
grid minor; title('OCV-SOC (2)'); legend('OCVdch', 'OCVchg', 'OCVavg'); 
xlabel('SOC'); ylabel('OCV (V)'), xlim([0,1]); ylim([2.4,4]); 

% Export OCV-SOC relationship
filename = append(folder, 'OCV_SOC.mat');
save(filename, 'OCV_SOC', 'eta', 'Q', 'dch', 'chg');
