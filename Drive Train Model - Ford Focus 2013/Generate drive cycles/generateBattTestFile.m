clear
% close all

%% Load drive cycles

% Load drive cycles & down sample to 0.1 sec
deltaT = 0.1;
n_sample = deltaT/1e-3;

load UDDS
UDDS_current = drivecycle_current(1:n_sample:end, :);
UDDS_power = drivecycle_power(1:n_sample:end, :);
UDDS_CRate = drivecycle_CRate(1:n_sample:end, :);
load HWFET
HWFET_current = drivecycle_current(1:n_sample:end, :);
HWFET_power = drivecycle_power(1:n_sample:end, :);
HWFET_CRate = drivecycle_CRate(1:n_sample:end, :);
load US06
US06_current = drivecycle_current(1:n_sample:end, :);
US06_power = drivecycle_power(1:n_sample:end, :);
US06_CRate = drivecycle_CRate(1:n_sample:end, :);
load NEDC
NEDC_current = drivecycle_current(1:n_sample:end, :);
NEDC_power = drivecycle_power(1:n_sample:end, :);
NEDC_CRate = drivecycle_CRate(1:n_sample:end, :);

%% Only take negative (discharge) currents
% for i = 1:height(UDDS_current)
%     if UDDS_current(i, 2) <= 0
%         UDDS_current_neg(i, 1) = UDDS_current(i, 2);
%         UDDS_power_neg(i, 1) = UDDS_power(i, 2);
%         UDDS_CRate_neg(i, 1) = UDDS_CRate(i, 2);
%     end
% end
% UDDS_current_neg = [UDDS_current(:, 1) UDDS_current_neg]; 
% UDDS_power_neg = [UDDS_power(:, 1) UDDS_power_neg]; 
% UDDS_CRate_neg = [UDDS_CRate(:, 1) UDDS_CRate_neg];
% 
% for i = 1:height(HWFET_current)
%     if HWFET_current(i, 2) <= 0
%         HWFET_current_neg(i, 1) = HWFET_current(i, 2);
%         HWFET_power_neg(i, 1) = HWFET_power(i, 2);
%         HWFET_CRate_neg(i, 1) = HWFET_CRate(i, 2);
%     end
% end
% HWFET_current_neg = [HWFET_current(:, 1) HWFET_current_neg]; 
% HWFET_power_neg = [HWFET_power(:, 1) HWFET_power_neg]; 
% HWFET_CRate_neg = [HWFET_CRate(:, 1) HWFET_CRate_neg];
% 
% for i = 1:height(US06_current)
%     if US06_current(i, 2) <= 0
%         US06_current_neg(i, 1) = US06_current(i, 2);
%         US06_power_neg(i, 1) = US06_power(i, 2);
%         US06_CRate_neg(i, 1) = US06_CRate(i, 2);
%     end
% end
% US06_current_neg = [US06_current(:, 1) US06_current_neg]; 
% US06_power_neg = [US06_power(:, 1) US06_power_neg]; 
% US06_CRate_neg = [US06_CRate(:, 1) US06_CRate_neg];

%% Scale drive cycles
% maxC = 0.9;
% 
% % UDDS
% maxC_UDDS = max(abs(UDDS_CRate(:, 2))); 
% scaleFactor = maxC / maxC_UDDS; 
% UDDS_current = UDDS_current * scaleFactor; 
% UDDS_power = UDDS_power * scaleFactor; 
% UDDS_CRate = UDDS_CRate * scaleFactor; 

%% Plot drive cycles
% CRate
figure; 
subplot(4, 1, 1); plot(UDDS_CRate(:, 1), UDDS_CRate(:, 2)); title('UDDS CRate'); grid on
subplot(4, 1, 2); plot(HWFET_CRate(:, 1), HWFET_CRate(:, 2)); title('HWFET CRate'); grid on
subplot(4, 1, 3); plot(US06_CRate(:, 1), US06_CRate(:, 2)); title('US06 CRate'); grid on
subplot(4, 1, 4); plot(NEDC_CRate(:, 1), NEDC_CRate(:, 2)); title('NEDC CRate'); grid on

% % CRate_neg
% figure; plot(UDDS_CRate_neg(:, 1), UDDS_CRate_neg(:, 2)); title('UDDS CRate (neg)'); grid on
% figure; plot(HWFET_CRate_neg(:, 1), HWFET_CRate_neg(:, 2)); title('HWFET CRate (neg)'); grid on
% figure; plot(US06_CRate_neg(:, 1), US06_CRate_neg(:, 2)); title('US06 CRate (neg)'); grid on
% 
% % power
% figure; 
% subplot(4, 1, 1); plot(UDDS_power(:, 1), UDDS_power(:, 2)); title('UDDS power'); grid on
% subplot(4, 1, 2); plot(HWFET_power(:, 1), HWFET_power(:, 2)); title('HWFET power'); grid on
% subplot(4, 1, 3); plot(US06_power(:, 1), US06_power(:, 2)); title('US06 power'); grid on
% subplot(4, 1, 4); plot(NEDC_power(:, 1), NEDC_power(:, 2)); title('NEDC power'); grid on
% 
% power_neg
% figure; plot(UDDS_power_neg(:, 1), UDDS_power_neg(:, 2)); title('UDDS power (neg)'); grid on
% figure; plot(HWFET_power_neg(:, 1), HWFET_power_neg(:, 2)); title('HWFET power (neg)'); grid on
% figure; plot(US06_power_neg(:, 1), US06_power_neg(:, 2)); title('US06 power (neg)'); grid on

%% Select segment for C rate testing
% Define segment
% Case = 6;
% switch Case
%     case 1 % UDDS Low C rate 
%         cycle_current = UDDS_current_neg;
%         cycle_power = UDDS_power_neg;
%         cycle_CRate = UDDS_CRate_neg;
%         
%         startTime = 472;
%         endTime = 491.1;
%         
%         iStart = find(cycle_power(:, 1) == startTime, 1); 
%         iEnd = find(cycle_power(:, 1) == endTime, 1); 
%         
%         seg_current = cycle_current(iStart:iEnd, :);
%         seg_power = cycle_power(iStart:iEnd, :);
%         seg_CRate = cycle_CRate(iStart:iEnd, :);
%         
%         drivecycle_name = 'UDDSLowC';
%         
%     case 2 % UDDS Mid C rate 
%         cycle_current = UDDS_current_neg;
%         cycle_power = UDDS_power_neg;
%         cycle_CRate = UDDS_CRate_neg;
%         
%         startTime = 1267;
%         endTime = 1280.2;
%         
%         iStart = find(cycle_power(:, 1) == startTime, 1); 
%         iEnd = find(cycle_power(:, 1) == endTime, 1); 
%         
%         seg_current = cycle_current(iStart:iEnd, :);
%         seg_power = cycle_power(iStart:iEnd, :);
%         seg_CRate = cycle_CRate(iStart:iEnd, :);
% 
%         drivecycle_name = 'UDDSMidC';
%         
%     case 3 % UDDS High C rate 
%         cycle_current = UDDS_current_neg;
%         cycle_power = UDDS_power_neg;
%         cycle_CRate = UDDS_CRate_neg;
%         
%         startTime = 187;
%         endTime = 207;
%         
%         iStart = find(cycle_power(:, 1) == startTime, 1); 
%         iEnd = find(cycle_power(:, 1) == endTime, 1); 
%         
%         seg_current = cycle_current(iStart:iEnd, :);
%         seg_power = cycle_power(iStart:iEnd, :);
%         seg_CRate = cycle_CRate(iStart:iEnd, :);
% 
%         
%         drivecycle_name = 'UDDSHighC';
%     
%     case 4 % US06 Low C rate 
%         cycle_current = US06_current_neg;
%         cycle_power = US06_power_neg;
%         cycle_CRate = US06_CRate_neg;
%         
%         startTime = 425;
%         endTime = 469.1;
%         
%         iStart = find(cycle_power(:, 1) == startTime, 1); 
%         iEnd = find(cycle_power(:, 1) == endTime, 1); 
%         
%         seg_current = cycle_current(iStart:iEnd, :);
%         seg_power = cycle_power(iStart:iEnd, :);
%         seg_CRate = cycle_CRate(iStart:iEnd, :);
%         
%         drivecycle_name = 'US06LowC';
%         
%     case 5 % US06 Mid C rate 
%         cycle_current = US06_current_neg;
%         cycle_power = US06_power_neg;
%         cycle_CRate = US06_CRate_neg;
%         
%         startTime = 49;
%         endTime = 70;
%         
%         iStart = find(cycle_power(:, 1) == startTime, 1); 
%         iEnd = find(cycle_power(:, 1) == endTime, 1); 
%         
%         seg_current = cycle_current(iStart:iEnd, :);
%         seg_power = cycle_power(iStart:iEnd, :);
%         seg_CRate = cycle_CRate(iStart:iEnd, :);
% 
%         drivecycle_name = 'US06MidC';
%         
%     case 6 % US06 High C rate 
%         cycle_current = US06_current_neg;
%         cycle_power = US06_power_neg;
%         cycle_CRate = US06_CRate_neg;
%         
%         startTime = 568;
%         endTime = 579.5;
%         
%         iStart = find(cycle_power(:, 1) == startTime, 1); 
%         iEnd = find(cycle_power(:, 1) == endTime, 1); 
%         
%         seg_current = cycle_current(iStart:iEnd, :);
%         seg_power = cycle_power(iStart:iEnd, :);
%         seg_CRate = cycle_CRate(iStart:iEnd, :);
% 
%         
%         drivecycle_name = 'US06HighC';
% end

% Plot segment
% figure; 
% hold on;
% plot(cycle_CRate(:, 1), cycle_CRate(:, 2));
% plot(seg_CRate(:, 1), seg_CRate(:, 2));
% hold off; 
% title('CRate'); grid on
% figure; 
% hold on;
% plot(cycle_power(:, 1), cycle_power(:, 2));
% plot(seg_power(:, 1), seg_power(:, 2));
% hold off; 
% title('power'); grid on
% 
% % Repeat segment to cover desired SOC range
% SOCRange = 1-0.85;
% seg_totalPower = sum(seg_power(:, 2))*deltaT;         % total power consumed by the drive cycle [W] 
% seg_length = seg_power(end, 1) - seg_power(1, 1);     % length of the drive cycle [s]
% Vnom = 3.2;                                           % nominal voltage of the battery [V]
% batt_power = Q * 3600 * Vnom;                         % total power of the battery [W*s]
% seg_n =  (batt_power*SOCRange)/ (-seg_totalPower);
% seg_n = fix(seg_n) + 1;                               % number of segments needed
% 
% % Construct profile (test_powerProfile) for testing
% rest = zeros(5/deltaT, 1);                                  % rest for 5 seconds
% test_powerProfile = [seg_power(:, 2); rest];
% test_length = (length(test_powerProfile) - 1) * deltaT;     % length of a single segment with rest [s]
% test_time = transpose(0:deltaT:test_length); 
% % figure; plot(test_time, test_powerProfile); grid on;
% 
% % Total length of the test
% test_totalLength = test_length*seg_n;                       % total length of the test [s]
% test_totalLength_hr = test_totalLength/3600;                % total length of the test [h]
% 
% % Write command for Digatron
% power = test_powerProfile;
% command = string(deltaT) + ' sec;;' + string(power) +  ';;';
% folder = "C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\5 - Drive cycle\";
% filename = folder + "EVE280_" + string(drivecycle_name) + "_repeat" + string(seg_n) + "_totalLength" + string(test_totalLength_hr)+ ".txt";
% writematrix(command, filename);

%% Drive cycles

% Select drive cycle
drivecycle = 4;
switch drivecycle
    case 1 % UDDS
        drivecycle_current = UDDS_current;
        drivecycle_power = UDDS_power; 
        drivecycle_name = 'UDDS';
    case 2 % HWFET
        drivecycle_current = HWFET_current;
        drivecycle_power = HWFET_power; 
        drivecycle_name = 'HWFET';
    case 3 % US06
        drivecycle_current = US06_current;
        drivecycle_power = US06_power; 
        drivecycle_name = 'US06';
    case 4 % NEDC 
        drivecycle_current = NEDC_current;
        drivecycle_power = NEDC_power; 
        drivecycle_name = 'NEDC';
end
test_powerProfile = drivecycle_power(:, 2);

% Repeat drive cylce to cover entire SOC range
SOCRange = 1-0;
drivecycle_totalPower = sum(drivecycle_power(:, 2))*deltaT;             % total power consumed by the drive cycle [W] 
drivecycle_length = drivecycle_power(end, 1) - drivecycle_power(1, 1);  % length of the drive cycle [s]
Vnom = 3.2;                                                             % nominal voltage of the battery [V]
batt_power = Q * 3600 * Vnom;                                           % total power of the battery [W*s]
drivecycle_n =  (batt_power*SOCRange)/ (-drivecycle_totalPower);
drivecycle_n = fix(drivecycle_n) + 1;                                   % number of cycle needed

% Length of a drive cycle & total length of the test
test_length = (length(test_powerProfile) - 1) * deltaT;     % length of a single segment with rest [s]
test_time = transpose(0:deltaT:test_length); 
test_totalLength = test_length*drivecycle_n;                % total length of the test [s]
test_totalLength_hr = test_totalLength/3600;                % total length of the test [h]

% Export to Digatron command
command = string(deltaT) + ' sec;;' + string(test_powerProfile) +  ';;';
folder = "C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\6 - Drive cycle\Power Profiles\";
filename = folder + "EVE280_" + string(drivecycle_name) + "_repeat" + string(drivecycle_n) + "_totalLength" + string(test_totalLength_hr)+ ".txt";
% writematrix(command, filename);
