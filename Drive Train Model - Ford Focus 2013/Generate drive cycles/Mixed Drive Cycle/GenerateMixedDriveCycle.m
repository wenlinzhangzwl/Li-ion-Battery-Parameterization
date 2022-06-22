clear
close all

%% Import data
addpath(cd)
folder = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\Scripts\Drive Train Model - Ford Focus 2013\Generate drive cycles\Mixed Drive Cycle\";
addpath(folder)

% data1 = readtable(folder + "SOCScan_CRate2_20140211201611_CustRec.csv");
% time1 = data1.RecordingTime;
% current1 = data1.PEC_Measured_Current;
% voltage1 = data1.PEC_Measured_Voltage;
% figure; plot(time1, current1)

data2 = readtable(folder + "SOCScan_DOD2_20140211124413_CustRec.csv");
time2 = data2.RecordingTime;
current2 = data2.PEC_Measured_Current;
voltage2 = data2.PEC_Measured_Voltage;
figure; plot(time2, current2)

%% Single cycle
cyclelength = 7184.9-1379.1;
ind = find(time2 >= cyclelength, 1);
time = time2(1:ind);
current = current2(1:ind);
figure; plot(time, current);

%% Scale current for EVE280Ah cell at 1C

% Battery specs
batt = "EVE280";
V_nom = 3.2; 
deltaT = 0.1;  
Q = 280; 
SOCRange = 1-0;

% Scale test
maxC = "0p8";
maxI = 0.8 * Q;
ratio = Q/max(abs(current));
current = current * ratio;
Ah = -sum(current) * deltaT / 3600; 
figure; plot(time, current); grid on; xlabel("Time (s)"); ylabel("Current (A)"); title("Current Profile of the Mixed Drive Cycle")

% % Power profile & total power
% power = current * V_nom;            % power profile
% power_total = sum(power)*deltaT;    % total power consumed by the drive cycle [W] 

% Length of drive cycle
n =  Q/(Ah*SOCRange);
n = fix(n) + 1;                       % number of cycle needed to cover SOCRange
test_length = cyclelength * n/3600;   % total length of the test

%% Scale current for Samsung 35E at 1C
% 
% % Battery specs
% batt = "Sam35E";
% Q = 3.4;        % Capacity [Ah]
% V_nom = 3.6;    % Nominal voltage [V]
% deltaT = 0.1;   % Sampling time [s]
% SOCRange = 1-0; % SOC range to cover
% 
% CRate = Q * 3; 
% k = CRate/max(abs(current));
% current = current * k;
% figure; plot(time, current/Q); xlabel("Time"); ylabel("SOC")
% figure; plot(time, current); xlabel("Time"); ylabel("Current [A]"); grid on
% 
% % Power profile & total power
% power = current * V_nom;            % power profile
% power_total = sum(power)*deltaT;    % total power consumed by the drive cycle [W] 
% 
% % Length of drive cycle
% power_batt = Q * 3600 * V_nom;
% n =  (power_batt*SOCRange)/ (-power_total);
% n = fix(n) + 1;                       % number of cycle needed to cover SOCRange
% test_length = cyclelength * n/3600;   % total length of the test

%% Export to Digatron command (-ve = discharge, +ve = charge)
% command = string(deltaT) + ' sec;;' + string(power) +  ';;';
command = string(deltaT) + ' sec;' + string(current) +  ';;;';
inc = fix(height(command)/5);
validateattributes(inc, {'double'}, {'integer'})

command1 = command(1:inc);
command2 = command(inc+1:inc*2);
command3 = command(inc*2+1:inc*3);
command4 = command(inc*3+1:inc*4);
command5 = command(inc*4+1:end);

folder = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\11 - Drive Cycle (0p8C)\Input Power Profiles\";
filename = folder + batt + "_Mixed_" + maxC + "C_repeat" + string(n) + "_totalLength" + string(test_length);

% writematrix(command1, filename + "_part1" + ".txt");
% writematrix(command2, filename + "_part2" + ".txt");
% writematrix(command3, filename + "_part3" + ".txt");
% writematrix(command4, filename + "_part4" + ".txt");
% writematrix(command5, filename + "_part5" + ".txt");