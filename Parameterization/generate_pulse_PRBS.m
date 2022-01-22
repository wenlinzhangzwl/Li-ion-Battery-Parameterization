clearvars -except meas
% close all

currentFolder = cd; 
addpath(currentFolder)
deltaT = 0.1; % Sampling time for all data
PRBSfreq = "10";  % 0.2Hz or 10 Hz
param.Q = 274; % Capacity set to lower so the entire signal can be covered

%% Generate MLBS signal

% Number of bits (n) -------see [Fairweather 2011]
param.n = 8; 

% Source clock period (t_clk)
switch PRBSfreq % n = 50 for f = 0.2Hz; n = 1 for f = 10 Hz
    case "0.2"
        n = 50;  
    case "10"
        n = 1; 
end
param.t_clk = n*0.1;              % Multiple of sampling time (0.1s)
param.f_clk = 1/param.t_clk;      % Broadest components for power spectrum are below 2 Hz, freq<=10Hz

% Amplitude (a)
param.a = 0.25 * param.Q; %max(abs(CI95_pdf)); 

% PRBS signal with amplitude param.a
mlbs = prbs(param.n, 2^param.n-1);
mlbs_signal = (mlbs' - 0.5) * 2 * param.a; % Convert to symmetric signal with amplitude a

% Upsample to cycler sampling time
mlbs_signal = upsample(mlbs_signal, n); 
ind = find(mlbs_signal ~= 0); 
for i = 1:height(ind)
    if i ~= height(ind)
        k = 1; 
        while ind(i)+k < ind(i+1)
            mlbs_signal(ind(i)+k) = mlbs_signal(ind(i)); 
            k = k + 1; 
        end
    else
        k = 1; 
        while ind(i) + k <= height(mlbs_signal)
            mlbs_signal(ind(i)+k) = mlbs_signal(ind(i)); 
            k = k + 1; 
        end
    end
end

% Time
time_total = (height(mlbs_signal)-1) * deltaT; % Time of one full period of the prbs signal [s]
time_mlbs = [0:deltaT:time_total]';

% Write to 'signals'
signals.mlbs = [time_mlbs, mlbs_signal];  % Maximum length binary sequence

%% Generate IRBS from MLBS

% mlbs = [1 0 1 1 0 0 1 0 0 0 1 1 1 1 0]; % Example given in [Zhu 2019]

% Generate irbs signal
irbs = [mlbs mlbs]; % Double the mlbs
ind = [1:2:length(irbs)];
irbs_inv = ~irbs;
irbs(ind) = irbs_inv(ind);

irbs_signal = (irbs' - 0.5) * 2 * param.a; % convert to symmetric signal with amplitude a

% Upsample to cycler sampling time
irbs_signal = upsample(irbs_signal, n); 
ind = find(irbs_signal ~= 0); 
for i = 1:height(ind)
    if i ~= height(ind)
        k = 1; 
        while ind(i)+k < ind(i+1)
            irbs_signal(ind(i)+k) = irbs_signal(ind(i)); 
            k = k + 1; 
        end
    else
        k = 1; 
        while ind(i) + k <= height(irbs_signal)
            irbs_signal(ind(i)+k) = irbs_signal(ind(i)); 
            k = k + 1; 
        end
    end
end

% Time
time_total = (height(irbs_signal)-1) * deltaT; % Time of one full period of the prbs signal [s]
time_irbs = [0:deltaT:time_total]';

% Write to 'signals'
signals.irbs = [time_irbs, irbs_signal];  % Maximum length binary sequence

% Plot generated mlbs & irbs signals
% figure; hold on
% ax1 = subplot(2, 1, 1); plot(time_mlbs, mlbs_signal, '.-'); grid on; title("MLBS"); xlabel("Time [s]"); ylabel("Current [A]")
% ax2 = subplot(2, 1, 2); plot(time_irbs, irbs_signal, '.-'); grid on; title("IRBS"); xlabel("Time [s]"); ylabel("Current [A]")
% linkaxes([ax1 ax2], 'x')

% Max consecutive elements in the signals
mlbs_max = max(diff(find(diff([NaN mlbs_signal' NaN])))); 
param.mlbs_SOCVariation = mlbs_max * deltaT * param.a / (param.Q * 3600) * 100; 

irbs_max = max(diff(find(diff([NaN irbs_signal' NaN])))); 
param.irbs_SOCVariation = irbs_max * deltaT * param.a / (param.Q * 3600) * 100; 

%% Generate pulse signals

% Pulse magnitude (I)
CRate = 0.8; 
param.I = CRate * param.Q;

% Short pulse duration
deltaSOC1 = 0.01;
param.t_pulse1 = param.Q * deltaSOC1 * 3600 / param.I;

% Create pulse1 signal
time = [0:deltaT:param.t_pulse1]';
pulse1 = ones(height(time), 1) * param.I;
signals.pulse1 = [time pulse1];

% Long pulse duration
deltaSOC2 = 0.05;
param.t_pulse2 = param.Q * deltaSOC2 * 3600 / param.I;

% Create pulse2 signal
time = [0:deltaT:param.t_pulse2]';
pulse2 = ones(height(time), 1) * param.I;
signals.pulse2 = [time pulse2];
        
%% Generate rest signal
param.t_relax = 3600; 
time = [0:deltaT:param.t_relax]';
rest = zeros(height(time), 1); 
signals.rest = [time rest];

%% Combinations of signals

% Combine pulse1 + rest for 1% SOC drop
sigToCombine = {signals.pulse1(:, 2); signals.rest(:, 2)}; 
signals.pulse1_mlbs_rest = combineSignals(sigToCombine, deltaT); 

% Combine pulse2 + rest for 5% SOC drop
sigToCombine = {signals.pulse2(:, 2); signals.rest(:, 2)}; 
signals.pulse1_mlbs_rest = combineSignals(sigToCombine, deltaT); 

% Combine pulse1 + MLBS + rest (1% SOC)
sigToCombine = {signals.pulse1(:, 2); signals.mlbs(:, 2); signals.rest(:, 2)}; 
signals.pulse1_mlbs_rest = combineSignals(sigToCombine, deltaT); 

% Combine pulse2 + MLBS + rest (5% SOC)
sigToCombine = {signals.pulse2(:, 2); signals.mlbs(:, 2); signals.rest(:, 2)}; 
signals.pulse2_mlbs_rest = combineSignals(sigToCombine, deltaT); 

% Combine pulse-MLBS signals to cover entire SOC range
time = signals.pulse1_mlbs_rest(:, 1); 
current = signals.pulse1_mlbs_rest(:, 2);
for i = 2:36
    % Long pulse or short pulse
    if i >=2 && i <= 10
        profile = signals.pulse1_mlbs_rest;
    elseif i >= 11 && i <= 26
        profile = signals.pulse2_mlbs_rest;
    elseif i >= 27 && i <= 36
        profile = signals.pulse1_mlbs_rest;
    end
    
    % Add to existing
    profile(:, 1) = time(end) + deltaT + profile(:, 1); 
    time = [time; profile(:, 1)]; 
    current = [current; profile(:, 2)]; 
end
signals.pulse_mlbs_rest_entireRange = [time current]; 

% Combine pulse1 + IRBS + rest (1% SOC)
sigToCombine = {signals.pulse1(:, 2); signals.irbs(:, 2); signals.rest(:, 2)}; 
signals.pulse1_irbs_rest = combineSignals(sigToCombine, deltaT); 

% Combine pulse2 + IRBS + rest (5% SOC)
sigToCombine = {signals.pulse2(:, 2); signals.irbs(:, 2); signals.rest(:, 2)}; 
signals.pulse2_irbs_rest = combineSignals(sigToCombine, deltaT); 

% Combine pulse-IRBS signals to cover entire SOC range
time = signals.pulse1_irbs_rest(:, 1); 
current = signals.pulse1_irbs_rest(:, 2);
for i = 2:36
    % Long pulse or short pulse
    if i >=2 && i <= 10
        profile = signals.pulse1_irbs_rest;
    elseif i >= 11 && i <= 26
        profile = signals.pulse2_irbs_rest;
    elseif i >= 27 && i <= 36
        profile = signals.pulse1_irbs_rest;
    end
    
    % Add to existing
    profile(:, 1) = time(end) + deltaT + profile(:, 1); 
    time = [time; profile(:, 1)]; 
    current = [current; profile(:, 2)]; 
end
signals.pulse_irbs_rest_entireRange = [time current]; 

%% Plotting & save all generated signals
figure
ax1 = subplot(6, 1, 1); plot(signals.pulse1_mlbs_rest(:, 1), signals.pulse1_mlbs_rest(:, 2), '.-'); grid on
ax2 = subplot(6, 1, 2); plot(signals.pulse2_mlbs_rest(:, 1), signals.pulse2_mlbs_rest(:, 2), '.-'); grid on
ax3 = subplot(6, 1, 3); plot(signals.pulse_mlbs_rest_entireRange(:, 1), signals.pulse_mlbs_rest_entireRange(:, 2), '.-'); grid on
ax4 = subplot(6, 1, 4); plot(signals.pulse1_irbs_rest(:, 1), signals.pulse1_irbs_rest(:, 2), '.-'); grid on
ax5 = subplot(6, 1, 5); plot(signals.pulse2_irbs_rest(:, 1), signals.pulse2_irbs_rest(:, 2), '.-'); grid on
ax6 = subplot(6, 1, 6); plot(signals.pulse_irbs_rest_entireRange(:, 1), signals.pulse_irbs_rest_entireRange(:, 2), '.-'); grid on
linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x')

figure
ax1 = subplot(2, 1, 1); plot(signals.pulse_mlbs_rest_entireRange(:, 1)/3600, signals.pulse_mlbs_rest_entireRange(:, 2), '.-'); 
grid on; title("MLBS"); xlabel("Time [h]"); ylabel("Current [A]")
ax2 = subplot(2, 1, 2); plot(signals.pulse_irbs_rest_entireRange(:, 1)/3600, signals.pulse_irbs_rest_entireRange(:, 2), '.-'); 
grid on; title("IRBS"); xlabel("Time [h]"); ylabel("Current [A]")
linkaxes([ax1 ax2], 'x')

folder = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\9 - Pulse-PRBS Tests\Profiles\"; 
filename = folder + "PulsePRBS.mat"; 
save(filename, 'signals', 'param', '-mat')

%% Format tests for Digatron

% Pulse1 + MLBS
sigToCombine = {signals.pulse1(:, 2); signals.mlbs(:, 2)}; 
signals.pulse1_mlbs = combineSignals(sigToCombine, deltaT); 
% figure; plot(signals.pulse1_mlbs(:, 1), signals.pulse1_mlbs(:, 2))
time = 0.1 * ones(height(signals.pulse1_mlbs), 1);
command = string(time) + ' sec;' + string(-signals.pulse1_mlbs(:, 2)) +  ';;;';
folder = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\9 - Pulse-PRBS Tests\Input\";
filename = folder + "EVE280_pulse1_mlbs.txt";
% writematrix(command, filename);

% % Pulse2 + MLBS
sigToCombine = {signals.pulse2(:, 2); signals.mlbs(:, 2)}; 
signals.pulse2_mlbs = combineSignals(sigToCombine, deltaT); 
% figure; plot(signals.pulse2_mlbs(:, 1), signals.pulse2_mlbs(:, 2))
time = 0.1 * ones(height(signals.pulse2_mlbs), 1);
command = string(time) + ' sec;' + string(-signals.pulse2_mlbs(:, 2)) +  ';;;';
folder = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\9 - Pulse-PRBS Tests\Input\";
filename = folder + "EVE280_pulse2_mlbs.txt";
% writematrix(command, filename);

%% Functions
function out = combineSignals(signals, deltaT)
% Combines singals assuming the same sampling time

    % Combine signals
    combined = []; 
    for i = 1:height(signals)
        combined = [combined; signals{i}];
    end

    % Time
    time_total = (height(combined)-1) * deltaT; % Time of one full period of the prbs signal [s]
    time = [0:deltaT:time_total]';

    out = [time combined];
end

function [P1, f] = fft_fcn(x, Fs)
% Compute power spectrum of input signal x with sample frequency Fs [Hz]
    
    T = 1/Fs;               % Sampling interval
    L = height(x);          % Number of time points

    Y = fft(x); 

    P2 = abs(Y)/L;

    P1 = P2(1:(fix(L/2)+1))';
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = Fs*(0:(L/2))'/L;
    f2 = ((0:1/L:1-1/L)*Fs).';

end