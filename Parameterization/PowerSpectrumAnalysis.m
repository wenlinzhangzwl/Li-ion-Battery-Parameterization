% clearvars -except data deltaT dataType CRate
% data = data(1:8, :); 

clear
close all

currentFolder = cd; 
addpath(currentFolder)
deltaT = 0.1; % Sampling time for all data

%% Generate pulse & rest signals

% Pulse magnitude (I)
CRate = 0.8; 
cellCapacity = 280; 
param.I = CRate * cellCapacity; 

% Pulse1 duration (t_pulse)
deltaSOC1 = 0.01;
param.t_pulse1 = (cellCapacity*3600) * deltaSOC1 / param.I;

% Create pulse1 signal
time = [0:deltaT:param.t_pulse1]';
pulse1 = ones(height(time), 1) * param.I;
signals.pulse1 = [time pulse1];

% Pulse2 duration (t_pulse)
deltaSOC2 = 0.05;
param.t_pulse2 = (cellCapacity*3600) * deltaSOC2 / param.I;

% Create pulse2 signal
time = [0:deltaT:param.t_pulse2]';
pulse2 = ones(height(time), 1) * param.I;
signals.pulse2 = [time pulse2];
        
% Create rest signal
param.t_relax = 3600; 
time = [0:deltaT:param.t_relax]';
rest = zeros(height(time), 1); 
signals.rest = [time rest];
        
%% Import data

folder_project = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\';

dataType = "entireTest";

switch dataType
    case "singlePeriod"
        %% Load dive cycles
        folder_data = append(folder_project, 'Scripts\Drive Train Model - Ford Focus 2013\Generate drive cycles\');
        addpath(folder_data)

        n_sample = deltaT/1e-3;
        
        data.HWFET = load('HWFET_1000Hz.mat', 'Q', 'drivecycle_current');
        data.HWFET.Time = data.HWFET.drivecycle_current(1:n_sample:end, 1);
        data.HWFET.Current = data.HWFET.drivecycle_current(1:n_sample:end, 2);
        data.HWFET.Fs = 1/deltaT;
        
        data.NEDC = load('NEDC_1000Hz.mat', 'Q', 'drivecycle_current');
        data.NEDC.Time = data.NEDC.drivecycle_current(1:n_sample:end, 1);
        data.NEDC.Current = data.NEDC.drivecycle_current(1:n_sample:end, 2);
        data.NEDC.Fs = 1/deltaT;
        
        data.UDDS = load('UDDS_1000Hz.mat', 'Q', 'drivecycle_current');
        data.UDDS.Time = data.UDDS.drivecycle_current(1:n_sample:end, 1);
        data.UDDS.Current = data.UDDS.drivecycle_current(1:n_sample:end, 2);
        data.UDDS.Fs = 1/deltaT;

        data.US06 = load('US06_1000Hz.mat', 'Q', 'drivecycle_current');
        data.US06.Time = data.US06.drivecycle_current(1:n_sample:end, 1);
        data.US06.Current = data.US06.drivecycle_current(1:n_sample:end, 2);
        data.US06.Fs = 1/deltaT; 
        
        %% Generate pulse + rest for 1% SOC drop

        % Combine pulse1 & rest
        time = signals.pulse1(end, 1) + deltaT + signals.rest(:, 1); 
        signals.pulse1_rest = [signals.pulse1; time signals.rest(:, 2)]; 

        % Load to 'data'
        data.pulse1.Time = signals.pulse1_rest(:, 1);
        data.pulse1.Current = signals.pulse1_rest(:, 2);
        data.pulse1.Fs = 1/deltaT; 

        % Combine pulse2 & rest
        time = signals.pulse2(end, 1) + deltaT + signals.rest(:, 1); 
        signals.pulse2_rest = [signals.pulse2; time signals.rest(:, 2)]; 
        
        % Load to 'data'
        data.pulse2.Time = signals.pulse2_rest(:, 1);
        data.pulse2.Current = signals.pulse2_rest(:, 2);
        data.pulse2.Fs = 1/deltaT; 
        
        % Plot two pulses
%         figure
%         ax1 = subplot(2, 1, 1); plot(data.pulse1.Time, data.pulse1.Current, '.-'); grid on
%         ax2 = subplot(2, 1, 2); plot(data.pulse2.Time, data.pulse2.Current, '.-'); grid on
%         linkaxes([ax1 ax2], 'x')

    case "entireTest"
        %% Load drive cycle tests
        folder_data = append(folder_project, '\Test Results\6 - Drive cycle\');
        addpath(folder_data)
        
        load('HWFET25.mat');
        data.HWFET.Time = meas.Time;
        data.HWFET.Current = -meas.Current;
        data.HWFET.Fs = 10;
        
        load('NEDC23.mat');
        data.NEDC.Time = meas.Time;
        data.NEDC.Current = -meas.Current;
        data.NEDC.Fs = 10;
        
        load('UDDS25.mat');
        data.UDDS.Time = meas.Time;
        data.UDDS.Current = -meas.Current;
        data.UDDS.Fs = 10;
        
        load('US0625.mat');
        data.US06.Time = meas.Time;
        data.US06.Current = -meas.Current;
        data.US06.Fs = 10;
        
        load('MIX25.mat');
        data.Mixed.Time = meas.Time;
        data.Mixed.Current = -meas.Current;
        data.Mixed.Fs = 10;

        %% Load pulse tests
        load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\9 - Modified pulse test\ModifiedPulse.mat');
        data.ModifiedPulse.Time = meas.Time;
        data.ModifiedPulse.Current = meas.Current;
        data.ModifiedPulse.Fs = 10;
        
        folder_data = append(folder_project, 'Test Results\7 - Pulse Test');
        addpath(folder_data); 

        load('PROCESSED_PUL_0p1_discharge.mat'); 
        data.Pulse0p1.Time = meas.Time;
        data.Pulse0p1.Current = meas.Current;
        data.Pulse0p1.Fs = 10;

        Pulse0p8 = load('PROCESSED_PUL25dch.mat'); 
        data.Pulse0p8.Time = meas.Time;
        data.Pulse0p8.Current = meas.Current;
        data.Pulse0p8.Fs = 10;
        
end

%% Generate power spectrums

profiles = fieldnames(data); 
data = struct2cell(data); 
data = [profiles data];

for i = 1:height(data)
    
    data{i, 2}.P1 = [];
    data{i, 2}.f = [];
    
    x = data{i, 2}.Current; 
    t = data{i, 2}.Time; 
    
    ind = find(isnan(x)); 
    if ~isempty(ind)
        x(ind) = []; 
        t(ind) = [];
        data{i, 2}.Current = x; 
        data{i, 2}.Time = t;
    end

    %% Compute fft   
    
    Fs = data{i, 2}.Fs; % Sampling frequency [Hz]
    L = height(t);      % Length of time signal
     
    % magnitude
    Y = fft(x);                     % fft of x 
    P2 = abs(Y)/L;                  % magnitude in [A] normalized by signal length(L) 
    P1 = P2(1:fix(L/2)+1);          % single-sided fft
    P1(2:end-1) = 2*P1(2:end-1);    % magnitude multiplied by 2 since negative side is ignored
    
    % frequency
    f = Fs*(0:(L/2))'/L;
    
    % write to 'data'
    data{i, 2}.f = f; 
    data{i, 2}.P1 = P1; 
end

%% Histogram & PDF of a selected drive cycle

datanum = 3; 
profile = data{datanum, 2}.Current; 
profileName = data{datanum, 1};

% Plot power spectrum
% figure; 
% plot(data{datanum, 2}.f, data{datanum, 2}.P1);%, 'LineWidth',0.1); 
% grid on; xlabel('f [Hz]'); ylabel('|Current(f)| [A]');
% title(profileName)
% set(gca, 'XScale', 'log');

% Plot histogram
figure_histogram = figure; 
H = histogram(profile, 'Normalization','probability'); 
grid on; title(profileName); xlabel('Drive Cycle Current [A]'); ylabel('Relative Frequency')

% Find & plot pdf
[mu, sigma] = normfit(profile); 
pdf_normal = pdf('Normal', profile, mu, sigma);
figure_histogram; hold on; 
plot(profile, pdf_normal, 'linewidth', 1);

% 95% interval of pdf
ts_pdf = tinv([0.05  0.95], length(profile)-1);      % T-Score
CI95_pdf = mu + ts_pdf * sigma;
y = [0:0.001:max(pdf_normal)]';
x = ones(height(y), 1) .* [CI95_pdf mu];
figure_histogram; hold on; 
plot(x(:, 1), y, 'linewidth', 2);
plot(x(:, 2), y, 'linewidth', 2);
plot(x(:, 3), y, 'linewidth', 2);

%% Generate MLBS signal

% Amplitude (a)
param.a = CRate * 280; %max(abs(CI95_pdf)); 

% Source clock period (t_clk)
f_clk = 0.2; % Broadest components for power spectrum are below 2 Hz, freq<=10Hz
param.t_clk = 1/f_clk; 

% Number of bits (n) -------see [Fairweather 2011]
param.n = 8; 

% PRBS signal with amplitude param.a
mlbs = prbs(param.n, 2^param.n-1);
mlbs_signal = (mlbs' - 0.5) * 2 * param.a; % convert to symmetric signal with amplitude a

% Upsample to cycler sampling time
mlbs_signal = upsample(mlbs_signal, param.t_clk/deltaT); 
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
irbs_signal = upsample(irbs_signal, param.t_clk/deltaT); 
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
% plot(time_mlbs, mlbs_signal, '.-');
% plot(time_irbs, irbs_signal, '.-'); 
% grid on; legend('mlbs', 'irbs')


%% Combine pulse1 + MLBS + rest (1% SOC)

% % Combine signals
sigToCombine = {signals.pulse1(:, 2); signals.mlbs(:, 2); signals.rest(:, 2)}; 
signals.pulse1_mlbs_rest = combineSignals(sigToCombine, deltaT); 

% Write to 'data' only if single periods of signals are compared
if dataType == "singlePeriod"
    % Write to 'data' 
    ind = height(data)+1;
    data{ind, 1} = "PULSE-MLBS 1C 1% SOC";
    data{ind, 2}.Time = signals.pulse1_mlbs_rest(:, 1);
    data{ind, 2}.Current = signals.pulse1_mlbs_rest(:, 2);
    data{ind, 2}.Fs = 1/deltaT; 
    
    % Compute power spectrum
    [data{ind, 2}.P1, data{ind, 2}.f] = fft_fcn(signals.pulse1_mlbs_rest(:, 2), 1/deltaT); 
end

%% Combine pulse2 + MLBS + rest (5% SOC)

% % Combine signals
sigToCombine = {signals.pulse2(:, 2); signals.mlbs(:, 2); signals.rest(:, 2)}; 
signals.pulse2_mlbs_rest = combineSignals(sigToCombine, deltaT); 

% Write to 'data' only if single periods of signals are compared
if dataType == "singlePeriod"
    % Write to 'data' 
    ind = height(data)+1;
    data{ind, 1} = "PULSE-MLBS 1C 5% SOC";
    data{ind, 2}.Time = signals.pulse2_mlbs_rest(:, 1);
    data{ind, 2}.Current = signals.pulse2_mlbs_rest(:, 2);
    data{ind, 2}.Fs = 1/deltaT; 
    
    % Compute power spectrum
    [data{ind, 2}.P1, data{ind, 2}.f] = fft_fcn(signals.pulse2_mlbs_rest(:, 2), 1/deltaT); 
end

%% Combine pulse-MLBS signals to cover entire SOC range

if dataType == "entireTest"
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
    
    % Write to 'data' 
    ind = height(data)+1; 
    data{ind, 1} = "Pulse-MLBS";
    data{ind, 2}.Time = signals.pulse_mlbs_rest_entireRange(:, 1);
    data{ind, 2}.Current = signals.pulse_mlbs_rest_entireRange(:, 2);
    data{ind, 2}.Fs = 1/deltaT; 
    
    % Compute power spectrum
    [data{ind, 2}.P1, data{ind, 2}.f] = fft_fcn(signals.pulse_mlbs_rest_entireRange(:, 2), 1/deltaT); 
    
end


%% Combine pulse1 + IRBS + rest (1% SOC)

% Combine signals
sigToCombine = {signals.pulse1(:, 2); signals.irbs(:, 2); signals.rest(:, 2)}; 
signals.pulse1_irbs_rest = combineSignals(sigToCombine, deltaT); 

% Write to 'data' only if single periods of signals are compared
if dataType == "singlePeriod"
    % Write to 'data' 
    ind = height(data)+1;
    data{ind, 1} = "PULSE-IRBS 1C 1% SOC";
    data{ind, 2}.Time = signals.pulse1_irbs_rest(:, 1);
    data{ind, 2}.Current = signals.pulse1_irbs_rest(:, 2);
    data{ind, 2}.Fs = 1/deltaT; 
    
    % Compute power spectrum
    [data{ind, 2}.P1, data{ind, 2}.f] = fft_fcn(signals.pulse1_irbs_rest(:, 2), 1/deltaT); 
end

%% Combine pulse2 + IRBS + rest (5% SOC)

% % Combine signals
sigToCombine = {signals.pulse2(:, 2); signals.irbs(:, 2); signals.rest(:, 2)}; 
signals.pulse2_irbs_rest = combineSignals(sigToCombine, deltaT); 

% Write to 'data' only if single periods of signals are compared
if dataType == "singlePeriod"
    % Write to 'data' 
    ind = height(data)+1;
    data{ind, 1} = "PULSE-IRBS 1C 5% SOC";
    data{ind, 2}.Time = signals.pulse2_irbs_rest(:, 1);
    data{ind, 2}.Current = signals.pulse2_irbs_rest(:, 2);
    data{ind, 2}.Fs = 1/deltaT; 
    
    % Compute power spectrum
    [data{ind, 2}.P1, data{ind, 2}.f] = fft_fcn(signals.pulse2_irbs_rest(:, 2), 1/deltaT); 
end

%% Combine pulse-IRBS signals to cover entire SOC range

if dataType == "entireTest"
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
    
    % Write to 'data' 
    ind = height(data)+1; 
    data{ind, 1} = "Pulse-IRBS";
    data{ind, 2}.Time = signals.pulse_irbs_rest_entireRange(:, 1);
    data{ind, 2}.Current = signals.pulse_irbs_rest_entireRange(:, 2);
    data{ind, 2}.Fs = 1/deltaT; 
    
    % Compute power spectrum
    [data{ind, 2}.P1, data{ind, 2}.f] = fft_fcn(signals.pulse_irbs_rest_entireRange(:, 2), 1/deltaT); 
    
end

%% Plot all generated signals
figure
ax1 = subplot(6, 1, 1); plot(signals.pulse1_mlbs_rest(:, 1), signals.pulse1_mlbs_rest(:, 2), '.-'); grid on
ax2 = subplot(6, 1, 2); plot(signals.pulse2_mlbs_rest(:, 1), signals.pulse2_mlbs_rest(:, 2), '.-'); grid on
ax3 = subplot(6, 1, 3); plot(signals.pulse_mlbs_rest_entireRange(:, 1), signals.pulse_mlbs_rest_entireRange(:, 2), '.-'); grid on
ax4 = subplot(6, 1, 4); plot(signals.pulse1_irbs_rest(:, 1), signals.pulse1_irbs_rest(:, 2), '.-'); grid on
ax5 = subplot(6, 1, 5); plot(signals.pulse2_irbs_rest(:, 1), signals.pulse2_irbs_rest(:, 2), '.-'); grid on
ax6 = subplot(6, 1, 6); plot(signals.pulse_irbs_rest_entireRange(:, 1), signals.pulse_irbs_rest_entireRange(:, 2), '.-'); grid on
linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x')

%% Power spectrums

figure; hold on; 
plot(data{4, 2}.f, data{4, 2}.P1, 'LineWidth',0.1); 
plot(data{9, 2}.f, data{9, 2}.P1, 'LineWidth',0.1); 
plot(data{10, 2}.f, data{10, 2}.P1, 'LineWidth',0.1); 
plot(data{3, 2}.f, data{3, 2}.P1, 'LineWidth',0.1); 
grid on; xlabel('f [Hz]'); ylabel('|Current(f)| [A]'); legend("US06", 'Pulse-MLBS', 'Pulse-IRBS', "UDDS")
title("Power Spectrums of Test Profiles")
% xlim([0, 0.4]); ylim([0, 20]); 
set(gca, 'XScale', 'log');


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