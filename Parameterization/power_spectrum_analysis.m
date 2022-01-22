% clearvars -except data deltaT dataType CRate
% data = data(1:8, :); 

clear
close all

currentFolder = cd; 
addpath(currentFolder)
deltaT = 0.1; % Sampling time for all data
        
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
        load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\9 - Pulse-PRBS Tests\ModifiedPulse.mat');
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

%% Histogram, pdf & power percentage of drive cycle current

% figure_histogram = figure; 
% figure_powerSpectrum = figure;
% figure_powerPercentage = figure; 

edges = [-280:5:280];
freqrange = [0.01:0.001:0.5];

distribution.mu = [];
distribution.sigma = [];
distribution.CI95 = []; 
distribution.f_cutoff = []; 
distribution.totalPower = []; 

for datanum = 1:5
    profile = data{datanum, 2}.Current; 
    profileName = data{datanum, 1};
    
    % Plot histogram
%     figure(figure_histogram)
%     subplot(5, 1, datanum); 
%     H = histogram(profile, edges, 'Normalization','pdf'); 
%     xlim([-150, 150]);
%     grid on; title(profileName); xlabel('I [A]'); ylabel('P(I)')
    
    % Find & plot pdf
    [mu, sigma] = normfit(profile);
    current_range = [-250:0.1:250];
    pdf_normal = pdf('Normal', current_range, mu, sigma);

%     figure(figure_histogram); hold on; 
%     plot(current_range, pdf_normal, 'color', 'r');

    distribution.mu = [distribution.mu; mu];
    distribution.sigma = [distribution.sigma; sigma]; 

    % 95% interval of pdf
    ts_pdf = tinv([0.05  0.95], length(profile)-1);     % T-Score
    CI95_pdf = mu + ts_pdf * sigma;                     % 95% confidence interval
    y = [0:0.001:max(pdf_normal)]';
    x = ones(height(y), 1) .* [CI95_pdf mu];

%     figure(figure_histogram); hold on; 
%     plot(x(:, 1), y, 'color', 'r');
%     plot(x(:, 2), y, 'color', 'r');
%     plot(x(:, 3), y, 'color', 'r');

    distribution.CI95 = [distribution.CI95; CI95_pdf]; 

    % Power percentage
    P1 = data{datanum, 2}.P1; 
    f = data{datanum, 2}.f; 
%     totalPower = sum(P1);
%     perc = cumsum(P1)/totalPower; 
%     ind = find(perc>= 0.9, 1); 
%     f_cutoff = f(ind); % Frequency at which 90% of power is covered
    perc = 0; 
    totalPower = bandpower(profile,10,[0 max(f)]); 
    for i = 1:length(freqrange)
        power = bandpower(profile,10,[0 freqrange(i)]); 
        perc = power/totalPower; 
        if perc >= 0.9
            break
        end
    end
    f_cutoff = freqrange(i);

    distribution.f_cutoff = [distribution.f_cutoff; f_cutoff]; 
    distribution.totalPower = [distribution.totalPower; totalPower]; 
    
    % Plot power percentage curves
%     figure(figure_powerPercentage);
%     subplot(5, 1, datanum); hold on;
%     plot(f, perc)
%     x = [0:0.1:5];
%     y = ones(length(x)) * 0.9; 
%     plot(x, y)
%     grid on; title(profileName)

    % Plot power spectrums
%     figure(figure_powerSpectrum);
%     subplot(5, 1, datanum); hold on;
%     plot(data{datanum, 2}.f, data{datanum, 2}.P1);
%     grid on; xlabel('f [Hz]'); ylabel('|Current(f)| [A]');
%     title(profileName)
%     set(gca, 'XScale', 'log');
end


%% Power spectrums
figure; hold on; 
plot(data{9, 2}.f, data{9, 2}.P1, 'LineWidth',0.1, 'color', '#EDB120'); % Pulse-MLBS
plot(data{10, 2}.f, data{10, 2}.P1, 'LineWidth',0.1, 'color', '#4DBEEE'); % Pulse-IRBS
grid on; xlabel('f [Hz]'); ylabel('|Current(f)| [A]'); legend("Pulse-MLBS", 'Pulse-IRBS')
ylim([0, 15])
title("Power Spectrums of Test Profiles")
set(gca, 'XScale', 'log');


figure; hold on; 
% plot(data{4, 2}.f, data{4, 2}.P1, 'LineWidth',0.1, 'color', '#0072BD');  % US06
plot(data{5, 2}.f, data{5, 2}.P1, 'LineWidth',0.1, 'color', '#0072BD');  % Mixed
plot(data{3, 2}.f, data{3, 2}.P1, 'LineWidth',0.1, 'color', '#7E2F8E'); % UDDS
plot(data{9, 2}.f, data{9, 2}.P1, 'LineWidth',0.1, 'color', '#EDB120'); % Pulse-MLBS
plot(data{8, 2}.f, data{8, 2}.P1, 'LineWidth',0.1, 'color', '#77AC30'); % Pulse 0.8C
grid on; xlabel('f [Hz]'); ylabel('|Current(f)| [A]'); legend("Mixed", "UDDS", 'Pulse-MLBS', 'Pulse 0.8C')
ylim([0, 15])
title("Power Spectrums of Test Profiles")
set(gca, 'XScale', 'log');

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
