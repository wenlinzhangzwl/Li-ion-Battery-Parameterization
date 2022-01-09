% clearvars -except data
clear

currentFolder = cd; 
addpath(currentFolder)

%% Import data

% Pulse tests
dataFolder = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\';
addpath(dataFolder); 

Pulse0p1 = load('PROCESSED_PUL_0p1_discharge.mat'); 
data.Pulse0p1 = table2struct(Pulse0p1.meas_t, 'ToScalar', true);
data.Pulse0p1.Fs = 10; % Sampling rate
clear Pulse0p1

Pulse0p8 = load('PROCESSED_PUL25dch.mat'); 
data.Pulse0p8 = table2struct(Pulse0p8.meas_t, 'ToScalar', true);
data.Pulse0p8.Fs = 10; % Sampling rate
clear Pulse0p8

% Drive cycles
DCtype = "experimental";
switch DCtype
    case "simulation"
        dataFolder = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Model\Drive Train Model - Ford Focus 2013\Generate drive cycles\';
        addpath(dataFolder)

        data.HWFET = load('HWFET.mat');
        data.HWFET.Fs = 1000;
        data.NEDC = load('NEDC.mat');
        data.NEDC.Fs = 1000;
        data.UDDS = load('UDDS.mat');
        data.UDDS.Fs = 1000;
        data.US06 = load('US06.mat');
        data.US06.Fs = 1000;
    case "experimental"
        dataFolder = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\6 - Drive cycle\';
        addpath(dataFolder)
        
        data.HWFET = load('HWFET25.mat');
        data.HWFET = data.HWFET.meas;
        data.HWFET.Fs = 10;
        data.NEDC = load('NEDC23.mat');
        data.NEDC = data.NEDC.meas;
        data.NEDC.Fs = 10;
        data.UDDS = load('UDDS25.mat');
        data.UDDS = data.UDDS.meas;
        data.UDDS.Fs = 10;
        data.US06 = load('US0625.mat');
        data.US06 = data.US06.meas;
        data.US06.Fs = 10;
        data.Mixed = load('MIX25.mat');
        data.Mixed = data.Mixed.meas;
        data.Mixed.Fs = 10;
        data.ModifiedPulse = load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\9 - Modified pulse test\ModifiedPulse.mat');
        data.ModifiedPulse = data.ModifiedPulse.meas;
        data.ModifiedPulse.Fs = 10;
end


% Generate & plot power spectrum

data_struct = data; 
data = struct2cell(data); 
if DCtype == "experimental"
    profiles = {"PULSE 0.1C"; "PULSE 0.8C"; "HWFET"; "NEDC"; "UDDS"; "US06"; "Mixed"; "ModifiedPulse"};  %#ok<CLARRSTR>
elseif DCtype == "simulation"
    profiles = {"PULSE 0.1C"; "PULSE 0.8C"; "HWFET"; "NEDC"; "UDDS"; "US06"};  %#ok<CLARRSTR>
end
data = [profiles data];


for i = 1:height(data)
    data{i, 2}.Y = [];
    data{i, 2}.P1 = [];
    data{i, 2}.f = [];
    data{i, 2}.f2 = [];
    
    if i <=2
        x = data{i, 2}.Current; 
        t = data{i, 2}.Time; 
    else
        if DCtype == "experimental"
            x = data{i, 2}.Current; 
            t = data{i, 2}.Time; 
            ind = find(isnan(x)); 
            if ~isempty(ind)
                x(ind) = []; 
                t(ind) = [];
            end
        elseif DCtype == "simulation"
            x = data{i, 2}.drivecycle_current(:, 2);
            t = data{i, 2}.drivecycle_current(:, 1);
        end
    end

    Fs = data{i, 2}.Fs; % Sampling frequency
    T = 1/Fs;           % Sampling interval
    L = height(x);      % Number of time points
    Y = fft(x); 
    P2 = abs(Y/L);
    P1 = P2(1:(fix(L/2)+1));
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))'/L;
    f2 = ((0:1/L:1-1/L)*Fs).';

    data{i, 2}.Y = Y; 
    data{i, 2}.P1 = P1; 
    data{i, 2}.f = f; 
    data{i, 2}.f2 = f2; 
end


%%
% fig_fft = figure; 
% figure(fig_fft); hold on; 
% plot(data{6, 2}.f, data{6, 2}.P1); 
% plot(data{5, 2}.f, data{5, 2}.P1); 
% plot(data{4, 2}.f, data{4, 2}.P1); 
% grid on; xlabel('f [Hz]'); ylabel('|P1(f)|'); legend("US06", "UDDS", "NEDC"); title("Power Spectrums of Test Profiles")
% % set(gca, 'YScale', 'log', 'XScale', 'log');
%  set(gca, 'XScale', 'log');

% %%
% fig_fft = figure; 
% figure(fig_fft); hold on; 
% plot(data{6, 2}.f, data{6, 2}.P1); 
% plot(data{7, 2}.f, data{7, 2}.P1); 
% plot(data{5, 2}.f, data{5, 2}.P1); 
% grid on; xlabel('f [Hz]'); ylabel('|P1(f)|'); legend("US06", "Mixed", "UDDS")
% % set(gca, 'YScale', 'log', 'XScale', 'log');
%  set(gca, 'XScale', 'log');
% 

%%
% fig_fft = figure; 
% figure(fig_fft); hold on; 
% plot(data{5, 2}.f, data{5, 2}.P1, 'LineWidth',0.1); 
% plot(data{6, 2}.f, data{6, 2}.P1, 'LineWidth',0.1); 
% plot(data{2, 2}.f, data{2, 2}.P1, 'LineWidth',0.1); 
% plot(data{1, 2}.f, data{1, 2}.P1, 'LineWidth',0.1); 
% grid on; xlabel('f [Hz]'); ylabel('|P1(f)|'); legend("UDDS", "US06", "Pulse 0.8C", "Pulse 0.1C")
% % set(gca, 'YScale', 'log', 'XScale', 'log');
%  set(gca, 'XScale', 'log');
% 
%%
fig_fft = figure; 
figure(fig_fft); hold on; 
plot(data{6, 2}.f, data{6, 2}.P1, 'LineWidth',0.1); 
plot(data{4, 2}.f, data{4, 2}.P1, 'LineWidth',0.1); 
plot(data{5, 2}.f, data{5, 2}.P1, 'LineWidth',0.1);
plot(data{2, 2}.f, data{2, 2}.P1, 'LineWidth',0.1); 
plot(data{1, 2}.f, data{1, 2}.P1, 'LineWidth',0.1); 
plot(data{7, 2}.f, data{7, 2}.P1, 'LineWidth',0.1); 
grid on; xlabel('f [Hz]'); ylabel('|P1(f)|'); legend("US06", "NEDC", "UDDS", "Pulse 0.8C", "Pulse 0.1C", "ModifiedPulse")
title("Power Spectrums of Test Profiles")
% set(gca, 'YScale', 'log', 'XScale', 'log');
 set(gca, 'XScale', 'log');
 
fig_fft = figure; 
figure(fig_fft); hold on; 
plot(data{6, 2}.f, data{6, 2}.P1, 'LineWidth',0.1); 
plot(data{4, 2}.f, data{4, 2}.P1, 'LineWidth',0.1); 
plot(data{5, 2}.f, data{5, 2}.P1, 'LineWidth',0.1);
plot(data{2, 2}.f, data{2, 2}.P1, 'LineWidth',0.1); 
plot(data{1, 2}.f, data{1, 2}.P1, 'LineWidth',0.1); 
grid on; xlabel('f [Hz]'); ylabel('|P1(f)|'); legend("US06", "NEDC", "UDDS", "Pulse 0.8C", "Pulse 0.1C")
title("Power Spectrums of Test Profiles")
% set(gca, 'YScale', 'log', 'XScale', 'log');
set(gca, 'XScale', 'log');

fig_fft = figure; 
figure(fig_fft); hold on; 
plot(data{6, 2}.f, data{6, 2}.P1, 'color', '#0072BD', 'LineWidth',0.1); 
plot(data{5, 2}.f, data{5, 2}.P1, 'color', '#EDB120', 'LineWidth',0.1);
plot(data{2, 2}.f, data{2, 2}.P1, 'color', '#7E2F8E', 'LineWidth',0.1); 
plot(data{7, 2}.f, data{7, 2}.P1, 'color', '#77AC30', 'LineWidth',0.1); 
grid on; xlabel('f [Hz]'); ylabel('|P1(f)|'); legend("US06", "UDDS", "Pulse 0.8C", "ModifiedPulse")
title("Power Spectrums of Test Profiles")
% set(gca, 'YScale', 'log', 'XScale', 'log');
 set(gca, 'XScale', 'log');

figure; hold on
ax1 = subplot(2, 1, 1); plot(data{8, 2}.Time, data{8, 2}.Current); 
grid on; title("Current Profile of Modified Pulse Test"); xlabel("Time [s]"); ylabel("Current [A]")
ax2 = subplot(2, 1, 2); plot(data{8, 2}.Time, data{8, 2}.Current); 
grid on; title("Current Profile of Modified Pulse Test (Zoomed-in View)"); xlabel("Time [s]"); ylabel("Current [A]")
ylim([0,250]); 

%% Bode Plot
% clearvars -except data
% 
% for i = [1, 2, 5, 6]
%     Fs = 10; 
%     y = data{i, 2}.Current; 
% 
%     ind = find(isnan(y)); 
%     if ~isempty(ind)
%         y(ind) = []; 
%     end
% 
%     NFFT = length(y); 
%     Y = fft(y, NFFT); 
%     F = Fs*(0:NFFT)'/NFFT;
% 
%     Ymag = abs(Y);              % Magnitude of the FFT
%     Yangle = unwrap(angle(Y));  % Phase of the FFT
% 
%     figure(1)
%     ind = fix(NFFT/2) + 1; 
%     subplot(2,1,1); 
%     hold on; plot(F(1:ind), 20*log10(Ymag(1:ind)));
%     title('Magnitude response of the audio signal'); xlabel('Frequency in Hz'); ylabel('dB'); grid on; axis tight 
%     subplot(2,1,2); 
%     hold on; plot(F(1:ind),Yangle(1:ind));
%     title('Magnitude response of the audio signal'); xlabel('Frequency in Hz'); ylabel('radians'); grid on;axis tight
% end
% subplot(2,1,1); set(gca, 'XScale', 'log', 'YScale', 'linear'); legend("Pulse 0.1", "Pulse 0.8", "UDDS", "US06")
% subplot(2,1,2); set(gca, 'XScale', 'linear', 'YScale', 'linear');