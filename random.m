%% Plot single file
% clear
% Q = 277.7845;

% UDDS = load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\6 - Drive cycle\UDDS25.mat');
% UDDS.meas_t = struct2table(UDDS.meas);
% ind = find(UDDS.meas.Voltage <=2.5, 1);
% UDDS.meas_t = UDDS.meas_t(1:ind, :);
% 
% UDDS.meas_t.Current = -UDDS.meas_t.Current;
% UDDS.meas_t.Ah = UDDS.meas_t.Ah - UDDS.meas_t.Ah(1);
% UDDS.meas_t.Time = UDDS.meas_t.Time - UDDS.meas_t.Time(1);
% UDDS.meas_t.SOC = 1 - (-UDDS.meas_t.Ah)/Q;
% 
% HWFET = load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\6 - Drive cycle\HWFET25.mat');
% HWFET.meas_t = struct2table(HWFET.meas);
% ind = find(HWFET.meas.Voltage <=2.5, 1);
% HWFET.meas_t = HWFET.meas_t(1:ind, :);
% 
% HWFET.meas_t.Current = -HWFET.meas_t.Current;
% HWFET.meas_t.Ah = HWFET.meas_t.Ah - HWFET.meas_t.Ah(1);
% HWFET.meas_t.Time = HWFET.meas_t.Time - HWFET.meas_t.Time(1);
% HWFET.meas_t.SOC = 1 - (-HWFET.meas_t.Ah)/Q;
% 
% US06 = load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\6 - Drive cycle\US0625.mat');
% US06.meas_t = struct2table(US06.meas);
% ind = find(US06.meas.Voltage <=2.5, 1);
% US06.meas_t = US06.meas_t(1:ind, :);
% 
% US06.meas_t.Current = -US06.meas_t.Current;
% US06.meas_t.Ah = US06.meas_t.Ah - US06.meas_t.Ah(1);
% US06.meas_t.Time = US06.meas_t.Time - US06.meas_t.Time(1);
% US06.meas_t.SOC = 1 - (-US06.meas_t.Ah)/Q;
% US06_1 = US06; 
% 
% 
% US06 = load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\8 - Capacity test & US06\US0625.mat');
% US06.meas_t = struct2table(US06.meas);
% ind = find(US06.meas.Voltage <=2.5, 1);
% US06.meas_t = US06.meas_t(1:ind, :);
% 
% US06.meas_t.Current = -US06.meas_t.Current;
% US06.meas_t.Ah = US06.meas_t.Ah - US06.meas_t.Ah(1);
% US06.meas_t.Time = US06.meas_t.Time - US06.meas_t.Time(1);
% US06.meas_t.SOC = 1 - (-US06.meas_t.Ah)/Q;
% US06_2 = US06; 
% 
% figure; hold on
% plot(US06_1.meas_t.Time, US06_1.meas_t.Voltage)
% plot(US06_2.meas_t.Time, US06_2.meas_t.Voltage)
% plot(UDDS.meas_t.Time, UDDS.meas_t.Voltage)
% plot(HWFET.meas_t.Time, HWFET.meas_t.Voltage)
% legend('US06_1', 'US06_2', 'UDDS', 'HWFET'); grid on

% figure; 
% ax1 = subplot(3,1,1); 
% plot(UDDS.meas_t.Time, UDDS.meas_t.Voltage); grid on; title('voltage'); xlabel('SOC'); ylabel('Voltage (V)')
% ax2 = subplot(3,1,2); 
% plot(HWFET.meas_t.Time, HWFET.meas_t.Voltage); hold on; grid on; title('voltage'); xlabel('SOC'); ylabel('Voltage (V)')
% ax3 = subplot(3,1,3); 
% plot(US06.meas_t.Time, US06.meas_t.Voltage); hold on; grid on; title('voltage'); xlabel('SOC'); ylabel('Voltage (V)')
% linkaxes([ax1, ax2, ax3], 'x')
% 
% figure; 
% ax1 = subplot(3,1,1); 
% plot(UDDS.meas_t.Time, UDDS.meas_t.Current); grid on; title('Current'); xlabel('SOC'); ylabel('Current (A)')
% ax2 = subplot(3,1,2); 
% plot(HWFET.meas_t.Time, HWFET.meas_t.Current); hold on; grid on; title('Current'); xlabel('SOC'); ylabel('Current (A)')
% ax3 = subplot(3,1,3); 
% plot(US06.meas_t.Time, US06.meas_t.Current); hold on; grid on; title('Current'); xlabel('SOC'); ylabel('Current (A)')
% linkaxes([ax1, ax2, ax3], 'x')


% % Docked figures
% figure('WindowStyle', 'docked')
% plot(UDDS.meas.Time, UDDS.meas.Voltage);
% figure('WindowStyle', 'docked')
% plot(meas.Time, meas.Current);


% % Trim data so it starts at desired point
% stepStart = 24;
% stepEnd = 33;
% meas_t = struct2table(meas);
% meas_t.SOC = 1 - (-meas_t.Ah)/Q;
% meas_t = meas_t(and(meas.Step >= stepStart, meas.Step <= stepEnd), :);

% figure; 
% subplot(2,1,1); 
% plot(meas_t.Time, meas_t.Current); grid on; title('current'); xlabel('Time'); ylabel('Current (A)') %ylim([-250,0]);
% subplot(2,1,2); 
% plot(meas_t.Time, meas_t.Voltage); hold on; grid on; title('voltage'); xlabel('Time'); ylabel('Voltage (V)')

%% Process capacity tests
% 
% folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\2 - Capacity test';
% cd(folder);
% addpath(folder);   
% files = dir('*.mat');
% 
% for i = 1:length(files)
%     filename = append(folder, '\', files(i).name);
%     load(filename);
%     files(i).meas = meas;
% %     figure; 
% %     subplot(3,1,1); 
% %     plot(meas.Time/3600, meas.Current); grid minor; title('current');
% %     subplot(3,1,2); 
% %     plot(meas.Time/3600, meas.Voltage); hold on; grid minor; title('voltage')
% %     subplot(3,1,3); 
% %     plot(meas.Time/3600, meas.Ah); hold on; grid minor; title('Ah')
%     files(i).Q = abs(meas.Ah(end) - meas.Ah(1));
% end

%% Compare SOCs
% file1 = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\1 - Characterization test 1\OCV_SOC';
% file2 = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\3 - Characterization test 2\OCV_SOC';
% file3 = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\4 - Characterization test 3\OCV_SOC';
% OCV1 = load(file1);
% OCV2 = load(file2);
% OCV3 = load(file3);
% 
% figure
% plot(OCV1.OCV_SOC.SOC*100, OCV1.OCV_SOC.OCV);
% hold on; plot(OCV2.OCV_SOC.SOC*100, OCV2.OCV_SOC.OCV);
% hold on; plot(OCV3.OCV_SOC.SOC*100, OCV3.OCV_SOC.OCV);
% grid minor; legend('OCV1', 'OCV2', 'OCV3');

%% Plot
% clear
% currentfolder = cd; 
% addpath(currentfolder); 
% 
% cd('C:\Users\Wenlin\OneDrive\SCHOOL\Samsung INR18650-35E Cell\Testing')
% % C20CCCV = load('C20CCCV');
% % CAP25 = load('CAP25');
% % dch_OCV = load('dch_OCV');
% PUL25dch = load('PUL25dch');
% PUL25chg = load('PUL25chg');
% 
% figure; 
% ax1 = subplot(2, 1, 1); plot(PUL25dch.meas.Time, PUL25dch.meas.Current); title('current'); grid on
% ax2 = subplot(2, 1, 2); plot(PUL25dch.meas.Time, PUL25dch.meas.Voltage); title('voltage'); grid on
% linkaxes([ax1, ax2], 'x');
% 
% figure; 
% ax1 = subplot(2, 1, 1); hold on; plot(PUL25chg.meas.Time, PUL25chg.meas.Current); title('current'); grid on
% ax2 = subplot(2, 1, 2); hold on; plot(PUL25chg.meas.Time, PUL25chg.meas.Voltage); title('voltage'); grid on
% linkaxes([ax1, ax2], 'x')

%% Plot pulse test
% clear
% currentfolder = cd; 
% addpath(currentfolder); 

% file = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\10degC\MIX_10degC_1p0.mat"; 
% load(file);

% meas = rmfield(meas, {'TimeStamp', 'StepTime', 'Procedure', 'Wh', 'Power'});

Q = abs(max(meas.Ah) - min(meas.Ah));
C_max = max(abs(meas.Current))/280; 

figure("WindowStyle","docked"); 
ax(1) = subplot(4, 1, 1); hold on; plot(meas.Time, meas.Voltage, '.-'); %plot(meas.Time(ind), meas.Voltage(ind), '*');
title("Voltage"); xlabel("Time [s]"); ylabel("Voltage [V]"); grid on

ax(2) = subplot(4, 1, 2); hold on; plot(meas.Time, meas.Current, '.-'); %plot(meas.Time(ind), meas.Current(ind), '*');
title("Current"); xlabel("Time [s]"); ylabel("Current [A]"); grid on

ax(3) = subplot(4, 1, 3); hold on; plot(meas.Time, meas.Ah, '.-'); %plot(meas.Time(ind), meas.Ah(ind), '*');
title("Ah"); xlabel("Time [s]"); ylabel("Ah"); grid on
linkaxes([ax], 'x')

% subplot(4, 1, 4)
% plot(meas.Time, meas.Battery_Temp_degC); 
% title("Temperature"); xlabel("Time [s]"); ylabel("Ah"); grid on
subplot(4, 1, 4)
plot(meas.Time, meas.Step); 
title("Step"); xlabel("Time [s]"); ylabel("Ah"); grid on

%%
load("C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\10 - Characterization\PUL_25degC_0p8_30min.mat")
figure; 
ax(1) = subplot(2, 1, 1); hold on; plot(meas.Time, meas.Voltage);
title("Voltage"); xlabel("Time [s]"); ylabel("Voltage [V]"); grid on
ax(2) = subplot(2, 1, 2); hold on; plot(meas.Time, meas.Current);
title("Current"); xlabel("Time [s]"); ylabel("Current [A]"); grid on
print("mixed_drive_cycle", "-dpng")

%% FFT
% clear
% Fs = 1000;
% T = 1/Fs; 
% L = 1500; 
% t = (0:L-1)'*T; % Time vector
% 
% S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
% X = S + 2*randn(size(t)); % Signal to take fft of
% 
% Y = fft(X);
% 
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1); % single sided magnitude
% P1(2:end-1) = 2*P1(2:end-1);
% 
% figure; plot(t, P2); grid on
% figure; plot(t(1:height(P1)), P1); grid on

