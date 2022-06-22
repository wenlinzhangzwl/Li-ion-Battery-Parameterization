clear
close all
clc

filename = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\2 - Samsung 35E Cylindrical Battery Pack\Test Results\25 degC\PUL_25degC_0p2C_60min.mat";
load(filename);
meas = struct2table(meas); 

figure("WindowStyle","docked"); 
ax(1) = subplot(3, 1, 1); hold on; plot(meas.Time, meas.Voltage, '.-'); %plot(meas.Time(ind), meas.Voltage(ind), '*');
title("Voltage"); xlabel("Time [s]"); ylabel("Voltage [V]"); grid on
ax(2) = subplot(3, 1, 2); hold on; plot(meas.Time, meas.Current, '.-'); %plot(meas.Time(ind), meas.Current(ind), '*');
title("Current"); xlabel("Time [s]"); ylabel("Current [A]"); grid on
ax(3) = subplot(3, 1, 3); hold on; plot(meas.Time, meas.Ah, '.-'); %plot(meas.Time(ind), meas.Ah(ind), '*');
title("Ah"); xlabel("Time [s]"); ylabel("Ah"); grid on
linkaxes([ax], 'x')

% Change these for each test: 
I_threshold_low = 0.03;
I_threshold_high = 0.75;

%% Remove current offset
ind = find(abs(meas.Current) <= I_threshold_low);
offset = mean(meas.Current(ind));
meas.Current = meas.Current - offset; 

figure("WindowStyle","docked"); 
ax(1) = subplot(3, 1, 1); hold on; plot(meas.Time, meas.Voltage, '.-'); %plot(meas.Time(ind), meas.Voltage(ind), '*');
title("Voltage"); xlabel("Time [s]"); ylabel("Voltage [V]"); grid on
ax(2) = subplot(3, 1, 2); hold on; plot(meas.Time, meas.Current, '.-'); %plot(meas.Time(ind), meas.Current(ind), '*');
title("Current"); xlabel("Time [s]"); ylabel("Current [A]"); grid on
ax(3) = subplot(3, 1, 3); hold on; plot(meas.Time, meas.Ah, '.-'); %plot(meas.Time(ind), meas.Ah(ind), '*');
title("Ah"); xlabel("Time [s]"); ylabel("Ah"); grid on
linkaxes([ax], 'x')

%% Set current to 0 when it's below threshold
ind = find(abs(meas.Current) <= I_threshold_low);
meas.Current(ind) = 0; 

figure("WindowStyle","docked"); 
ax(1) = subplot(3, 1, 1); hold on; plot(meas.Time, meas.Voltage, '.-'); %plot(meas.Time(ind), meas.Voltage(ind), '*');
title("Voltage"); xlabel("Time [s]"); ylabel("Voltage [V]"); grid on
ax(2) = subplot(3, 1, 2); hold on; plot(meas.Time, meas.Current, '.-'); %plot(meas.Time(ind), meas.Current(ind), '*');
title("Current"); xlabel("Time [s]"); ylabel("Current [A]"); grid on
ax(3) = subplot(3, 1, 3); hold on; plot(meas.Time, meas.Ah, '.-'); %plot(meas.Time(ind), meas.Ah(ind), '*');
title("Ah"); xlabel("Time [s]"); ylabel("Ah"); grid on
linkaxes([ax], 'x')

%% Remove noise
current_abs = abs(meas.Current);
ind = find(and(current_abs>=I_threshold_low, current_abs<=I_threshold_high));
meas(ind, :) = []; 
ind = find(current_abs>=I_threshold_high+0.01);
meas(ind, :) = []; 

figure("WindowStyle","docked"); 
ax(1) = subplot(3, 1, 1); hold on; plot(meas.Time, meas.Voltage, '.-'); %plot(meas.Time(ind), meas.Voltage(ind), '*');
title("Voltage"); xlabel("Time [s]"); ylabel("Voltage [V]"); grid on
ax(2) = subplot(3, 1, 2); hold on; plot(meas.Time, meas.Current, '.-'); %plot(meas.Time(ind), meas.Current(ind), '*');
title("Current"); xlabel("Time [s]"); ylabel("Current [A]"); grid on
ax(3) = subplot(3, 1, 3); hold on; plot(meas.Time, meas.Ah, '.-'); %plot(meas.Time(ind), meas.Ah(ind), '*');
title("Ah"); xlabel("Time [s]"); ylabel("Ah"); grid on
linkaxes([ax], 'x')

%% Redo Coulomb counting
deltaT = meas.Time(2:end) - meas.Time(1:end-1);
deltaT = [0.1; deltaT];
Ah = deltaT .* meas.Current; 
meas.Ah = cumsum(Ah);

figure("WindowStyle","docked"); 
ax(1) = subplot(3, 1, 1); hold on; plot(meas.Time, meas.Voltage, '.-'); %plot(meas.Time(ind), meas.Voltage(ind), '*');
title("Voltage"); xlabel("Time [s]"); ylabel("Voltage [V]"); grid on
ax(2) = subplot(3, 1, 2); hold on; plot(meas.Time, meas.Current, '.-'); %plot(meas.Time(ind), meas.Current(ind), '*');
title("Current"); xlabel("Time [s]"); ylabel("Current [A]"); grid on
ax(3) = subplot(3, 1, 3); hold on; plot(meas.Time, meas.Ah, '.-'); %plot(meas.Time(ind), meas.Ah(ind), '*');
title("Ah"); xlabel("Time [s]"); ylabel("Ah"); grid on
linkaxes([ax], 'x')
