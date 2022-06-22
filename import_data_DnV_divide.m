clear
close all

folder_current = cd; 
folder_data = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\2 - Samsung 35E Cylindrical Battery Pack\Test Results\25 degC\";
addpath(folder_current, folder_data)

filename = "March14_35E Characterization - WZ - All Data";
data_raw = readtable(folder_data + filename + ".csv");

% Convert time stamp to double
Time = data_raw.Time - data_raw.Time(1); 
Time.Format = 's'; 
Time = seconds(Time); 

% Plot data
figure; 
ax1 = subplot(3, 1, 1); plot(Time, data_raw.Voltage, '.-'); ylabel("Voltage (V)"); grid on
ax2 = subplot(3, 1, 2); plot(Time, data_raw.Current, '.-'); ylabel("Current (A)"); grid on
ax3 = subplot(3, 1, 3); plot(Time, data_raw.Ah, '.-'); ylabel("Ah"); grid on
linkaxes([ax1, ax2, ax3], 'x')

%% Find capacity test
time1 = 12283; 
time2 = 3e4; 
ind1 = find(Time>=time1, 1);
ind2 = find(Time>=time2, 1);
ind_max = find(data_raw.Voltage(ind1:ind2) == max(data_raw.Voltage(ind1:ind2))) + ind1;
ind_min = find(data_raw.Voltage(ind1:ind2) == min(data_raw.Voltage(ind1:ind2))) + ind_max-2;

data = data_raw(ind_max:ind_min, :); 

% Plot data
figure; 
ax1 = subplot(3, 1, 1); plot(data.Time, data.Voltage, '.-'); ylabel("Voltage (V)"); grid on
ax2 = subplot(3, 1, 2); plot(data.Time, data.Current, '.-'); ylabel("Current (A)"); grid on
ax3 = subplot(3, 1, 3); plot(data.Time, data.Ah, '.-'); ylabel("Ah"); grid on
linkaxes([ax1, ax2, ax3], 'x')

writetable(data, folder_data+"CAP25.csv")

%% Find OCV test (discharge)
time1 = 44618; 
time2 = 99240; 
ind1 = find(Time>=time1, 1);
ind2 = find(Time>=time2, 1);
ind_max = ind1+1;
ind_min = ind2-1;

data = data_raw(ind_max:ind_min, :); 

% Plot data
figure; 
ax1 = subplot(3, 1, 1); plot(data.Time, data.Voltage, '.-'); ylabel("Voltage (V)"); grid on
ax2 = subplot(3, 1, 2); plot(data.Time, data.Current, '.-'); ylabel("Current (A)"); grid on
ax3 = subplot(3, 1, 3); plot(data.Time, data.Ah, '.-'); ylabel("Ah"); grid on
linkaxes([ax1, ax2, ax3], 'x')

writetable(data, folder_data+"OCVdch25.csv")

%% Find OCV test (charge)
time1 = 102842; 
time2 = 207134; 
ind1 = find(Time>=time1, 1);
ind2 = find(Time>=time2, 1);
ind_max = ind1+1;
ind_min = ind2-1;

data = data_raw(ind_max:ind_min, :); 

% Plot data
figure; 
ax1 = subplot(3, 1, 1); plot(data.Time, data.Voltage, '.-'); ylabel("Voltage (V)"); grid on
ax2 = subplot(3, 1, 2); plot(data.Time, data.Current, '.-'); ylabel("Current (A)"); grid on
ax3 = subplot(3, 1, 3); plot(data.Time, data.Ah, '.-'); ylabel("Ah"); grid on
linkaxes([ax1, ax2, ax3], 'x')

writetable(data, folder_data+"OCVchg25.csv")

%% Find pulse test
time1 = 207134; 
time2 = 360704; 
ind1 = find(Time>=time1, 1);
ind2 = find(Time>=time2, 1);
ind_max = ind1;
ind_min = ind2-1;

data = data_raw(ind_max:ind_min, :); 

% Plot data
figure; 
ax1 = subplot(3, 1, 1); plot(data.Time, data.Voltage, '.-'); ylabel("Voltage (V)"); grid on
ax2 = subplot(3, 1, 2); plot(data.Time, data.Current, '.-'); ylabel("Current (A)"); grid on
ax3 = subplot(3, 1, 3); plot(data.Time, data.Ah, '.-'); ylabel("Ah"); grid on
linkaxes([ax1, ax2, ax3], 'x')

writetable(data, folder_data+"PUL25_0p.csv")


%% Find drive cycle test
time1 = 369512; 
ind1 = find(Time>=time1, 1);
ind_max = ind1+1;

data = data_raw(ind_max:end, :); 

% Plot data
figure; 
ax1 = subplot(3, 1, 1); plot(data.Time, data.Voltage, '.-'); ylabel("Voltage (V)"); grid on
ax2 = subplot(3, 1, 2); plot(data.Time, data.Current, '.-'); ylabel("Current (A)"); grid on
ax3 = subplot(3, 1, 3); plot(data.Time, data.Ah, '.-'); ylabel("Ah"); grid on
linkaxes([ax1, ax2, ax3], 'x')

writetable(data, folder_data+"MIX25.csv")
