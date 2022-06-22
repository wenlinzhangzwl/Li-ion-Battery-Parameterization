clear

folder_current = cd; 
folder_data = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\2 - Samsung 35E Cylindrical Battery Pack\Test Results\25 degC\Back up\";
addpath(folder_current, folder_data)

% Load single file
filename = "April 5_35E Characterization 10 degC - WZ - Mixed drive cycle";
data = readtable(folder_data + filename + ".csv");

% Test parameters
deltaT = 0.1; 
% testType = "pulse"; % OCV, pulse, drivecycle, capacity

%%
% Convert time stamp to double
Time = data.Time - data.Time(1); 
Time.Format = 's'; 
Time = seconds(Time); 
if deltaT == 1
    meas.Time = Time; 
elseif deltaT == 0.1
    meas.Time = [0:0.1:height(Time)/10-0.1]';
end

meas.Voltage = data.Voltage; 
meas.Current = data.Current; 
meas.Ah = data.Ah;
meas.Battery_Temp_degC = data.Temperature; 

Q = abs(meas.Ah(end) - meas.Ah(1));

% % Set steps
% if testType == "pulse"
%     
% end

%%
% Plot data
figure; 
ax1 = subplot(3, 1, 1); plot(meas.Time, meas.Voltage, '.-'); ylabel("Voltage (V)"); grid on
ax2 = subplot(3, 1, 2); plot(meas.Time, meas.Current, '.-'); ylabel("Current (A)"); grid on
ax3 = subplot(3, 1, 3); plot(meas.Time, meas.Ah, '.-'); ylabel("Ah"); grid on
linkaxes([ax1, ax2, ax3], 'x')

%%
% Save data
save(folder_data + filename + ".mat", "meas", '-mat')