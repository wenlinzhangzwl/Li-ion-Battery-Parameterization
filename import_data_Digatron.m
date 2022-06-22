clear

folder = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\10 - Characterization\Back up\";
file = "PUL_25degC_0p8_1min"; 
data = readtable(folder + file + ".csv");

save_var = 1; 
Q_nom = 280; 

% Fields to convert
meas.Time = (data.ProgTime - data.ProgTime(1));
meas.Time = seconds(meas.Time);
meas.Step = data.Step; 
meas.Voltage = data.Voltage; 
meas.Current = data.Current; 
meas.Ah = data.Capacity; 
meas.Battery_Temp_degC = data.Temperature;

% Save as table
meas = struct2table(meas); 

% Delete beginning & end if necessary
% i_beg = 12; 
% i_end = 16; 
% ind = find(or(meas.Step == i_beg, meas.Step == i_end)); 
% meas(ind:end, :) = []; 

meas.Ah = meas.Ah - meas.Ah(1); 


%%
figure; 
subplot(5, 1, 1); plot(meas.Time, meas.Voltage, '.-'); title("Voltage"); xlabel("Time [s]"); ylabel("Voltage [V]"); grid on
subplot(5, 1, 2); plot(meas.Time, meas.Current, '.-'); title("Current"); xlabel("Time [s]"); ylabel("Current [A]"); grid on
subplot(5, 1, 3); plot(meas.Time, meas.Ah, '.-'); title("Ah"); xlabel("Time [s]"); ylabel("Ah"); grid on
subplot(5, 1, 4); plot(meas.Time, meas.Battery_Temp_degC, '.-'); title("Temperature"); xlabel("Time [s]"); ylabel("Ah"); grid on
subplot(5, 1, 5); plot(meas.Time, meas.Step, '.-'); title("Step"); xlabel("Time [s]"); ylabel("Ah"); grid on

Q = abs(max(meas.Ah) - min(meas.Ah));
C_max = max(abs(meas.Current))/Q_nom; 

%%
if save_var == 1
    filename = folder+file+".mat"; 
    save(filename, "meas");
end