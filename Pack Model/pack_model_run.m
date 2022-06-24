clear

folder_current = cd; 
folder_scripts = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\Scripts\Pack Model";
addpath(folder_current, genpath(folder_scripts))

% Input
temp_ambient = 25; 
[Q_Ah, cell_parameters] = load_parameters("CylindricalCell", 80);

% Initial conditions
batt.Q_Ah = Q_Ah; 
batt.SOC_init = 1;
batt.efficiency = 1; 

% Load & scale current input
load("UDDS_10Hz.mat", "drivecycle_current"); 
drivecycle_current(:, 2) = -drivecycle_current(:, 2); % -ve is charge, +ve is discharge
C_rate = 1; % Maximum C rate
C_max = max(abs(drivecycle_current(:, 2)));
drivecycle_current(:, 2) = drivecycle_current(:, 2)*(C_rate * batt.Q_Ah)/C_max; 

% Model settings
sim_duration = drivecycle_current(end, 1);

% Run model
sim("pack_model.slx")