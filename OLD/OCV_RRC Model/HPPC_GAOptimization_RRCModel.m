clear

%% Get test data & conditions
% Load files
currentFolder = cd; 
addpath(currentFolder);
folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Prismatic Cell\Test Results\1 - Characterization test 1\';
cd(folder);
load('12-11-20_08.01 1336_Charge3_HPPC.mat');
load('OCV_SOC.mat');
load('HPPC_GetResistance');

% Calculate model parameters 
Q = abs(meas.Ah(end) - meas.Ah(1));     % Total capacity
meas.SOC = 1 - (-meas.Ah)/Q;            % SOC

%% Define GA Optimization Options
global Batt
Batt = struct('Time', meas.Time, ...
              'Current', -meas.Current,...
              'SOC_exp', meas.SOC,...
              'Vt_exp', meas.Voltage,...
              'Q', Q, ...
              'OCVSOC', OCV_SOC);
Batt.R0 = Resistance;

% Lower and upper bounds of varialbes
R1_lb = 0; R1_ub = 0.001; 
R2_lb = 0; R2_ub = 0.001; 
C2_lb = 0; C2_ub = 4;
C1_lb = 0; C1_ub = 4;

lb = [R1_lb, R2_lb, C1_lb, C2_lb];
ub = [R2_ub, R2_ub, C1_ub, C2_ub];

clearvars -except Batt lb ub

% Other options
nvars = 4; 
PopInitRange_Data = [lb; ub];
PopulationSize_Data = 50; 
Generations_Data = 2000; 

%% Runs Optimizer
R1_init = 2e-4;
R2_init = 2e-4; 
C2_init = 1; 
C1_init = 1; 

InitialPopulation_Data = [R1_init, R2_init, C1_init, C2_init];

[results, fval, exitflag, output, population, score] = Run_GA(nvars, lb, ub, PopInitRange_Data, PopulationSize_Data, Generations_Data, InitialPopulation_Data);

%% Functions

function [x, fval, exitflag, output, population, score] = ...
    Run_GA(nvars, lb, ub, PopInitRange_Data, PopulationSize_Data, Generations_Data, InitialPopulation_Data)
% Runs the GA optimization function ga() of MATLAB

    % Start with default options
    options = optimoptions('ga'); 
    
    % Modify default options
    options = optimoptions(options, 'PopInitRange', PopInitRange_Data,...
                                    'PopulationSize', PopulationSize_Data,...
                                    'Generations', Generations_Data,...
                                    'InitialPopulation', InitialPopulation_Data,...
                                    'Display', 'diagnose',...
                                    'PlotFcns', {@gaplotbestf @gaplotbestindiv},...
                                    'Vectorized', 'off',...
                                    'UseParallel', 'never',...
                                    'FitnessLimit', 1);
    
    % Run ga() with objective function "Main_Combined_Optimization"
    [x, fval, exitflag, output, population, score] = ...
        ga(@RunOCVRRC_Optimization, nvars, [], [], [], [], lb, ub, [], [], options);

end

function [OBJ] = RunOCVRRC_Optimization(Opt_Param)
% Objective function for optimization
    
    %% Run Model with Optimization Parameters
    global Batt;   
    t = Batt.Time;
    I = Batt.Current; 
    SOC_exp = Batt.SOC_exp; 
    Vt_exp = Batt.Vt_exp;
    
    %% Call the Optimizer
    SOCinit = SOC_exp(1);  % Initial SOC of the battery
    
    [Vt_sim, SOC_sim] = OCVRRCModel_Optimization(Opt_Param);
    error_Vt = Vt_exp - Vt_sim;
    OBJ = sum(trapz(t, error_Vt.^2));

    %% Calculate the Terminal Voltage RMSE
    RMSE_Vt  = sqrt((sum((Vt_exp(1:end)) - Vt_sim(1:length(Vt_sim))).^2)) / ((length(Vt_exp)-1));
    RMSE_SOC = sqrt((sum((SOC_exp(1:end) - SOC_sim).^2)) / ((length(SOC_exp)-1)));
end

function [Vt, SOC] = OCVRRCModel_Optimization(param)
% Simulates battery with parameters from optimization

    % Optimization Parameters
    R1 = param(1);
    R2 = param(2);
    C1 = param(3); 
    C2 = param(4); 
    
    % Other parameters
    tau1        = R1 * C1;
    tau2        = R2 * C2 ;
    Irc1_old    = 0;       % Current through the RC branch at t = k-1
    Irc1        = 0;       % Current through the RC branch at t = k
    Irc2_old    = 0;       % Current through the RC branch at t = k-1
    Irc2        = 0;       % Current through the RC branch at t = k

    % Battery Fixed/Known Parameters 
    global Batt
    Q           = Batt.Q*3600;           % Capacity [Amp*s]
    I           = Batt.Current;          % Current {A]
    
    % Create R0 look up tables
    R0 = Batt.R0;
    Size = size(R0);
    for i = 1:Size(1)
        R0SOC(i,1) = R0{i, 1};
        R0chg(i,1) = R0{i, 3};
        R0dch(i,1) = R0{i, 8};
    end
    
    %% Run Actual Battery Model
    SOC = [Batt.SOC_exp(1)];
    Vt  = [Batt.Vt_exp(1)];

for k      = 2:length(I)
    
   oSOC    = SOC(k-1);
   oOCV    = pchip(Batt.OCVSOC.SOC, Batt.OCVSOC.OCV, oSOC);
   deltaT  = Batt.Time(k) - Batt.Time(k-1);
   
   % Current through RC branchs
    Irc1 =  ( 1 - exp(-deltaT/tau1) ) * I(k) + exp(-deltaT/tau1) * Irc1_old;
    Irc1_old   = Irc1; 
    Irc2 =  ( 1 - exp(-deltaT/tau2) ) * I(k) + exp(-deltaT/tau2) * Irc2_old;
    Irc2_old   = Irc1; 
   
    % Terminal voltage 
    if I(k) >= 0      % Discharging
        R0 = pchip(R0SOC, R0dch, oSOC);
        oVt   = oOCV - (R0  * I(k)) - R1 * Irc1 - R2 * Irc2;      
    elseif I(k)<0     % Charging
        R0 = pchip(R0SOC, R0chg, oSOC);
        oVt   = oOCV - (R0  * I(k)) - R1 * Irc1 - R2 * Irc2; 
    end
    
    % SOC Update
    oSOC  = oSOC - (deltaT / Q)* I(k);                  

    Vt  = [Vt; oVt];
    SOC = [SOC; oSOC];

end

%     figure; plot(Batt.Time, I); grid on; title('current')
%     figure; plot(Batt.Time, Vt);grid on; title('voltage')
%     figure; plot(Batt.Time, SOC);grid on; title('soc')
%     hold on; plot(Batt.Time, Batt.SOC_exp, 'r');
end