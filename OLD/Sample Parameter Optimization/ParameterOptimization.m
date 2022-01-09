
clear

%% Get test data
% folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Prismatic Cell\Test Results\1 - Characterization test 1\';
% file = append(folder, '12-11-20_08.01 1336_Charge3_HPPC.mat');
% 
% load(file);
% 
% global Batt
% Batt = struct('RecordingTime', meas.Time, ...
%               'I', meas.Current,...
%               'SOC_Actual', meas.SOC,...
%               'V_Actual', meas.Voltage);

folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Prismatic Cell\MATLAB\';
file = append(folder, 'SampleData.mat');
load(file);

global Batt;

%% Define GA Optimization Options
nvars = 7; 

% Lower and upper bounds of varialbes
Rchg_lb = 0; Rchg_ub = 0.01; 
Rdch_lb = 0; Rdch_ub = 0.01; 
K0_lb = 2; K0_ub = 4; 
K1_lb = 0; K1_ub = 0.1; 
K2_lb = 0; K2_ub = 0.1; 
K3_lb = 0; K3_ub = 0.1; 
K4_lb = 0; K4_ub = 0.1;

lb = [Rchg_lb, Rdch_lb, K0_lb, K1_lb, K2_lb, K3_lb, K4_lb];
ub = [Rchg_ub, Rdch_ub, K0_ub, K1_ub, K2_ub, K3_ub, K4_ub];

% Other options
PopInitRange_Data = [lb; ub];
PopulationSize_Data = 10; 
Generations_Data = 20; 

%% Runs Optimizer
Rchg_init = 0.0006408;
Rdch_init = 0.009884;
K0_init = 0.0001; 
K1_init = 0.0001; 
K2_init = 0.0001; 
K3_init = 0.0001; 
K4_init = 0.0001; 

InitialPopulation_Data = [Rchg_init, Rdch_init, K0_init, K1_init, K2_init, K3_init, K4_init];

[x, fval, exitflag, output, population, score] = ...
    Run_GA(nvars, lb, ub, PopInitRange_Data, PopulationSize_Data, Generations_Data, InitialPopulation_Data);

%% Functions

function [x, fval, exitflag, output, population, score] = ...
    Run_GA(nvars, lb, ub, PopInitRange_Data, PopulationSize_Data, Generations_Data, InitialPopulation_Data)
% Runs the GA optimization function ga() of MATLAB

    % Start with default options
    options = optimoptions('ga'); 
    
    % Modify default options
    options = optimoptions(options, 'PopInitRange', PopInitRange_Data); 
    options = optimoptions(options, 'PopulationSize', PopulationSize_Data); 
    options = optimoptions(options, 'Generations', Generations_Data); 
    options = optimoptions(options, 'InitialPopulation', InitialPopulation_Data);
    options = optimoptions(options, 'Display', 'diagnose');
    options = optimoptions(options, 'PlotFcns', {@gaplotbestf @gaplotbestindiv});
    options = optimoptions(options, 'Vectorized', 'off');
    options = optimoptions(options, 'UseParallel', 'never');
    
    % Run ga() with objective function "Main_Combined_Optimization"
    [x, fval, exitflag, output, population, score] = ...
        ga(@Main_Combined_Optimization, nvars, [], [], [], [], lb, ub, [], [], options);

end

function [VTerminal_Optim, SOC_Optim] = Combined_Optimization(Current, param, SOC_init)
% Simulates battery with parameters from optimization

    %% Define Parameters
    % Optimization Parameters
    Rchg = param(1);     % Charging resistance [ohm]
    Rdch = param(2);     % Discharging resistance [ohm]
    K0 = param(3); 
    K1 = param(4); 
    K2 = param(5); 
    K3 = param(6); 
    K4 = param(7); 
    
    % *************************************************************************************************
    % Battery Fixed/Known Parameters    
    DeltaT = 0.1;   % Sampling interval [s]
    Cn = 5.4*3600;  % Capacity [Amp*s]
    X = SOC_init;   % Battery initial SOC [N/A]
    eta = 1;        % Efficiency [N/A]
    % *************************************************************************************************
    
    %% Run Actual Battery Model
    SOC_Optim = [];
    VTerminal_Optim = [];
    SOC = X; 

    for k = 1:length(Current)
        % Update SOC
        U = Current(k);                 % Current at time k 
        CoffB3 = - (eta * DeltaT / Cn); % Coefficient
        SOC = SOC + (CoffB3 * U);       % New SOC

        % Run Updated Model
        if U >= 0
            VTerminal = K0 - Rchg*U - K1/SOC - K2*SOC + K3*log(SOC) + K4*log(1-SOC);
        else
            VTerminal = K0 - Rdch*U - K1/SOC - K2*SOC + K3*log(SOC) + K4*log(1-SOC);
        end
        
        % Export VTerminal_Optim & SOC_Optim
        if ~isnan(VTerminal) && isreal(VTerminal)
            VTerminal_Optim = [VTerminal_Optim; VTerminal];
        else
            VTerminal_Optim = [VTerminal_Optim; 1e10];
        end
        SOC_Optim = [SOC_Optim; SOC];
    end
    
    
end

function [OBJ] = Main_Combined_Optimization(Opt_Param)
% Objective function for optimization
    
    %% Run Model with Optimization Parameters
    global Batt; 
    format long
    
    t = Batt.RecordingTime;
    I = Batt.I; 
    SOC = Batt.SOC_Actual; 
    V = Batt.V_Actual;
    
    %% Call the Optimizer
    SOCinit = 0.9;  % Initial SOC of the battery
    
    [VTerminal_Optim, SOC_Optim] = Combined_Optimization(I, Opt_Param, SOCinit);
    
    error_VT = V -VTerminal_Optim;
    OBJ = sum(trapz(t, error_VT.^2));

    %% Plot Functions
%     figure
%     plot(V);
%     grid on
% 
%     hold all
%     plot(TerminalVoltage_Optimization(1:length(TerminalVoltage_Optimization)));
%     ylabel('Terminal Voltage (V)')
%     xlabel('Time (sec)')
%     legend('Experimental Battery', 'Optimized Model')
%     grid on
% 
%     magnifyOnFigure
    % figure
    % plot(SOC*100)
    % hold all
    % plot(SOC_Optimization(1:length(SOC_Optimization))*100)
    % grid on 

    %% Calculate the Terminal Voltage RMSE
    RMSE_VT  = sqrt((sum((V(1:end)) - VTerminal_Optim(1:length(VTerminal_Optim))).^2))...
        /((length(V)-1));

    RMSE_SOC = sqrt((sum((SOC(1:end) - SOC_Optim).^2))...
        /((length(SOC)-1)));
end