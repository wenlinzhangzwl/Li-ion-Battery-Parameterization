clear

%% Get test data & conditions
folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Prismatic Cell\Test Results\1 - Characterization test 1\';
file1 = append(folder, '12-11-20_08.01 1336_Charge3_HPPC.mat');
file2 = append(folder, 'OCV_SOC.mat');
load(file1);
load(file2);
samplingTime = 10; 

% Calculate other model parameters
Q = abs(meas.Ah(end) - meas.Ah(1));                                         % Total capacity
meas.SOC = 1 - (-meas.Ah)/Q;                                                % SOC
meas.OCV = interp1(OCV_SOC.SOC, OCV_SOC.OCV, meas.SOC, 'linear', 'extrap'); % OCV
meas.Vdyn = meas.Voltage - meas.OCV;                                        % Dynamic part of voltage

%% Define GA Optimization Options
global Batt
Batt = struct('RecordingTime', meas.Time(4000:36000), ...
              'I', -meas.Current(4000:36000),...
              'SOC_Actual', meas.SOC(4000:36000),...
              'V_Actual', meas.Voltage(4000:36000),...
              'samplingTime', samplingTime,...
              'Q', Q,...
              'eta', eta25);

clearvars -except Batt
nvars = 7; 

% Lower and upper bounds of varialbes
Rchg_lb = 0; Rchg_ub = 0.2; 
Rdch_lb = 0; Rdch_ub = 0.2; 
K0_lb = 3; K0_ub = 4; 
K1_lb = 0; K1_ub = 0.2; 
K2_lb = 0; K2_ub = 0.2; 
K3_lb = 0; K3_ub = 0.2; 
K4_lb = 0; K4_ub = 0.2;

lb = [Rchg_lb, Rdch_lb, K0_lb, K1_lb, K2_lb, K3_lb, K4_lb];
ub = [Rchg_ub, Rdch_ub, K0_ub, K1_ub, K2_ub, K3_ub, K4_ub];

% Other options
PopInitRange_Data = [lb; ub];
PopulationSize_Data = 30; 
Generations_Data = 500; 

%% Runs Optimizer
Rchg_init = 8.000000000000000e-04;
Rdch_init = 6.520480632781982e-04;
K0_init =  3.401679851698876; 
K1_init = 9.207639694213884e-05; 
K2_init = 0.070935855841637; 
K3_init = 0.093063778877259; 
K4_init = 1.798133850097663e-05; 

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
    SOCinit = Batt.SOC_Actual(1);  % Initial SOC of the battery
    
    [VTerminal_Optim, SOC_Optim] = Combined_Optimization(Opt_Param);
    
    error_VT = V - VTerminal_Optim;
    OBJ = sum(trapz(t, error_VT.^2));

    %% Calculate the Terminal Voltage RMSE
    RMSE_VT  = sqrt((sum((V(1:end)) - VTerminal_Optim(1:length(VTerminal_Optim))).^2))...
        /((length(V)-1));

    RMSE_SOC = sqrt((sum((SOC(1:end) - SOC_Optim).^2))...
        /((length(SOC)-1)));
end

function [VTerminal_Optim, SOC_Optim] = Combined_Optimization(param)
% Simulates battery with parameters from optimization

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
    global Batt
    deltaT = Batt.samplingTime;     % Sampling interval [s]
    Cn = Batt.Q*3600;               % Capacity [Amp*s]
    SOC = Batt.SOC_Actual(1);       % Battery initial SOC [N/A]
    eta = Batt.eta;                 % Efficiency [N/A]
	Current = Batt.I;

    % *************************************************************************************************
    
    %% Run Actual Battery Model
    SOC_Optim = [];
    VTerminal_Optim = [];

    for k = 1:length(Current)
        % Update SOC
        U = Current(k);                 % Current at time k 
        CoffB3 = (eta * deltaT) / Cn;
        SOC = SOC - (CoffB3 * U);       % New SOC

        % Run Updated Model
        if Current(k) >= 0
            VTerminal = K0 - Rchg*Current(k) - K1/SOC - K2*SOC + K3*log(SOC) + K4*log(1-SOC);
        else
            VTerminal = K0 - Rdch*Current(k) - K1/SOC - K2*SOC + K3*log(SOC) + K4*log(1-SOC);
        end
        
        % Export VTerminal_Optim & SOC_Optim
        if ~isnan(VTerminal) && isreal(VTerminal)
            VTerminal_Optim = [VTerminal_Optim; VTerminal];
        else
            VTerminal_Optim = [VTerminal_Optim; 1e100];
        end
%         VTerminal_Optim = [VTerminal_Optim; VTerminal];
        SOC_Optim = [SOC_Optim; SOC];
    end
    
%     figure; plot(Batt.RecordingTime, Current); grid on; title('current')
%     figure; plot(Batt.RecordingTime, VTerminal_Optim);grid on; title('voltage')
%     figure; plot(Batt.RecordingTime, SOC_Optim);grid on; title('soc')
%     hold on; plot(Batt.RecordingTime, Batt.SOC_Actual, 'r');
end