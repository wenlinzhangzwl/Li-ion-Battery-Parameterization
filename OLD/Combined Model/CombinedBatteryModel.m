function [SOCActual, VTerminalActual] = CombinedBatteryModel(BatteryParameters)
    % Simulates battery with a current input
    % Based on http://hevpdd.ca/module-3-battery-pack-modeling-and-control/

    %% Define Battery Parameters
    Current = BatteryParameters.I; 
    DeltaT = BatteryParameters.samplingTime;    % Sampling interval [s]
    Cn = BatteryParameters.Q*3600;              % Capacity [Amp*s]
    eta = BatteryParameters.eta;                % Efficiency [-]
    SOC = BatteryParameters.SOC_Actual(1);         % Initial SOC

    %% Initialize Model Parameters
    Rchg = 7.668896198272705e-04;     % Charging resistance [ohm]
    Rdch = 6.520629644393921e-04;     % Discharging resistance [ohm]
    K0 = 3.398194410490990;
    K1 = 8.961915969850354e-07;
    K2 = 0.067967350816727;
    K3 = 0.091428391337395;
    K4 = 1.798133850097663e-05;

    %% Run Actual Battery Model
    SOCActual = [];
    VTerminalActual = [];

    for k = 1:length(Current)
        % SOC Update
        U = Current(k); % Current at time k 
        CoffB3 = - (eta * DeltaT / Cn);
        SOC = SOC + (CoffB3 * U);

        %% Run Updated Model
        if U >= 0
            VTerminal = K0 - Rchg*U - K1/SOC - K2*SOC + K3*log(SOC) + K4*log(1-SOC);
        else
            VTerminal = K0 - Rdch*U - K1/SOC - K2*SOC + K3*log(SOC) + K4*log(1-SOC);
        end

        VTerminalActual = [VTerminalActual; VTerminal];
        SOCActual = [SOCActual; SOC];
    end

end