function [SOCActual, VTerminalActual] = CombinedBatteryModel(I, Time)
    % Simulates battery with a current input
    % Based on http://hevpdd.ca/module-3-battery-pack-modeling-and-control/

    %% Define Battery Parameters
    Current = I; 
    DeltaT = 0.1;    % Sampling interval [s]
    Cn = 5.4*3600;  % Capacity [Amp*s]
    eta = 1;        % Efficiency [-]
    SOC = 0.9;      % Initial SOC

    %% Initialize Model Parameters
    Rchg = 0.1;     % Charging resistance [ohm]
    Rdch = 0.1;     % Discharging resistance [ohm]
    K0 = 3;
    K1 = 0.01;
    K2 = 0.01;
    K3 = 0.01;
    K4 = 0.01;

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