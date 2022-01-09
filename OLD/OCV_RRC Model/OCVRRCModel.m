function [Vt, SOC] = OCVRRCModel(Batt)

%% Define Battery Parameters
% Model parameters
% results = [2.5e-4, 9e-5, 9e-5, 9e-5, 8e4, 8e4, 8e4];
results = [6.30476136720992e-06,0.000359823357828901,2.52306589774566e-05,0.000242966882953506,319629.003587578,253503.943186590,304976.715514364];
R0 = results(1);
R1 = results(2);
R2 = results(3);
R3 = results(4);
C1 = results(5);
C2 = results(6);
C3 = results(7);
tau1 = R1 * C1;
tau2 = R2 * C2;
tau3 = R3 * C3;


%% Simulate the Model

% Battery Fixed/Known Parameters 
Q = Batt.Q*3600;           % Capacity [Amp*s]
Current = Batt.Current;    % Current {A]

SOC = [Batt.SOC_exp(1)];
Vt  = [Batt.Voltage_exp(1)];
Irc1 = [0];
Irc2 = [0];
Irc3 = [0];

for k = 2:length(Current)
    
   oSOC = SOC(k-1);
   OCV = Batt.Voltage_exp(end);%pchip(Batt.OCVSOC.SOC, Batt.OCVSOC.OCV, oSOC);
   deltaT = Batt.Time(k) - Batt.Time(k-1);
   
   % Current through RC branchs
    oIrc1 =  ( 1 - exp(-deltaT/tau1) ) * Current(k) + exp(-deltaT/tau1) * Irc1(k-1);
    oIrc2 =  ( 1 - exp(-deltaT/tau2) ) * Current(k) + exp(-deltaT/tau2) * Irc2(k-1);
    oIrc3 =  ( 1 - exp(-deltaT/tau3) ) * Current(k) + exp(-deltaT/tau3) * Irc3(k-1);
    
    % Calculate current Vt & SOC
    oVt   = OCV + (R0  * Current(k)) + R1 * oIrc1 + R2 * oIrc2 + R3 * oIrc3;
    oSOC  = oSOC + (deltaT / Q)* Current(k);
    
    % Update parameters
    Irc1 = [Irc1; oIrc1];
    Irc2 = [Irc2; oIrc2];
    Irc3 = [Irc2; oIrc3];
    Vt  = [Vt; oVt];
    SOC = [SOC; oSOC];
end
