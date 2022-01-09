function [SOC, Vt] = OCVRRCModel(Batt)

%% Define Battery Parameters
% Model parameters
R0      = 0.0096;
R1      = 0.0049;
C1      = 3860.14;
tau1    = C1 * R1;
R2      = 0.0049;
C2      = 3860.14;
tau2    = C2 * R2;
 
Irc1_old = 0;   % Current through the RC branch at t = k-1
Irc1 = 0;       % Current through the RC branch at t = k
Irc2_old = 0;   % Current through the RC branch at t = k-1
Irc2 = 0;       % Current through the RC branch at t = k
 

%% Define Battery Fixed/known Parameteres 
deltaT      = 10;
Q           = Batt.Q * 3600;
Current     = Batt.Current;
eta         = 0.99; 

%% Simulate the Model
SOC = [Batt.SOC_exp(1)];
Vt  = [Batt.V_exp(1)];

for k      = 2 : 1 : length(Current)
    
   oSOC    = SOC(k-1);
   oOCV     = pchip(Batt.OCVSOC.SOC, Batt.OCVSOC.OCV, oSOC);
   
   % Current through RC branchs
    Irc1 =  ( 1 - exp(-deltaT/tau1) ) * Current(k) + exp(-deltaT/tau1) * Irc1_old;
    Irc1_old   = Irc1; 
    Irc2 =  ( 1 - exp(-deltaT/tau2) ) * Current(k) + exp(-deltaT/tau2) * Irc2_old;
    Irc2_old   = Irc1; 
   
    oVt   = oOCV - (R0  * Current(k)) - R1 * Irc1 - R2 * Irc2;      % Terminal voltage 
    oSOC  = oSOC - (eta * deltaT / Q)* Current(k);                  % SOC Update

    Vt  = [Vt; oVt];
    SOC = [SOC; oSOC];

end

