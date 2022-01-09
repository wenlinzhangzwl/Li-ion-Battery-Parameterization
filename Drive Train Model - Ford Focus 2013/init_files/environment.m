%% ENVIRONMENT PARAMETERS

enviro.temp                 = 20;
enviro.pressure             = 1;
enviro.air_cap				= 1009;         % J/kgK  ave cp of air
enviro.gravity 			    = 9.81;         % kg/m-s^2
enviro.dens_air 			= 1.23;
enviro.boltzmann_cnst       = 1.38e-23;     %Boltzmann's constant in J/K 
enviro.air_mol_wght         = 28.97;        %g/mol, kg/kmol

env.init.percent_humidity_ambient = 20;% percentage relaticve humidity
Tambient_K = 295;
Pambient_Pa = 101325;
rho_amb = 0.001*(0.01*env.init.percent_humidity_ambient*18+(1-0.01*env.init.percent_humidity_ambient)*29)*Pambient_Pa/(8.31447*Tambient_K); %kg/m^3
g_per_mole_air = 29;
g_per_mole_H2O = 18;
UGC_Joules_per_mole_per_Kelvin = 8.31447;
