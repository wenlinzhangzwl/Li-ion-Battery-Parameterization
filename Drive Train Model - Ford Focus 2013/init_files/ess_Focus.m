%% LOAD DATA

% Focus config is 86S5P
% total energy is 75 Ah
% total capacity for EV is 23kWh

load focus_rchg_pack_v2
load focus_rdis_pack_v2
load focus_sococv_pack_v2

%% BATTERY PARAMETERS

nser = 86;
npar = 5;
ess.soc_min             = 0.15;
ess.soc_max 			= 0.90;
ess.eff                 = 0.99;

% soc vs ocv
ess.soc_idx = [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]/100;
ess.idx_soc = ess.soc_idx;
ess.ocv_cell = focus_sococv_pack(2,:)/nser;

% resistance
ess.rdis.idx_soc = ess.soc_idx;
ess.rdis_cell= focus_rdis_pack_v2(2,:)/nser/npar;
ess.rchg.idx_soc = ess.soc_idx;
ess.rchg_cell= focus_rchg_pack_v2(2,:)/nser/npar;

% capacity for the battery
ess.cap = 75;
