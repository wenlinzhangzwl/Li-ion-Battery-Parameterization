%% windsor motor 
% based on data from windsor motor testing
% data tested at 400V

%{ 
3/18/2019
changed the motor efficiency map to only include 400V (testing voltage),
one sheet of map instead of 2
changed respective voltage inputs to the efficiency map and trq-spd curve
increased trq-spd curve to max 150 Nm x 20.6 scale = 3090 Nm peak
decreased final drive ratio to 1:1, or 3:1
%}

% import files
load wmot_alum_mot_map
load wmot_copp_mot_map
load wmot_idx_spd
load wmot_idx_trq

%% constants

rpm_rads_conv= 1/60*2*pi;
motor.inertia               = 0.0226;% kg-m^2 
mot.plant.init.coeff_regen  = 0.99;

% change depending on vehicle *
mult_trq = 3090/156; % 250Nm for focus, peak 156Nm for windsor testing
motor.prated = 107000; % 107kW for focus, may have to change

%% define motor type, efficiency maps

switch motor_type
    case 1 % copper windings
        wmot_copp_map = [zeros(1,19); wmot_copp_map]; % add an extra row for trq=0
        wmot_copp_map(1,:) = wmot_copp_map(2,:);
        motor.eff_map(:,:,1) = wmot_copp_map'/100; % efficiency @ 400V
        motor.spd_ind = [0:(1275/18):1275] * rpm_rads_conv;
        motor.trq_ind = [0 wmot_idx_trq] * mult_trq; %add low torque row
%         motor.vol_ind = [400]; % testing voltage
    case 2 % aluminum windings
        wmot_alum_map = [zeros(1,19); wmot_alum_map]; % add an extra row for trq=0
        wmot_alum_map(1,:) = wmot_alum_map(2,:);
        motor.eff_map(:,:,1) = wmot_alum_map'/100; 
        motor.spd_ind = [0:(1275/18):1275] * rpm_rads_conv;
        motor.trq_ind = [0 wmot_idx_trq] * mult_trq;
%         motor.vol_ind = [400]; 
end

%% torque-speed curve

motor.torquespeed_RPMidx=[0:(1275/8):1275] * rpm_rads_conv;
% motor.torquespeed_Vidx = [400];
motor.torquespeed_curve400=[150,150,150,150,100,60,55,45,0] * mult_trq;

% no voltage dependencies
motor.torquespeed_curveTOT=[motor.torquespeed_curve400];

%% braking

motor.brake_idx=[0,10,20,30,40,50,60,70,80,90,100]/3.6; 
motor.brake_map=[0 0.5 0.65 0.7 0.75 0.85 0.9 0.95 1 1 1];