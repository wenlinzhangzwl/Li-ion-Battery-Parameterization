%% run drivecycle for focus

clear all
% close all

%{
3/18/2019
added the driven trq-spd curve overlaid over the motor efficiency map
changed the vehicle to direct drive with 1:1 fd
%}

%% DEFINE CONSTANTS

mph_to_ms = 0.4470;

%% LOAD AND DEFINE TEST CYCLES

load UDDS
sch_cycle_city = sch_cycle;
sch_grade_city = sch_grade;
load HWFET
sch_cycle_hwy = sch_cycle;
sch_grade_hwy = sch_grade;
load US06
sch_cycle_us06 = sch_cycle;
sch_grade_us06 = sch_grade;
load NEDC
sch_cycle_NEDC = sch_cycle;
sch_grade_NEDC = sch_grade;

sim_architecture = 'FordFocus2013_DirectDrive';

% CHOOSE TEST TO RUN
test = 3; %****************************************************change
%{
test1 = UDDS, city test
test2 = HWFET, hwy test
test3 = US06, aggressive
test4 = NEDC, New European Drive Cycle
%}

switch test
    case 1
        drivecycle_speed = [sch_cycle_city(:,1) sch_cycle_city(:,2)]; 
        drivecycle_grade = sch_grade_city;
        drivecycle_name = "UDDS";
    case 2
        drivecycle_speed = [sch_cycle_hwy(:,1) sch_cycle_hwy(:,2)]; 
        drivecycle_grade = sch_grade_hwy;
        drivecycle_name = "HWFET";
    case 3
        drivecycle_speed = [sch_cycle_us06(:,1) sch_cycle_us06(:,2)]; 
        drivecycle_grade = sch_grade_us06;
        drivecycle_name = "US06";
    case 4
        drivecycle_speed = [sch_cycle_NEDC(:,1) sch_cycle_NEDC(:,2)]; 
        drivecycle_grade = sch_grade_NEDC;
        drivecycle_name = "NEDC";
        
%         nedc = readtable('nedc.csv');
%         nedc.start_velocity = nedc.start_velocity * 1000 / 3600;
%         nedc.end_velocity = nedc.end_velocity * 1000 / 3600;
%         speed = [];
%         time = [];
%         deltaT = 1; 
%         for i = 1:height(nedc)
%             
%             if i == 1
%                 t = [0];
%             else
%                 t = [sum(nedc.duration(1:i-1)) + deltaT];
%             end
%             
%             v_start = nedc.start_velocity(i); 
%             v_end = nedc.end_velocity(i); 
%             a = nedc.acceleration(i);
%             
%             v = [v_start];
%             ind = 1;
%             while t(end) < sum(nedc.duration(1:i))-0.01
%                 if v < v_end
%                     v_inst = v(ind) + a * deltaT;
%                 else
%                     v_inst = v(ind);
%                 end
%                 ind = ind +1;
%                 v = [v; v_inst];
%                 t = [t; t(end) + deltaT];        
%             end
%             
%             speed = [speed; v]; % km/h
%             time = [time; t];   % sec
%         end
%         
% %         speed_kmh = speed /1000 * 3600; 
% %         figure; plot(time, speed_kmh)
% %         figure; plot(time, speed)
% 
%         drivecycle_speed = [time speed];
%         drivecycle_grade = [0 0; height(drivecycle_speed)-1 0;];
%         drivecycle_name = 'NEDC';
end

%% Battery pack parameters
n_series = 15*8;
n_parallel = 1;
Q = 280;
V_nom = 3.2;
% k_battery = 1; % Amount of power provided by the designed battery pack

%% LOAD INIT FILES / PARAMETERS

cd('init_files')

% environmental, driver, power electronics
driver_ctrl;
environment;
accessory_Focus;

% ford focus plant/drivetrain
finalDrive_Focus;
wheel_Focus;
chassis_Focus;
ess_Focus;

% motor files
% choose motor type, 1 = copper, 2 = aluminum
motor_type = 2; % *************************************************change
motor_Windsor;

cd ..

%% SIMULATIONS PARAMETERS

% set parameters, decimations, load system
logvars = 4;
dec = 10;
len = max(drivecycle_speed(:,1));
load_system(sim_architecture)
set_param(sim_architecture, 'StopTime', 'len')

%% SETUP 

ser1 = 86; % ford pack is 86S5P
par1 = 5;

% finalize ess maps
ess.ocv_map = ess.ocv_cell * ser1;
ess.rdis_map = ess.rdis_cell * ser1 * par1;
ess.rchg_map = ess.rchg_cell * ser1 * par1;

%% RUN SIMULATION

ess.soc_start = 0.90; % ANL data is good for 15 - 90%

sim(sim_architecture)

%% DATA

% for mpge and efficiency
total_km = chas_plant_distance_out_simu(end)/1000;
total_mi = chas_plant_distance_out_simu(end)/1000/1.609;
ess_wh = ess1_Ws(end)/3600; 
ess_wh_mi = ess_wh / total_mi; %Em
ess_kwh = ess_wh/1000;

% for motor 

ave_mot_eff = sum(motor_eff)/length(motor_eff);
peak_t = max(motor_t);
peak_w = max(motor_w)/rpm_rads_conv;
% around 75% efficient in the city, 61% efficient on the highway

% % plot the torque speed maximum curves and efficiency map
% figure(1)
% hold on
% plot(motor_w/rpm_rads_conv,abs(motor_t),'k');
% plot(motor.torquespeed_RPMidx/rpm_rads_conv,motor.torquespeed_curve400,'r');
% % plot(motor.torquespeed_RPMidx,motor.torquespeed_curve225,'b');
% title('Torque-Speed Travelled vs. Limits @ 400V');
% xlabel('Speed (RPM)');ylabel('Torque (Nm)');
% legend('Actual','400V Limit','225V Limit');
% axis tight
% 
% figure(2)
% hold on
% surf(motor.spd_ind/rpm_rads_conv,motor.trq_ind,motor.eff_map(:,:,1)');
% plot3(motor_w/rpm_rads_conv,abs(motor_t),ones(length(motor_w),1),'k');
% title('Motor Efficiency Map @400V');
% xlabel('Speed (RPM)');ylabel('Torque (Nm)');
% axis tight
% colormap jet
% colorbar

%% Battery testing data
drivecycle_power = [drivecycle_power.time -drivecycle_power.signals.values];
% figure; plot(drivecycle_power(:, 1), drivecycle_power(:, 2)); title('power'); grid on

drivecycle_current = [drivecycle_current.time -drivecycle_current.signals.values];
% figure; plot(drivecycle_current(:, 1), drivecycle_current(:, 2)); title('current'); grid on

drivecycle_CRate = [drivecycle_current(:, 1) drivecycle_current(:, 2)/Q];
figure; 
ax1 = subplot(2, 1, 1); plot(drivecycle_speed(:, 1), drivecycle_speed(:, 2)); title(drivecycle_name)
ax2 = subplot(2, 1, 2); plot(drivecycle_CRate(:, 1), drivecycle_CRate(:, 2)); title('C rate'); grid on
linkaxes([ax1, ax2], 'x')

folder = "C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Model\Drive Train Model - Ford Focus 2013\Generate drive cycles\";
% save(folder + drivecycle_name + ".mat", 'drivecycle_power', 'drivecycle_current', 'drivecycle_CRate', 'Q')
