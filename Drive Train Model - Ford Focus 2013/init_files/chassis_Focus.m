%% Parameters

rpm_rads_conv= 1/60*2*pi;

chassis.area     = 2.42424; %Focus: 1.82 m wide, 1.48 m high, x 0.9
chassis.mass     = 1791; %from ANL sheet % prev 1643+80; %Curb weight + 80kg driver
chassis.cd       = 0.295; % need to confirm, on focus help page

veh_init_speed = 0;

%max wheel toque from dyno
vpc.prop.init.whl_trq_max.idx1_eng_spd= rpm_rads_conv*[1500 2000 2500 3000 3500 4000 4250];
vpc.prop.init.whl_trq_max.map= [240 340 410 450 485 500 435];
