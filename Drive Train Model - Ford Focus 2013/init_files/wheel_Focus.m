%% PARAMETERS

% tire size: P225/5017R17 for cmax
% size is the same for the ford focus electric 2013

mph_to_ms = 0.44704;

wheel.num_wheels = 2;
wheel.effFraction = 1;
wheel.inert_per_wheel   = 1.0;				% kg-m^2 (standard)
wheel.inertia = wheel.inert_per_wheel*wheel.num_wheels;

wheel.rated_radius = 0.4318;   
wheel.coeffFriction = 0.95;     
wheel.radius = wheel.rated_radius * wheel.coeffFriction;

wheel.speed_thresh = 1.00;

wheel.u1 = 0.006;
wheel.u2 = 0.0001;
