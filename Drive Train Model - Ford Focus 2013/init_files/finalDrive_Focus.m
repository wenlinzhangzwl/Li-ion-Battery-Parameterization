%% FINAL DRIVE PARAMETERS

%{
3/18/2019
changed the final drive ratio to 1:1 to reflect direct drive application
changed the efficiency of the final drive to 100%
%}

mph_to_ms = 0.4470;

fd.ratio   			 = 1; % focus drive ratio from ANL is 7.82, direct drive is 1
fd.inert     		 = 0;
fd.thresh   		 = 10;

fd.torque      = [51.40,52.40,104.7,157.1,209.4,261.8,314.2,366.5,418.9,471.2,523.6];
fd.speed       = [0.500,6.000,33.90,67.80,101.7,135.6,169.5,203.4,237.3,271.2,305.1,339];
% fd.effmap      = ones(size(fd.torque,2),size(fd.speed,2)).* 0.987;
fd.effmap      = ones(size(fd.torque,2),size(fd.speed,2));
