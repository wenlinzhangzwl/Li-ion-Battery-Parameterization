%% IMM_SVSF_VBL STRATEGY with SoC bias Estimation
% Sara Rahimifard, rahimifs@mcmaster.ca April 2022
% For more details regarding the algorithm read the following paper: 
% paper: Interacting Multiple Model Strategy for Electric Vehicle Batteries State of Charge/Health/ Power Estimation, Access 2021
% URL: https://ieeexplore.ieee.org/abstract/document/9507489

clear all
close all
clc
% A1 = [-1369.24347228066,7421.71822250840,-17648.3540656459,24203.0063080228,-21104.6830331269,12111.8742123879,-4555.43282767973,1081.32513507328,-150.890441725417,11.6990766992578,3.17153839630177];
%% week1 
% parameters of a third-order RC model with SoH = 100%, obtained by pulse discharge test. 
    param.model1.R0 = [0.0075 0.0089285 0.0083931 0.0085407 0.0076668 0.0078668 0.0085128 0.0073136 0.0074673 0.0081496 0.0072858 0.0078491 0.0077757 0.0079393 0.0089019];
    param.model1.R1 = [0.0010469 0.0021494 0.0022383 0.0022167 0.0024672 0.0019318 0.001029 0.0019691 0.002171 0.0010662 0.002096 0.0020532 0.0018474 0.0022544 0.0013452];
    param.model1.R2 = [0.0015 0.0018662 0.0026024 0.0026655 0.0025175 0.0026477 0.0028206 0.0027098 0.0027313 0.0025242 0.0029337 0.002909 0.002787 0.0029233 0.0016296];
    param.model1.R3 = [0.0085 0.0065716 0.0075182 0.0059005 0.0072186 0.0053778 0.0039391 0.0048358 0.0075815 0.0065545 0.0028382 0.01 0.0081938 0.001809 0.0075096];
    param.model1.tau1 = [1 0.84817 0.33065 0.60917 0.49397 0.59296 0.76472 0.56026 0.63189 0.7807 0.71081 0.64937 0.55428 0.54859 0.54865];
    param.model1.tau2 = [10 8.4758 7.3269 13.181 11.033 11.865 11.814 11.814 13.559 10.747 12.376 13.118 11.929 12.437 7.7433];
    param.model1.tau3 = [100 97.952 97.155 88.827 108.58 90.41 76.25 89.734 130 86.683 70.506 117.07 111.57 92.252 122.45];
   param.model1.Cn = 5.4*ones(1,15);
    A1 = [-2092.5 9351.9 -17282 16961 -9387.9 2860.7 -492.24 136.88 -72.802 20.041 1.6105];
    %% week30
% parameters of a third-order RC model with SoH = 90%, obtained by pulse discharge test. 
 
    param.model2.R0 = [0.0091707 0.010066 0.013172 0.013308 0.0093895 0.011691 0.010873 0.011875 0.0096341 0.010043 0.013011 0.010459 0.011152 0.0098861 0.0071416];
    param.model2.R1 = [0.018107 0.014002 0.0076243 0.0057001 0.00826 0.0052229 0.0052406 0.0038002 0.0056543 0.0049584 0.0019843 0.0039876 0.0034365 0.0045449 0.0078978];
    param.model2.R2 = [0.0078349 0.0056479 0.004503 0.0038334 0.0036244 0.0049028 0.0041169 0.0036426 0.0037729 0.0038389 0.0047267 0.0031594 0.0034328 0.0040039 0.00038265];
    param.model2.R3 = [0.00014214 0.010313 0.0040451 0.0070796 0.0070506 0.004583 0.01028 0.0091275 0.0090881 0.0080761 0.0015705 0.0036376 0.017164 1.9848e-05 0.027572];
    param.model2.tau1 = [0.35803 0.29585 0.33169 0.32434 0.16853 0.19787 0.18721 0.21347 0.1621 0.15698 0.34798 0.16291 0.20571 0.1562 0.14404];
    param.model2.tau2 = [10 12.199 11.753 12.01 10.937 16.946 14.624 12.47 13.336 12.972 18.977 11.665 12.811 12.704 10.097];
    param.model2.tau3 = [240.54 171.33 227.92 250 150.01 232.77 249.5 249.91 249.99 249.52 150.36 150.05 249.86 249.99 249.99];
    param.model2.Cn = 4.77*ones(1,15);
      A2 = [67006.6412715162	-355951.237408991	830816.525021041	-1120201.57641690	964579.013854285	-553286.067441687	213743.979431011	-54830.4962158105	8927.77493524045	-832.058541135094	37.1987933639492];
 %% week45
% parameters of a third-order RC model with SoH = 80%, obtained by pulse discharge test. 
  
    param.model3.R0 = [0.017917 0.017601 0.015062 0.015 0.015002 0.015 0.015004 0.015458 0.015386 0.0153 0.015328 0.0153 0.015243 0.01521 0.015282];
    param.model3.R1 = [0.024718 0.011197 0.011262 0.0082635 0.0067136 0.0056299 0.0046007 0.0038367 0.0035903 0.0031969 0.0029377 0.0028532 0.0026446 0.0026623 0.0026902];
    param.model3.R2= [0.0026418 0.0068817 0.005351 0.0047763 0.0049192 0.0049489 0.0044063 0.0042652 0.0041161 0.0045352 0.0038498 0.0050979 0.0029469 0.0030304 0.0020633];
    param.model3.R3 = [0.088737 8.9907e-06 0.0030686 0.0097389 0.0077232 0.0060569 0.0098072 0.010402 0.0097264 0.009793 0.010335 0.0049609 0.0023403 0.018746 0.013109];
    param.model3.tau1 = [1 0.74855 0.48457 0.33822 0.30246 0.31428 0.27296 0.53285 0.596 0.39839 0.53206 0.67917 0.7449 0.43681 0.53374];
    param.model3.tau2 = [15 15.569 15.737 15.014 15.006 17.866 15.714 18.234 18.166 17.513 17.967 26.451 15.06 17.264 16.853];
    param.model3.tau3 = [350 350 349.7 350 250.26 250.3 281.29 349.98 350 350 349.88 250.39 305.8 348.1 260.92];
    param.model3.Cn = 4.33*ones(1,15);
    A3 = [7867.61770046033,-52053.6516069609,148124.836069431,-239805.009792072,245253.317990107,-165834.200826406,75160.6883321978,-22563.9811919106,4296.71317091153,-468.314195927266,25.6920380353340];
%% model preperation for IMM algorithm
% initial parameters of the models 
SOC_LUT = (0.2:0.05:.9); %SoC look-up table 
Eta = 1; % Cell Coulombic Efficiency
DeltaT = 0.1; % sample time
% x_init = [0;0;0;0.9]; 
     
imm_svsf.model.mdl1 = param.model1;
imm_svsf.model.mdl2 = param.model2;
imm_svsf.model.mdl3 = param.model3;


 imm_svsf.sim.mdl1 = [imm_svsf.model.mdl1.R1; imm_svsf.model.mdl1.R2; imm_svsf.model.mdl1.R3; imm_svsf.model.mdl1.tau1; imm_svsf.model.mdl1.tau2; imm_svsf.model.mdl1.tau3; imm_svsf.model.mdl1.Cn;param.model1.R0];
 imm_svsf.sim.mdl2 = [imm_svsf.model.mdl2.R1; imm_svsf.model.mdl2.R2; imm_svsf.model.mdl2.R3; imm_svsf.model.mdl2.tau1; imm_svsf.model.mdl2.tau2; imm_svsf.model.mdl2.tau3; imm_svsf.model.mdl2.Cn;param.model2.R0];
 imm_svsf.sim.mdl3 = [imm_svsf.model.mdl3.R1; imm_svsf.model.mdl3.R2; imm_svsf.model.mdl3.R3; imm_svsf.model.mdl3.tau1; imm_svsf.model.mdl3.tau2; imm_svsf.model.mdl3.tau3; imm_svsf.model.mdl3.Cn;param.model3.R0];
%% initialization for the algorithm, IMM
x0 = [0; 0;0;.85;0;0.01]; %initial states for all the models [v1,v2,v3,soc,soc bias, internal resistance]
R = [.5e-2 0;0  9.5e-2]; % Measurement noise covariance matrix
Q = diag([1e-8,1e-8,1e-8,1e-10,1e-10,1e-10]); %System noise covariance matrix
P = diag([1e1,1e1,1e1,1e1,1e1,1e1]);%State error covariance matrix
Gamma = .38; %SVSF ‘‘convergence’’ or memory parameter
Psi =3.2; %SVSF smoothing boundary layer width
delta = 0.1;% sample time
bias = 0; %added bias to the calculated SoC from CC, this can be kept zero for real-time as there is already bias in the data
Wsoh = [1, .9 ,.8]; %weightning vector for SoH calculation using the mode probability 
alpha = 0.5; %factor for the estimated SoH based on the internal resistance and mode probability
%% input data 
load SOCScan_DOD2_20131129044058_CustRec_week1.mat %validation data for a battery with SoH = 100% 
imm_svsf.data.t = [week1(:,1)]; %time 
imm_svsf.data.I1 = [week1(:,1),week1(:,2)]; %input, current
imm_svsf.data.V1 = [week1(:,1),week1(:,3)]; % measurement 1, terminal voltage
imm_svsf.data.soc1 = [week1(:,1),week1(:,5)];% measurement 2, calculated SoC using coulomb counting method


%% run the simulink

sim('IMM_SVSF_VBL.slx'); 


 %% plotting 

% figure
% plot(Verror.time,Verror.signals(1).values,'LineWidth',1)
% hold on    
% plot(Verror.time, Verror.signals(2).values,'LineWidth',1)
% hold on 
% plot(Verror.time, Verror.signals(3).values,'LineWidth',1)
% legend('100\% capacity','90\% capacity','80\% capacity','Interpreter','latex')
% ylabel("Terminal Voltage Error(v)",'interpreter','latex')
% xlabel("time(seconds)")
% ylim([-0.1 0.1])
% grid on 
% magnifyOnFigure
% 
%% RMSE calculation  
% RMSE_SOC100    = sqrt((sum((SOCin.signals(1).values - SOCes.signals(1).values ).^2))/((length(SOCin.signals(1).values))))
% RMSE_V100    = sqrt((sum((Vac.signals(1).values(20:end) - Ves.signals(1).values(20:end) ).^2))/((length(Vac.signals(1).values(20:end)))))
