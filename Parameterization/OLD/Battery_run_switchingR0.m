clear

% Folders
folder_current = cd; 
folder_data_pulseTest = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\'; % where experimental data is saved
folder_data_OCVtest = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\';
folder_result = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\'; % where results from this script is saved
folder_functions = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Model\Electric Model\Parameterization\functions\'; % where functions used in this script is saved
addpath(folder_current);
addpath(folder_data_pulseTest);
addpath(folder_data_OCVtest);
addpath(folder_result);
addpath(folder_functions);

% Load charging R0
load('PROCESSED_PUL25chg_R0.mat'); 
param_chg(:, 2) = -param_chg(:, 2); 

% Set test profile & parameter set 
profileSet = ["US06"; "HWFET"; "UDDS"; "MIX"; "NEDC"];
paramSet = ["param_init"];
% paramSet = ["US06_OCVpoly9_c"; "US06_OCVpoly10_c"; "US06_OCVpoly11_c";"US06_OCVpoly15_c";"US06_OCVpoly18_c";"US06_OCVpoly19_c";"US06_OCVpoly21_c"; ];

% RMSE for each case
rmse = {}; 

for paramNum = 1:height(paramSet)
    
    figurename = append('Vt_', string(paramSet(paramNum)));
    figure('WindowStyle', 'docked', 'Name', figurename);

    % Parameter sets
    inputParameter = paramSet(paramNum); 
    switch inputParameter
        case 'param_pulse'
            useOCVpolynomial = 0;
            load('param_pulse.mat');
            load('C20CCCV_OCVSOC11.mat');

            SOC = param.SOC;
            OCV = param.OCV_init;
            R0 = [param.R0_init param_chg(:, 2)];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;
        case 'param_init'
            useOCVpolynomial = 0;
            load('param_init.mat');
            load('C20CCCV_OCVSOC11.mat');

            SOC = param.SOC;
            OCV = param.OCV_init;
            R0 = [param.R0_init param_chg(:, 2)];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;
        case 'US06_36_OCVest_uc'
            useOCVpolynomial = 0;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_unconstrained_US06_36.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_OCVest_c'
            useOCVpolynomial = 0;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_US06_OCVest_36.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'pulse_OCVest_c'
            useOCVpolynomial = 0;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_pulse_OCVest_36.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC11.mat');
            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'NEDC_OCVest_c'
            useOCVpolynomial = 0;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_NEDC_OCVest_36.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC11.mat');
            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'MIX_OCVest_c'
            useOCVpolynomial = 0;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_MIX_OCVest_35.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC11.mat');
            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-35;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'UDDS_OCVest_c'
            useOCVpolynomial = 0;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_UDDS_OCVest_36.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC11.mat');
            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_OCVpoly9_c'
            useOCVpolynomial = 1;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC9.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_US06_OCVpoly9_36.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_OCVpoly10_c'
            useOCVpolynomial = 1;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC10.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_US06_OCVpoly10_36.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_OCVpoly10_c_new'
            useOCVpolynomial = 1;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC10_new.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_US06_OCVpoly10_36_new.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_OCVpoly11_c'
            useOCVpolynomial = 1;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC11.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_US06_OCVpoly11_36.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_OCVpoly15_c'
            useOCVpolynomial = 1;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC15.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_US06_OCVpoly15_36.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_OCVpoly18_c'
            useOCVpolynomial = 1;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC18.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_US06_OCVpoly18_36.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_OCVpoly19_c'
            useOCVpolynomial = 1;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC19.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_US06_OCVpoly19_36.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_OCVpoly21_c'
            useOCVpolynomial = 1;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC21.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_US06_OCVpoly21_36.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_OCVest_c_withPolynomial_10'
            useOCVpolynomial = 1;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_US06_OCVest_36.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC10.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_OCVest_c_withPolynomial_11'
            useOCVpolynomial = 1;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_US06_OCVest_36.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_OCVest_c_withPolynomial_15'
            useOCVpolynomial = 1;
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Test\Results\param_opt_constrained_US06_OCVest_36.mat');
            load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\4 - Characterization test 3\Results\C20CCCV_OCVSOC15.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = [param.R0_init param.R0_init_chg];
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = [param.R0_opt(ind1:ind2) param.R0_opt_chg(ind1:ind2)];
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
    end
    validateattributes(SOC, {'double'}, {'increasing'})
        
    for runNum = 1:height(profileSet)
        %% Load files
        clear meas_t meas sim_result I SOC_exp t Vt_exp sim_SOC sim_t sim_Vt sim_Vt_exp Time

        % Input current profiles
        inputProfile = profileSet(runNum);
        switch inputProfile
            case 'pulse'
                data_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\7 - Pulse Discharge\';
                cd(data_folder);
                load('PUL25dch_meas_t.mat');    % Raw data
                Q = abs(meas_t.Ah(end) - meas_t.Ah(1));
                meas_t.Current = -meas_t.Current;
            case 'US06'
                filename = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\6 - Drive cycle\US0625';
                load(append(filename, '.mat'))

                meas_t = struct2table(meas);
                meas_t.Current = -meas_t.Current;
                Q = abs(meas_t.Ah(end) - meas_t.Ah(1));

                ind = find(meas.Voltage <=2.5, 1);
                if ~isempty(ind)
                    meas_t = meas_t(1:ind, :);
                    Q = abs(meas_t.Ah(end)-meas_t.Ah(1));
                end
                meas_t.Ah = meas_t.Ah - meas_t.Ah(1);
                meas_t.Time = meas_t.Time - meas_t.Time(1);
                meas_t.SOC = 1 - (-meas_t.Ah)/Q;
            case 'HWFET'
                current_folder = cd;
                data_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\6 - Drive cycle\';
                cd(data_folder);
                filename = 'HWFET25';
                file = append(data_folder, filename, '.mat');
                load(file)
                cd(current_folder)

                meas_t = struct2table(meas);
                meas_t.Current = -meas_t.Current;
                Q = abs(meas_t.Ah(end) - meas_t.Ah(1));

                ind = find(meas.Voltage <=2.5, 1);
                if ~isempty(ind)
                    meas_t = meas_t(1:ind, :);
                    Q = abs(meas_t.Ah(end)-meas_t.Ah(1));
                end
                meas_t.Ah = meas_t.Ah - meas_t.Ah(1);
                meas_t.Time = meas_t.Time - meas_t.Time(1);
                meas_t.SOC = 1 - (-meas_t.Ah)/Q;
            case 'UDDS'
                current_folder = cd;
                data_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\6 - Drive cycle\';
                cd(data_folder);
                filename = 'UDDS25';
                file = append(data_folder, filename, '.mat');
                load(file)
                cd(current_folder)

                meas_t = struct2table(meas);
                meas_t.Current = -meas_t.Current;
                Q = abs(meas_t.Ah(end) - meas_t.Ah(1));

                ind = find(meas.Voltage <=2.5, 1);
                if ~isempty(ind)
                    meas_t = meas_t(1:ind, :);
                    Q = abs(meas_t.Ah(end)-meas_t.Ah(1));
                end
                meas_t.Ah = meas_t.Ah - meas_t.Ah(1);
                meas_t.Time = meas_t.Time - meas_t.Time(1);
                meas_t.SOC = 1 - (-meas_t.Ah)/Q;

                % Delete inf & NaN
                meas_t = meas_t( ~any( isnan( meas_t.Current ) | isinf( meas_t.Current ), 2 ),: );
            case 'NEDC'
                current_folder = cd;
                data_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\6 - Drive cycle\';
                cd(data_folder);
                filename = 'NEDC23';
                file = append(data_folder, filename, '.mat');
                load(file)
                cd(current_folder)

                meas_t = struct2table(meas);
                meas_t.Current = -meas_t.Current;
                Q = abs(meas_t.Ah(end) - meas_t.Ah(1));

                ind = find(meas.Voltage <=2.5, 1);
                if ~isempty(ind)
                    meas_t = meas_t(1:ind, :);
                    Q = abs(meas_t.Ah(end)-meas_t.Ah(1));
                end
                meas_t.Ah = meas_t.Ah - meas_t.Ah(1);
                meas_t.Time = meas_t.Time - meas_t.Time(1);
                meas_t.SOC = 1 - (-meas_t.Ah)/Q;

                % Delete inf & NaN
                meas_t = meas_t( ~any( isnan( meas_t.Current ) | isinf( meas_t.Current ), 2 ),: );
            case 'MIX'
                current_folder = cd;
                data_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\6 - Drive cycle\';
                cd(data_folder);
                filename = 'MIX25';
                file = append(data_folder, filename, '.mat');
                load(file)
                cd(current_folder)

                meas_t = struct2table(meas);
                meas_t.Current = -meas_t.Current;
                Q = abs(meas_t.Ah(end) - meas_t.Ah(1));

                ind = find(meas.Voltage <=2.5, 1);
                if ~isempty(ind)
                    meas_t = meas_t(1:ind, :);
                    Q = abs(meas_t.Ah(end)-meas_t.Ah(1));
                end
                meas_t.Ah = meas_t.Ah - meas_t.Ah(1);
                meas_t.Time = meas_t.Time - meas_t.Time(1);
                meas_t.SOC = 1 - (-meas_t.Ah)/Q;

                % Delete inf & NaN
                meas_t = meas_t( ~any( isnan( meas_t.Current ) | isinf( meas_t.Current ), 2 ),: );
        end

        % Model to run
        model = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Model\Electric Model\Parameterization\Battery_switchingR0.slx';

        %% Simulink model parameters

        % SOCRange: 
        time_begin = 0;
        time_end = meas_t.Time(end); %**************************************Change
        I1_init = 0;
        I2_init = 0;
        I3_init = 0;

        % Beginning & end of the segment for optimization
        i_begin = find(meas_t.Time >= time_begin, 1);
        i_end = find(meas_t.Time >= time_end, 1);
        Time = meas_t.Time(i_begin:i_end);
        Time = Time - Time(1);
        runTime = Time(end);

        % Input data
        stepSize = 0.1;
        t = [Time Time];
        I = [Time meas_t.Current(i_begin:i_end)];
        SOC_init = meas_t.SOC(1);
        SOC_exp = [Time meas_t.SOC(i_begin:i_end)];

        % Input voltage
        Vt_exp = [Time meas_t.Voltage(i_begin:i_end)];

        %% Run simulation
        sim_result = sim(model);

        % Results
        sim_t = sim_result.Vt.time;
        sim_Vt_exp = sim_result.Vt.signals(2).values;
        sim_Vt = sim_result.Vt.signals(1).values;
        sim_SOC = sim_result.SOC.signals(2).values;

        % Overall RMSE
        rmse_Vt = sqrt(mean((sim_Vt - sim_Vt_exp).^2));
        
        % RMSE for each segment
        SOC_flip = flip(SOC);
        rmse_seg = [];
        for i = 1:length(SOC_flip)-1
            ind_start = find(sim_SOC <= SOC_flip(i), 1);
            ind_end = find(sim_SOC <= SOC_flip(i+1), 1);
            t_sim_i = sim_t(ind_start:ind_end);
            Vt_sim_i = sim_Vt_exp(ind_start:ind_end);
            Vt_exp2_i = sim_Vt(ind_start:ind_end);
            SOC_sim_i = sim_SOC(ind_start:ind_end);
            rmse_Vt_i = sqrt(mean((Vt_exp2_i - Vt_sim_i).^2));
            rmse_seg(i, runNum) = rmse_Vt_i;
        end
        
        % Store rmse for this case
        if isempty(rmse)
            rmse = {inputParameter inputProfile rmse_Vt rmse_seg};
        else
            caseNum = height(rmse) + 1; 
            rmse{caseNum, 1} = inputParameter;
            rmse{caseNum, 2} = inputProfile;
            rmse{caseNum, 3} = rmse_Vt;
            rmse{caseNum, 4} = rmse_seg;
        end
        
        %% Plot voltage curves

    %     % Vt vs Time
    %     figurename = append('Vt ', '(', string(input), ', ', string(paramset), ')');
    %     figure('WindowStyle', 'docked', 'Name', figurename); hold on
    %     plot(sim_t, sim_Vt_exp, '-');
    %     plot(sim_t, sim_Vt, '--');
    %     hold off; grid on; title(figurename, 'interpreter', 'none'); legend(['exp'; 'sim'])
    %     annotation('textbox', [0.14, 0.82, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt))

    %     % Vt_err vs Time
    %     figurename = append('Vt_err', '(', string(input), ', ', string(paramset), ')');
    %     figure('WindowStyle', 'docked', 'Name', figurename); hold on
    %     plot(sim_t, sim_Vt - sim_Vt_exp); 
    %     grid on; title(figurename, 'Interpreter', 'none');
    %     annotation('textbox', [0.14, 0.82, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt))

        % Vt vs SOC
        axes{runNum, 1} = subplot(height(profileSet), 1, runNum);
        hold on
        plot(sim_SOC, sim_Vt);
        plot(sim_SOC, sim_Vt_exp, '--');
        hold off; grid on; 
        xlim([-0.01, 1.01]); 
        legend(['sim'; 'exp']); title(inputProfile)
%         annotation('textbox', [0.14, 1.12-0.3*runNum, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt)) % 3 plots
        annotation('textbox', [0.13, 0.82-0.172*(runNum-1), 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt)) % 5 plots

    end

    %% Plot parameters
    
    % % Plot parameters (before & after)
    % figure('WindowStyle', 'docked', 'Name', 'Params');
    % ax1 = subplot(2, 5, 1); hold on; plot(SOC, param.R0_init, '.-'); plot(SOC, R0(:, 1), '.-'); title('R0_dch'); grid on; xlabel('SOC'); ylabel('R [Ohm]'); legend('before', 'after', 'Interpreter', 'none')
    % ax2 = subplot(2, 5, 2); hold on; plot(SOC, param.R1_init, '.-'); plot(SOC, R1, '.-'); title('R1'); grid on; xlabel('SOC'); ylabel('R [Ohm]'); legend('before', 'after')
    % ax3 = subplot(2, 5, 3); hold on; plot(SOC, param.R2_init, '.-'); plot(SOC, R2, '.-'); title('R2'); grid on; xlabel('SOC'); ylabel('R [Ohm]'); legend('before', 'after')
    % ax4 = subplot(2, 5, 4); hold on; plot(SOC, param.R3_init, '.-'); plot(SOC, R3, '.-'); title('R3'); grid on; xlabel('SOC'); ylabel('R [Ohm]'); legend('before', 'after')
    % ax5 = subplot(2, 5, 6); hold on; plot(SOC, param.R0_init_chg, '.-'); plot(SOC, R0(:, 2), '.-'); title('R0_chg'); grid on; xlabel('SOC'); ylabel('R [Ohm]'); legend('before', 'after', 'Interpreter', 'none')
    % ax6 = subplot(2, 5, 7); hold on; plot(SOC, param.tau1_init, '.-'); plot(SOC, tau1, '.-'); title('tau1'); grid on; xlabel('SOC'); ylabel('tau [s]'); legend('before', 'after')
    % ax7 = subplot(2, 5, 8); hold on; plot(SOC, param.tau2_init, '.-'); plot(SOC, tau2, '.-'); title('tau2'); grid on; xlabel('SOC'); ylabel('tau [s]'); legend('before', 'after')
    % ax8 = subplot(2, 5, 9); hold on; plot(SOC, param.tau3_init, '.-'); plot(SOC, tau3, '.-'); title('tau3'); grid on; xlabel('SOC'); ylabel('tau [s]'); legend('before', 'after')
    % ax9 = subplot(2, 5, 5); hold on; plot(SOC, param.OCV_init, '.-'); plot(SOC, OCV, '.-'); title('OCV'); grid on; xlabel('SOC'); ylabel('OCV [V]'); legend('before', 'after')
    % linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8, ax9], 'x')

    % Find slope of OCV curve wrt SOC
    OCV_diff = diff(OCV);
    % figure;
    % ax1 = subplot(2, 1, 1); plot(SOC, OCV); title('OCV'); grid on
    % ax2 = subplot(2, 1, 2); plot(SOC(1:length(OCV_diff)), OCV_diff); title('slope'); grid on
    % linkaxes([ax1, ax2], 'x')
    
    % Plot parameters (after)
    figure('WindowStyle', 'docked', 'Name', inputParameter);
    ax1 = subplot(2, 5, 1); hold on; plot(SOC, R0(:, 1), '.-'); title('R0_dch', 'interpreter', 'none'); grid on; xlabel('SOC'); ylabel('R [Ohm]');
    ax2 = subplot(2, 5, 2); hold on; plot(SOC, R1, '.-'); title('R1'); grid on; xlabel('SOC'); ylabel('R [Ohm]');
    ax3 = subplot(2, 5, 3); hold on; plot(SOC, R2, '.-'); title('R2'); grid on; xlabel('SOC'); ylabel('R [Ohm]');
    ax4 = subplot(2, 5, 4); hold on; plot(SOC, R3, '.-'); title('R3'); grid on; xlabel('SOC'); ylabel('R [Ohm]');
    ax5 = subplot(2, 5, 6); hold on; plot(SOC, R0(:, 2), '.-'); title('R0_chg', 'interpreter', 'none'); grid on; xlabel('SOC'); ylabel('R [Ohm]');
    ax6 = subplot(2, 5, 7); hold on; plot(SOC, tau1, '.-'); title('tau1'); grid on; xlabel('SOC'); ylabel('tau [s]');
    ax7 = subplot(2, 5, 8); hold on; plot(SOC, tau2, '.-'); title('tau2'); grid on; xlabel('SOC'); ylabel('tau [s]');
    ax8 = subplot(2, 5, 9); hold on; plot(SOC, tau3, '.-'); title('tau3'); grid on; xlabel('SOC'); ylabel('tau [s]');
    ax9 = subplot(2, 5, 5); hold on; plot(SOC, OCV, '.-'); title('OCV'); grid on; xlabel('SOC'); ylabel('OCV [V]');
    ax10 = subplot(2, 5, 10); hold on; plot(SOC(1:length(OCV_diff)), OCV_diff, '.-'); title('OCV slope'); grid minor; grid on
    linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10], 'x')

    figure('WindowStyle', 'docked', 'Name', inputParameter+" RMSE"); hold on; grid on
    for i = 1:width(rmse{end, 4})
        plot(SOC_flip(2:end), rmse{end+1-i, 4}(:, width(rmse{end, 4})+1-i)', '.-')
    end
    legend(rmse{end+1-width(rmse{end, 4}):end, 2}, 'interpreter', 'none'); xlim([-0.05, 1.05])


end


% %% Plot RMSE for different cases
% addpath('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\Projects\Meeting Slides\Optimization Results')
% load('CaseAll_RMSE.mat')
% optProfile_init = cell2mat(rmse(1:5, 3))'*1000;
% 
% %% Constrained vs unconstrained
% % optProfile_US06_uc = cell2mat(rmse(6:10, 3))'*1000;
% % optProfile_US06_c = cell2mat(rmse(11:15, 3))'*1000;
% % 
% % x = categorical({'Initial','Unconstrained','Constrained'});
% % x = reordercats(x,{'Initial','Unconstrained','Constrained'});
% % y = [optProfile_init; optProfile_US06_uc; optProfile_US06_c];
% % figure;
% % b = bar(x, y);
% % for i = 1:5
% %     xtips = b(i).XEndPoints;
% %     ytips = b(i).YEndPoints;
% %     labels = string(round(b(i).YData, 1));
% %     text(xtips,ytips,labels,'HorizontalAlignment','center',...
% %         'VerticalAlignment','bottom')
% % end
% % grid on
% % title('RMSE of Constrained & Unconstrained Optimization')
% % legend('US06', 'HWFET', 'UDDS', 'MIXED', 'NEDC')
% % xlabel("Profile Optimized With"); ylabel("RMSE [mV]"); 
% 
% %% Different drive cycle optimization
% optProfile_US06 = cell2mat(rmse(11:15, 3))'*1000;
% optProfile_pulse = cell2mat(rmse(16:20, 3))'*1000;
% optProfile_NEDC = cell2mat(rmse(21:25, 3))'*1000;
% optProfile_Mixed = cell2mat(rmse(26:30, 3))'*1000;
% optProfile_UDDS = cell2mat(rmse(31:35, 3))'*1000;
% 
% x = categorical({'Initial','Pulse','UDDS','US06','NEDC','Mixed'});
% x = reordercats(x,{'Initial','Pulse','UDDS','US06','NEDC','Mixed'});
% y = [optProfile_init; optProfile_pulse; optProfile_UDDS; optProfile_US06; optProfile_NEDC; optProfile_Mixed];
% figure;
% b = bar(x, y);
% for i = 1:5
%     xtips = b(i).XEndPoints;
%     ytips = b(i).YEndPoints;
%     labels = string(round(b(i).YData, 1));
%     text(xtips,ytips,labels,'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom')
% end
% grid on
% title('RMSE with Different Drive Cycles')
% legend('US06', 'HWFET', 'UDDS', 'MIXED', 'NEDC')
% xlabel("Profile Optimized With"); ylabel("RMSE [mV]"); 
% % 
% %% Different polynomial degrees
% optProfile_poly9 = cell2mat(rmse(36:40, 3))'*1000;
% optProfile_poly10 = cell2mat(rmse(41:45, 3))'*1000;
% optProfile_poly11 = cell2mat(rmse(46:50, 3))'*1000;
% optProfile_poly15 = cell2mat(rmse(51:55, 3))'*1000;
% optProfile_poly18 = cell2mat(rmse(56:60, 3))'*1000;
% optProfile_poly19 = cell2mat(rmse(61:65, 3))'*1000;
% optProfile_poly21 = cell2mat(rmse(66:70, 3))'*1000;
% optProfile_estOCV_poly10 = cell2mat(rmse(76:80, 3))'*1000;
% optProfile_estOCV_poly11 = cell2mat(rmse(81:85, 3))'*1000;
% optProfile_estOCV_poly15 = cell2mat(rmse(86:90, 3))'*1000;
% 
% % x = categorical({'PULSE', 'US06 opt', '9','10','11','15','18','19','21'});
% % x = reordercats(x,{'PULSE', 'US06 opt', '9','10','11','15','18','19','21'});
% % y = [optProfile_init; optProfile_US06; optProfile_poly9; optProfile_poly10; optProfile_poly11; optProfile_poly15; optProfile_poly18; optProfile_poly19; optProfile_poly21];
% % figure;
% % b = bar(x, y);
% % for i = 1:5
% %     xtips = b(i).XEndPoints;
% %     ytips = b(i).YEndPoints;
% %     labels = string(round(b(i).YData, 1));
% %     text(xtips,ytips,labels,'HorizontalAlignment','center',...
% %         'VerticalAlignment','bottom')
% % end
% % grid on
% % title('RMSE with Different Polynomial Degrees')
% % legend('US06', 'HWFET', 'UDDS', 'MIXED', 'NEDC')
% % xlabel("Polynomial Fit Degree"); ylabel("RMSE [mV]"); 
% 
% % x = categorical(  {'PULSE','US06 est', '10 est','10 poly','11 est','11 poly','15 est','15 poly'});
% % x = reordercats(x,{'PULSE','US06 est', '10 est','10 poly','11 est','11 poly','15 est','15 poly'});
% % y = [optProfile_init; optProfile_US06; optProfile_estOCV_poly10; optProfile_poly10; ...
% %                                        optProfile_estOCV_poly11; optProfile_poly11; ...
% %                                        optProfile_estOCV_poly15; optProfile_poly15;];
% % figure;
% % b = bar(x, y);
% % for i = 1:5
% %     xtips = b(i).XEndPoints;
% %     ytips = b(i).YEndPoints;
% %     labels = string(round(b(i).YData, 1));
% %     text(xtips,ytips,labels,'HorizontalAlignment','center',...
% %         'VerticalAlignment','bottom')
% % end
% % grid on
% % title('RMSE with Different Polynomial Degrees')
% % legend('US06', 'HWFET', 'UDDS', 'MIXED', 'NEDC')
% % xlabel("Polynomial Fit Degree"); ylabel("RMSE [mV]"); 