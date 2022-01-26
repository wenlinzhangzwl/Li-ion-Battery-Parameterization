clear

% Folders
folder_current = cd; 
folder_project = 'C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design';
folder_functions = append(folder_project, '\Scripts\Parameterization\functions\'); 
folder_data = append(folder_project, '\Test Results'); % where experimental data is saved
folder_result = append(folder_project, '\Test Results\7 - Pulse Test\Results\'); % where results from this script is saved
addpath(folder_current);
addpath(folder_functions);
addpath(genpath(folder_data));
addpath(folder_result);

% Set test profile & parameter set 
% profileSet = ["Pulse"; "US06"; "HWFET"; "UDDS"; "MIX"; "NEDC"];
% paramSet = ["US06_Poly7_36";"US06_Poly8_36";"US06_Poly9_36";"US06_Poly10_36";"US06_Poly11_36";"US06_Poly15_36";"US06_Poly18_36";"US06_Poly19_36";"US06_Poly21_36"];

profileSet = ["US06"; "HWFET"; "UDDS"; "MIX"; "NEDC"];
paramSet = ["MLBS0p2_LS"]; % "MLBS0p2_PSO", "MLBS0p2_LS"

% RMSE for each case
rmse = zeros(height(profileSet)*height(paramSet), 5); 
rmse = num2cell(rmse);

for paramNum = 1:height(paramSet)
    
    figurename = append('Vt_', string(paramSet(paramNum)));
    figure('WindowStyle', 'docked', 'Name', figurename);

    % Load parameter set
    inputParameter = paramSet(paramNum); 
    switch inputParameter
        case 'param_Pulse0p1'
            useOCVpolynomial = 0;
            load('param_Pulse0p1.mat');
            load('C20CCCV_OCVSOC11.mat');

            SOC = param.SOC;
            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;
        case 'param_Pulse0p8'
            useOCVpolynomial = 0;
            load('param_Pulse0p8.mat');
            load('C20CCCV_OCVSOC11.mat');

            SOC = param.SOC;
            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;
        case 'Pulse0p1_Poly0_36'
            useOCVpolynomial = 0;
            load('param_opt_Pulse0p1_Poly0_36.mat');
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'Pulse0p8_Poly0_36'
            useOCVpolynomial = 0;
            load('param_opt_Pulse0p8_Poly0_36.mat');
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'UDDS_Poly0_36'
            useOCVpolynomial = 0;
            load('param_opt_UDDS_Poly0_36.mat');
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_Poly0_36'
            useOCVpolynomial = 0;
            load('param_opt_US06_Poly0_36.mat');
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'HWFET_Poly0_36'
            useOCVpolynomial = 0;
            load('param_opt_HWFET_Poly0_36.mat');
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'MIXED_Poly0_36'
            useOCVpolynomial = 0;
            load('param_opt_MIXED_Poly0_36.mat');
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'NEDC_Poly0_36'
            useOCVpolynomial = 0;
            load('param_opt_NEDC_Poly0_36.mat');
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_Poly7_36'
            useOCVpolynomial = 1;
            polyDeg = 7;
            filename = "param_opt_US06_Poly" + polyDeg + "_36.mat"; 
            load(filename);
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_Poly8_36'
            useOCVpolynomial = 1;
            polyDeg = 8;
            filename = "param_opt_US06_Poly" + polyDeg + "_36.mat"; 
            load(filename);
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_Poly9_36'
            useOCVpolynomial = 1;
            polyDeg = 9;
            filename = "param_opt_US06_Poly" + polyDeg + "_36.mat"; 
            load(filename);
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_Poly10_36'
            useOCVpolynomial = 1;
            polyDeg = 10;
            filename = "param_opt_US06_Poly" + polyDeg + "_36.mat"; 
            load(filename);
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_Poly11_36'
            useOCVpolynomial = 1;
            polyDeg = 11;
            filename = "param_opt_US06_Poly" + polyDeg + "_36.mat"; 
            load(filename);
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_Poly15_36'
            useOCVpolynomial = 1;
            polyDeg = 15;
            filename = "param_opt_US06_Poly" + polyDeg + "_36.mat"; 
            load(filename);
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_Poly18_36'
            useOCVpolynomial = 1;
            polyDeg = 18;
            filename = "param_opt_US06_Poly" + polyDeg + "_36.mat"; 
            load(filename);
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_Poly19_36'
            useOCVpolynomial = 1;
            polyDeg = 19;
            filename = "param_opt_US06_Poly" + polyDeg + "_36.mat"; 
            load(filename);
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case 'US06_Poly21_36'
            useOCVpolynomial = 1;
            polyDeg = 21;
            filename = "param_opt_US06_Poly" + polyDeg + "_36.mat"; 
            load(filename);
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case "MLBS0p2_LS"
            useOCVpolynomial = 0;
            load('LS_param_opt_MLBS_0p2Hz_Poly0_36.mat');
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
        case "MLBS0p2_PSO"
            useOCVpolynomial = 0;
            load('PSO_param_opt_MLBS_0p2Hz_Poly0_36.mat');
            load('C20CCCV_OCVSOC11.mat');

            param = flip(param);
            SOC = param.SOC;

            OCV = param.OCV_init;
            R0 = param.R0_init;
            R1 = param.R1_init;
            R2 = param.R2_init;
            R3 = param.R3_init;
            tau1 = param.tau1_init;
            tau2 = param.tau2_init;
            tau3 = param.tau3_init;

            ind1 = 37-36;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2, :) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
    end
    validateattributes(SOC, {'double'}, {'increasing'})
        
    % Variable 'parameters' as an input to the model
    [parameters.SOC, parameters.R0, parameters.R1, parameters.R2, parameters.R3, parameters.T1, parameters.T2, parameters.T3, parameters.OCV] = ...
        deal(SOC, R0, R1, R2, R3, tau1, tau2, tau3, OCV);

    % Run model for the current case
    parfor profileNum = 1:height(profileSet)
        %% Load test profiles
%         clear meas_t meas sim_result I SOC_exp t Vt_exp sim_SOC sim_t sim_Vt sim_Vt_exp Time

        inputProfile = profileSet(profileNum);
        switch inputProfile
            case 'Pulse'
                data = load('PROCESSED_PUL25dch.mat', 'meas_t');
                meas_t = data.meas_t;
%                 clear meas_resampled pulse
            case 'US06'
                data = load('US0625.mat');
                meas = data.meas; 
%                 clear data
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
                data_folder = append(folder_project, '\Test Results\6 - Drive cycle\');
                cd(data_folder);
                filename = 'HWFET25';
                file = append(data_folder, filename, '.mat');
                data = load(file);
                meas = data.meas; 
%                 clear data
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
                data_folder = append(folder_project, '\Test Results\6 - Drive cycle\');
                cd(data_folder);
                filename = 'UDDS25';
                file = append(data_folder, filename, '.mat');
                data = load(file);
                meas = data.meas; 
%                 clear data
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
                data_folder = append(folder_project, '\Test Results\6 - Drive cycle\');
                cd(data_folder);
                filename = 'NEDC23';
                file = append(data_folder, filename, '.mat');
                data = load(file);
                meas = data.meas; 
%                 clear data
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
                data_folder = append(folder_project, '\Test Results\6 - Drive cycle\');
                cd(data_folder);
                filename = 'MIX25';
                file = append(data_folder, filename, '.mat');
                data = load(file);
                meas = data.meas; 
%                 clear data
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
            case "MLBS0p2"
                current_folder = cd;
                data_folder = append(folder_project, '\Test Results\\9 - Pulse-PRBS Tests\');
                cd(data_folder);
                filename = 'MLBS_0p2Hz_RT_280Ah';
                file = append(data_folder, filename, '.mat');
                data = load(file);
                meas = data.meas; 
%                 clear data
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
        meas = data.meas; 
        
        %% Run model

        % Range to simulation 
        time_begin = 0;
        time_end = meas_t.Time(end);
        i_begin = find(meas_t.Time >= time_begin, 1);
        i_end = find(meas_t.Time >= time_end, 1);

        % Initial states
        initStates.SOC = meas_t.SOC(i_begin); 
        initStates.I_RC = [0 0 0]; 
        initStates.Vt = meas_t.Voltage(i_begin);
        
        deltaT = 0.1;
        curr_dmd = meas_t.Current(i_begin:i_end);

        [Vt, ~, ~] = Model_3RC(curr_dmd, initStates, parameters, Q, deltaT);
        initStates = [];  % Overwrite var so it can be recognized as a temporary variable by parfor

        %% Calculate simulation results

        % Calculate error
        Vt_err = Vt - meas_t.Voltage;
        Vt_rmse = sqrt(mean((meas_t.Voltage(1:height(Vt)) - Vt).^2));
        Vt_max_err = max(abs(Vt_err)); 

        % Write result to "rmse_temp"
        rmse_temp{profileNum, :} = [{inputParameter}, {inputProfile}, {Vt_rmse}, {Vt_max_err}, {[meas_t.Time, meas_t.SOC, meas_t.Voltage, Vt]}]; 
    end

    % Write 'rmse_temp' to 'rmse'
    ind = (paramNum-1)*height(profileSet) + 1; 
    for i = 1:height(rmse_temp)
        [rmse{ind, 1}, rmse{ind, 2}, rmse{ind, 3}, rmse{ind, 4}, rmse{ind, 5}] = ...
            deal(rmse_temp{i,1}{1, 1}, rmse_temp{i,1}{1, 2}, rmse_temp{i,1}{1, 3}, rmse_temp{i,1}{1, 4}, rmse_temp{i,1}{1, 5}); 
        ind = ind + 1; 
    end

    %% Plot results for the current set of parameters 

    % Plot Vt_err
    figurename = append('Vt_err (', string(paramSet), ')');
    figure('WindowStyle', 'docked', 'Name', figurename); hold on
    legend_txt = []; 
    for i = 1:height(rmse_temp)
        SOC_i = rmse_temp{i, 1}{1, 5}(:, 2); 
        Vt_err_i = rmse_temp{i, 1}{1, 5}(:, 4) - rmse_temp{i, 1}{1, 5}(:, 3); 
        plot(SOC_i, Vt_err_i); 

        legend_txt = [legend_txt; rmse_temp{i, 1}{1, 2}];
    end
    grid on; title(figurename, 'Interpreter', 'none');
    legend(legend_txt, 'Interpreter', 'none');

    % Plot Vt
    figurename = append('Vt (', string(paramSet), ')');
    figure('WindowStyle', 'docked', 'Name', figurename);
    for i = 1:height(rmse_temp)
        axes{i, 1} = subplot(height(rmse_temp), 1, i); hold on

        SOC_i = rmse_temp{i, 1}{1, 5}(:, 2); 
        Voltage_i = rmse_temp{i, 1}{1, 5}(:, 3);
        Vt_sim_i = rmse_temp{i, 1}{1, 5}(:, 4); 
        inputProfile_i = rmse_temp{i, 1}{1, 2}; 
        rmse_i = rmse_temp{i, 1}{1, 3}; 

        plot(SOC_i, Voltage_i);
        plot(SOC_i, Vt_sim_i);

        grid on; legend(['exp'; 'sim']); title(inputProfile_i)
        xlim([-0.01, 1.01]); 
    %         annotation('textbox', [0.14, 1.12-0.3*runNum, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt*1000)) % 3 plots
        annotation('textbox', [0.13, 0.82-0.172*(i-1), 0.1, 0.1], 'String', "RMSE = " + string(rmse_i*1000)) % 5 plots
    %         annotation('textbox', [0.13, 0.82-0.145*(runNum-1), 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt*1000)) % 6 plots
    end

    % Plot parameters (before & after)
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
    % OCV_diff = diff(OCV);
    % figure;
    % ax1 = subplot(2, 1, 1); plot(SOC, OCV); title('OCV'); grid on
    % ax2 = subplot(2, 1, 2); plot(SOC(1:length(OCV_diff)), OCV_diff); title('slope'); grid on
    % linkaxes([ax1, ax2], 'x')
    
    % Plot parameters (after)
    figure('WindowStyle', 'docked', 'Name', inputParameter);
    ax1 = subplot(2, 4, 1); hold on; plot(SOC, R0(:, 1), '.-'); title('R0_dch', 'interpreter', 'none'); grid on; xlabel('SOC'); ylabel('R [Ohm]');
    ax2 = subplot(2, 4, 2); hold on; plot(SOC, R1, '.-'); title('R1'); grid on; xlabel('SOC'); ylabel('R [Ohm]');
    ax3 = subplot(2, 4, 3); hold on; plot(SOC, R2, '.-'); title('R2'); grid on; xlabel('SOC'); ylabel('R [Ohm]');
    ax4 = subplot(2, 4, 4); hold on; plot(SOC, R3, '.-'); title('R3'); grid on; xlabel('SOC'); ylabel('R [Ohm]');
    ax5 = subplot(2, 4, 5); hold on; plot(SOC, tau1, '.-'); title('tau1'); grid on; xlabel('SOC'); ylabel('tau [s]');
    ax6 = subplot(2, 4, 6); hold on; plot(SOC, tau2, '.-'); title('tau2'); grid on; xlabel('SOC'); ylabel('tau [s]');
    ax7 = subplot(2, 4, 7); hold on; plot(SOC, tau3, '.-'); title('tau3'); grid on; xlabel('SOC'); ylabel('tau [s]');
    ax8 = subplot(2, 4, 8); hold on; plot(SOC, OCV, '.-'); title('OCV'); grid on; xlabel('SOC'); ylabel('OCV [V]');
    linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8], 'x')

end










%****************************************************************************************************
%****************************************************************************************************

%% Plot RMSE for different cases
% load('C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Model\Electric Model\Parameterization\Optimization Results\RMSE_All.mat')
% 
% %% Different drive cycle optimization
% optProfile_pulse0p1 = cell2mat(rmse(1:6, 3))'*1000;
% optProfile_pulse0p8 = cell2mat(rmse(7:12, 3))'*1000;
% optProfile_UDDS = cell2mat(rmse(13:18, 3))'*1000;
% optProfile_US06 = cell2mat(rmse(19:24, 3))'*1000;
% optProfile_HWFET = cell2mat(rmse(25:30, 3))'*1000;
% optProfile_Mixed = cell2mat(rmse(31:36, 3))'*1000;
% optProfile_NEDC = cell2mat(rmse(37:42, 3))'*1000;
% 
% x = categorical({'Pulse(0.1C)','Pulse(0.8C)','UDDS','US06','HWFET','Mixed','NEDC'});
% x = reordercats(x,{'Pulse(0.1C)','Pulse(0.8C)','UDDS','US06','HWFET','Mixed','NEDC'});
% y = [optProfile_pulse0p1(2:end); optProfile_pulse0p8(2:end); optProfile_UDDS(2:end); optProfile_US06(2:end); optProfile_HWFET(2:end); optProfile_Mixed(2:end); optProfile_NEDC(2:end)];
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


%% Different polynomial degrees
% 
% optProfile_poly7 = cell2mat(rmse(43:48, 3))'*1000;
% optProfile_poly8 = cell2mat(rmse(49:54, 3))'*1000;
% optProfile_poly9 = cell2mat(rmse(55:60, 3))'*1000;
% optProfile_poly10 = cell2mat(rmse(61:66, 3))'*1000;
% optProfile_poly11 = cell2mat(rmse(67:72, 3))'*1000;
% optProfile_poly15 = cell2mat(rmse(73:78, 3))'*1000;
% optProfile_poly18 = cell2mat(rmse(79:84, 3))'*1000;
% optProfile_poly19 = cell2mat(rmse(85:90, 3))'*1000;
% optProfile_poly21 = cell2mat(rmse(91:96, 3))'*1000;
% 
% x = categorical({'Estimated OCV','7','8','9','10','11','15','18','19','21'});
% x = reordercats(x,{'Estimated OCV','7','8','9','10','11','15','18','19','21'});
% y = [optProfile_US06(3:end); optProfile_poly7(3:end); optProfile_poly8(3:end); optProfile_poly9(3:end); optProfile_poly10(3:end);...
%      optProfile_poly11(3:end); optProfile_poly15(3:end); optProfile_poly18(3:end); optProfile_poly19(3:end); optProfile_poly21(3:end)];
% figure;
% b = bar(x, y);
% for i = 1:4
%     xtips = b(i).XEndPoints;
%     ytips = b(i).YEndPoints;
%     labels = string(round(b(i).YData, 1));
%     text(xtips,ytips,labels,'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom')
% end
% grid on
% title('RMSE with Different Polynomial Degrees')
% legend('HWFET', 'UDDS', 'MIXED', 'NEDC')
% xlabel("Polynomial Fit Degree"); ylabel("RMSE [mV]"); ylim([6,20]);
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