clear

% Folders
folder_current = cd; 
folder_functions = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\Scripts\Parameterization\functions"; 
folder_profile = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\11 - Drive Cycle (0p8C)\"; % where profiles are saved
folder_result = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\Results\"; % where results from this script is saved
addpath(folder_current, folder_functions, folder_profile, genpath(folder_result));

% Test settings
cell = "EVE280";
deltaT = 0.1;
switch cell
    case "EVE280"
        Vmax = 3.65; 
        Vmin = 2.5;
end

save_rmse_of_each_parameter_set = 0;
save_rmse_of_all_parameter_sets = 0; 

% Set test profile & parameter set 
profileSet = ["US06"; "HWFET"];
% profileSet = ["US06"; "HWFET"; "UDDS"; "MIX"; "NEDC"];

paramSet = ["Parameters_Optimized_US06_25degC_1p0_PSO_Poly0"];
% paramSet = ["Parameters_Layered_PUL_25degC_0p8_30min";...
%             "Parameters_Layered_PUL_25degC_0p8_60min";...
%             "Parameters_Layered_PUL_25degC_0p8_120min";...
%             "Parameters_Optimized_PUL_25degC_0p8_30min_PSO_Poly0";...
%             "Parameters_Optimized_PUL_25degC_0p8_60min_PSO_Poly0";...
%             "Parameters_Optimized_PUL_25degC_0p8_120min_PSO_Poly0";...
%             "Parameters_Optimized_HWFET_25degC_0p8_PSO_Poly0";...
%             "Parameters_Optimized_MIX_25degC_0p8_PSO_Poly0";...
%             "Parameters_Optimized_NEDC_25degC_0p8_PSO_Poly0";...
%             "Parameters_Optimized_UDDS_25degC_0p8_PSO_Poly0";...
%             "Parameters_Optimized_US06_25degC_0p8_PSO_Poly0";
%             "Parameters_Optimized_HWFET_25degC_0p47_PSO_Poly0";...
%             "Parameters_Optimized_MIX_25degC_1p0_PSO_Poly0";...
%             "Parameters_Optimized_NEDC_23degC_0p43_PSO_Poly0";...
%             "Parameters_Optimized_UDDS_25degC_0p45_PSO_Poly0";...
%             "Parameters_Optimized_US06_25degC_1p0_PSO_Poly0";];
% paramSet = ["Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly7_36";...
%             "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly8_36";
%             "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly9_36";...
%             "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly10_36";...
%             "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly11_36";...
%             "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly12_36";...
%             "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly13_36";...
%             "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly14_36";];

% RMSE for each case
rmse = zeros(height(profileSet)*height(paramSet), 5); 
rmse = num2cell(rmse);

for paramNum = 1:height(paramSet)

    % Load parameter set
    inputParameter = paramSet(paramNum); 
    load(inputParameter + ".mat", "param")

    % Extract polynomial degree
    str = extractAfter(inputParameter, "Poly"); 
    ind = strfind(str, "_");
    if ismissing(str)
        polyDeg = "0";
    elseif isempty(ind)
        polyDeg = str; 
    else
        polyDeg = extractBefore(str, "_");
    end

    % Load OCV polynomial coefficients if necessary
    if polyDeg ~= "0"
        load("OCVpoly_" + polyDeg + "_25degC.mat", "coeff")
    else
        coeff = 0; 
    end

    % Make sure SOC is increasing
    if param.SOC(1) > param.SOC(2)
        param = flip(param);
    end
    validateattributes(param.SOC, {'double'}, {'increasing'})

    % Set variable 'parameters' which is an input to the model. Use optimized parameters if exist. 
    isColumn = any(strcmp('R0_opt', param.Properties.VariableNames)); 
    if ~isColumn
        % Use initial parameters if no optimized parameters available
        [parameters.SOC, parameters.R0, parameters.R1, parameters.R2, parameters.R3, parameters.T1, parameters.T2, parameters.T3, parameters.OCV] = ...
            deal(param.SOC, param.R0_init, param.R1_init, param.R2_init, param.R3_init, param.tau1_init, param.tau2_init, param.tau3_init, param.OCV_init);
    else
        % Use optimized parameters if exist
        [parameters.SOC, parameters.R0, parameters.R1, parameters.R2, parameters.R3, parameters.T1, parameters.T2, parameters.T3, parameters.OCV] = ...
            deal(param.SOC, param.R0_opt, param.R1_opt, param.R2_opt, param.R3_opt, param.tau1_opt, param.tau2_opt, param.tau3_opt, param.OCV_opt);
    end
    
    % Check data formats
    if size(parameters.SOC) ~= size(parameters.R0)
        error("Check here")
    end

    % Run model for the current profile
    for profileNum = 1:height(profileSet) %parfor profileNum = 1:height(profileSet)
        %% Load test profiles

        inputProfile = profileSet(profileNum);

        % Load data
        switch inputProfile
            case 'US06'
                load("US06_25degC_0p8.mat", "meas");
            case 'HWFET'
                load("HWFET_25degC_0p8.mat", "meas");
            case 'UDDS'
                load("UDDS_25degC_0p8.mat", "meas");
            case 'NEDC'
                load("NEDC_25degC_0p8.mat", "meas");
            case 'MIX'
                load("MIX_25degC_0p8.mat", "meas");
        end

        % Convert struct to table
        if inputProfile ~= "Pulse"
            meas_t = struct2table(meas);
            meas_t.Current = -meas_t.Current; % Assume test data is in the opposite convension
        end

        % Formatting & calculate SOC
        Q = abs(meas_t.Ah(end) - meas_t.Ah(1));
        meas_t.Ah = meas_t.Ah - meas_t.Ah(1);
        meas_t.Time = meas_t.Time - meas_t.Time(1);
        meas_t.SOC = 1 - (-meas_t.Ah)/Q;
        meas_t = meas_t( ~any( isnan( meas_t.Current ) | isinf( meas_t.Current ), 2 ),: ); % Delete inf & NaN

        % Delete anything after Vmin was reached
        ind = find(meas_t.Voltage <= Vmin, 1);
        if ~isempty(ind)
            meas_t = meas_t(1:ind, :);
            Q = abs(meas_t.Ah(end)-meas_t.Ah(1));
        end
        
        %% Run model

        % Initial states
        initStates.SOC = meas_t.SOC(1); 
        initStates.I_RC = [0 0 0]; 
        initStates.Vt = meas_t.Voltage(1);
        
        [Vt, ~, ~] = Model_3RC(meas_t.Current, initStates, parameters, Q, deltaT, coeff);

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
        % Each parameter set run
        [rmse_param{i, 1}, rmse_param{i, 2}, rmse_param{i, 3}, rmse_param{i, 4}, rmse_param{i, 5}] = ...
            deal(rmse_temp{i,1}{1, 1}, rmse_temp{i,1}{1, 2}, rmse_temp{i,1}{1, 3}, rmse_temp{i,1}{1, 4}, rmse_temp{i,1}{1, 5}); 
        
        % All runs
        [rmse{ind, 1}, rmse{ind, 2}, rmse{ind, 3}, rmse{ind, 4}, rmse{ind, 5}] = ...
            deal(rmse_temp{i,1}{1, 1}, rmse_temp{i,1}{1, 2}, rmse_temp{i,1}{1, 3}, rmse_temp{i,1}{1, 4}, rmse_temp{i,1}{1, 5}); 
        ind = ind + 1; 
    end

    %% Plot results for the current set of parameters 

%     % Plot Vt_err
%     figurename = append('Vt_err (', string(inputParameter), ')');
%     figure('WindowStyle', 'docked', 'Name', figurename); hold on
%     legend_txt = []; 
%     for i = 1:height(rmse_temp)
%         SOC_i = rmse_temp{i, 1}{1, 5}(:, 2); 
%         Vt_err_i = rmse_temp{i, 1}{1, 5}(:, 4) - rmse_temp{i, 1}{1, 5}(:, 3); 
%         plot(SOC_i, Vt_err_i); 
% 
%         legend_txt = [legend_txt; rmse_temp{i, 1}{1, 2}];
%     end
%     grid on; title(figurename, 'Interpreter', 'none');
%     legend(legend_txt, 'Interpreter', 'none');

    % Plot Vt
    figurename = append('Vt (', string(inputParameter), ')');
    figure('WindowStyle', 'docked', 'Name', figurename);
    for i = 1:height(rmse_temp)
        axes{i, 1} = subplot(height(rmse_temp), 1, i); hold on

        SOC_i = rmse_temp{i, 1}{1, 5}(:, 2); 
        Voltage_i = rmse_temp{i, 1}{1, 5}(:, 3);
        Vt_sim_i = rmse_temp{i, 1}{1, 5}(:, 4); 
        inputProfile_i = rmse_temp{i, 1}{1, 2}; 
        rmse_i = round(rmse_temp{i, 1}{1, 3}*1000, 2); 
        max_err_i = round(rmse_temp{i, 1}{1, 4}*1000, 2); 

        plot(SOC_i, Voltage_i);
        plot(SOC_i, Vt_sim_i);

        grid on; legend(['exp'; 'sim']); 
        title_txt = inputProfile_i + " (RMSE:" + string(rmse_i) + "mV, Max error:" + string(max_err_i) + "mV)";
        title(title_txt)
        xlim([-0.01, 1.01]); 
    end

%     % Plot parameters (after)
%     figure('WindowStyle', 'docked', 'Name', inputParameter);
%     ax1 = subplot(2, 4, 1); hold on; plot(parameters.SOC, parameters.R0, '.-'); title('R0_dch', 'interpreter', 'none'); grid on; xlabel('SOC'); ylabel('R [Ohm]');
%     ax2 = subplot(2, 4, 2); hold on; plot(parameters.SOC, parameters.R1, '.-'); title('R1'); grid on; xlabel('SOC'); ylabel('R [Ohm]');
%     ax3 = subplot(2, 4, 3); hold on; plot(parameters.SOC, parameters.R2, '.-'); title('R2'); grid on; xlabel('SOC'); ylabel('R [Ohm]');
%     ax4 = subplot(2, 4, 4); hold on; plot(parameters.SOC, parameters.R3, '.-'); title('R3'); grid on; xlabel('SOC'); ylabel('R [Ohm]');
%     ax5 = subplot(2, 4, 5); hold on; plot(parameters.SOC, parameters.T1, '.-'); title('tau1'); grid on; xlabel('SOC'); ylabel('tau [s]');
%     ax6 = subplot(2, 4, 6); hold on; plot(parameters.SOC, parameters.T2, '.-'); title('tau2'); grid on; xlabel('SOC'); ylabel('tau [s]');
%     ax7 = subplot(2, 4, 7); hold on; plot(parameters.SOC, parameters.T3, '.-'); title('tau3'); grid on; xlabel('SOC'); ylabel('tau [s]');
%     ax8 = subplot(2, 4, 8); hold on; plot(parameters.SOC, parameters.OCV, '.-'); title('OCV'); grid on; xlabel('SOC'); ylabel('OCV [V]');
%     linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8], 'x')

    % Save results
    filename = folder_result + "RMSE_" + inputParameter; 
    rmse_param_noData = rmse_param(:, 1:4); 
    
    if save_rmse_of_each_parameter_set == 1
        save(filename, "rmse_param", 'rmse_param_noData');
    end

end

filename = folder_result + "RMSE"; 
rmse_noData = rmse(:, 1:4); 
if save_rmse_of_all_parameter_sets == 1
    save(filename, "rmse", 'rmse_noData');
end