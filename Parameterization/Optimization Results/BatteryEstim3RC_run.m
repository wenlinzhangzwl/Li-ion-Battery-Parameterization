clear
% clearvars -except meas meas_t

%% Load and process data from the experiment

currentFolder = cd; 
addpath(currentFolder);

% Cell parameters
Q = 277.758;
CapacityAh = Q;
BattCapInit = CapacityAh; 
InitialCapVoltage = [0 0 0];

runCases = {'US06',    'param';... 
            'US06',    'param_opt';...
            'HWFET',     'param';... 
            'HWFET',     'param_opt';...
            'UDDS',     'param';... 
            'UDDS',     'param_opt';};

% runCases = {'pulse',    'param'};

%% Model parameters
for runNum = 1:height(runCases)
    
    clear meas meas_t param
    
    % Input current profiles
    input = runCases{runNum, 1};
    switch input
        case 'pulse'
            data_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Discharge\';
            cd(data_folder);
            load('PUL25dch_meas_t.mat');    % Raw data
        case 'US06'
            filename = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\6 - Drive cycle\US0625';
            load(append(filename, '.mat'))

            meas_t = struct2table(meas);
            ind = find(meas.Voltage <=2.5, 1);
            meas_t = meas_t(1:ind, :);
            meas_t.Ah = meas_t.Ah - meas_t.Ah(1);
            meas_t.Time = meas_t.Time - meas_t.Time(1);
            meas_t.SOC = 1 - (-meas_t.Ah)/Q;
        case 'HWFET'
            current_folder = cd;
            data_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\6 - Drive cycle\';
            cd(data_folder);
            filename = 'HWFET25';
            file = append(data_folder, filename, '.mat');
            load(file)

            meas_t = struct2table(meas);
            ind = find(meas.Voltage <=2.5, 1);
            meas_t = meas_t(1:ind, :);
            meas_t.Ah = meas_t.Ah - meas_t.Ah(1);
            meas_t.Time = meas_t.Time - meas_t.Time(1);
            meas_t.SOC = 1 - (-meas_t.Ah)/Q;
        case 'UDDS'
            current_folder = cd;
            data_folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\6 - Drive cycle\';
            cd(data_folder);
            filename = 'UDDS25';
            file = append(data_folder, filename, '.mat');
            load(file)

            meas_t = struct2table(meas);
            ind = find(meas.Voltage <=2.5, 1);
            meas_t = meas_t(1:ind, :);
            meas_t.Ah = meas_t.Ah - meas_t.Ah(1);
            meas_t.Time = meas_t.Time - meas_t.Time(1);
            meas_t.SOC = 1 - (-meas_t.Ah)/Q;
    end
    
    % Parameter sets
    paramset = runCases{runNum, 2};
    switch paramset
        case 'param'
            load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Discharge\Results\param_v6_opt_34.mat');

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
        case 'param_opt'
            load('C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\7 - Pulse Discharge\Results\param_v6_opt_34.mat');

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

            ind1 = 37-34;
            ind2 = 37;
            OCV(ind1:ind2) = param.OCV_opt(ind1:ind2);
            R0(ind1:ind2) = param.R0_opt(ind1:ind2);
            R1(ind1:ind2) = param.R1_opt(ind1:ind2);
            R2(ind1:ind2) = param.R2_opt(ind1:ind2);
            R3(ind1:ind2) = param.R3_opt(ind1:ind2);
            tau1(ind1:ind2) = param.tau1_opt(ind1:ind2);
            tau2(ind1:ind2) = param.tau2_opt(ind1:ind2);
            tau3(ind1:ind2) = param.tau3_opt(ind1:ind2);
    end
    
    % Model to run
    model = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Model\Electric Model\Parameterization\BatteryEstim3RC_PTBS_edit.slx';

    SOC_LUT = SOC;
    R0 = R0;
    Rx = [R1'; R2'; R3';];
    Tx = [tau1'; tau2'; tau3';];
    Em = OCV;


    %% Inputs
    stepSize = 0.05;
    input_time = meas_t.Time;
    input_current = -meas_t.Current;
    input_Vt = meas_t.Voltage;

    ind = find(isnan(input_current) | isinf(input_current));
    input_time(ind,:) = [];
    input_current(ind,:) = [];
    input_Vt(ind,:) = [];
    input_Vt = timeseries(input_Vt, input_time);

    runtime = input_time(end);

    %% Run model
    sim(model);

    sim_t = sim_Vt.time;
    sim_Vt_exp = sim_Vt.signals(2).values;
    sim_Vt = sim_Vt.signals(1).values;
    sim_SOC = sim_SOC.signals.values;

    rmse_Vt = sqrt(mean( (sim_Vt_exp - sim_Vt).^2 ));

    %% Plot results

    % Vt vs Time
    figurename = append('Vt ', '(', string(input), ', ', string(paramset), ')');
    figure('WindowStyle', 'docked', 'Name', figurename); 
    hold on
    plot(sim_t, sim_Vt_exp, '-');
    plot(sim_t, sim_Vt, '--');
    hold off; grid on; title(figurename, 'interpreter', 'none'); legend(['exp'; 'sim'])
    annotation('textbox', [0.14, 0.82, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt))

%     % Vt_err vs Time
%     figurename = append('Vt_err ', '(', string(input), ', ', string(paramset), ')');
%     figure('WindowStyle', 'docked', 'Name', figurename); 
%     plot(sim_t, sim_Vt - sim_Vt_exp); 
%     grid on; title(figurename, 'interpreter', 'none');
%     annotation('textbox', [0.14, 0.82, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt))
%     
%     % Vt vs SOC
%     figurename = append('Vt', '(', string(input), ', ', string(paramset), ')');
%     figure('WindowStyle', 'docked', 'Name', figurename);
%     hold on
%     plot(sim_SOC, sim_Vt_exp);
%     plot(sim_SOC, sim_Vt, '--');
%     hold off; grid on; title(figurename, 'interpreter', 'none'); legend(['exp'; 'sim'])
%     annotation('textbox', [0.14, 0.82, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt))

    rmse = zeros(35, 1);
%     SOC_flip = flip(SOC);
    SOC_flip = [1, 0.98, 0.04, min(sim_SOC)];
    for i = 1:length(SOC_flip)-1
        ind_start = find(sim_SOC <= SOC_flip(i), 1);
        ind_end = find(sim_SOC <= SOC_flip(i+1), 1);
        t_sim_i = sim_t(ind_start:ind_end);
        Vt_sim_i = sim_Vt_exp(ind_start:ind_end);
        Vt_exp2_i = sim_Vt(ind_start:ind_end);
        SOC_sim_i = sim_SOC(ind_start:ind_end);
        rmse_Vt_i = sqrt(mean((Vt_exp2_i - Vt_sim_i).^2));
        
        figurename = 'pulse' + string(i) + ' (' + string(round(SOC_flip(i), 2)) + '-' + string(round(SOC_flip(i+1), 2)) + ')';
        figure('WindowStyle', 'docked', 'Name', figurename); hold on
        plot(SOC_sim_i, Vt_exp2_i);
        plot(SOC_sim_i, Vt_sim_i, '--');
        hold off; grid on; legend(['exp'; 'sim'])
        title(append('Vt ', '(', string(input), ',', string(paramset), ')'), 'interpreter', 'none');
        annotation('textbox', [0.14, 0.82, 0.1, 0.1], 'String', "RMSE = " + string(rmse_Vt_i))
        
        rmse(i) = rmse_Vt_i;
    end

% %     Plot parameters
%     figure;
%     ax1 = subplot(2, 4, 1); hold on; plot(SOC, param.R0_init, '.-'); plot(SOC, R0, '.-'); title('R0'); grid on; xlabel('SOC'); ylabel('R [Ohm]'); legend('before', 'after')
%     ax2 = subplot(2, 4, 2); hold on; plot(SOC, param.R1_init, '.-'); plot(SOC, R1, '.-'); title('R1'); grid on; xlabel('SOC'); ylabel('R [Ohm]'); legend('before', 'after')
%     ax3 = subplot(2, 4, 3); hold on; plot(SOC, param.R2_init, '.-'); plot(SOC, R2, '.-'); title('R2'); grid on; xlabel('SOC'); ylabel('R [Ohm]'); legend('before', 'after')
%     ax4 = subplot(2, 4, 4); hold on; plot(SOC, param.R3_init, '.-'); plot(SOC, R3, '.-'); title('R3'); grid on; xlabel('SOC'); ylabel('R [Ohm]'); legend('before', 'after')
%     ax5 = subplot(2, 4, 5); hold on; plot(SOC, param.tau1_init, '.-'); plot(SOC, tau1, '.-'); title('tau1'); grid on; xlabel('SOC'); ylabel('tau [s]'); legend('before', 'after')
%     ax6 = subplot(2, 4, 6); hold on; plot(SOC, param.tau2_init, '.-'); plot(SOC, tau2, '.-'); title('tau2'); grid on; xlabel('SOC'); ylabel('tau [s]'); legend('before', 'after')
%     ax7 = subplot(2, 4, 7); hold on; plot(SOC, param.tau3_init, '.-'); plot(SOC, tau3, '.-'); title('tau3'); grid on; xlabel('SOC'); ylabel('tau [s]'); legend('before', 'after')
%     ax8 = subplot(2, 4, 8); hold on; plot(SOC, param.OCV_init, '.-'); plot(SOC, OCV, '.-'); title('OCV'); grid on; xlabel('SOC'); ylabel('OCV [V]'); legend('before', 'after')
%     linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8], 'x')
% 
%     figure; hold on
%     plot(SOC, param.R0_init, '.-');
%     plot(SOC, param.R1_init, '.-');
%     plot(SOC, param.R2_init, '.-');
%     plot(SOC, param.R3_init, '.-');
%     hold off; grid on
%     legend('R0', 'R1', 'R2', 'R3')
end
