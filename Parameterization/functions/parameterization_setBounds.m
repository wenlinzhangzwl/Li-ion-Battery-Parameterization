function param = parameterization_setBounds(initialGuess)

    % Set upper & lower bounds
    if isempty(initialGuess)
        % If no initial guess is available
        R_ub = 1e-3; 
        R_lb = 1e-5; 

        param.OCV_init = 3.65; 
        param.OCV_ub =  3.75;
        param.OCV_lb =  2.4;
        param.R0_ub =  R_ub;
        param.R0_lb =  R_lb;
        param.R1_ub =  R_ub;
        param.R1_lb =  R_lb;
        param.R2_ub =  R_ub;
        param.R2_lb =  R_lb;
        param.R3_ub =  R_ub;
        param.R3_lb =  R_lb;

        param.tau1_init = 2;
        param.tau1_ub = 10;
        param.tau1_lb = 0.1;

        param.tau2_init = 50;
        param.tau2_ub = 500;
        param.tau2_lb = 10;

        param.tau3_ub = 2000;
        param.tau3_lb = 500;
    else
        % Use bounds from intial guess
        load(initialGuess, 'param')

        % Update bounds based on optimization results
        param.OCV_ub =  param.OCV_init*1.01;
        param.OCV_lb =  param.OCV_init*0.99;
        param.R0_ub =  param.R0_init*1.1;
        param.R0_lb =  param.R0_init*0.9;
        param.R1_ub =  param.R1_init*1.1;
        param.R1_lb =  param.R1_init*0.9;
        param.R2_ub =  param.R2_init*1.1;
        param.R2_lb =  param.R2_init*0.9;
        param.R3_ub =  param.R3_init*1.1;
        param.R3_lb =  param.R3_init*0.9;
        param.tau1_ub = param.tau1_init*1.1;
        param.tau1_lb = param.tau1_init*0.9;
        param.tau2_ub = param.tau2_init*1.1;
        param.tau2_lb = param.tau2_init*0.9;
        param.tau3_ub = param.tau3_init*1.1;
        param.tau3_lb = param.tau3_init*0.9;
        
    end
    
    %% Plot the parameters & their bounds
    figure('WindowStyle', 'docked');

    ax1 = subplot(2, 4, 1); hold on; 
    plot(param.SOC, param.R0_init, '.-'); plot(param.SOC, param.R0_lb, '.-'); plot(param.SOC, param.R0_ub, '.-'); 
    title('R0'); grid on; xlabel('SOC'); ylabel('R [Ohm]');

    ax2 = subplot(2, 4, 2); hold on; 
    plot(param.SOC, param.R1_init, '.-'); plot(param.SOC, param.R1_lb, '.-'); plot(param.SOC, param.R1_ub, '.-'); 
    title('R1'); grid on; xlabel('SOC'); ylabel('R [Ohm]');

    ax3 = subplot(2, 4, 3); hold on; 
    plot(param.SOC, param.R2_init, '.-'); plot(param.SOC, param.R2_lb, '.-'); plot(param.SOC, param.R2_ub, '.-'); 
    title('R2'); grid on; xlabel('SOC'); ylabel('R [Ohm]');

    ax4 = subplot(2, 4, 4); hold on; 
    plot(param.SOC, param.R3_init, '.-'); plot(param.SOC, param.R3_lb, '.-'); plot(param.SOC, param.R3_ub, '.-'); 
    title('R3'); grid on; xlabel('SOC'); ylabel('R [Ohm]');

    ax5 = subplot(2, 4, 5); hold on; 
    plot(param.SOC, param.tau1_init, '.-'); plot(param.SOC, param.tau1_lb, '.-'); plot(param.SOC, param.tau1_ub, '.-'); 
    title('tau1'); grid on; xlabel('SOC'); ylabel('tau [s]');

    ax6 = subplot(2, 4, 6); hold on; 
    plot(param.SOC, param.tau2_init, '.-'); plot(param.SOC, param.tau2_lb, '.-'); plot(param.SOC, param.tau2_ub, '.-'); 
    title('tau2'); grid on; xlabel('SOC'); ylabel('tau [s]');

    ax7 = subplot(2, 4, 7); hold on; 
    plot(param.SOC, param.tau3_init, '.-'); plot(param.SOC, param.tau3_lb, '.-'); plot(param.SOC, param.tau3_ub, '.-');  
    title('tau3'); grid on; xlabel('SOC'); ylabel('tau [s]');
    
    ax8 = subplot(2, 4, 8); hold on; 
    plot(param.SOC, param.OCV_init, '.-'); plot(param.SOC, param.OCV_lb, '.-'); plot(param.SOC, param.OCV_ub, '.-'); 
    title('OCV'); grid on; xlabel('SOC'); ylabel('OCV [V]');

    linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8], 'x')

    %% Validate OCV ub is greater than OCV lb 
    OCVdiff = param.OCV_ub - param.OCV_lb; 
    errorInd = find(OCVdiff<=0, 1); 
    if ~isempty(errorInd)
        error('Check ub & lb of OCV (see errorInd)');
    end

end

