function [param] = getInitialcondition(meas_t, pulse, dch)
% Returns a table ('param') containing the initial guess of each parameters & their upper/lower bounds

param = table('Size', [height(pulse)+2, 8*3],...
              'VariableTypes', ["double", "double", "double", "double", "double", "double", "double", "double",...
                                "double", "double", "double", "double", "double", "double", "double", "double",...
                                "double", "double", "double", "double", "double", "double", "double", "double",],...
              'VariableNames',["R0_init", "R1_init", "R2_init", "R3_init", "tau1_init", "tau2_init", "tau3_init", "OCV_init",...
                               "R0_lb", "R1_lb", "R2_lb", "R3_lb", "tau1_lb", "tau2_lb", "tau3_lb", "OCV_lb",... 
                               "R0_ub", "R1_ub", "R2_ub", "R3_ub", "tau1_ub", "tau2_ub", "tau3_ub", "OCV_ub",]);


% Parameters at the extreme ends
if dch == 1
    param.OCV_init(1) = 3.65;
    param.OCV_init(end) = 2.5;
    param.SOC = [1; pulse.SOC; 0];
elseif dch == 0
    param.OCV_init(1) = 2.5;
    param.OCV_init(end) = 3.65;
    param.SOC = [0; pulse.SOC; 1];
end

for i = 1:height(pulse)
    
    ind1 = pulse.load{i}(1,1) - 1;
    ind2 = pulse.load{i}(1,1);
    ind3 = pulse.load{i}(1,2);
    ind4 = ind3 + 1;
    
    if i ~= height(pulse)
        ind5 = pulse.relaxation{i+1}(1,2);
        ind6 = ind5 + 1;
    else
        ind5 = height(meas_t);
        ind6 = NaN;
    end
    
    %% Calculate OCV
    param.OCV_init(i+1) = meas_t.Voltage(ind5);
    
    %% Calculate R0
    param.R0_init(i+1) = getR0(ind3, ind4, ind5, ind6, meas_t);
    if i == 1
        param.R0_init(i) = param.R0_init(i+1);
    end
    
%     %% Calculate tau   
%     Vfcn = fittype(...
%                     @(tau1, tau2, tau3, V1i, V2i, V3i, t)...
%                     (V1i*exp(-t/tau1) + V2i*exp(-t/tau2) + V3i*exp(-t/tau3) - V1i - V2i - V3i),...
%                     'independent', 't'...
%                    );
%                
%     Vt_relax = meas_t.Voltage(ind3:ind5);
%     Vt_relax = Vt_relax - Vt_relax(1);  % Zero Vexp since curve fitting starts at 0 
%     
%     t_relax = meas_t.Time(ind3:ind5);
%     t_relax = t_relax - t_relax(1);
%     
%     coeff = fit(t_relax, Vt_relax, Vfcn, 'Lower', [0,10,100], 'Upper', [100, 1000, 10000], 'StartPoint', [1, 50, 1000, 0, 0, 0] );
%     
% %     % Plot the the fit for validation
% %     figure; 
% %     plot(coeff, '*', t, Vt_relax_exp, '.')
% 
%     param.tau1_init(i+1) = coeff.tau1;
%     param.tau2_init(i+1) = coeff.tau2;
%     param.tau3_init(i+1) = coeff.tau3;
%     
%     if i == 1
%         param.tau1_init(i) = param.tau1_init(i+1);
%         param.tau2_init(i) = param.tau2_init(i+1);
%         param.tau3_init(i) = param.tau3_init(i+1);
%     end
%     
%     %% Calculate Rx
%     
%     % Load segment
%     segment = meas_t(ind2:ind3, :);
%     
%     % Calcualte deltaT
%     deltaT = zeros(height(segment), 1);
%     deltaT(1) = meas_t.Time(ind1) - meas_t.Time(ind1 - 1);
%     for j = 2:height(segment)
%         deltaT(j, 1) = segment.Time(j) - segment.Time(j-1);
%     end
%     
%     % Calculate I1, I2, I3 (assuming they start at 0)
%     tau1 = param.tau1_init(i);
%     tau2 = param.tau2_init(i);
%     tau3 = param.tau3_init(i);
%     for j = 1:height(segment)
%         if j ~= 1
%             segment.I1(j) = ( 1 - exp(-deltaT(j) / tau1) ) * segment.Current(j) + exp(-deltaT(j) / tau1) * segment.I1(j-1);
%             segment.I2(j) = ( 1 - exp(-deltaT(j) / tau2) ) * segment.Current(j) + exp(-deltaT(j) / tau2) * segment.I2(j-1);
%             segment.I3(j) = ( 1 - exp(-deltaT(j) / tau3) ) * segment.Current(j) + exp(-deltaT(j) / tau3) * segment.I3(j-1);
%         else
%             segment.I1(j) = ( 1 - exp(-deltaT(j) / tau1) ) * segment.Current(j);
%             segment.I2(j) = ( 1 - exp(-deltaT(j) / tau2) ) * segment.Current(j);
%             segment.I3(j) = ( 1 - exp(-deltaT(j) / tau3) ) * segment.Current(j);
%         end
%     end
% %     figName = 'Pulse ' + string(i);
% %     figure('Name', figName, 'WindowStyle', 'docked');
% %     hold on
% %     ax1 = subplot(4, 1, 1); plot(segment.Time, segment.Current, '.-'); title('I'); grid on
% %     ax2 = subplot(4, 1, 2); plot(segment.Time, segment.I1, '.-'); title('I1'); grid on
% %     ax3 = subplot(4, 1, 3); plot(segment.Time, segment.I2, '.-'); title('I2'); grid on
% %     ax4 = subplot(4, 1, 4); plot(segment.Time, segment.I3, '.-'); title('I3'); grid on
% %     linkaxes([ax1,ax2,ax3,ax4], 'x')
% 
%     % Calculate C
%     C = zeros(height(segment.Voltage), height(SOC)*3);
%     
%     SOC1 = SOC(i);    % previous SOC breakpoint
%     SOC2 = SOC(i+1);  % next SOC breakpoint
%     
%     column_R1 = i;
%     column_R2 = column_R1 + 38; 
%     column_R3 = column_R2 + 38; 
%     
%     interp1 = (SOC2 - segment.SOC)/(SOC2 - SOC1);
%     interp2 = (segment.SOC - SOC1)/(SOC2 - SOC1);
%     
%     C(:, column_R1) = interp1 .* segment.I1;
%     C(:, column_R1+1) = interp2 .* segment.I1;
%     C(:, column_R2) = interp1 .* segment.I2;
%     C(:, column_R2+1) = interp2 .* segment.I2;
%     C(:, column_R3) = interp1 .* segment.I3;
%     C(:, column_R3+1) = interp2 .* segment.I3;
%     
%     % Solve for Rx
%     Vt_load = segment.Voltage - segment.Voltage(1);
%     lb = zeros(height(Vt_load), 1);
%     [lb(column_R1), lb(column_R1+1), lb(column_R2), lb(column_R2+1), lb(column_R3), lb(column_R3+1)] = deal(1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5);
%     ub = zeros(height(Vt_load), 1);
%     [ub(column_R1), ub(column_R1+1), ub(column_R2), ub(column_R2+1), ub(column_R3), ub(column_R3+1)] = deal(0.2, 0.2, 0.2, 0.2, 0.2, 0.2);
%     x_opt = lsqlin(C, Vt_load, [], [], [], [], lb, ub);
%     
%     % Plot results for validation
%     figName = 'Pulse ' + string(i);
%     figure('Name', figName, 'WindowStyle', 'docked');
%     hold on
%     plot(segment.Time, Vt_load, '.-')
%     Vt_fit = C * x_opt; 
%     plot(segment.Time, Vt_fit, '.-')
%     legend('data', 'fit'); grid on
%     
%     %% Write Rx estimations to 'param'
%     if i ~= 1
%         param.R1_init(i+1) = x_opt(column_R1+1);
%         param.R2_init(i+1) = x_opt(column_R2+1);
%         param.R3_init(i+1) = x_opt(column_R3+1);
%     else
%         param.R1_init(i) = x_opt(column_R1);
%         param.R2_init(i) = x_opt(column_R2);
%         param.R3_init(i) = x_opt(column_R3);
%         
%         param.R1_init(i+1) = x_opt(column_R1+1);
%         param.R2_init(i+1) = x_opt(column_R2+1);
%         param.R3_init(i+1) = x_opt(column_R3+1);
%     end
end

param(end, 1:7) = param(end-1, 1:7);
end

function [R0] = getR0(ind3, ind4, ind5, ind6, data)
    if isenum(ind6)
        R0a = (data.Voltage(ind4)-data.Voltage(ind3))/(data.Current(ind4)-data.Current(ind3));
        R0b = (data.Voltage(ind6)-data.Voltage(ind5))/(data.Current(ind6)-data.Current(ind5));
        R0 = (R0a + R0b)/2;
    else
        R0 = (data.Voltage(ind4)-data.Voltage(ind3))/(data.Current(ind4)-data.Current(ind3));
    end
end