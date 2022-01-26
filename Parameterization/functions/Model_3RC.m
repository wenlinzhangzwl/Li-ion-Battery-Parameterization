function [Vt, SOC, endStates] = Model_3RC(current_dmd, initStates, parameters, cellCapacity, deltaT)

    % Convert capacity to A*s 
    Q = cellCapacity * 3600; 

    % Calculate SOC at each time step
    SOC_discharged = current_dmd * deltaT/Q; 
    SOC = initStates.SOC * ones(height(current_dmd), 1) - cumsum(SOC_discharged); 

    % Calculate parameters at each time step
    SOC_max = max(parameters.SOC);
    SOC_min = min(parameters.SOC); 
    SOC_interp = max(SOC_min, min(SOC, SOC_max));
    OCV = interp1(parameters.SOC, parameters.OCV, SOC_interp, 'linear'); 
    R0 = interp1(parameters.SOC, parameters.R0, SOC_interp, 'linear'); 
    R1 = interp1(parameters.SOC, parameters.R1, SOC_interp, 'linear'); 
    R2 = interp1(parameters.SOC, parameters.R2, SOC_interp, 'linear'); 
    R3 = interp1(parameters.SOC, parameters.R3, SOC_interp, 'linear'); 
    T1 = interp1(parameters.SOC, parameters.T1, SOC_interp, 'linear'); 
    T2 = interp1(parameters.SOC, parameters.T2, SOC_interp, 'linear'); 
    T3 = interp1(parameters.SOC, parameters.T3, SOC_interp, 'linear'); 
    R = [R1, R2, R3]; 
    T = [T1, T2, T3]; 

%     figure;
%     time = 1:height(current_dmd); 
%     plot(time, current_dmd); grid on
    
    % Initial states
    I_RC = initStates.I_RC;
    Vt = zeros(height(current_dmd), 1); 
    Vt(1) = initStates.Vt;

    V = zeros(height(current_dmd), 2);  % V = [V_R0, V_RC]
    V(1, :) = [0 0];

    % Calculate voltage at each time step
    for i = 2:height(current_dmd)

        V_R0 = R0(i)*current_dmd(i); 

        I_RC(i, :) = (1 - exp(-deltaT./T(i, :))) .* current_dmd(i-1) + exp(-deltaT./T(i, :)) .* I_RC(i-1, :);
        I_RC(i, :) = max(I_RC(i, :), 1e-100*ones(1, 3)); 
        V_RC = R(i, :) * I_RC(i, :)';

        Vt(i) = OCV(i) - V_R0 - V_RC; 
        V(i, 1) = V_R0; 
        V(i, 2) = V_RC; 
    end

    endStates.I_RC = I_RC(end, :);
    endStates.Vt = Vt(end);
    endStates.SOC = SOC(end);

    errInd = find(isinf(Vt), 1); 
    if ~isempty(errInd)
        error("Check errInd")
    end
end

