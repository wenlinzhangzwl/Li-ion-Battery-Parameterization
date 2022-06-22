function obj = parameterization_Objective(param_opt, param_prev, param_SOC, BatteryModel, data, initStates, cellCapacity, deltaT, optimMethod, OCVcoeff)
    
    parameters.SOC = param_SOC; 

    if isempty(param_prev)
        parameters.OCV = [param_opt(1); param_opt(9)];
        parameters.R0 = [param_opt(2); param_opt(10)];
        parameters.R1 = [param_opt(3); param_opt(11)];
        parameters.R2 = [param_opt(4); param_opt(12)];
        parameters.R3 = [param_opt(5); param_opt(13)];
        parameters.T1 = [param_opt(6); param_opt(14)];
        parameters.T2 = [param_opt(7); param_opt(15)];
        parameters.T3 = [param_opt(8); param_opt(16)];
    else
        parameters.OCV = [param_prev(1); param_opt(1)];
        parameters.R0 = [param_prev(2); param_opt(2)];
        parameters.R1 = [param_prev(3); param_opt(3)];
        parameters.R2 = [param_prev(4); param_opt(4)];
        parameters.R3 = [param_prev(5); param_opt(5)];
        parameters.T1 = [param_prev(6); param_opt(6)];
        parameters.T2 = [param_prev(7); param_opt(7)];
        parameters.T3 = [param_prev(8); param_opt(8)];
    end

    current_dmd = data.Current; 
    [Vt, ~, ~] = BatteryModel(current_dmd, initStates, parameters, cellCapacity, deltaT, OCVcoeff); 

    error_Vt = data.Voltage - Vt; 

    if optimMethod == "LS" % Return vector obj function
        obj = error_Vt; 
    else % Return scalar obj function

%         % Weight dynamic & rest portions proportional to the # of data points
%         indDyn = find(data.Current ~= 0);
%         indRest = find(data.Current == 0);
%         numDyn = height(indDyn);
%         numRest = height(indRest);
%         ratio = numRest/numDyn; 
%         error_Vt(indDyn) = error_Vt(indDyn) * ratio; 

        obj = sum(error_Vt.^2);
    end
    
end