function [Q_Ah, model_parameters] = load_parameters(cell, SOH)
    switch cell
        case "EVE280"
            Q_Ah = 280; % Nominal capacity [Ah]

            load("Parameters_EVE280_Layered_PUL_25degC_0p8_1min.mat", "param")
            model_parameters.SOC = param.SOC; 
            model_parameters.OCV = param.OCV_init; 
            model_parameters.R0 = param.R0_init;
            model_parameters.RX = [param.R1_init, param.R2_init, param.R3_init;];
            model_parameters.TX = [param.tau1_init, param.tau2_init, param.tau3_init;];
        case "CylindricalCell"
            
            load("Parameters_CylindricalCell.mat") % Cell data from "IMM_SVSF_VBL_init.m" from Sara
            model_parameters.SOC = SOC_LUT'; 

            switch SOH
                case 100
                    % week1, parameters of a third-order RC model with SoH = 100%, obtained by pulse discharge test. 
                    Q_Ah = 5.4; 
                    model_parameters.OCV = polyval(A1, model_parameters.SOC)'; 
                    model_parameters.R0 = param.model1.R0';
                    model_parameters.RX = [param.model1.R1', param.model1.R2', param.model1.R3'];
                    model_parameters.TX = [param.model1.tau1', param.model1.tau2', param.model1.tau3'];
                case 90
                    % week30, parameters of a third-order RC model with SoH = 90%, obtained by pulse discharge test. 
                    Q_Ah = 4.77;
                    model_parameters.OCV = polyval(A2, model_parameters.SOC)'; 
                    model_parameters.R0 = param.model2.R0';
                    model_parameters.RX = [param.model2.R1', param.model2.R2', param.model2.R3'];
                    model_parameters.TX = [param.model2.tau1', param.model2.tau2', param.model2.tau3'];
                case 80
                    % week45, parameters of a third-order RC model with SoH = 80%, obtained by pulse discharge test. 
                    Q_Ah = 4.33;
                    model_parameters.OCV = polyval(A3, model_parameters.SOC)'; 
                    model_parameters.R0 = param.model3.R0';
                    model_parameters.RX = [param.model3.R1', param.model3.R2', param.model3.R3'];
                    model_parameters.TX = [param.model3.tau1', param.model3.tau2', param.model3.tau3'];
            end
    end
end

