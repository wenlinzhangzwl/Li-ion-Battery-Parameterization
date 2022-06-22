
clear

file = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\1 - Prismatic 280Ah Battery Pack Design\Test Results\25degC\Results\RMSE_ALL.mat";

load(file, "rmse_noData")

paramName = string(rmse_noData(:, 1)); 
profileName = string(rmse_noData(:, 2)); 
rmse = cell2mat(rmse_noData(:, 3)) * 1000;
maxErr = cell2mat(rmse_noData(:, 4)) * 1000; 

%% Pulse (Layered)
paramName_pulses = ["Parameters_Layered_PUL_25degC_0p8_1min";...
                    "Parameters_Layered_PUL_25degC_0p8_5min";...
                    "Parameters_Layered_PUL_25degC_0p8_15min";...
                    "Parameters_Layered_PUL_25degC_0p8_30min";...
                    "Parameters_Layered_PUL_25degC_0p8_60min";...
                    "Parameters_Layered_PUL_25degC_0p8_120min";];  % Make sure order is the same as x

y_rmse = zeros(height(paramName_pulses), 5); 
y_rmse_avg = zeros(height(paramName_pulses), 1);
y_maxErr = zeros(height(paramName_pulses), 5); 
y_maxErr_avg = zeros(height(paramName_pulses), 1);

for i = 1:height(paramName_pulses)
    ind = find(paramName == paramName_pulses(i));

    y_rmse(i, :) = transpose(rmse(ind));
    y_rmse_avg(i) = mean(y_rmse(i, :));

    y_maxErr(i, :) = transpose(maxErr(ind)); 
    y_maxErr_avg(i) = mean(y_maxErr(i, :));
end

cases = ["1 min", "5 min", "15 min", "30 min", "60 min", "120 min"]; 
x = categorical(cases); % Profile optimized with
x = reordercats(x, cases); % Profile validated on

% Plot RMSE
figure("WindowStyle", "docked"); hold on
b = bar(x, y_rmse);
% for i = 1:5
%     xtips = b(i).XEndPoints;
%     ytips = b(i).YEndPoints;
%     labels = string(round(b(i).YData, 1));
%     text(xtips,ytips,labels,'HorizontalAlignment','center', 'VerticalAlignment','bottom')
% end
plot(x, y_rmse_avg, '.-', 'MarkerSize',20)
grid on
ylim([0, 35]);
title("RMSE with Different Pulse Lengths (Analytical Method)")
legend('US06', 'HWFET', 'UDDS', 'MIXED', 'NEDC', 'Average')
xlabel("Profile Optimized With"); ylabel("RMSE [mV]"); 

% Plot max error
figure("WindowStyle", "docked"); hold on
b = bar(x, y_maxErr);
% for i = 1:5
%     xtips = b(i).XEndPoints;
%     ytips = b(i).YEndPoints;
%     labels = string(round(b(i).YData, 1));
%     text(xtips,ytips,labels,'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom')
% end
plot(x, y_maxErr_avg, '.-', 'MarkerSize',20);
grid on
ylim([0, 400]);
title('Max Error with Different Pulse Lengths (Analytical Method)')
legend('US06', 'HWFET', 'UDDS', 'MIXED', 'NEDC', 'Average')
xlabel("Profile Optimized With"); ylabel("Max Error [mV]"); 

%% Pulse (optimized directly)
paramName_pulses = ["Parameters_Optimized_PUL_25degC_0p8_1min_PSO_Poly0";...
                    "Parameters_Optimized_PUL_25degC_0p8_5min_PSO_Poly0";...
                    "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly0";...
                    "Parameters_Optimized_PUL_25degC_0p8_30min_PSO_Poly0";...
                    "Parameters_Optimized_PUL_25degC_0p8_60min_PSO_Poly0";...
                    "Parameters_Optimized_PUL_25degC_0p8_120min_PSO_Poly0";];  % Make sure order is the same as x

y_rmse = zeros(height(paramName_pulses), 5); 
y_rmse_avg = zeros(height(paramName_pulses), 1);
y_maxErr = zeros(height(paramName_pulses), 5); 
y_maxErr_avg = zeros(height(paramName_pulses), 1);

for i = 1:height(paramName_pulses)
    ind = find(paramName == paramName_pulses(i));

    y_rmse(i, :) = transpose(rmse(ind));
    y_rmse_avg(i) = mean(y_rmse(i, :));

    y_maxErr(i, :) = transpose(maxErr(ind)); 
    y_maxErr_avg(i) = mean(y_maxErr(i, :));
end

cases = ["1 min", "5 min", "15 min", "30 min", "60 min", "120 min"]; 
x = categorical(cases); % Profile optimized with
x = reordercats(x, cases); % Profile validated on

% Plot RMSE
figure("WindowStyle", "docked"); hold on
b = bar(x, y_rmse);
% for i = 1:5
%     xtips = b(i).XEndPoints;
%     ytips = b(i).YEndPoints;
%     labels = string(round(b(i).YData, 1));
%     text(xtips,ytips,labels,'HorizontalAlignment','center', 'VerticalAlignment','bottom')
% end
plot(x, y_rmse_avg, '.-', 'MarkerSize',20)
grid on
ylim([0, 35]);
title("RMSE with Different Pulse Lengths (Direct Optimization Method)")
legend('US06', 'HWFET', 'UDDS', 'MIXED', 'NEDC', 'Average')
xlabel("Profile Optimized With"); ylabel("RMSE [mV]"); 

% Plot max error
figure("WindowStyle", "docked"); hold on
b = bar(x, y_maxErr);
% for i = 1:5
%     xtips = b(i).XEndPoints;
%     ytips = b(i).YEndPoints;
%     labels = string(round(b(i).YData, 1));
%     text(xtips,ytips,labels,'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom')
% end
plot(x, y_maxErr_avg, '.-', 'MarkerSize',20);
grid on
ylim([0, 400]);
title('Max Error with Different Pulse Lengths (Direct Optimization Method)')
legend('US06', 'HWFET', 'UDDS', 'MIXED', 'NEDC', 'Average')
xlabel("Profile Optimized With"); ylabel("Max Error [mV]"); 

%% Different drive cycle optimization

paramName_pulses = ["Parameters_Optimized_US06_25degC_0p8_PSO_Poly0";...
                    "Parameters_Optimized_HWFET_25degC_0p8_PSO_Poly0";...
                    "Parameters_Optimized_UDDS_25degC_0p8_PSO_Poly0";...
                    "Parameters_Optimized_MIX_25degC_0p8_PSO_Poly0";...
                    "Parameters_Optimized_NEDC_25degC_0p8_PSO_Poly0";];  % Make sure order is the same as x

y_rmse = zeros(height(paramName_pulses), 5); 
y_rmse_avg = zeros(height(paramName_pulses), 1);
y_maxErr = zeros(height(paramName_pulses), 5); 
y_maxErr_avg = zeros(height(paramName_pulses), 1);

for i = 1:height(paramName_pulses)
    ind = find(paramName == paramName_pulses(i));

    y_rmse(i, :) = transpose(rmse(ind));
    y_rmse_avg(i) = mean(y_rmse(i, :));

    y_maxErr(i, :) = transpose(maxErr(ind)); 
    y_maxErr_avg(i) = mean(y_maxErr(i, :));
end

cases = ["US06", "HWFET", "UDDS", "MIXED", "NEDC" ]; 
x = categorical(cases); % Profile optimized with
x = reordercats(x, cases); % Profile validated on

% Plot RMSE
figure("WindowStyle", "docked"); hold on; 
b = bar(x, y_rmse);
% for i = 1:5
%     xtips = b(i).XEndPoints;
%     ytips = b(i).YEndPoints;
%     labels = string(round(b(i).YData, 1));
%     text(xtips,ytips,labels,'HorizontalAlignment','center', 'VerticalAlignment','bottom')
% end
plot(x, y_rmse_avg, '.-', 'MarkerSize',20)
grid on
ylim([0, 35]);
title('RMSE with Different Drive Cycles')
legend('US06', 'HWFET', 'UDDS', 'MIXED', 'NEDC', 'Average')
xlabel("Profile Optimized With"); ylabel("RMSE [mV]"); 

% Plot max error
figure("WindowStyle", "docked"); hold on; 
b = bar(x, y_maxErr);
% for i = 1:5
%     xtips = b(i).XEndPoints;
%     ytips = b(i).YEndPoints;
%     labels = string(round(b(i).YData, 1));
%     text(xtips,ytips,labels,'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom')
% end
plot(x, y_maxErr_avg, '.-', 'MarkerSize',20);
grid on
ylim([0, 400]);
title('Max Error with Different Drive Cycles')
legend('US06', 'HWFET', 'UDDS', 'MIXED', 'NEDC', 'Average')
xlabel("Profile Optimized With"); ylabel("Max Error [mV]"); 

%% Different polynomial degrees
paramName_poly = ["Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly0";...
                "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly5";...
                "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly6";...
                "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly7";...
                "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly8";...
                "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly9";...
                "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly10";...
                "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly11";...
                "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly12";...
                "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly13";...
                "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly14";...
                "Parameters_Optimized_PUL_25degC_0p8_15min_PSO_Poly15";];  % Make sure order is the same as x

y_rmse = zeros(height(paramName_poly), 5); 
y_rmse_avg = zeros(height(paramName_poly), 1);
y_maxErr = zeros(height(paramName_poly), 5); 
y_maxErr_avg = zeros(height(paramName_poly), 1);

for i = 1:height(paramName_poly)
    ind = find(paramName == paramName_poly(i));

    y_rmse(i, :) = transpose(rmse(ind));
    y_rmse_avg(i) = mean(y_rmse(i, :));

    y_maxErr(i, :) = transpose(maxErr(ind)); 
    y_maxErr_avg(i) = mean(y_maxErr(i, :));
end

cases = ["0", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"]; 
x = categorical(cases); % Profile optimized with
x = reordercats(x, cases); % Profile validated on

% Plot RMSE
figure("WindowStyle", "docked"); hold on
b = bar(x, y_rmse);
% for i = 1:5
%     xtips = b(i).XEndPoints;
%     ytips = b(i).YEndPoints;
%     labels = string(round(b(i).YData, 1));
%     text(xtips,ytips,labels,'HorizontalAlignment','center', 'VerticalAlignment','bottom')
% end
plot(x, y_rmse_avg, '.-', 'MarkerSize',20)
grid on
ylim([0, 35]);
title("RMSE with Different Polynomial Degrees")
legend('US06', 'HWFET', 'UDDS', 'MIXED', 'NEDC', 'Average')
xlabel("OCV Polynomial Degree"); ylabel("RMSE [mV]"); 

% Plot max error
figure("WindowStyle", "docked"); hold on
b = bar(x, y_maxErr);
% for i = 1:5
%     xtips = b(i).XEndPoints;
%     ytips = b(i).YEndPoints;
%     labels = string(round(b(i).YData, 1));
%     text(xtips,ytips,labels,'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom')
% end
plot(x, y_maxErr_avg, '.-', 'MarkerSize',20);
grid on
ylim([0, 400]);
title('Max Error with Different Polynomial Degrees')
legend('US06', 'HWFET', 'UDDS', 'MIXED', 'NEDC', 'Average')
xlabel("OCV Polynomial Degree"); ylabel("Max Error [mV]"); 