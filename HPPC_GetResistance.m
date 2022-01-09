clear

%% Load file
current_folder=cd;  
folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Prismatic Cell\Test Results\1 - Characterization test 1';
cd(folder);
file_array = dir('*.mat');
addpath(current_folder);
file = append(folder, '\12-11-20_08.01 1336_Charge3_HPPC.mat');
load(file)

%% Test Conditions
Q = 285;
% current = [0.25,0.5,0.75,1,0.25,0.5,0.75,1];
% maxvolt = 3.6;
% minvolt = 2.55;
% 
% HPPC={[],'Vdis[V]',[],[],[],'Vchg[V]',[],[],[],[],'Rdis[ohm]',[],[],[],'Rchg[ohm]',[],[],[];'SOC',...
%       current(1),current(2),current(3),current(4),current(5),current(6),current(7),current(8),[],...
%       current(1),current(2),current(3),current(4),current(5),current(6),current(7),current(8)};         % Stores R0 extracted

%% Find the beginning & end all segments in a step

i025_1 = 2;     % First pulse at 0.25C
i025_2 = 6;     % Third pulse at 0.25C
i050_1 = 12;    % First pulse at 0.50C
i050_2 = 16;    % Third pulse at 0.50C
i075_1 = 22;    % First pulse at 0.75C
i075_2 = 26;    % Third pulse at 0.75C
i100_1 = 32;    % First pulse at 1.00C
i100_2 = 36;    % Third pulse at 1.00C

steps = [i025_1 i025_2 i050_1 i050_2 i075_1 i075_2 i100_1 i100_2];

for i = 1:length(steps)
    ind{i,1} = steps(i);
    ind{i,2} = find(meas.Step == steps(i));
    indices = ind{i,2};             % Indices of data points with this step number
    segments = [indices(1) 0];      % Start & end points of each segment with this start number
    k = 1;                          % Counter
    
    for j = 1:length(indices)
        if j == length(indices)
            segments(k,2) = indices(j);
            break
        elseif indices(j+1) ~= indices(j)+1
            segments(k,2) = indices(j);
            segments(k+1,1) = indices(j+1);
            k = k+1;
        end
    end
    
    ind{i,3} = segments; 
end

%% Calculate the SOCs tested
SOCPt = ind{1,3}(:, 1);     % Point at which SOC is calculated
SOC = 1 + meas.Ah(SOCPt)/Q;

%% Calculate charging/discharging resistances at each SOC at each current rate

for i = 1:length(SOC)
    Resistance1{i,1} = SOC(i);      % First pulse
    Resistance2{i,1} = SOC(i);      % Third pulse
    Resistance{i,1} = SOC(i);       % Average resistance
end

% Calculate resistance for each pulse
for i = 1:length(steps)
    for j = 1:length(SOC)
        % Start & end points of the segment
        Pt1 = ind{i,3}(j,1);
        Pt2 = ind{i,3}(j,2);
        
        % Charging & discharging resistances
        R = [(meas.Voltage(Pt1+1) - meas.Voltage(Pt1))/meas.Current(Pt1+1),...
             (meas.Voltage(Pt2+1) - meas.Voltage(Pt2))/meas.Current(Pt2)];
        Rchg = R(R>=0); 
        Rdch = R(R<0); 
        
        % Fill the calculated resistances into Resistance1 & Resistance2
        step = ind{i,1};
        switch step
            case 2                          % 0.25C
                Resistance1{j, 3} = Rchg;    % Charging resistance at 0.25C
                Resistance1{j, 8} = -Rdch;    % Discharging resistance at 0.25C
            case 6                          % 0.25C
                Resistance2{j, 3} = Rchg;
                Resistance2{j, 8} = -Rdch;
            case 12                         % 0.5C
                Resistance1{j, 4} = Rchg;
                Resistance1{j, 9} = -Rdch;
            case 16                         % 0.5C
                Resistance2{j, 4} = Rchg;
                Resistance2{j, 9} = -Rdch;
            case 22                         % 0.75C
                Resistance1{j, 5} = Rchg;
                Resistance1{j, 10} = -Rdch;
            case 26                         % 0.75C
                Resistance2{j, 5} = Rchg;
                Resistance2{j, 10} = -Rdch;
            case 32                         % 1C
                Resistance1{j, 6} = Rchg;
                Resistance1{j, 11} = -Rdch;
            case 36                         % 1C
                Resistance2{j, 6} = Rchg;
                Resistance2{j, 11} = -Rdch;
        end
    end
end

% Calculate the average of Resistance1 & Resistance2
size = size(Resistance1); 
for i = 1:size(1)
    for j = 3:11
        R1 = Resistance1{i,j}; 
        R2 = Resistance2{i,j};
        if isnumeric(R1) & ~isinf(R1) & ~isnan(R1)
            if isnumeric(R2) & ~isinf(R2) & ~isnan(R2)
                R = (R1+R2)/2;
            else
                R = R1;
            end
        elseif isnumeric(R2) & ~isinf(R2) & ~isnan(R2)
            R = R2;
        else
            R = NaN;
        end
        Resistance{i,j} = R;
    end
end

% Export Resistances as excel and a look up table
filename = append(folder, '\HPPC_GetResistance.mat');
save(filename, 'Resistance');

%% Plot Resistances
Rchg025 = [];
Rchg050 = [];
Rchg075 = [];
Rchg100 = [];
Rdch025 = [];
Rdch050 = [];
Rdch075 = [];
Rdch100 = [];

for i = 1:length(SOC)
    Rchg025 = [Rchg025, Resistance{i,3}];
    Rchg050 = [Rchg050, Resistance{i,4}];
    Rchg075 = [Rchg075, Resistance{i,5}];
    Rchg100 = [Rchg100, Resistance{i,6}];
    Rdch025 = [Rdch025, Resistance{i,8}];
    Rdch050 = [Rdch050, Resistance{i,9}];
    Rdch075 = [Rdch075, Resistance{i,10}];
    Rdch100 = [Rdch100, Resistance{i,11}];
end

figure
hold on
plot(SOC, Rchg025)
plot(SOC, Rchg050)
plot(SOC, Rchg075)
plot(SOC, Rchg100)
hold off
grid minor
legend('0.25', '0.5', '0.75', '1')
title('Rchg')

figure
hold on
plot(SOC, Rdch025)
plot(SOC, Rdch050)
plot(SOC, Rdch075)
plot(SOC, Rdch100)
hold off
grid minor
legend('0.25', '0.5', '0.75', '1')
title('Rdch')
