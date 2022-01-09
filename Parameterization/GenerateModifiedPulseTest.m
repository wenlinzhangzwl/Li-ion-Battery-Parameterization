clear

Q_Ah = 280;  % Capacity
Q_As = Q_Ah * 3600; 
SOC = [0:0.01:0.1 0.15:0.05:0.9 0.91:0.01:1]';

% Current magnitude
CRate = 0.8;
currentMag = Q_Ah * CRate; 

% Pulse length for 1% & 5% SOC
pulseTime_1 = 0.01 * Q_As / currentMag ; 
pulseTime_5 = 0.05 * Q_As / currentMag; 
restTime_short = 2;
restTime_long = 2 * 3600;

% Generate current profile
n = 5; 
current = [];
deltaT = 0.1; 

for i = 2:height(SOC)
    if SOC(i)-SOC(i-1) == 0.01
        pulseTime = pulseTime_1; 
    else
        pulseTime = pulseTime_5; 
    end
    
    pulse = ones(pulseTime/n/deltaT, 1); 
    rest_short = zeros(restTime_short/deltaT, 1); 
    rest_long = zeros(restTime_long/deltaT, 1); 
    
    for k = 1:n
        current = [current; currentMag*pulse; rest_short];
    end
    
    current = [current; rest_long]; 
end

time = 0:0.1:(height(current)-1)*0.1;
time = time'; 

figure; 
plot(time, current); 
grid on; 

meas.Time = time; 
meas.Current = current; 

filename = "C:\Users\Wenlin\OneDrive\SCHOOL\Projects\48V Battery Pack Design\Test Results\9 - Modified pulse test\ModifiedPulse"; 
save(filename, "meas")