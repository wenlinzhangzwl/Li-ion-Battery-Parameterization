function [pulse] = getPulse(meas_t, ind)
% Returns a table ('pulse') that contains:
% (refer to 'Pulse Identification Reference' for references of point number)
%     1. load:          Vector containing the start & end points of a load event
%     2. relaxation:    Vector containing the start & end points of a relaxation event
%     3. segment1_t:    All points in a relaxation segment
%     4. segment1:      Resampled points in a relaxation segment
%     5. segment2_t:    All points in a load-relaxation segment
%     6. segment2:      Resampled points in a load-relaxation segment
%     7. SOC:           SOC after each load event (aka SOC during relaxation)
%     8. param:         Initial guess of the parameters, will be filled later

% pulse: m x 8 table where m is the number of available load-relax segments
% meas_t: Table containing raw data (meas_t)
% ind: Vector containing the start & end index of each pulse & relaxation segment in meas_t

pulse = table('Size', [1, 9],...
               'VariableTypes', ["cell", "cell", "cell", "cell", "cell", "cell", "double", "cell", "cell"],...
               'VariableNames',["load", "relaxation", "segment1_t", "segment1", "segment2_t", "segment2", "SOC", "param", "Optim"]);    

% Start & end indices of all load (ind1) & relaxation (ind2) periods
ind1 = [];
ind2 = [];
for i = 1:2:(height(ind)-2)
    ind2 = [ind2; ind.segment{i}];    % Relaxation
    ind1 = [ind1; ind.segment{i+1}];  % Load
end
ind2 = [ind2; ind.segment{height(ind)-1}; ind.segment{height(ind)}]; % Add the last relaxation pulse

% % Validate the indices found
% figure; 
% ax1 = subplot(2, 1, 1); hold on
% plot(meas_t.Time, meas_t.Voltage, '-')
% plot(meas_t.Time(ind1(:,1)), meas_t.Voltage(ind1(:,1)), '*r')
% plot(meas_t.Time(ind1(:,2)), meas_t.Voltage(ind1(:,2)), '*r')
% plot(meas_t.Time(ind2(:,1)), meas_t.Voltage(ind2(:,1)), '*g')
% plot(meas_t.Time(ind2(:,2)), meas_t.Voltage(ind2(:,2)), '*g')
% ax2 = subplot(2, 1, 2); hold on
% plot(meas_t.Time, meas_t.Current, '-')
% plot(meas_t.Time(ind1(:,1)), meas_t.Current(ind1(:,1)), '*r')
% plot(meas_t.Time(ind1(:,2)), meas_t.Current(ind1(:,2)), '*r')
% plot(meas_t.Time(ind2(:,1)), meas_t.Current(ind2(:,1)), '*g')
% plot(meas_t.Time(ind2(:,2)), meas_t.Current(ind2(:,2)), '*g')
% linkaxes([ax1,ax2], 'x')

% Divide data into segments
num = length(ind1);
for i = 1:num
    
    % Get segment1 for time constant curve fitting
    [segment1_t, segment1] = getSeg1(meas_t, ind1, ind2, i);
    pulse.segment1_t(i) = {segment1_t};
    pulse.segment1(i) = {segment1};
    
    % Get segment2 for parameter optimization
    [segment2_t, segment2] = getSeg2(meas_t, ind1, ind2, i, num);
    pulse.segment2_t(i) = {segment2_t};
    pulse.segment2(i) = {segment2};
    
    % Record pulse info 
    pulse.load(i) = {ind1(i,:)};
    pulse.relaxation(i) = {ind2(i, :)};
    pulse.SOC(i) = segment1_t.SOC(2);
    
end

end


function [segment1_t, segment1] = getSeg1(meas_t, ind1, ind2, i)
% Divide into relaxation segments & resample for curve fitting (segment_t1, pulse.segment1)

iBeg = ind1(i, 2);
iEnd = ind2(i+1, 2);
segment1_t = meas_t(iBeg:iEnd, :);

div1 = 100;
div2 = 200;
div3 = 250;
inc1 = 1;
inc2 = 40;
inc3 = fix((height(segment1_t) - div1*inc1 - div2*inc2) / div3);
div = div1 + div2 + div3 + 2; % Total number of points sampled
p = 1;
for k = 1:div
    if k <= div1
        segment1(k,:) = segment1_t(p,:); %#ok<*AGROW>
        p = p + inc1;
    elseif k>div1 && k<=div1+div2 
        segment1(k,:) = segment1_t(p,:);
        p = p + inc2;
    elseif k>div1+div2 && k <= div-2
        segment1(k,:) = segment1_t(p,:);
        p = p + inc3; 
    elseif k == div-1
        segment1(k,:) = segment1_t(end-1,:);
    elseif k == div
        segment1(k,:) = segment1_t(end,:);
    end
end
    
end

function [segment2_t, segment2] = getSeg2(meas_t, ind1, ind2, i, num)
% Divide into load-relax segments & resample for optimization (segment_t2, pulse.segment2)

    % Entire pulse segment
    iBeg = ind1(i, 1);
    iEnd = ind2(i+1, 2);
    segment2_t = meas_t(iBeg:iEnd, :);
    
    % Segment layered with previous segment
    segment2_overlap1 = meas_t(ind1(i,1)-20:1:ind1(i,1)-1, :);
    
    % Segment layered with next segment
    if i ~= height(ind2)-1
        segment2_overlap2 = meas_t(ind2(i+1,2)+1:1:ind2(i+1,2)+20, :);
    else
        segment2_overlap2 = [];
    end
    
    % Load segment
    segment2_load = meas_t(ind1(i,1):ind1(i,2), :);

    % Resample relexation segment
    segment2_relax_t = meas_t(ind2(i+1,1):ind2(i+1, 2), :);
    div1 = 100;
    div2 = 1000;
    div3 = 5000;
    inc1 = 1;
    inc2 = 10;
    inc3 = fix((height(segment2_relax_t) - div1*inc1 - div2*inc2) / div3);
    div3 = fix((height(segment2_relax_t) - div1*inc1 - div2*inc2) / inc3);
    div = div1 + div2 + div3 + 1; % Total number of points sampled
    p = 1;
    segment2_relax(1, :) = segment2_relax_t(p,:); 
    p = p + inc1; 
    for k = 2:div       
        if k <= div1
            segment2_relax(k,:) = segment2_relax_t(p,:);
            p = p + inc1;
        elseif k>div1 && k<=div1+div2
            segment2_relax(k,:) = segment2_relax_t(p,:);
            p = p + inc2;
        elseif k>div1+div2 && k<=div-1 
            segment2_relax(k,:) = segment2_relax_t(p,:);
            p = p + inc3;
        elseif k == div
            segment2_relax(k,:) = segment2_relax_t(end,:);
        end
    end
    
% %     Plot fit for validation
%     figure; hold on
%     plot(segment2_overlap1.Time, segment2_overlap1.Voltage, '.-');
%     plot(segment2_overlap2.Time, segment2_overlap2.Voltage, '.-');
%     plot(segment2_load.Time, segment2_load.Voltage, '.-');
%     plot(segment2_relax.Time, segment2_relax.Voltage, '.-');
%     hold off
    
    segment2 = [segment2_overlap1; segment2_load; segment2_relax; segment2_overlap2];
    pulse.segment2(i) = {segment2};
    
    Time = segment2.Time;
    for i = 1:length(Time)-1
        if Time(i+1) <= Time(i)
            break
        end
    end
    validateattributes(Time, {'double'},{'increasing'})
end