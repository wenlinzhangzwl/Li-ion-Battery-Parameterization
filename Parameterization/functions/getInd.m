function [ind] = getInd(Data, StepNum)
% Returns a table ('ind') that contains the indices of points of each step number ('indices' column) & 
% the start & end points of all the segments of a step number ('segment' column)
% ind: m x 3 table where m is the number of step numbers searched 
% Data: Table containing raw data (meas_t)
% StepNum: Vector containing step numbers of interest


% Find indices of each step in meas stored as table 'ind': 
ind = table('Size', [length(StepNum), 3],...
            'VariableTypes', ["double","cell", "cell"],...
            'VariableNames',["stepNum", "indices", "segment"]);
for i = 1:length(StepNum)
    ind.stepNum(i) = StepNum(i);
    ind.indices(i) = {find(Data.Step == StepNum(i))};
    indices = ind.indices{i};       % Indices of data points with this step number
    segments = [indices(1) 0];      % Start & end points of each segment with this start number
    p = 1;                          % Counter
    
    for j = 1:length(indices)
        if j == length(indices)
            segments(p,2) = indices(j);
        elseif indices(j+1) ~= indices(j)+1
            segments(p,2) = indices(j);
            segments(p+1,1) = indices(j+1);
            p = p + 1;
        end
    end
    ind.segment(i) = {segments}; 
end