% Generates MATLAB data from Digatron .csv data

clear

folder = 'C:\Users\Wenlin\OneDrive\SCHOOL\48V Battery Pack Design\Test Results\3 - Characterization test 2';
cd(folder);
files = dir('*.csv');
addpath(folder);    

for q = 1:length(files)
    filename = files(q).name;
    meas = importfile(filename,30);
    filename2=[datestr(meas.TimeStamp(1),'mm-dd-yy_HH.MM '),filename(1:(end-4)),'.mat'];
    save(filename2,'meas')
end

clear ALL
function [meas] = importfile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [BATTERYNAME,SBLIMOTIVE52AH_DELPHI_SN090,VARNAME3,VARNAME4,VARNAME5,VARNAME6,VARNAME7,VARNAME8,VARNAME9,VARNAME10,VARNAME11,VARNAME12,VARNAME13,VARNAME14,VARNAME15,VARNAME16]
%   = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [BATTERYNAME,SBLIMOTIVE52AH_DELPHI_SN090,VARNAME3,VARNAME4,VARNAME5,VARNAME6,VARNAME7,VARNAME8,VARNAME9,VARNAME10,VARNAME11,VARNAME12,VARNAME13,VARNAME14,VARNAME15,VARNAME16]
%   = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [BatteryName,SBLimotive52Ah_Delphi_SN090,VarName3,VarName4,VarName5,VarName6,VarName7,VarName8,VarName9,VarName10,VarName11,VarName12,VarName13,VarName14,VarName15,VarName16] = importfile('10degC_SBLimotive.csv',28, 34);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2017/12/04 20:11:31

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 28;
    endRow = inf;
end

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,6,7,9,10,11,12,13,14,15,16]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch 
        end
    end
end

% Convert the contents of columns with dates to MATLAB datetimes using date
% format string.
try
    dates{4} = datenum(dataArray{4}, 'HH:MM:SS.FFF');
catch
    try
        % Handle dates surrounded by quotes
        dataArray{4} = cellfun(@(x) x(2:end-1), dataArray{4}, 'UniformOutput', false);
        dates{4} = datenum(dataArray{4}, 'HH:MM:SS.FFF');
    catch
        dates{4} = repmat(datetime([NaN NaN NaN]), size(dataArray{4}));
    end
end

%anyBlankDates = cellfun(@isempty, dataArray{4});
%anyInvalidDates = isnan(dates{4}.Hour) - anyBlankDates;
dates = dates(:,4);

%% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [2,6,7,9,10,11,12,13,14,15,16]);
rawCellColumns = raw(:, [1,3,5,8]);


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Allocate imported array to column variable names
 VarName1 = rawCellColumns(:, 1);
 VarName2 = cell2mat(rawNumericColumns(:, 1));
 VarName3 = rawCellColumns(:, 2);
 VarName4 = dates{:, 1};    
 VarName5 = rawCellColumns(:, 3);
 VarName6 = cell2mat(rawNumericColumns(:, 2));
 VarName7 = cell2mat(rawNumericColumns(:, 3));
 VarName8 = rawCellColumns(:, 4);
 VarName9 = cell2mat(rawNumericColumns(:, 4));
 VarName10 = cell2mat(rawNumericColumns(:, 5));
 VarName11 = cell2mat(rawNumericColumns(:, 6));
 VarName12 = cell2mat(rawNumericColumns(:, 7));
 VarName13 = cell2mat(rawNumericColumns(:, 8));
 VarName14 = cell2mat(rawNumericColumns(:, 9));
 VarName15 = cell2mat(rawNumericColumns(:, 10));
 VarName16 = cell2mat(rawNumericColumns(:, 11));

TimeStamp = VarName1;
Step = VarName2;
Status = VarName3;
ProgTime = VarName4;
StepTime = VarName5;
Cycle = VarName6;
CycleLevel = VarName7;
Procedure = VarName8;
Voltage = VarName9;
Current = VarName10;
Batt_Temp = VarName11;
AhAccu = VarName12;
Energy = VarName14;

%% Assign variables to measured data array and create time and power variables

meas.Time=(ProgTime-ProgTime(1))*24*60*60; %Convert to seconds
% for i= 1:length(ProgTime)
%     integ=floor(meas.Time(i));
%     fract=meas.Time(i)-integ;
%     meas.TimeStamp{i}=datestr(datenum(TimeStamp{i})+fract/86400,'dd/mm/yyyy HH:MM:SS.FFF PM');%add progtime to date
% end
meas.TimeStamp=TimeStamp;
meas.Step = Step;
meas.StepTime = StepTime;
meas.Procedure = Procedure;
meas.Voltage=Voltage;
meas.Current=Current;
meas.Ah=AhAccu;
meas.Wh=Energy;
meas.Power=meas.Voltage.*meas.Current;
meas.Battery_Temp_degC=Batt_Temp;
%meas.Chamber_Temp_degC=ChamberT;

beep
pause(0.3)
beep
pause(0.3)
beep

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% VarName4=datenum(VarName4);

end