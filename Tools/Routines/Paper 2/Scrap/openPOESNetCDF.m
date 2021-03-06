%% Test: opening netcdf
poesCDFPath = 'G:\My Drive\Research\Projects\Paper 2\Data\NOAA17';
poesnCDFName = 'POES_combinedSpectrum_n17_90_20080326.nc';
source = [poesCDFPath,filesep,poesnCDFName];
ncid = netcdf.open(source);
finfo = ncinfo(source);

%% Test: Opening the text file
poesnTxtName = 'n17_20080326.txt';
sourcetxt = [poesCDFPath,filesep,poesnTxtName];
noaa17 = importfile(sourcetxt);


function n1 = importfile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   N1 = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   N1 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   n1 = importfile('n17_20080326.txt', 2, 43201);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/04/19 16:15:07

%% Initialize variables.
delimiter = ' ';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
n1 = table;
n1.year = cell2mat(raw(:, 1));
n1.doy = cell2mat(raw(:, 2));
n1.hr = cell2mat(raw(:, 3));
n1.min = cell2mat(raw(:, 4));
n1.sec = cell2mat(raw(:, 5));
n1.lat = cell2mat(raw(:, 6));
n1.lon = cell2mat(raw(:, 7));
n1.altkm = cell2mat(raw(:, 8));
n1.orbit = cell2mat(raw(:, 9));
n1.mep0P1 = cell2mat(raw(:, 10));
n1.mep0P2 = cell2mat(raw(:, 11));
n1.mep0P3 = cell2mat(raw(:, 12));
n1.mep0P4 = cell2mat(raw(:, 13));
n1.mep0P5 = cell2mat(raw(:, 14));
n1.mep0P6 = cell2mat(raw(:, 15));
n1.mep0E1 = cell2mat(raw(:, 16));
n1.mep0E2 = cell2mat(raw(:, 17));
n1.mep0E3 = cell2mat(raw(:, 18));
n1.mep90P1 = cell2mat(raw(:, 19));
n1.mep90P2 = cell2mat(raw(:, 20));
n1.mep90P3 = cell2mat(raw(:, 21));
n1.mep90P4 = cell2mat(raw(:, 22));
n1.mep90P5 = cell2mat(raw(:, 23));
n1.mep90P6 = cell2mat(raw(:, 24));
n1.mep90E1 = cell2mat(raw(:, 25));
n1.mep90E2 = cell2mat(raw(:, 26));
n1.mep90E3 = cell2mat(raw(:, 27));
n1.mepOmniP6 = cell2mat(raw(:, 28));
n1.mepOmniP7 = cell2mat(raw(:, 29));
n1.mepOmniP8 = cell2mat(raw(:, 30));
n1.mepOmniP9 = cell2mat(raw(:, 31));
n1.ted01 = cell2mat(raw(:, 32));
n1.ted02 = cell2mat(raw(:, 33));
n1.ted03 = cell2mat(raw(:, 34));
n1.ted04 = cell2mat(raw(:, 35));
n1.ted05 = cell2mat(raw(:, 36));
n1.ted06 = cell2mat(raw(:, 37));
n1.ted07 = cell2mat(raw(:, 38));
n1.ted08 = cell2mat(raw(:, 39));
n1.ted301 = cell2mat(raw(:, 40));
n1.ted302 = cell2mat(raw(:, 41));
n1.ted303 = cell2mat(raw(:, 42));
n1.ted304 = cell2mat(raw(:, 43));
n1.ted305 = cell2mat(raw(:, 44));
n1.ted306 = cell2mat(raw(:, 45));
n1.ted307 = cell2mat(raw(:, 46));
n1.ted308 = cell2mat(raw(:, 47));
n1.ted0s1 = cell2mat(raw(:, 48));
n1.ted0s2 = cell2mat(raw(:, 49));
n1.ted0s3 = cell2mat(raw(:, 50));
n1.ted0s4 = cell2mat(raw(:, 51));
n1.ted0s5 = cell2mat(raw(:, 52));
n1.ted0s6 = cell2mat(raw(:, 53));
n1.ted0s7 = cell2mat(raw(:, 54));
n1.ted0s8 = cell2mat(raw(:, 55));
n1.ted30s1 = cell2mat(raw(:, 56));
n1.ted30s2 = cell2mat(raw(:, 57));
n1.ted30s3 = cell2mat(raw(:, 58));
n1.ted30s4 = cell2mat(raw(:, 59));
n1.ted30s5 = cell2mat(raw(:, 60));
n1.ted30s6 = cell2mat(raw(:, 61));
n1.ted30s7 = cell2mat(raw(:, 62));
n1.ted30s8 = cell2mat(raw(:, 63));
n1.tedback11 = cell2mat(raw(:, 64));
n1.tedback12 = cell2mat(raw(:, 65));
n1.tedback13 = cell2mat(raw(:, 66));
n1.tedback14 = cell2mat(raw(:, 67));
n1.tedback21 = cell2mat(raw(:, 68));
n1.tedback22 = cell2mat(raw(:, 69));
n1.tedback23 = cell2mat(raw(:, 70));
n1.tedback24 = cell2mat(raw(:, 71));
n1.tedfx1 = cell2mat(raw(:, 72));
n1.tedfx2 = cell2mat(raw(:, 73));
n1.tedfx3 = cell2mat(raw(:, 74));
n1.tedfx4 = cell2mat(raw(:, 75));
n1.tedfx5 = cell2mat(raw(:, 76));
n1.tedfx6 = cell2mat(raw(:, 77));
n1.tedfx7 = cell2mat(raw(:, 78));
end

