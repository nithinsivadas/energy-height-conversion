function [data, info] = parse_superMag_data(fileStr)
%parse_superMag_data Parse the superMag data file that contains time, sml
%                    and smu, downloaded from superMag indices website.
% Input
%           filestr - SuperMag indices text file
%                     Ex: 'G:\My Drive\Research\Projects\Data\20201217-07-37-supermag.txt';
% Ouput
%           data.time - posix/unix time
%           data.SML  - Supermag electrojet index, lower (nT)
%           data.SMU  - Supermag electroject index, upper (nT)
%           info      - character vector with comments from superMag
%
% Created by: Nithin Sivadas
% Creation date: 17 Dec 2020


% Read the content of the file
content = fileread(fileStr);

% Skip the intial comments, and store it in info
temp = strfind(content, '>');
infoEnd = temp(end)+1;
info = content(1:infoEnd);

% Parse the data, starting after the initial commentary
cellarr=textscan(content(infoEnd:1:end),'%4u %2u %2u %2u %2u %2u %6.2f %6.2f');

% Store it in a data structure, convert time into posix/unix units
data.time = posixtime(datetime(cellarr{1},cellarr{2},cellarr{3},cellarr{4},cellarr{5},cellarr{6}));
data.SML = cellarr{7};
data.SMU = cellarr{8}; 

end

