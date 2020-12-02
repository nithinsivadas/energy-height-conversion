function [Numlist, Wavelength] = fitsfiletimestamp(file_list)
% fitsfiletimestamp.m
% Numlist = fitsfiletimestamp(file_list)
% by John Swoboda
% This function will give you an array of time stamps (same size and shape 
% as file_list) in the MATLAB datenum format for a cell array of FITS files
% from the Poker Flat allsky camera.  It is assumed that file_list is a
% cell array of strings.
%% Inputs
% file_list - A Nx1 or 1xN cell array with the file names.
%% Outputs
% Numlist - A Nx1 or 1xN array of datenum time stamps.
% Wavelength - A Nx1 or 1xN array of wavelength in [nm] that the image
% corresponds to. 
%% Main
Wavelength = string();
if ischar(file_list) || isstring(file_list)
    temp_list = replace(file_list,'__','_');
    N  = size(temp_list,1);
    file_list = {};
    
    i = 1;
    for k = 1:N
        if length(deblank(char(temp_list(k,:))))<=38 && length(deblank(char(temp_list(k,:))))>=36
            file_list{i} = temp_list(k,:);
            i=i+1;
        end
    end
end
N = length(file_list);
Numlist = zeros(size(file_list));

makewarn = false;

% if N<k
%     warning([num2str(k-N),' of the FITS file names have a different format.',...
%         ' They have string lengths > 38']);
% end

if N==0
    error('No FITS files have the right file name length = 38 or 37');
end

for k = 1:N
    temp_str = file_list{k};
    split_cell = regexp(temp_str,'\.','split');
    
    time_split = regexp(split_cell{1},'\_','split');
    if length(time_split{end-1})==6
        yearstr = ['20',time_split{end-1}];
%         makewarn =true;
    elseif length(time_split{end-1})==8
        yearstr = time_split{end-1};
    end
        
    time_str = [yearstr,'-',time_split{end},'.',split_cell{2}];
    
    Numlist(k) = datenum(time_str,'yyyymmdd-HHMMSS.FFF');
    Wavelength(1,k) = string(time_split{3});
end

if makewarn
    warning('Some of the strings year strings were only 2 digits');
end