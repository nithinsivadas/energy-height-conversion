function red_file_list = fitslistparce(file_list,timeextremes,varargin)
% fitslistparce.m
% red_file_list = fitslistparce(file_list,timeextreams)
% by John Swoboda
% This function will parce through a cellarray of strings that conatin
% names of FITS files that are time stamped according to University of
% Alaska's naming convention  and will find all of the files that are in a
% desired time period and if specified, the wavelength
%% Inputs
% file_list - A Nx1 or 1xN cell array with the file names.
% timeextreams - A 1x2 cell array containing the times in datestr format or
% a 1x2 numerical array holding datenum format numbers for the extremes of
% the window one wants to look at.
% wl - A scalar that is the wavelength in nm. (optional)
%% Outputs
% red_file_list - A reduced version of file_list with only the desired
% files.
%%
p = inputParser;
addOptional(p,'wl',[])
p.parse(varargin{:})
U = p.Results;
%%
if iscell(timeextremes)
    lowerextstr = timeextremes{1};
    upperextstr = timeextremes{2};
    lowe = datenum(lowerextstr);
    uppe = datenum(upperextstr);

elseif isnumeric(timeextremes)
    lowe = timeextremes(1);
    uppe = timeextremes(2);

else
    error('timeextreams need to be either a string or a datenum number');
end


% for wavelength
if isfinite(U.wl)
    str_2find =sprintf('%0.4d',varargin{1});
    keep_wl = false(size(file_list));

    for ifit = 1:length(file_list)
        fits_parts = regexp(file_list{ifit},'\_','split');
        keep_wl(ifit) = strcmp(fits_parts{3},str_2find);
    end
else
    keep_wl = true(size(file_list));
end

Numlist = fitsfiletimestamp(file_list);
keep_nums = (Numlist>=lowe & Numlist<=uppe) & keep_wl;
red_file_list = file_list(keep_nums);
