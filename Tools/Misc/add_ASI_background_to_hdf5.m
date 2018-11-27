function background = add_ASI_background_to_hdf5(siteName,outputH5FileStr)
%add_ASI_background_to_hdf5.m Adds background pixel intensity to hdf5 file
% Calculates the night-time background pixel intensity of the
% All-Sky-Imager referred to by 'siteName' in the HDF5
%file, and writes it back into the hdfile under the dataset 'background'
%
% Syntax: [background] = add_ASI_background_to_hdf5(siteName, outputH5FileStr)
%
% Inputs:
%    siteName - Name of All-sky-imager site, already in the HDF5 file, e.g. 'gako'
%    outputH5FileStr - H5 file in which you want to add background pixel intensity
%
% Outputs:
%    background - Background pixel intensity
%    This will be written on to the HDF5 file e.g., /GAKO/background
%
% Example:
%    background = add_ASI_background_to_hdf5('gako','test.h5')
%
% Other m-files required: none
% Subfunctions: ish5dataset
% MAT-files required: none
%
% See also:
% Author: Nithin Sivadas, Ph.D. Candidate,
% Center for Space Physics, Boston University
% email address: nithin@bu.edu
% Website:
% November 2018; Last revision: -
%------------- BEGIN CODE --------------

if strcmp(siteName,'pokerFlat')
    siteName = 'dasc';
end
groupName = ['/',upper(siteName),'/'];

[datasetExists, datasetInfo]=ish5dataset(outputH5FileStr,groupName);
if ~datasetExists
    warning(['Dataset ',groupName,' does not exist']);
    background = nan;
else
    tempData = (h5read(outputH5FileStr,[groupName,'ASI']));
    tempData(tempData>10000) = nan; % Removing all bright objects (like the moon)

    background = nanmedian(tempData(:)); % Taking the median of pixels in all time
                                         % and since most of the time, there is no aurora
                                         % you can expect the result to be the background
                                         % night-time pixel intensity.
                                         % (Or at least that is what this assumes.)

    [datasetExists, datasetInfo]=ish5dataset(outputH5FileStr,[groupName,'background']);
    if ~datasetExists
        h5create(outputH5FileStr,[groupName,'background'],1);
        h5writeatt(outputH5FileStr,[groupName,'background'],...
                'Dimensions','1 - estimated pixel intensity value of background');
        h5writeatt(outputH5FileStr,[groupName,'background'],...
                'Units','1 - estimated pixel intensity value of background');
        h5writeatt(outputH5FileStr,[groupName,'background'],...
                'updated_date',datestr(now));
    end
    h5write(outputH5FileStr,[groupName,'background'],background);
end

end
%------------- END CODE --------------
