function background = add_ASI_background_to_hdf5(siteName,outputH5FileStr)
%add_ASI_background_to_hdf5.m Calculates the night-time background pixel
%intensity of the All-Sky-Imager referred to by 'siteName' in the HDF5
%file, and writes it into the appropriate database

if strcmp(siteName,'pokerFlat')
    siteName = 'dasc';
end

groupName = ['/',upper(siteName),'/'];
[datasetExists, datasetInfo]=ish5dataset(outputH5FileStr,groupName);
if ~datasetExists
    error(['Dataset ',groupName,' does not exist']);
end
tempData = (h5read(outputH5FileStr,[groupName,'ASI']));
tempData(tempData>10000) = nan; % Removing all bright objects
background = nanmedian(tempData(:));

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

