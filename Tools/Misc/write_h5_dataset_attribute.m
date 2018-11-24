function write_h5_dataset_attribute(h5FileStr,...
    datasetPath,Description,Dimensions,Units,preset)
%WRITE_H5_DATASET_ATTRIBUTE Add attributes to a dataset.
% Need to check if there already exists an attribute. 

if nargin < 6
    preset = 'none';
end

if strcmp(preset,'none')
    if isempty(Description)||isempty(Dimensions)||isempty(Units)
        error('Not allowed to skip any attributes without a valid preset');
    end
elseif strcmp(preset,'time')
    Description = 'Time';
    Dimensions = 'nTime x 1';
    Units = '[Posix time]';
else
    error('Invalid preset');
end

h5writeatt(h5FileStr,datasetPath,'Description',Description);
h5writeatt(h5FileStr,datasetPath,'Dimensions',Dimensions);
h5writeatt(h5FileStr,datasetPath,'Units',Units);
h5writeatt(h5FileStr,datasetPath,'Last_updated',datestr(now));
end

