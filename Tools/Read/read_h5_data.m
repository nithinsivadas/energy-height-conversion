function dataset = read_h5_data (h5FilePath, groupPath)
%READ_H5_DATA outputs a table containing hdf5 file data within specified group 
%
% SYNOPSIS: data = read_h5_data(h5FilePath, groupPath)
%
% INPUT h5FilePath - The hdf5 file path including file name
%		
%		groupPath  - The path of the group inside h5 file
%		             default: '/'
%
% OUTPUT data     - Table containing hdf5 file data in its format
%
% REMARKS
%
% created with MATLAB ver.: 9.3.0.713579 (R2017b) on Microsoft Windows 10 Enterprise Version 10.0 (Build 17134)
%
% created by: Nithin Sivadas
% DATE: 23-Nov-2018
% UPDATE: 25-Nov-2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin<2
        groupPath = '/';
    end
    
    datasetPaths = read_h5_dataset_paths(h5FilePath,groupPath);
    nDatasets = length(datasetPaths);
    
    for thisDataset = 1:1:nDatasets
        info = h5info(h5FilePath, char(datasetPaths(thisDataset)));
        groupNames = strsplit(datasetPaths(thisDataset),'/');
        nGroups = length(groupNames);
        dataset{thisDataset,1} = datasetPaths(thisDataset); % Path
        dataset{thisDataset,2} = groupNames(nGroups); % Name
        dataset{thisDataset,3} = h5read(h5FilePath, char(datasetPaths(thisDataset))); % Data
        dataset{thisDataset,4} = info.Dataspace.Size; % Size
        dataset{thisDataset,5} = info.Dataspace.MaxSize; % maxSize
        dataset{thisDataset,6} = info.ChunkSize; % chunkSize
        if ~isempty(info.Attributes)
            dataset{thisDataset,7} = string({info.Attributes.Name}); % Attribute Names
            dataset{thisDataset,8} = string({info.Attributes.Value}); % Attribute Values
        else
            dataset{thisDataset,7} = " "; % Attribute Names
            dataset{thisDataset,8} = " "; % Attribute Values
        end
        % Permuting H5 data, to resolve the matlab/hdf5 conversion problem
        permuteVar = fliplr(1:1:length(size(dataset{thisDataset,3})));
        if ~isempty(permuteVar)
            dataset{thisDataset,3} = ipermute(dataset{thisDataset,3},permuteVar);
        end
        dataset{thisDataset,4} = fliplr(dataset{thisDataset,4});
        dataset{thisDataset,5} = fliplr(dataset{thisDataset,5});
        dataset{thisDataset,6} = fliplr(dataset{thisDataset,6});
        
    end
dataset = array2table(dataset,'VariableNames',...
    {'Path','Name','Data','Size','MaxSize','ChunkSize','Attribute_Names','Attribute_Values'});
end
