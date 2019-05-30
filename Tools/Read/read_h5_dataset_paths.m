function datasetPaths = read_h5_dataset_paths (h5FilePath, groupPath)
%READ_H5_DATASET_PATHS Outputs string array of the paths of all datasets in HDF5 file
%
% SYNOPSIS: datasetPaths = read_h5_dataset_paths (h5FilePath, groupPath)
%
% INPUT h5FilePath - The hdf5 file path including file name
%		
%		groupPath  - The path of the group inside h5 file
%		
%		                   in which you'd like the function to search
%		
%		                   for datasets                                   
%
% OUTPUT datasetPaths - String array of the paths of datasets
%			
%			                       within the h5 file                
%
% REMARKS
%
% created with MATLAB ver.: 9.3.0.713579 (R2017b) on Microsoft Windows 10 Enterprise Version 10.0 (Build 17134)
%
% created by: Nithin Sivadas
% DATE: 23-Nov-2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin<2
        groupPath = [];
    end
    if isempty(groupPath)
        info = h5info(h5FilePath);
    else
        info = h5info(h5FilePath,groupPath);
    end
    nGroups = length(info.Groups);
    datasetPaths = [];
    if nGroups>0
        for thisGroup = 1:1:nGroups
            nGroupsTemp = length(info.Groups(thisGroup).Groups);
            if ~isempty(info.Groups(thisGroup).Datasets)
                datasetPaths = [datasetPaths;...
                        strcat([info.Groups(thisGroup).Name,'/'],...
                        string({info.Groups(thisGroup).Datasets.Name})')];
            end
            if nGroupsTemp >=1
                tempDatasetPaths =...
                    read_h5_dataset_paths(h5FilePath,info.Groups(thisGroup).Name);
                datasetPaths = [datasetPaths; tempDatasetPaths];
            end
        end
    else
        datasetPaths = [datasetPaths;...
                        strcat([info.Name,'/'],...
                        string({info.Datasets.Name})')];
        if isempty(datasetPaths)
            warning(['No groups or datasets in file: ',info.Filename, '|| Group: ',info.Name]);
        end
    end
end