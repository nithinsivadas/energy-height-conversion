function [data,h5DataLoc]=h5read_data(h5FilePath)
% Reads all h5 data from an h5 file into a structure data.
% no longer needs to be in use. See read_h5_data.m
info=h5info(h5FilePath);
nGroups=length(info.Groups);
for thisGroup=1:1:nGroups
    for thisDataSet = 1:1:length(info.Groups(thisGroup).Datasets)
    
    thisGroupName=info.Groups(thisGroup).Name(2:end);
    thisDataSetName=info.Groups(thisGroup).Datasets(thisDataSet).Name;
    
    h5DataLoc{thisGroup,thisDataSet} = ...
        strcat('/',thisGroupName,'/',thisDataSetName);
    
    temp=h5read(h5FilePath,h5DataLoc{thisGroup,thisDataSet});
    if isstruct(temp)
        data.(thisGroupName).(thisDataSetName) = ...
            structfun(@transpose,...
            temp,...
            'UniformOutput',false);
    else
        data.(thisGroupName).(thisDataSetName) = ...
            transpose(temp);
    end
end
end