function [ varargout ] = read_h5( filename )
%% read_h5.m
% By John Swoboda
% This function will read in data from the specifed h5 format for the
% GeoData class.
% Inputs
% filename - A string for the file name that holds the data.
% Output
% Varargout - A cell array with the GeoData values. See GeoData
% Documentation for details.
% {'data','coordnames','dataloc','sensorloc','times'};

varnames = {'data','coordnames','dataloc','sensorloc','times'};

fileinfo = hdf5info(filename);
varargout = cell(1,length(varnames));
for k = 1:length(varnames)
    ivar = varnames{k};
    isstructvar = false;
    groupnames = {fileinfo.GroupHierarchy.Groups(:).Name};
    % for the structs
    for l = 1:length(groupnames)
        igro = groupnames{l};
        if strcmp(ivar,igro(2:end))
            
            isstructvar=true;
            tempstruct = struct();
            setnames = {fileinfo.GroupHierarchy.Groups.Datasets(:).Name};
            for m=1:length(setnames)
                iset = setnames{m};
                isetname = iset(length(ivar)+3:end);
                temparr = h5read(filename,iset);
                tempstruct.(isetname) = permute(temparr,ndims(temparr):-1:1);
            end
            varargout{k} = tempstruct;
        end
    end 
    if isstructvar
        continue;
    end
    % For everything else.
    setnames = {fileinfo.GroupHierarchy.Datasets(:).Name};
    for l = 1:length(setnames)
        iset = setnames{l};
        if strcmp(ivar,iset(2:end))
            temparr = h5read(filename,iset);
            if ~ischar(temparr)
                temparr=permute(temparr,ndims(temparr):-1:1);
            end
            varargout{k} = temparr;
        end
    end
end

