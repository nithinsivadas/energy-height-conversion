%Okay, Matlab's h5write(filename, dataset, data) function doesn't work for
%strings. It hasn't worked with strings for years. The forum post that
%comes up first in Google about it is from 2009. Yeah. This is terrible,
%and evidently it's not getting fixed. So, low level functions. Fun fun.
%
%What I've done here is adapt examples, one from the hdf group's website
%https://support.hdfgroup.org/HDF5/examples/api18-m.html called
%"Read / Write String Datatype (Dataset)", the other by Jason Kaeding.
%
%I added functionality to check whether the file exists and either create
%it anew or open it accordingly. I wanted to be able to likewise check the
%existence of a dataset, but it looks like this functionality doesn't exist
%in the API, so I'm doing a try-catch to achieve the same end. Note that it
%appears you can't just create a dataset or group deep in a heirarchy: You
%have to create each level. Since I wanted to accept dataset names in the
%same format as h5read(), in the event the dataset doesn't exist, I loop
%over the parts of the dataset's path and try to create all levels. If they
%already exist, then this action throws errors too; hence a second
%try-catch.
%
%I've made it more advanced than h5create()/h5write() in that it all
%happens in one call and can accept data inputs of variable size. I take
%care of updating the dataset's extent to accomodate changing data array
%sizes. This is important for applications like adding a new timestamp
%every time the file is modified.
%
%@author Pavel Komarov pavel@gatech.edu 941-545-7573
function h5createwritestr(filename, dataset, str)

    %"The class of input data must be cellstring instead of char when the
    %HDF5 class is VARIABLE LENGTH H5T_STRING.", but also I don't want to
    %force the user to put braces around single strings, so this.
    if ischar(str)
        str = {str};
    end

    %check whether the specified .h5 exists and either create or open
    %accordingly
    if ~exist(filename, 'file')
        file = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
    else
        file = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
    end

    %set variable length string type
    vlstr_type = H5T.copy('H5T_C_S1');
    H5T.set_size(vlstr_type,'H5T_VARIABLE');

    % There is no way to check whether a dataset exists, so just try to
    % open it, and if that fails, create it.
    try
        dset = H5D.open(file, dataset);
        H5D.set_extent(dset, fliplr(size(str)));
    catch
        %create the intermediate groups one at a time because evidently the
        %API's functions aren't smart enough to be able to do this themselves.
        slashes = strfind(dataset, '/');
        for i = 2:length(slashes)
            url = dataset(1:(slashes(i)-1));%pull out the url of the next level
            try
                H5G.create(file, url, 1024);%1024 "specifies the number of
            catch   %bytes to reserve for the names that will appear in the group"
            end
        end

        %create a dataspace for cellstr
        H5S_UNLIMITED = H5ML.get_constant_value('H5S_UNLIMITED');
        spacerank = max(1, sum(size(str) > 1));

        if ndims(str) == 2 && min(size(str)) == 1
            dspace = H5S.create_simple(spacerank, length(str), ones(1, spacerank)*H5S_UNLIMITED);        
        else
            dspace = H5S.create_simple(spacerank, fliplr(size(str)), ones(1, spacerank)*H5S_UNLIMITED);
        end
        %create a dataset plist for chunking. (A dataset can't be unlimited
        %unless the chunk size is defined.)
        plist = H5P.create('H5P_DATASET_CREATE');
        chunksize = ones(1, spacerank);
        chunksize(1) = 2;
        H5P.set_chunk(plist, chunksize);% 2 strings per chunk
        dset = H5D.create(file, dataset, vlstr_type, dspace, plist);

        %close things
        H5P.close(plist);
        H5S.close(dspace);
    end
    
   
    %write data
    H5D.write(dset, vlstr_type, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', str);

    %close file & resources
    H5T.close(vlstr_type);
    H5D.close(dset);
    H5F.close(file);
end