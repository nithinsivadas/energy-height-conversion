function write_sc_to_hdf5(h5OutputFile,probeName,time,XYZ_GEO,magFieldStr,NFoot,Lm)
    % write_sc_to_hdf5.m Writes spacecraft state data to hdf5, as well
    % footpoint and Lm data. 
    
    NFoot = NFoot(:,[2 3 1]); % Converting GDZ to lat,lon,alt format
    time  = posixtime(datetime(time,'ConvertFrom','datenum'));
    groupName = ['/SC/',upper(probeName),'/'];
    
    % Creating dataset if previously non-existent
    [datasetExists, ~]=ish5dataset(h5OutputFile,groupName);
    if ~datasetExists
        h5create(h5OutputFile,[groupName,'time'],...
                [1 Inf],'ChunkSize',[1 10],'Deflate',9);
        h5writeatt(h5OutputFile,[groupName,'time'],...
            'Dimensions','nTime x 1');
        h5writeatt(h5OutputFile,[groupName,'time'],...
            'Units','Posix time');
        
        h5create(h5OutputFile,[groupName,'XYZ_GEO'],...
                [3 Inf],'ChunkSize',[3 10],'Deflate',9);
        h5writeatt(h5OutputFile,[groupName,'XYZ_GEO'],...
            'Dimensions','nTime x 3');
        h5writeatt(h5OutputFile,[groupName,'time'],...
            'Units','RE (XYZ_GEO)');
    end
    magGroupName = [groupName,upper(magFieldStr),'/'];
    [datasetExists, ~]=ish5dataset(h5OutputFile, magGroupName);
    if ~datasetExists
        h5create(h5OutputFile,[magGroupName,'NFoot'],...
                [3 Inf],'ChunkSize',[3 10],'Deflate',9);
        h5writeatt(h5OutputFile,[magGroupName,'NFoot'],...
            'Dimensions','nTime x 3');
        h5writeatt(h5OutputFile,[magGroupName,'NFoot'],...
            'Units','RE (Lat,Lon,Alt_GDZ)');
        
        h5create(h5OutputFile,[magGroupName,'Lm'],...
                [1 Inf],'ChunkSize',[1 10],'Deflate',9);
        h5writeatt(h5OutputFile,[magGroupName,'Lm'],...
            'Dimensions','nTime x 1');
        h5writeatt(h5OutputFile,[magGroupName,'Lm'],...
            'Units','RE');
    end
    
    % Write time
    h5write(h5OutputFile,[groupName,'time'],(time'),[1 1],size(time'));
    h5write(h5OutputFile,[groupName,'XYZ_GEO'],(XYZ_GEO)',[1 1], size(XYZ_GEO'));
    h5write(h5OutputFile,[magGroupName,'NFoot'],(NFoot)',[1 1], size((NFoot)'));
    h5write(h5OutputFile,[magGroupName,'Lm'],(Lm)',[1 1], size((Lm)'));
        
    if ~datasetExists
    h5writeatt(h5OutputFile,magGroupName,...
        'creation_date',datestr(now));
    else
    h5writeatt(h5OutputFile,magGroupName,...
        'updated_date',datestr(now));
    end

    
end