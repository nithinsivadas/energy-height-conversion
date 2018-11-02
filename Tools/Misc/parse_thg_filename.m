function data=parse_thg_filename(fileName)
    % Parses the themis gbo filename to decipher information regarding the
    % site name, data type (asf,asc,ast,ask), level (l0,l1,l2),datetime
    
    split = strsplit(fileName,filesep);
    fileSplit = strsplit(split{end},'_');
    
    data.site = fileSplit{4};
    data.type = fileSplit{3};
    data.level = fileSplit{2};
    data.time = datenum(fileSplit{5},'yyyymmddHH');
end

