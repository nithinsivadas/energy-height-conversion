function thgdata = parse_thg_cdfData(dataFileStr,calFileStr)
% parse_thg_cdfData.m Extracts data from at thg_asf_ data file, with
% corresponding cal file
    fileinfo = parse_thg_filename(dataFileStr);
    dataInfo = spdfcdfinfo(dataFileStr);
    calInfo = spdfcdfinfo(calFileStr);
    thgdata.ASI = spdfcdfread(dataFileStr,'variables',dataInfo.Variables{1});
    thgdata.alt = spdfcdfread(calFileStr,'variables',calInfo.Variables{11});
    tempAz  = spdfcdfread(calFileStr,'variables',calInfo.Variables{5});
    thgdata.az = tempAz(:,:,1);
    tempEl  = spdfcdfread(calFileStr,'variables',calInfo.Variables{6});
    thgdata.el = tempEl(:,:,1);
    tempGlat    = cdfread(calFileStr,'variables',calInfo.Variables{8});
    thgdata.glat = tempGlat{1}(2:end,2:end,:);
    tempGlon = cdfread(calFileStr,'variables',calInfo.Variables{7});
    thgdata.glon = tempGlon{1}(2:end,2:end,:);
    tempMlat = cdfread(calFileStr,'variables',calInfo.Variables{10});
    thgdata.mlat = tempMlat{1}(2:end,2:end,:);
    tempMlon = cdfread(calFileStr,'variables',calInfo.Variables{9});
    thgdata.mlon = tempMlon{1}(2:end,2:end,:);
    thgdata.time = spdfcdfread(dataFileStr,'variables',dataInfo.Variables{4}); %Unix time   
end

