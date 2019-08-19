function [status] = download_DASC_FITS_for_storm(urls,wavelengthStr,stormTime,storeDir,h5FileStr)
    wavelengthStr = string(wavelengthStr);
    wavelengths = unique(wavelengthStr);
    for w = 1:1:length(wavelengths)
       tempStorePath = strcat(storeDir,'temp',filesep,wavelengths(w));
       folderStr = datestr(stormTime,'yyyymmdd');
       if ~isdir(tempStorePath)
           mkdir(tempStorePath);
       end
       urlFile = 'tempURL.txt';
       urlFilePath = strcat(tempStorePath,filesep,urlFile);
       fileID = fopen(urlFilePath,'w'); fprintf(fileID,'%s\r\n',urls');fclose(fileID);
       if isunix
       [status,cmdout]=unix(strcat('aria2c -V -c -j 50 ',' -d "',tempStorePath,'" -i "',urlFilePath,'"'));
       else
       [status,cmdout]=system(strcat('aria2c -V -c -j 50 ',' -d "',tempStorePath,'" -i "',urlFilePath,'"'));
       end
       read_all_local_FITS_and_create_HDF5(h5FileStr, wavelengths(w), folderStr, tempStorePath);
       % Need remove all the temp files
       rmdir(tempStorePath,'s');
       disp(w)
    end
end
