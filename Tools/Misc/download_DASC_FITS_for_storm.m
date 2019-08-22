function [status] = download_DASC_FITS_for_storm(urls,wavelengthStr,stormID,storeDir)
    wavelengthStr = string(wavelengthStr);
    wavelengths = unique(wavelengthStr);
    for w = 1:1:length(wavelengths)
       tempStorePath = strcat(storeDir,'temp',filesep,num2str(stormID),filesep,wavelengths(w));
       if ~isfolder(tempStorePath)
           mkdir(tempStorePath);
       end
       urlFile = 'tempURL.txt';
       urlFilePath = strcat(tempStorePath,filesep,urlFile);
       urlsThisW = urls(strcmp(wavelengthStr,wavelengths(w)));
       fileID = fopen(urlFilePath,'w'); fprintf(fileID,'%s\r\n',urlsThisW');fclose(fileID);
       if isunix
       [status,cmdout]=system(strcat('aria2c -V -c -j 50 ',' -d "',tempStorePath,'" -i "',urlFilePath,'"'));
       else
       [status,cmdout]=system(strcat('aria2c -V -c -j 50 ',' -d "',tempStorePath,'" -i "',urlFilePath,'"'));
       end
%        tempStr = strsplit(h5FileStr,filesep);
%        tempStr{end} = strcat(num2str(stormID),'_',stormDateStr,'_',tempStr{end});
%        h5FileStr = strjoin(tempStr,filesep);
%        read_all_local_FITS_and_create_HDF5(h5FileStr, wavelengths(w), tempStorePath);
%        % Remove All Temp Files
%        try
%         rmdir(tempStorePath,'s');
%         catch ME
%        end
%        if ~isempty(dir(strcat(tempStorePath,filesep,'*.FITS')))
%             error('The folder or its contents could not be erased');
%        end
    end
end

