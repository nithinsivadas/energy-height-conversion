function [daten, f107] = get_solar_flux_from_penticton(inputDateNum, localFolder)
    % Get the penticton solar irradiation flux table, average values of
    % f107 per day (at noon local time or 20:00UT). 
    if nargin<1
        localFolder = strcat(initialize_root_path,'LargeFiles',filesep,'KP_AP',filesep);
    end
    
    if ~exist([localFolder,'f107.txt'],'file')
        fid = fopen([localFolder,'f107_2.txt'],'a');
        fwrite(fid,'0')
        fclose(fid);
    end
    
    fid1= fopen([localFolder,'f107_2.txt']);
    update = fscanf(fid1,'%f');
    fclose(fid1);
    readfile = [localFolder 'fluxtable.txt'];
    
    % update if new data available   
    if inputDateNum > update
    
        f = ftp('ftp.seismo.nrcan.gc.ca');
        solar_dir = '/spaceweather/solar_flux/daily_flux_values/';
        cd(f,solar_dir);
        latest = dir(f);
        last = latest(6).datenum; %checking last updated date of fluxtable.txt
        
        if ~exist(readfile,'file') || last>(update)
        
        disp('downloading solar data...') 
        mget(f,'fluxtable.txt',localFolder);
        %     filename = [localFolder 'listing_drao_noontime-flux-observed_daily_2.txt'];
        
        fid = fopen([localFolder,'f107_2.txt'],'w');
            fprintf(fid,'%f',last);
        fclose(fid);
        
        end
    
        close(f);
    end
    
    % read existing file, and estimate the necessary ouputs
    fid = fopen(readfile, 'r');
    if fid == 1
        error('Could not open file');
    end
    dat = textscan(fid, '%8d %6d %f %f %f %f %f','HeaderLines',2);
    fclose(fid);
    
    daten = datenum(num2str(dat{1}(dat{2}==200000)),'yyyymmdd');
    f107 = dat{5}(dat{2}==200000);
%     f107a = dat{6}(dat{2}==200000);
    f107(f107==0)=NaN; % NaN empty fields
%     f107(f107a==0)=NaN; % NaN empty fields
    
end