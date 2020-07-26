omniH5FileStr = 'C:\Users\nithin\Documents\GitHub\LargeFiles\omni\omni.h5';
outputH5FileStr ='G:\My Drive\Research\Projects\Referee\Shen et al 2020\snkq_sheng_2020.h5'; 
minTimeStr = '2014-08-27 06:00:00';
maxTimeStr = '2014-08-27 06:30:00';
dateFormat = 'yyyy-mm-dd HH:MM:SS';

%%
create_thg_hdf5('snkq',outputH5FileStr,'minTimeStr',minTimeStr,'maxTimeStr',maxTimeStr,...
'localStorePath','G:\My Drive\Research\Projects\Referee\Shen et al 2020\');

%%
T=read_h5_data('G:\My Drive\Research\Projects\Referee\Shen et al 2020\snkq_sheng_2020.h5');
%% Fix Download themis ftp to web :/
add_thm_hdf5('tha',outputH5FileStr,omniH5FileStr,...
        'minTimeStr',minTimeStr,...
        'maxTimeStr',maxTimeStr,...
        'dateFormat',dateFormat,...
        'magneticFieldModel','TS96');

%%
figure; 
iDateStr='2014-08-27 06:20:30'; 
i=find_time(unixtime2matlab(T.Data{11}),iDateStr);  
plot_DASC_geodetic(squeeze(T.Data{1}(i,:,:)),unixtime2matlab(T.Data{11}(i)),T.Data{6},T.Data{7},256,[50,58],[-95,-65]); 
colorbar; 
set(gca,'ColorScale','log'); 
title(datestr(unixtime2matlab(T.Data{11}(i)))); 
caxis([3*10^3 10^4]); colormap(gray);