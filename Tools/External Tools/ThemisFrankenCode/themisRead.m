%% Read raw themis data
% created by J. Brent Parham 3/19/2018
% takes in date time object and prints out data for plotting
function camera=themisRead(timeImage)

%% FSIM
system(['/usr/local/bin/convert  ','asi_fulldata/fsim_themis08/ut',datestr(timeImage,'HH'),'/20180114_',datestr(timeImage,'HHMM'),'_fsim_themis08_full.pgm',' -compress none outfile.pgm']);
system('split -p P2 outfile.pgm');

display(timeImage)

index=1;
alt=[90;110;150];
Re=6378.1;
el_lim=12;
h=alt(index);
p_slant=(pi - el_lim*pi/180 -pi/2 - asin(Re*sind(90+el_lim)/(Re+h)))*Re;

caldata = cdfread('thg_l2_asc_fsim_19700101_v01.cdf','variables','thg_asf_fsim_glon');
lon=double(caldata{1,1}(2:end,2:end,index)-360);
caldata = cdfread('thg_l2_asc_fsim_19700101_v01.cdf','variables','thg_asf_fsim_glat');
lat=double(caldata{1,1}(2:end,2:end,index));

listIm=dir('x*');

% restrict the image data to quality center
cam_cent=[61.762,238.779-360];

index2=min(round(second(timeImage)/60*20)+1,length(listIm));
themImage=rot90((double(imread(listIm(index2).name))),2);
alph = sqrt((lat-cam_cent(1)).^2+(cosd(lat).*(lon-cam_cent(2))).^2);
themImage(alph>(p_slant./6378.1*180/pi))=nan;

% rescale for better plotting 
themImage=(themImage-nanmedian(themImage(:)))/(max(themImage(:))-nanmedian(themImage(:)));

camera(1).themImage=themImage;
camera(1).lat=lat;
camera(1).lon=lon;

system('rm x* outfile.pgm');

%% FSMI
system(['/usr/local/bin/convert  ','asi_fulldata/fsmi_themis10/ut',datestr(timeImage,'HH'),'/20180114_',datestr(timeImage,'HHMM'),'_fsmi_themis10_full.pgm',' -compress none outfile.pgm']);
system('split -p P2 outfile.pgm');


caldata = cdfread('thg_l2_asc_fsmi_19700101_v01.cdf','variables','thg_asf_fsmi_glon');
lon=double(caldata{1,1}(2:end,2:end,index)-360);
caldata = cdfread('thg_l2_asc_fsmi_19700101_v01.cdf','variables','thg_asf_fsmi_glat');
lat=double(caldata{1,1}(2:end,2:end,index));

listIm=dir('x*');

% restrict the image data to quality center
cam_cent=[59.984,248.158-360];

index2=min(round(second(timeImage)/60*20)+1,length(listIm));
themImage=rot90((double(imread(listIm(index2).name))),2);
alph = sqrt((lat-cam_cent(1)).^2+(cosd(lat).*(lon-cam_cent(2))).^2);
themImage(alph>(p_slant./6378.1*180/pi))=nan;

% rescale for better plotting 
themImage=(themImage-nanmedian(themImage(:)))/(max(themImage(:))-nanmedian(themImage(:)));

camera(2).themImage=themImage;
camera(2).lat=lat;
camera(2).lon=lon;

system('rm x* outfile.pgm');

% figure(9); clf;
%     worldmap([40 75], [-110 -75])
%     geoshow('landareas.shp','FaceColor',[.9 .9 .9])
%     hold on;
%     surfm(lat,lon,themImage,'edgecolor','none','facealpha',0.5); view([0 0 1]); axis tight
%     colormap bone
%     caxis([nanmedian(themImage(:)) 0.2*max(themImage(:))]); 
%     title(datestr(todatenum(imagedata{window,6})))
% 
%     [hcirc,circlelat,circlelon]=circlem(cam_cent(1),cam_cent(2),p_slant);
%     set(hcirc,'edgecolor','w')
% figure(2); clf
%     surf(double(imagedata{window(jj),1}),'edgecolor','none'); view([0 0 1]); axis tight
%     colormap bone
%     caxis([double(median(imagedata{window(jj),1}(:))) 0.2*double(max(imagedata{window(jj),1}(:)))])
%     title(datestr(todatenum(imagedata{window(jj),6})))