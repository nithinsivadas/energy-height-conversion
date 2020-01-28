%% plot poes data distribution
load('G:\My Drive\Nithin\Job Search\JEPostDoc\Data\poes15_orbit.mat');
%%
time = datenum(orbit.time);
[orbit.MLAT, orbit.MLON, orbit.MLT] =...
    get_magnetic_coordinates([orbit.alt,orbit.lat,convert_longitude(orbit.lon,'360to180')],time);

%%
% figure; 
% h = histogram2Polar(orbit.lon,90-abs(orbit.lat),5,'RTicks',90-[90,60,30,0]);
% h.RAxisLabels(1).String='90';
% h.RAxisLabels(2).String='60';
% h.RAxisLabels(3).String='30';
% h.RAxisLabels(4).String='0';
% %,'RTicks',90-[90,80,70,60,50]
% set(gca,'ColorScale','log');
% colormap(inferno);

%%
storeImageDir = 'G:\My Drive\Nithin\Job Search\JEPostDoc\Data\';
imageName = 'plot2';

figureHandle=figure;
resize_figure(figureHandle,150,150);p = panel();
p.pack(1,1);
p.select('all');
text(0.5,1.15,{'2004-2005 NOAA-15','Total Samples'},...
    'FontWeight','bold','Units','normalized','HorizontalAlignment','center');
p.margintop=15;
p.marginright=15;
set(gca,'XColor', 'none','YColor','none');
index = orbit.MLAT>0;
h = histogram2Polar(orbit.MLT(index)*360./24,90-(orbit.MLAT(index)),1,'RTicks',90-[90,60,30,0]);
h.RAxisLabels(1).String='90';
h.RAxisLabels(2).String='60';
h.RAxisLabels(3).String='30';
h.RAxisLabels(4).String='0';
h.ThetaAxisLabels(1).String='00';
h.ThetaAxisLabels(2).String='06';
h.ThetaAxisLabels(3).String='12';
h.ThetaAxisLabels(4).String='18';
h.ThetaZeroLocation='bottom';
h.Colorbar.Label.String='Samples per 1° grid';
set(gca,'ColorScale','log','CLim',[1 10^4]);
colormap(inferno);
% set(gca,'FontName','palatino linotype','FontSize',12);
export_fig(strcat(storeImageDir,imageName,'.pdf'),'-r300','-pdf','-nocrop');
close(figureHandle);


