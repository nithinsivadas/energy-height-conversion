%%  Routine to Plot PFISR based energy spectra with Error
% We use SIC model, and constant alpha model to calculate q
% And the maximum entrop method to invert the spectra

clear all;
initialize_geodata;

%Global Variables
computer=getenv('computername');
if strcmp(computer,'NITHIN-SURFACE')||strcmp(computer,'NITHIN-CARBON')
    pfisrDataFileNameStr =...
        'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\PFISR_Energy_Spectra\Data\DataFile_2008_1.h5';
    thmDataFileNameStr =...
        'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Data_Mar_08_Event\space\Espectra_thm.mat';
    sicDataFileNameStr='C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Data\Ground\summaryresQinversion.mat';
else 
    pfisrDataFileNameStr =...
        '/home/nithin/Documents/git-repos/energy-height-conversion/PFISR_Energy_Spectra/Data/DataFile_2008_1.h5';
    thmDataFileNameStr =...
        '/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/Espectra_thm.mat';
    sicDataFileNameStr='/home/nithin/Documents/git-repos/ISR_School_2016/Antti  SIC Model/summaryresQinversion.mat';
end

load (thmDataFileNameStr);

timeMin=datenum('26 Mar 2008 08:00');
timeMax=datenum('26 Mar 2008 13:00');
altMin = 60;
altMax = 150;


% data_thm -> E, ebin, time
thm.energyFlux = data_thm.E(:,19:1:end)'*10^4; %[eV m-2 sr-1 s-1 eV-1]
thm.energyBin  = data_thm.ebin(19:1:end)';
thm.time       = data_thm.time;
clear data_thm;
energyBin = thm.energyBin;

%%  Mode 1: From measurements averaged across the PFISR beams
[ePopDenAvg, altPop, timePop] = get_pfisr_variable(pfisrDataFileNameStr, 'popl', 0);
ePopDenAvg = ePopDenAvg';

% Tidying up Electron Density
[NeAvg, time] = crop_time(ePopDenAvg, timePop, timeMin, timeMax);
[NeAvg, alt] = crop_altitude(NeAvg', altPop,altMin, altMax);
NeAvg = interp_nans(NeAvg);
NeAvg(NeAvg<0)=10^9;

%%  Mode 1.1: Error in PFISR Measurements
[dePopDenAvg, altPop, timePop] = get_pfisr_variable(pfisrDataFileNameStr, 'dpopl', 0);
dePopDenAvg = dePopDenAvg';

% Tidying up Electron Density
[dNeAvg, time] = crop_time(dePopDenAvg, timePop, timeMin, timeMax);
[dNeAvg, alt] = crop_altitude(dNeAvg', altPop,altMin, altMax);
dNeAvg = interp_nans(dNeAvg);
dNeAvg(dNeAvg<0)=10^9;

% Generate Production Rate
% [qAvg,qAvgTime, alpha] = get_production_rate(NeAvg, alt, time, 2);
[dqAvg, qAvg, qAvgTime] = get_error_in_q(NeAvg, dNeAvg, alt, time, 2);

%% Invert production rates
A =get_energy_dep_matrix(alt,energyBin);
[dataAvg] = get_inverted_flux(qAvg, dqAvg, qAvgTime, alt, energyBin, A);



% Generate Production Rate
% [dqAvg0,dqAvgTime, alpha] = get_production_rate(dNeAvg, alt, time, 2);
% dqAvg = 2*alpha.*NeAvg.*dNeAvg; %Simon wedlund et al. 2013
% Invert production rates
% [dataAvgError0] = get_inverted_flux(dqAvg, dqAvgTime, alt, energyBin);
% dataAvgError = dataAvg;
% dataAvgError.energyFlux = (get_error_in_energyFlux(dqAvg', dataAvg.A, dataAvg.energyBin, dataAvg.energyFlux, dataAvg.time,-1))';

% [dE,dEFrac] = get_error_in_energyFlux(dq, A, dataAvg.energyBin, dataAvg.energyFlux, dataAvg.time);
% 882 is the number of singular values (A<0.0001) in the forward model
% dataAvgError.energyFlux = dataAvgError0.energyFlux-dataAvg.energyFlux;
% dataAvgError.flux = dataAvgError0.flux-dataAvg.flux;
%% Mode 2: From SIC Measurements

load (sicDataFileNameStr);
[ionrate_raw, sic.alt] = crop_altitude(ionrate_raw, h,altMin, altMax);
[ionrate_smoothed, sic.alt] = crop_altitude(ionrate_smoothed, h,altMin, altMax);
sic.qRaw = ionrate_raw*10^6;
sic.qSmooth = ionrate_smoothed*10^6;
sic.Ne = Nedata*10^6;
sic.NeSIC = NeSIC*10^6;
sic.time = t;
clear h; clear ionrate_raw; clear ionrate_smoothed; clear Nedata; clear NeSIC; clear t; clear UT;

% Invert production rates
[dataAvgSIC] = get_inverted_flux(sic.qSmooth',dqAvg(1:end-1,:), sic.time, sic.alt, energyBin, A);


%% Plotting 2-D Electron Density neasurements
% figure;
% 
% subplot (2,1,1)
% plot_2D_time_series(time, alt, log10(NeAvg), 0.5, 1);
% title('Uncorrected Electron Density: Averaged across all beams log_1_0 [m^-^3]');
% caxis([6, 12]);

%% Plotting 2-D Energy Spectra measurements
% figure;
% 
% subplot (3,1,1)
% plot_2D_time_series(dataAvg.time, energyBin, log10(dataAvg.energyFlux), 0.5, 3);
% title('Differential Energy Flux: Averaged across all beams using const \alpha_e_f_f log_1_0 [eV m^-^2 sr^-^1 s^-^1 eV^-^1]');
% caxis([6, 12]);
% 
% subplot (3,1,2)
% plot_2D_time_series(dataAvgSIC.time, energyBin, log10(dataAvgSIC.energyFlux), 0.5, 3);
% title('Differential Energy Flux: Averaged across all beams using SIC Model log_1_0 [eV m^-^2 sr^-^1 s^-^1 eV^-^1]');
% caxis([6, 12]);
% 
% subplot (3,1,3)
% plot_2D_time_series(thm.time, energyBin, log10(thm.energyFlux), 0.5, 3);
% title('Differential Energy Flux: THEMIS D Spacecraft log_1_0 [eV m^-^2 sr^-^1 s^-^1 eV^-^1]');
% caxis([6, 12]);

%%



thm.alt = dataAvg.alt;
thm.qForward = dataAvg.A*energy_to_num(thm.energyFlux, thm.time, thm.energyBin);
thm.NeForward = q_to_Ne(thm.qForward, thm.alt, thm.time);

%% Creating Panels
hFig=figure(2);
resize_figure(hFig);

clf
p=panel(hFig);
% p.pack(1,2);
panelSize = 80; %in mm
p.pack({{panelSize} {panelSize}},3);
p.marginleft=20;
p.marginright=10;
p.de.margin=4;
p(1,3).marginleft=20;
p(2,3).marginleft=20;
p(1).margintop=10;
p(2).margintop=6;
% p.fontsize=12;
p.select('all');
for iPanel=1
    for jPanel=1:1:3
        p(iPanel,jPanel).select();
        set(gca,'XTickLabel','');
    end;
end;
for iPanel=1:1:2
    for jPanel=2
        p(iPanel,jPanel).select();
        set(gca,'YTickLabel','');
    end;
end;

%% Plotting 1-D Validation Plots


plotThisTime(1,:)='26-Mar-2008 11:46';
plotThisTime(2,:)='26-Mar-2008 11:13';

color1='k'; lineStyle1='-'; lineMarker1 ='none'; lineWidth1=1;
color2='k'; lineStyle2='--'; lineMarker2 ='none'; lineWidth2=1;
color3='r'; lineStyle3='-'; lineMarker3 ='none'; lineWidth3=1;
color4='m'; lineStyle4='-.'; lineMarker4 ='none'; lineWidth4=1;
color5='b'; lineStyle5='-'; lineMarker5 ='none'; lineWidth5=1;

for iTime = 1:1:2

    p(iTime,1).select();
   
    q(1)=plot_1D_time_slice_with_error(dataAvg.time, dataAvg.alt, NeAvg, dNeAvg, plotThisTime(iTime,:), 0);
    q(1).Color = color1; q(1).LineStyle = lineStyle1; q(1).Marker = lineMarker1; q(1).LineWidth = lineWidth1;
    hold on;
    q(2)=plot_1D_time_slice(dataAvg.time, dataAvg.alt, dNeAvg, plotThisTime(iTime,:), 0);
    q(2).Color = color2; q(2).LineStyle = lineStyle2; q(2).Marker = lineMarker2; q(2).LineWidth = lineWidth2;
    hold on;

    q(3)=plot_1D_time_slice(dataAvg.time, dataAvg.alt, q_to_Ne(dataAvg.qInverted', dataAvg.alt, dataAvg.time), plotThisTime(iTime,:), 0);
    q(3).Color = color4; q(3).LineStyle = lineStyle4; q(3).Marker = lineMarker4; q(3).LineWidth = lineWidth4;
    hold on;
    q(4)=plot_1D_time_slice(thm.time, thm.alt, thm.NeForward, plotThisTime(iTime,:), 0);
    q(4).Color = color5; q(4).LineStyle = lineStyle5; q(4).Marker = lineMarker5; q(4).LineWidth = lineWidth5;
    hold on;
    
    x_lim=[10^9 10^12];
    y_lim=[60 150];
    xlim(x_lim);
    ylim(y_lim);
    
    if iTime==1
     text(0.025,0.975,'a)','Units','normalized','FontWeight','bold');
     patch([x_lim(1) x_lim(2) x_lim(2) x_lim(1)], [y_lim(1) y_lim(1) 85 85], 'k', 'FaceAlpha',0.2,'EdgeColor','none');
    else
     text(0.025,0.975,'b)','Units','normalized','FontWeight','bold');
     patch([x_lim(1) x_lim(2) x_lim(2) x_lim(1)], [y_lim(1) y_lim(1) 70 70], 'k', 'FaceAlpha',0.2,'EdgeColor','none');
    end
    
    set(gca,'XTick',10.^[9, 10, 11, 12]);
    if(iTime==2)
    set(gca,'XTickLabel',{'10^9', '10^1^0', '10^1^1', '10^1^2'});
    lgd=legend(char({'Measured','by PFISR'}),...
        char({'Noise in','PFISR measurement'}),...
        char({'Derived from', 'forward modelling','of inverted result','(Vickrey)'}),...
        char({'Derived from','THEMIS','measurement'}),'Location','southoutside');
        set(lgd,'FontSize',8);
    end;
    p(iTime,2).select();
    q(1)=plot_1D_time_slice_with_error(dataAvg.time, dataAvg.alt, qAvg', dqAvg', plotThisTime(iTime,:), 0);
    q(1).Color = color1; q(1).LineStyle = lineStyle1; q(1).Marker = lineMarker1; q(1).LineWidth = lineWidth1;
    hold on;

    q(2)=plot_1D_time_slice(dataAvg.time, dataAvg.alt, dqAvg', plotThisTime(iTime,:), 0);
     q(2).Color = color2; q(2).LineStyle = lineStyle2; q(2).Marker = lineMarker2; q(2).LineWidth = lineWidth2;
    hold on;

    q(2)=plot_1D_time_slice(sic.time, sic.alt, sic.qSmooth, plotThisTime(iTime,:), 0);
    q(2).Color = color3; q(2).LineStyle = lineStyle3; q(2).Marker = lineMarker3; q(2).LineWidth = lineWidth3;
    hold on;

    q(3)=plot_1D_time_slice(dataAvg.time, dataAvg.alt, dataAvg.qInverted', plotThisTime(iTime,:), 0);
    q(3).Color = color4; q(3).LineStyle = lineStyle4; q(3).Marker = lineMarker4; q(3).LineWidth = lineWidth4;
    hold on;

    q(4)=plot_1D_time_slice(thm.time, thm.alt, thm.qForward, plotThisTime(iTime,:), 0);
    q(4).Color = color5; q(4).LineStyle = lineStyle5; q(4).Marker = lineMarker5; q(4).LineWidth = lineWidth5;
    hold on;
    
    x_lim=[10^6 10^12];
    y_lim=[60 150];
    xlim(x_lim);
    ylim(y_lim);
    
    if iTime==1
     text(0.025,0.975,'c)','Units','normalized','FontWeight','bold');
     patch([x_lim(1) x_lim(2) x_lim(2) x_lim(1)], [y_lim(1) y_lim(1) 85 85], 'k', 'FaceAlpha',0.2,'EdgeColor','none');
    else
     text(0.025,0.975,'d)','Units','normalized','FontWeight','bold');
     patch([x_lim(1) x_lim(2) x_lim(2) x_lim(1)], [y_lim(1) y_lim(1) 70 70], 'k', 'FaceAlpha',0.2,'EdgeColor','none');
    end
    
     title(datestr(dataAvg.time(find_time(dataAvg.time, (plotThisTime(iTime,:)))),'DD-mmm-YYYY HH:MM'),'FontSize',12);
     set(gca,'XTick',10.^[7 8 9 10 11]);
     if(iTime==2)
%     lgd=legend('q_P_F_I_S_R','\Delta q_P_F_I_S_R ','q_P_F_I_S_R uing SIC',...
%         'q_P_F_I_S_R Inverted', 'q_T_H_E_M_I_S','Location','southoutside',...
%         'Orientation','horizontal');

    lgd=legend(char({'Derived from','PFISR & (Vickrey)'}),...
        char({'Noise derived from','PFISR & (Vickrey)'}),...
        char({'Derived from','PFISR & (SIC Model)'}),...
        char({'Derived from forward', 'modelling of inverted','result & (Vickrey)'}),...
        char({'Derived from THEMIS','measurement'}),'Location','southoutside');
    set(gca,'XTickLabel',{'10^7', '', '10^9' , '', '10^1^1'});
        set(lgd,'FontSize',8);
    end;


    p(iTime,3).select();
    q(1)=plot_1D_time_slice_with_error(dataAvg.time, dataAvg.energyBin, dataAvg.energyFlux, dataAvg.dEnergyFlux, plotThisTime(iTime,:), -1);
    q(1).Color = color1; q(1).LineStyle = lineStyle1; q(1).Marker = lineMarker1; q(1).LineWidth = lineWidth1;
    hold on;

%     q(2)=plot_1D_time_slice(dataAvg.time, dataAvg.energyBin, dataAvg.dEnergyFlux, plotThisTime(iTime,:), -1);
%     q(2).Color = color2; q(2).LineStyle = lineStyle2; q(2).Marker = lineMarker2; q(2).LineWidth = lineWidth2;
%     hold on;

%     q(2)=plot_1D_time_slice(dataAvgSIC.time, dataAvgSIC.energyBin, dataAvgSIC.energyFlux, plotThisTime(iTime,:), -1);
%     q(2).Color = color3; q(2).LineStyle = lineStyle3; q(2).Marker = lineMarker3; q(2).LineWidth = lineWidth3;
%     hold on;
    q(4)=plot_1D_time_slice(thm.time, thm.energyBin, thm.energyFlux, plotThisTime(iTime,:), -1);
    q(4).Color = color5; q(4).LineStyle = lineStyle5; q(4).Marker = lineMarker5; q(4).LineWidth = lineWidth5;
    hold on;

    
    x_lim=[10^3 10^6];
    y_lim=[10^6 2*10^12];
    xlim(x_lim);
    ylim(y_lim);
    
    if iTime==1
     text(0.025,0.975,'e)','Units','normalized','FontWeight','bold');
     patch([10^5 x_lim(2) x_lim(2) 10^5], [y_lim(1) y_lim(1) y_lim(2) y_lim(2)], 'k', 'FaceAlpha',0.2,'EdgeColor','none');
    else
     text(0.025,0.975,'f)','Units','normalized','FontWeight','bold');
     patch([4*10^5 x_lim(2) x_lim(2) 4*10^5], [y_lim(1) y_lim(1) y_lim(2) y_lim(2)], 'k', 'FaceAlpha',0.2,'EdgeColor','none');
    end
    
    set(gca,'XTick',10.^[3 4 5 5.477 6]);    
    if(iTime==2)
    set(gca,'XTickLabel',{'1', '10', '100', '300', ''});
        lgd=legend(char({'Inversion result','from PFISR','(Vickrey)'}),...
        char({'Noise from','PFISR measurement','(Vickrey)'}),...
        char({'Inversion result ','from PFISR','(SIC Model)'}),...
        char({'Measured by','THEMIS'}),'Location','southoutside');
    set(lgd,'FontSize',8);
    end;
end;
p.pack({{60} {60}},3);
p(1,2).select();
 title(datestr(dataAvg.time(find_time(dataAvg.time, (plotThisTime(1,:)))),'DD-mmm-YYYY HH:MM'),'FontWeight','bold');
p(2,2).select();
 title(datestr(dataAvg.time(find_time(dataAvg.time, (plotThisTime(2,:)))),'DD-mmm-YYYY HH:MM'),'FontWeight','bold');
p(1,1).select();
ylabel('Altitude [$km$]','Interpreter','latex');
set(gca,'YTick',[60, 70, 80, 90, 100, 110, 120, 130, 140, 150]);
p(2,1).select();
ylabel('Altitude [$km$]','Interpreter','latex');
set(gca,'YTick',[60, 70, 80, 90, 100, 110, 120, 130, 140, 150]);
xlabel('Electron density [$m^{2}$]','FontWeight','bold','Interpreter','latex');
p(2,2).select();
xlabel('Production rate [$m^{-3} s^{-1}$]','FontWeight','bold','Interpreter','latex');
p(1,3).select();
ylabel('Diff. Energy Flux [$eV m^{-2} sr^{-1} s^{-1} eV^{-1}$]','Interpreter','latex');
set(gca,'YTick',10.^[6, 7, 8, 9, 10, 11, 12]);
p(2,3).select();
xlabel('Energy [$keV$]','FontWeight','bold','Interpreter','latex');
ylabel('Diff. Energy Flux [$eV m^{-2} sr^{-1} s^{-1} eV^{-1}$]','Interpreter','latex');
set(gca,'YTick',10.^[6, 7, 8, 9, 10, 11, 12]);
%%
if computer=='NITHIN-SURFACE'
    export_fig 'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Plots\Draft 12 Plots\QualityOfInversion.pdf' -pdf 
    save('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Plots\Draft 12 Plots\QualityOfInversion.svg');
    export_fig 'C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Plots\Draft 12 Plots\QualityOfInversion.png' -r600 -png -nocrop
else
    export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/QualityOfInversion.pdf' -pdf 
    save('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/QualityOfInversion.svg');
    export_fig '/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Plots/Draft10_Plots/PDF/QualityOfInversion.png' -r600 -png -nocrop
end