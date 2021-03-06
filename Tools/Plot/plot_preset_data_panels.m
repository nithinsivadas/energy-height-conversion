function plot_preset_data_panels( dataChoice )
%% plot_preset_data_panels Plots data from the choices in panel form 
% Given below is the list of pane;ls to plot in the order that you wish in
% a letter-size page
%--------------------------------------------------------------------------
%Input
%-----
% dataChoice - Array of numbers - each number associated with a data set 
%              to be plotted in order
%----------------------------------------------------------------------------
% Modified: 24th Jan 2017 [needs more clarification]
% Created : 24th Jan 2017
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------


include_new_colormaps;
data=get_all_time_series_data();
totalDataNo=size(data,2);

if nargin<1
    display('Global Parameters');
    display('------------------');
    display('1: IMF Bz');
    display('2: AL & AU indices');
    display('                   ');
    display('Ground Measurements');
    display('--------------------');
    display('3: ASI Keogram: FYKN');
    display('4: Electron density from PFISR');
    display('5: Ionosphere Inverse Energy Spectra (Vickrey Model)');
    display('6: Ionosphere Inverse Energy Spectra (SIC Model)');
    display('                   ');
    display('Space Measurements');
    display('-------------------');
    display('7: B field Z and X from FGM ');
    display('8: Energy spectra from thd');
    display('9: Ion Anisotropy');
    display('10: Electron Anisotropy');
    display('11: Wave power 1-1000Hz from FBK');
    display('12: Wave power <1Hz from SCM along X');
    display('13: Wave power <1Hz from EFI along Y ');
    display('14: Wave power <1Hz from FGM along X');
    display('15: ASI Keogram: GAKO');
    display('                   ');

    prompt = 'Select plots: ';
    dataChoice=input(prompt);
end;

% dataChoice=[1,2,7,9,10,3];

totalPanelNo=length(dataChoice);

hFig=figure(1);


clf
p=panel();
p.pack(1);

panelSize = 30; %in mm
demargin = 4;
panelDefinition=cell(1,totalPanelNo);
for i=1:1:totalPanelNo
    panelDefinition(i)={{panelSize}}; 
end
p(1).pack(panelDefinition);

p.marginleft=35;
p.marginright=25;
p(1).de.margin=demargin;
% p.fontsize=12;
p.select('all');
timeMinStr='2008-03-26/08:05:00';
timeMaxStr='2008-03-26/12:55:00';
DateNumBeg=datenum(timeMinStr);
DateNumEnd=datenum(timeMaxStr);
timeTick=0.5;
  
resize_figure(hFig, panelSize*totalPanelNo+4*(totalPanelNo+1)+20);

q=p(1);
k=1;
for thisData=dataChoice
    q(k).select();
    text(-0.175,0.975,[char(96+k),')'],'Units','normalized','FontWeight','bold');
    hold on;
    if (thisData==3 || thisData==4 || thisData==5 || thisData==6 ||...
            thisData==8 || thisData==11 || thisData==12 || thisData==13 ||...
            thisData==14 || thisData==15 || thisData==16 ||thisData==17 ||...
            thisData==18)
        
        plot_2D_time_series(data(thisData).time',data(thisData).yAxis, log10(data(thisData).zValue),0,datestr(DateNumBeg),datestr(DateNumEnd));
        label_time_axis(data(thisData).time, false, timeTick, timeMinStr, timeMaxStr);
        % Setting Grid Lines
        grid off;

        % Making sure the color bar will be outside the axes
        axPos = get(gca, 'position');
        colormap(inferno);
        
        if (thisData==16 ||thisData==17 || thisData==18)        
            hold on;
            h_lgd=plot(data(thisData).time, data(thisData).yValue, 'Color', 'c');
            legend(h_lgd,'Median energy','Location','northwest','Orientation','vertical');
            hold off;
        end;
        
        if (thisData==13)        
            hold on;
            HI_lgd=plot(data(thisData).time, data(thisData).y.f_HI, '--y');
            OI_lgd=plot(data(thisData).time, data(thisData).y.f_OI, '-.g');
            HeI_lgd=plot(data(thisData).time, data(thisData).y.f_HeI, '-c');
            HeII_lgd=plot(data(thisData).time, data(thisData).y.f_HeII, ':c');
            legend([HI_lgd,HeII_lgd,HeI_lgd,OI_lgd],...
                'H^+','He^2^+','He^+','O^+','Location','southwest','Orientation','vertical');
%             legend(OI_lgd,'f_O_+','Location','southwest','Orientation','vertical');
            hold off;
        end;
        
        
        if(thisData~=3 || thisData~=15) % If data is not an ASI Keogram
            set(gca,'Yscale','log','YTick',data(thisData).y.tick,...
                'YTickLabel',data(thisData).y.tickLabel,...
                'YLim',data(thisData).y.lim);
            c=colorbar('eastoutside','YLim',log10(data(thisData).color.lim),...
                'YTick',log10(data(thisData).color.tick),...
                'YTickLabel', data(thisData).color.tickLabel);
            caxis(log10(data(thisData).color.lim));
            
        else
            c=colorbar('eastoutside','YLim',log10(data(thisData).color.lim),...
                'YTick',log10(data(thisData).color.tick),...
                'YTickLabel', data(thisData).color.tickLabel);
            caxis(log10(data(thisData).color.lim));
 
        end;
        ylabel(c,data(thisData).color.label,'FontSize',8);
        c.TickLength=0.05;
        set(gca,'position',axPos);
        
        % reducing color bar thickness
        cPos=get(c,'Position');
        cPos(3)=0.2*cPos(3);
        set(c, 'Position',cPos);
        
    elseif (thisData==1)
        plot(data(thisData).time, data(thisData).yAxis,'k');
        label_time_axis(data(thisData).time, false, timeTick, timeMinStr, timeMaxStr);
         set(gca,'YTick',data(thisData).y.tick,...
             'YTickLabel',data(thisData).y.tickLabel,...
             'YLim',data(thisData).y.lim);
         grid on;
    elseif (thisData==2||thisData==7||thisData==19)
        plot(data(thisData).time, data(thisData).yAxis,'k');
        hold on;
        plot(data(thisData).time, data(thisData).yAxis2,'--r');
        label_time_axis(data(thisData).time, false, timeTick, timeMinStr, timeMaxStr);
         set(gca,'YTick',data(thisData).y.tick,...
             'YTickLabel',data(thisData).y.tickLabel,...
             'YLim',data(thisData).y.lim);
         legend(data(thisData).legend,'Location','northwest','Orientation','horizontal');
         grid on;
    elseif (thisData==9||thisData==10)
        plot(data(thisData).time, data(thisData).yAxis,'k');
        hold on;
        plot(data(thisData).time2, data(thisData).yAxis2,'--r');
        label_time_axis(data(thisData).time2, false, timeTick, timeMinStr, timeMaxStr);
         set(gca,'YTick',data(thisData).y.tick,...
             'YTickLabel',data(thisData).y.tickLabel,...
             'YLim',data(thisData).y.lim);
         grid on;
         legend(data(thisData).legend,'Location','northwest','Orientation','horizontal');
    end;    
    q(k).fontsize=10;
    ylabel(data(thisData).label,'FontSize',8);
    
    %(data(thisData).label));
    hold on;
    k=k+1;
    box on;
end;
computer=getenv('computername');
if computer=='NITHIN-SURFACE'
    load ('C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\Projects\Paper 1\Data\Space\thd_data_26_Mar_2008.mat')
else
    load ('/home/nithin/Documents/git-repos/energy-height-conversion/Tools/Paper 1/Data/Space/thd_data_26_Mar_2008.mat')
end;
label_time_axis(data(thisData).time, true, timeTick, timeMinStr, timeMaxStr);
[TTick, DateNumBeg, DateNumEnd] = get_axes_time_tick_values( data(thisData).time, timeTick, timeMinStr, timeMaxStr );
TTickLim=[DateNumBeg DateNumEnd];
add_horizontal_axes( TTick, TTickLim, thd.state.time, thd.state.mlt, 'Thm-D MLT [Hr]', 2);
add_horizontal_axes( TTick, TTickLim, thd.state.time, thd.state.mlat, 'Thm-D MLAT [deg]', 3);
add_horizontal_axes( TTick, TTickLim, thd.state.time, thd.state.lShell, 'Thm-D LShell', 4);
set(gcf,'color','w'); 

% 
% q(1).ylabel('Pixels');
% % q(2).ylabel('Energy [eV]');
% q(3).ylabel('       Energy [eV]');
% q(4).ylabel('[\muT]');
% q(5).ylabel({'IMZ B_z','[nT]'});
% q(6).ylabel({'THEMIS-D', 'B-Field','[nT]'});
% q(7).ylabel({'Ion', 'Anisotropy'});
% q(8).ylabel({'Electron', 'Anisotropy'});
% 
% q(10).ylabel('Frequency [Hz]');
% q(6).ylabel('Frequency [Hz]');
    
% p(1).xlabel('data value (furlongs per fortnight)');
% p(1).ylabel('normalised frequency (%)');



end

