% Plot all the data

clear all;

load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/ground/asi_fy.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/ground/Espectra_pfisr.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/ground/electron_density_avg_pfisr.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/Espectra_thm.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/Global Measurements/aindex.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/Global Measurements/IMF_ACE.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/Bfield_thd.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/pa_anisotropy.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/Waves/fb_wv.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/Waves/FGM_xwv.mat')
load ('/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/space/Waves/EFI_ywv.mat')

%  Ionosphere ASI Keogram
i=1;
Data(i).time = asi_fy.time;
Data(i).ebin = asi_fy.pixels;
Data(i).E    = asi_fy.pow;

% Ionosphere Inversion (Energy Spectra)
i=2;
Data(i).time = data_pfisr.time;
Data(i).ebin = data_pfisr.ebin;
Data(i).E    = data_pfisr.E;

% Plasmasheet Measurement (Energy Spectra)
i=3;
Data(i).time = data_thm.time;
Data(i).ebin = data_thm.ebin;
Data(i).E    = data_thm.E;

%  AL Index
i=4;
Data(i).time = aindex.time;
Data(i).al   = aindex.al;
Data(i).au   = aindex.au;

%  Northward  Bz
i=5;
Data(i).time = IMF.time;
Data(i).Bz   = IMF.Bz;

% THEMIS Magnetic Field
i=6;
Data(i).time = thd_FGM.time;
Data(i).Bx   = thd_FGM.Bx;
Data(i).Bz   = thd_FGM.Bz;

% Ion Anisotropy
i=7;
Data(i).time_H = sst_ion.time;
Data(i).A_H    = sst_ion.A;
Data(i).time_L = esa_ion.time;
Data(i).A_L    = esa_ion.A;

% Electron Anisotropy
i=8;
Data(i).time_H = sst_e.time;
Data(i).A_H    = sst_e.A;
Data(i).time_L = esa_e.time;
Data(i).A_L    = esa_e.A;

%  Wave power (1-1000 Hz)
i=9;
Data(i).time = fb.time;
Data(i).ebin = fb.freq;
Data(i).E    = fb.pow;
Data(i).E(Data(i).E<=0)=10^-10;
Data(i).cmin = -2;
Data(i).cmax = 1;

%  Wave power (along y)
i=10;
Data(i).time = ywv.time;
Data(i).ebin = ywv.freq;
Data(i).E    = ywv.pow;
Data(i).E(Data(i).E<=0)=10^-10;
Data(i).cmin = -5;
Data(i).cmax = 2;

%  Wave power (along x)
i=11;
Data(i).time = xwvfgm.time;
Data(i).ebin = xwvfgm.freq;
Data(i).E    = xwvfgm.pow;
Data(i).E(Data(i).E<=0)=10^-10;
Data(i).cmin = -5;
Data(i).cmax = 2;

figure(1)
clf

p=panel();

p.pack(1);
p(1).pack(11);

p.marginleft=40;
p.marginright=20;
p(1).de.margin=4;
p.fontsize=12;
p.select('all');

DateNumBeg=datenum('2008-03-26/08:05:00');
DateNumEnd=datenum('2008-03-26/12:55:00');
dt=0.5;
TTick=[DateNumBeg,(floor(DateNumBeg)+(ceil((DateNumBeg-floor(DateNumBeg))/(dt/24))/(24/dt)):dt/24:DateNumEnd)];
  

q=p(1);
k=1;
for i=1:1:11
    q(k).select();
    if (i<=3)
        plot_2D_time_series(Data(i));
        set(gca,'XTickLabel',''); 
        if(i==1)
            caxis([log10(min(min(Data(i).E))) log10(max(max(Data(i).E)))-0.2]);
        end;
        
        k=k+1;
        if (i>1)
            set(gca,'YScale','log');
            set(gca,'YTick',[10^2 10^3 10^4 10^5],'YLim',[10^2 10^6]);
            
        end;
    elseif (i==4)
        plot(Data(i).time,Data(i).al/1000,'--black'); hold on;
        plot(Data(i).time,Data(i).au/1000,'-r');
        set(gca,'XTickLabel','','YLim',[-1.500 0.300]); 
        xlim([DateNumBeg DateNumEnd]);
        set(gca, 'Layer','top') 
        set(gca,'XTick',TTick,'ygrid','on');
        k=k+1;
        legend('AL','AU','Location','southwest','Orientation','horizontal');
   elseif (i==5)
        plot(Data(i).time,Data(i).Bz,'-black');
        set(gca,'XTickLabel','','YLim',[min(Data(i).Bz) max(Data(i).Bz)]); 
        xlim([DateNumBeg DateNumEnd]);
        set(gca, 'Layer','top') 
        set(gca,'XTick',TTick,'ygrid','on');
        k=k+1;
%         legend('IMF B_z','Location','southwest','Orientation','horizontal');
   elseif (i==6)
        plot(Data(i).time,Data(i).Bz,'--black'); hold on;
        plot(Data(i).time,Data(i).Bx,'-red');
        set(gca,'XTickLabel','','YLim',[min(Data(i).Bx)-5 max(Data(i).Bz)+5]); 
        xlim([DateNumBeg DateNumEnd]);
        set(gca, 'Layer','top') 
        set(gca,'XTick',TTick,'ygrid','on');
        k=k+1;
        legend('B_z','B_x','Location','northwest');  
   elseif (i==7)
        Data(i).time_H=(Data(i).time_H-floor(Data(i).time_H))+datenum('2008-03-26');
        Data(i).time_L=(Data(i).time_L-floor(Data(i).time_L))+datenum('2008-03-26');
        plot(Data(i).time_H,Data(i).A_H,'--black'); hold on;
        plot(Data(i).time_L,Data(i).A_L,'-red');
        set(gca,'XTickLabel','','YLim',[-1 1]); 
        xlim([DateNumBeg DateNumEnd]);
        set(gca, 'Layer','top') 
        set(gca,'XTick',TTick,'ygrid','on');
        k=k+1;
        legend('>25 keV','< 25 keV','Location','northwest','Orientation','horizontal');
    elseif (i==8)
        Data(i).time_H=(Data(i).time_H-floor(Data(i).time_H))+datenum('2008-03-26');
        Data(i).time_L=(Data(i).time_L-floor(Data(i).time_L))+datenum('2008-03-26');
        plot(Data(i).time_H,Data(i).A_H,'--black'); hold on;
        plot(Data(i).time_L,Data(i).A_L,'-red');
        set(gca,'XTickLabel','','YLim',[-1 1]); 
        xlim([DateNumBeg DateNumEnd]);
        set(gca, 'Layer','top') 
        set(gca,'XTick',TTick,'ygrid','on');
        k=k+1;     
        legend('>30 keV','<30 keV','Location','northwest','Orientation','horizontal');
    elseif (i>=9)
        plot_2D_time_series(Data(i));
        if(i~=11)
        set(gca,'XTickLabel',''); 
        end;
        caxis([Data(i).cmin Data(i).cmax]); 
        k=k+1;
        set(gca,'YScale','log');
        if(i==9)
            set(gca,'YTick',[10^1 10^2 10^3]);
        else
            set(gca,'YTick',[10^-2 10^-1 10^0]);
        end;
        if(i==11)

            set(gca,'XTickLabel',{datestr(TTick,'HH:MM')});
            set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
        end;
    end;
        %     else
%     plot(randn(10,2));
%     set(gca,'XTickLabel',[]);
%     end
%      
    hold on;
end;
%  q(6).select();
%  plot(randn(10,2));
 
q(1).ylabel('Pixels');
% q(2).ylabel('Energy [eV]');
q(3).ylabel('       Energy [eV]');
q(4).ylabel('[\muT]');
q(5).ylabel({'IMZ B_z','[nT]'});
q(6).ylabel({'THEMIS-D', 'B-Field','[nT]'});
q(7).ylabel({'Ion', 'Anisotropy'});
q(8).ylabel({'Electron', 'Anisotropy'});

q(10).ylabel('Frequency [Hz]');
% q(6).ylabel('Frequency [Hz]');
    
% p(1).xlabel('data value (furlongs per fortnight)');
% p(1).ylabel('normalised frequency (%)');

