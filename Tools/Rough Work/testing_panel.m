hFig=figure(2);
resize_figure(hFig);

clf
p=panel(hFig);
% p.pack(1,2);
panelSize = 30; %in mm
p.pack({{60} {60}},3);
p.marginleft=20;
p.marginright=10;
p.de.margin=4;
p(1,3).marginleft=15;
p(2,3).marginleft=15;
p(1).margintop=10;
p(2).margintop=10;
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
p(1,1).select();
ylabel('Altitude [km]');
p(2,1).select();
ylabel('Altitude [km]');
xlabel('Electron density [m^-^3]');
p(1,3).select();
ylabel('Energy [keV]');
p(2,3).select();
ylabel('Energy [keV]');
xlabel('Energy Flux [eV/m^2 sr s eV^-^1');


