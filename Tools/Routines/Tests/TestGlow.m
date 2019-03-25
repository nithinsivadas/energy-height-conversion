clear all;

time = datenum(2008,03,26,11,30,0);
glat = 65.1;
glon = -147.5;
[f107a, f107, aph] = f107_aph(time);
%%
Ap = aph(:,1);

%% solar radio flux [10-22 W m-2]
f107p = (f107a + f107)./2;
%% Particle flux [erg cm-2 == mW m-2]
Q= 0.5;
%% characteristic energy [eV]
Echar = 100e3;

% Energy Grid
Emin = 1;
Emax = 1e6;
Nbins = 250;
[Ebins, Phitop] = monoenergetic(Emin, Emax, Nbins, Echar); 
% Phitop = energetic electron differential number flux into top of atmosphere: cm-2 s-1 eV-1
Phitop = Phitop*100;
diffE = diff(Ebins);
diffE(end+1) = diffE(end);
iono0 = glowenergy(time, glat, glon, f107a, f107, f107p, Ap, Ebins, Phitop*0);
%% sergeinko ivanov production rates
alt = 70:1:200;
A = get_energy_dep_matrix(alt', Ebins', glat.*ones(length(alt),1),...
    glon.*ones(length(alt),1), time);
A = A/100;
A(isnan(A))=0;
%% Convert A from [m-1 eV] to [cm-1 eV]

%% glow model
% Axxxx wavelength in angstrom, intensity in Rayleigh 10^6 photons cm-2
% density cgs cm-3
figure;
ECharArr = [10,30,100,300];
colorArr = ['k','m','c','b'];
for i = 1:1:length(ECharArr)
    Echar = ECharArr(i).*1e3;
    [Ebins, Phitop] = monoenergetic(Emin, Emax, Nbins, Echar); 
    Phitop=Phitop*1000;
    iono = glowenergy(time, glat, glon, f107a, f107, f107p, Ap, Ebins, Phitop);
    q = A*Phitop';
    semilogx(q,alt,'Color',colorArr(i),'LineStyle','--');
    hold on;
    semilogx(iono.ionrate-iono0.ionrate,iono.altkm,'Color',colorArr(i));
    
end
ylim([60 200]);
xlim([1 10^6]);
legend('S&I:10 keV','G:10 keV','S&I:30 keV','G:30 keV','S&I:100 keV','G:100 keV','S&I:300 keV','G:300 keV');
%%

% function iono = run_glow(time, glat, glon, f107a, f107, f107p, Ap, Q, Echar)
% 
% currentPath = pwd;
% cd("C:\Users\nithin\Documents\GitHub\NCAR-GLOW\matlab");
% iono = glow(time, glat, glon, f107a, f107, f107p, Ap, Q, Echar);
% cd(currentPath);
% 
% end