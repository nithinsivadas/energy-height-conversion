%% Plot 1D energy spectra of GEEPs
clear all;
h5FilePath = 'G:\My Drive\Research\Projects\Paper 2\Data\Events List\'; 

h5FileStr{1} = [h5FilePath,'20080326.001_bc_2min-energyFlux.h5'];
h5InputFileStr{1} = [h5FilePath,'20080326.001_bc_2min-fitcal.h5'];
thisTimeStr{1} = '26 Mar 2008 11:24';
thisTimeStrQuiet{1} = '26 Mar 2008 11:00';

h5FileStr{2} = [h5FilePath,'20100528.001_bc_2min-energyFlux.h5'];
h5InputFileStr{2} = [h5FilePath,'20100528.001_bc_2min-Ne-cal.h5'];
thisTimeStr{2} = '28 May 2010 10:40';
thisTimeStrQuiet{2} = '28 May 2010 09:45';

h5FileStr{3} = [h5FilePath,'20101018.001_bc_2min-energyFlux.h5'];
h5InputFileStr{3} = [h5FilePath,'20101018.001_bc_2min-Ne-cal.h5'];
thisTimeStr{3} = '18 Oct 2010 08:17';
thisTimeStrQuiet{3} = '18 Oct 2010 08:00';

%%
figure;

for i=1:1:length(h5InputFileStr)
    amisrData = read_amisr(h5InputFileStr{i});
    energyFlux=permute(h5read(h5FileStr{i},'/energyFluxFromMaxEnt/energyFlux'),[ 3 2 1]);
    energyBin=h5read(h5FileStr{i},'/energyFluxFromMaxEnt/energyBin')';
    time=h5read(h5FileStr{i},'/energyFluxFromMaxEnt/time')';
    thisTimeIndx = find_time(time,thisTimeStr{i});
    energyFlux1D = squeeze(energyFlux(:,amisrData.magBeamNo,:));
    energyGreaterThan10keV(i) = integrate_energy(energyFlux1D(thisTimeIndx,:),energyBin,10*1000,900*1000);
    energyGreaterThan30keV(i) = integrate_energy(energyFlux1D(thisTimeIndx,:),energyBin,30*1000,900*1000);
    energyGreaterThan100keV(i) = integrate_energy(energyFlux1D(thisTimeIndx,:),energyBin,100*1000,900*1000);
    totalEnergy(i) = integrate_energy(energyFlux1D(thisTimeIndx,:),energyBin,1*1000,900*1000);
    p1{i}=plot_1D_time_slice(time,energyBin'/1000,energyFlux1D',thisTimeStr{i},-1);
    xlabel('Energy [keV]');
    ylabel('Differential energy Flux [eV/m^2 sr s eV]');
    set(gca,'XTick',[3 10 30 100 300],'YLim',[10^8 10^12]);
    hold on;
%     p2{i}=plot_1D_time_slice(time,energyBin'/1000,squeeze(energyFlux(:,amisrData.magBeamNo,:))',thisTimeStrQuiet{i},-1);
%     hold on;
    clearvars amisrData energyFlux energyBin time 
   
end
% legend(thisTimeStr{1}, thisTimeStrQuiet{1},...
%     thisTimeStr{2}, thisTimeStrQuiet{2},...
%     thisTimeStr{3}, thisTimeStrQuiet{3});
legend(thisTimeStr{1},...
    thisTimeStr{2},...
    thisTimeStr{3});
set(p1{1},'Color','k');
% set(p2{1},'Color','k','LineStyle','--');
set(p1{2},'Color','r');
% set(p2{2},'Color','r','LineStyle','--');
set(p1{3},'Color','m');
% set(p2{3},'Color','m','LineStyle','--');
grid off;

disp('Percentage of energy > 10 keV');
(energyGreaterThan10keV./totalEnergy)*100

disp('Percentage of energy > 30 keV');
(energyGreaterThan30keV./totalEnergy)*100

disp('Percentage of energy > 100 keV');
(energyGreaterThan100keV./totalEnergy)*100

%% function to calculate the % of energy > 10 keV
function result = integrate_energy(energyFlux,energyBin,minE,maxE)
    X = minE:0.1:maxE;
    Y = interp1(energyBin,energyFlux,X,'linear');
    result = trapz(X,Y);
end
