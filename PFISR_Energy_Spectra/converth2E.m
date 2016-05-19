function [data] = converth2E( data_thm, filename, DateNumBeg, DateNumEnd )
%CONVERTA2E_V2 Summary of this function goes here

% h: [N x 1]  - A column vector of the altitude for which the 
%               for which in mean electron density is required [in km]

% Forward and Backward Model to convert PFISR data to energy
% Requirements
% +Add Class GeoData to path: run /home/nithin/NithinBU/PFISR/GeoDataMATLAB/setup.m
% +
% """
% Models unit incident flux from the magnetosphere ionizing N2,O,O2
% vs. altitude and inital beam energy in the ionosphere.
% The electron production rate q is related to auroral optical intensity via
% known proportions.
% Variables:
% q: e- volume production rate [cm^-3 s^-1]
% Phi: "flux" [cm^-2 s^-1 eV^-1]
% A: energy deposition matrix [eV cm^-1]
% References:
% Rees 1989
% Wedlund et al "Electron Energy Spectra and Auroral Arcs" JGR 2013
% Sergienko and Ivanov 1993 "A new approach to calculate the excitation of atmospheric gases by auroral electron impact"
% """
% Output:
% data.z        : Altitude in [km] 
% data.Ebin     : Energy bin values in [eV]
% data.A        : The Kernal for transformation q = A*phi from the REES Model
% data.ne       : Electron density from PFISR Measurements
% data.q_PFISR  : Production rate derived from PFISR ne measurements
% data.eflux    : Electron flux derived by using the inversion procedure
% data.q_THEMIS : Production rates estimated back using the derived
%                 electron flux in the following way q = A*phi_derived
% data.max_iter : Maximum number of iteration of the mem_solve inversion
% data.time_ne  : Time steps in matlab time for electron density
% data.time_q   : Time steps in matlab tinme for production rates (it is 3
%                 steps lesser than time_ne, because it is a derivative of ne
% data.chi      : The chi^2 value of the fitting/ inversion

% Loading data from the HDF5 file created by a Glow run of precipitating electrons
% fileName  ='/home/nithin/Dropbox/aurora_data (1)/eigenprofiles/GLOWrates.h5';
% fileName1 ='q_2008_03_26_rees.h5';
% fileName1 ='/home/nithin/Documents/git-repos/aurora-rep/To be reviewed/Scrap/2Dec15/A thd E Normal Ext.h5'

% q 		  = hdf5read(fileName1,'prod/eigenprofile');%Production rates (Nreactionz x Nalt x NEnergy)
% q         = q+min(min(q(q>0)));

% z         = hdf5read(fileName1,'altitude');			%Altitude
% % diffflux  = hdf5read(fileName,'diffnumflux');		%Differential number flux (cm^-2 s^-1 eV^-1)
% %In Rees the diffnumber flux is 1 cm^-2 s^-1 eV^-1
% %Hence converting it into [cm^-2 s^-1]

% diffflux  = ones(size(Ebins,1),1).*eb_binsize';		
% 
% sensorloc = hdf5read(fileName1,'sensorloc'); 		%Unit in degrees, Geographical coordinates
% Edim0     = 1:1:length(Ebins);						%Energy available for computation

% [phi]     = diag(diffflux(Edim0));					%Total Input differemtial number flux [cm^-2 sec^-1]
% phi_inv   = inv(phi);
% 
% 
% [Z,X]=meshgrid(z(zdim0),Ebins(Edim0)/(1000));			%For surface plots
% 



% Downloading the PFISR Data
% fileNameStr=downloadPFISR('2008-03-26 06:00','2008-03-26 15:00','yyyy-mm-dd HH:MM');
% [ne,altitude,t]=pfisr_mag_e('DataFile_2008_1.h5',1);

[ne_temp,altitude,t_temp]  =pfisr_avg_e_v2(filename); % Using the default value of the altitude grid

% Cropping the time array and density matrix to that prescribed by the user
[c,TimeBegIndex] = min(abs(t_temp-DateNumBeg));
[c,TimeEndIndex] = min(abs(t_temp-DateNumEnd));
t=t_temp(TimeBegIndex+1:TimeEndIndex);
ne=ne_temp(:,TimeBegIndex+1:TimeEndIndex);

zdim0     = 1:1:length(altitude);	%Altitude available for computation
Ebin_0    = 9;
Ebins	  = data_thm.ebin(Ebin_0:end)';	%Energy bin values

A =generate_A(altitude(zdim0),data_thm.ebin(Ebin_0:end)');
data.A=A;
%Plotting A
% figure; surf(X,Z,log10(A((zdim0-zdim0(1)+1),Edim0))','EdgeColor','none');
% set(get(gca,'XLabel'),'String','Energy [eV]');
% set(get(gca,'YLabel'),'String','Altitude [km]');
% title('A Matrix');
% view([0 90]);
% set(gca,'XScale','log');
% xlim([min(min(X)) max(max(X))]);
% ylim([min(min(Z)) max(max(Z))]);






ebin_yi  = interp1(1:1:size(Ebins,1),Ebins,0.5:1:(size(Ebins,1)+0.5),'linear','extrap');
ebin_high= ebin_yi(2:end);
ebin_low = ebin_yi(1:end-1);
eb_binsize=ebin_high-ebin_low;
% ne is in m^-3 so converting to cm^-3
ne=ne*10^-6;

% [ne1,altitude1,t1]=pfisr_avg_e('DataFile_2008_1.h5');
% min_ne=min(min(ne));
% % Removing NANs
% [jj ii] = find(isnan(ne.'));
% % ne(ii,jj) = (ne(ii-1,jj)+ne(ii+1,jj))/2;
% for i=1:1:size(ii)
%     ne(ii(i),jj(i))=min_ne;
% end;

% Removing nans
x=altitude(zdim0);
for ty=1:1:size(ne,2)
    y=ne(:,ty);
    xi=x(find(~isnan(y)));
    yi=y(find(~isnan(y)));
    ne(:,ty)=interp1(xi,yi,x,'linear','extrap');
%     ne(:,ty)= smooth(x,ne(:,ty),0.1,'moving');
end;
 
%% Section 2
% Implementing the Maximum Entropy Method
% Requires:
% zdim0    : Range of interested altitude
% altitude : An array of altitude values
% t        : Time array
% ne       : electron density in [cm^-3]
% Ebins    : Electron bin values
% A        : Requires matrix A

alpha   = 2.5*10^-12*exp(-altitude(zdim0)./51.2); % [m^3/s] %z is in [km]
% Converting to cm^3/s
alpha   = alpha*10^6;
t1      = t(1:1:end-2);
t2      = t(3:1:end);
dt      = 86400*(t2-t1);
dn      = ne(zdim0,3:1:end)-ne(zdim0,1:1:end-2);
[DT,DNE]= meshgrid(dt,ne(zdim0,1));
dn_dt   = dn./(2*DT);

% phi_new = 10^13*ones(size(diag(phi)',2),1);   %[cm^-2 s^-1]
kb = 1.38*10^-23;
eV = 1.602E-19;
T_w = 10*data_thm.ebin(25)*eV/kb;
phi_new = kappa(data_thm.ebin(Ebin_0:end),data_thm.ebin(25),T_w,160);
phi_new = 10^-3*phi_new*sum(data_thm.E(2,Ebin_0:end).*diff(data_thm.E(2,Ebin_0-1:end))); 
% phi_new = 10^-1*(-(Ebins-max(Ebins))).^2;
beta    = 20;

% for i = 300:1:size(t1)
tdim    = 1:1:size(t1);
q_1_1   = dn_dt(:,tdim);
q_1_2   = diag(alpha)*(ne(zdim0,tdim).^2);
q_1     = q_1_1 + q_1_2;                         %[cm^-3 s^-1]
var_q   =var(q_1);

data.q_PFISR = q_1; %Collecting production rate

W=ones(size(data.A,1),1);
W = gaussmf((1:1:size(W,1)),[20 60])'+1;
W=W./sum(W);

% poolobj= parpool(2);
% figure;
% parfor i=tdim
for i=tdim
      [phi_new1(:,i),q_new(:,i),chi2(i),max_iter(i)] =mem_solve(q_1(:,i),A,beta,phi_new',sqrt(var_q(i)),5000,W);
      display([num2str(i/max(tdim)*100), ' %',num2str(chi2(i))]);
%       subplot(2,1,1);
%       plot(phi_new1(:,i));hold on; plot(phi_new,'g');
%       subplot(2,1,2);
%       plot(A*phi_new1(:,i));hold on; plot(q_1(:,i),'g');
      
end;

[T1,X]        = meshgrid(t2,Ebins./(2*pi));

data.eflux    = phi_new1.*X; %Converting diff number flux into diff energy flux
data.q_THEMIS = q_new;
data.chi      = chi2;
data.max_iter = max_iter;
data.time_ne  = t;
data.time_q   = t(2:1:end-1); % Due to differentiation the time for production rates is different

data.ne=ne;             % Electron density in 
data.h=altitude;
data.ebin=data_thm.ebin(Ebin_0:end);
data.time_E=data_thm.time;

% delete(poolobj);

% %% Figure Plot phi_new1
% % figure; plot(z_2,q_1(:,303)); hold on; plot(z_2,q_new(:,303));
% figure; plot(t2,chi2);
% set(get(gca,'XLabel'),'String','Time');datetick('x', 'HH:MM');
% set(get(gca,'YLabel'),'String','Chi^2');
% 
% figure; plot(t2,max_iter);datetick('x', 'HH:MM');
% set(get(gca,'XLabel'),'String','Time');datetick('x', 'HH:MM');
% set(get(gca,'YLabel'),'String','Max_Iterations');
% 
% [T1,X]=meshgrid(data.time_q,data.Ebin);	
% figure; surf(T1,X,log10(data.eflux),'EdgeColor','none'); % Energy Flux [eV (cm^-2 sec^-1 eV^-1)]
% set(get(gca,'XLabel'),'String','Time');datetick('x', 'HH:MM');
% set(get(gca,'YLabel'),'String','Energy [keV]');
% title('Energy Spectrum derived from PFISR Electron Densities');
% view([0 90]);
% set(gca,'YScale','log');
% xlim([min(min(T1)) max(max(T1))]);
% ylim([min(min(X)) max(max(X))]);
% % 
% %% Figure Plot q_new
% [T1,Z]=meshgrid(data.time_q,data.z);	
% figure; surf(T1,Z,real(log10(data.q_THEMIS)),'EdgeColor','none');
% set(get(gca,'XLabel'),'String','Time');datetick('x', 'HH:MM');
% set(get(gca,'YLabel'),'String','Altitude [km]');
% title('Production rates derived from the Energy Spectrum');
% view([0 90]);
% set(gca,'YScale','log');
% xlim([min(min(T1)) max(max(T1))]);
% ylim([min(min(Z)) max(max(Z))]);
% % 
% %% Figure Plot q_old
% [T1,Z]=meshgrid(data.time_q,data.z);	
% figure; surf(T1,Z,real(log10(data.q_PFISR)),'EdgeColor','none');
% set(get(gca,'XLabel'),'String','Time');datetick('x', 'HH:MM');
% set(get(gca,'YLabel'),'String','Altitude [km]');
% title('Production rates derived from the Energy Spectrum');
% view([0 90]);
% set(gca,'YScale','log');
% xlim([min(min(T1)) max(max(T1))]);
% ylim([min(min(Z)) max(max(Z))]);
% 
% %% Figure Plot ne
% 
% [T,ALT]=meshgrid(t2,z_2);	
% figure; surf(T,ALT,real(log10(q_1)),'EdgeColor','none');
% set(get(gca,'XLabel'),'String','Time');datetick('x', 'HH:MM');
% set(get(gca,'YLabel'),'String','Altitude [km]');
% title('Electron Density measured from PFISR');
% view([0 90]);
% set(gca,'YScale','log');
% xlim([min(min(T)) max(max(T))]);
% ylim([min(min(ALT)) max(max(ALT))]);

end

