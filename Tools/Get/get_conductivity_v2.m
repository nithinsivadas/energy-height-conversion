function [ data, input ] = get_conductivity_v2( alt, electronDensity,...
    latitude, longitude, time, mode, outputs, setPlotConductivity,...
    iriData, msisData)
%get_conductivity_v2.m returns the hall and pederson conductivity given an
%input electron density value, along with some other output parameters.
%Uses msis and iri models that run locally, instead on online servers.
%-------------------------------------------------------------------------
% Input
%------
%   alt                      - Altitude array [nHx1] [km]
%   electronDensity          - Electron denisity input [nHx1] [m^-3]
%   latitude                 - Latitude of the measurment station [1x1]
%   longitude                - Longitude of the measurement station [1x1]
%   time                     - Date-time measurments in matlab units
%	  mode                   - 0: uses electronDensity as input
%                              1: uses IRI2016 Ne as input
%-------------------------------------------------------------------------
% Output
%---------
%   outputArguments.sigma_P  - Pederson conductivity [S/m]
%   outputArguments.sigma_H  - Hall conductivity [S/m]
%   outputArguments.alt      - Altitude [km]
%   and  a lot of other things
%----------------------------------------------------------------------------
% Modified: 20th Nov 2018,12th Jan 2018
% Created : 27th Sep 2017
% Author  : Nithin Sivadas
% Ref     :
% Comments: 1. The formulae needs to be looked at carefully. [Pending]
%           2. Speed up the code using matrix multiplication instead of
%              for-loops
%           3. Inputs only one lat, and lon, check if this is acceptable assumption. 
%----------------------------------------------------------------------------

% Initilization
    if nargin <10
        msisData = [];
    end
    
    if nargin < 9
        iriData = [];
    end
    
    if nargin<8
        setPlotConductivity = false;
    end

    if nargin<7 || isempty(outputs)
        outputs = {'default'}; % uses input electronDensity to estimate conductivity
    end

	if nargin<6 || isempty(mode)
        mode = 0; % uses input electronDensity to estimate conductivity
    end

	if nargin<5 || isempty(time)
	    time        = datenum([2008 03 26 10 00 00]);
	end

	if nargin<4 || isempty(longitude)
	    longitude   = -147.5; % Degrees.
	end

	if nargin<3 || isempty(latitude)
	    latitude    = 65; % Degrees.
    end

% IRI2016
    if isempty(iriData)
        iriData = iri2016f90(time, alt, latitude, longitude);
    end
    Ti = iriData.Ti';
    Te = iriData.Te';

% Concentration of ion species
    C = [iriData.nOI',... %C(:,1) = [O+]
        iriData.nNOI',... %C(:,2) = [NO+]
        iriData.nO2I',... %C(:,3) = [O2+]
        iriData.nNI',...  %C(:,4) = [N+]
        iriData.nHI',...  %C(:,5) = [H+]
        iriData.nHeI']... %C(:,6) = [He+]
        ./iriData.totalIonDensity';

% Correcting IRI [If concentrations are negative]
    C(C<0)=0;


% MSIS
	alt1=alt(:)';
    lat1=latitude*ones(size(alt1));
    lon1=longitude*ones(size(alt1));
	
    if isempty(msisData)
        msisData = msis_irbem(time, [alt1',lat1',lon1']);
    end
    
    Nn(:,1)     = msisData.O; %O
	Nn(:,2)     = msisData.N2; %N2
	Nn(:,3)     = msisData.O2; %O2

  Tn = msisData.AltTemp;

% B field
    [Bx, By, Bz] = igrf(time, latitude, longitude, alt1, 'geodetic');
    B = ((Bx(:).^2 + By(:).^2 + Bz(:).^2).^0.5)*10^-9;

% Coefficients of ion-neutral interactions % Shunk and Negy, Table 4.4
           %   O       N2      O2
    C_in = [   0   ,   6.84,   6.64;   % OI
               2.44,   4.34,   4.27;   % NOI
               2.31,   4.13,   0   ;   % O2I
               4.42,   7.47,   7.25;   % NI
               0   ,   33.6,   32.0;]; % HI

% Converting to SI units
    C_in = C_in*10^-16;

% Necessary Constants
    konst = define_universal_constants();
    m_e = konst.me; % Electron mass
    m_p = konst.mp; % Proton mass
    Z_n = [16; 28; 32];
    m_n = Z_n*m_p; % Mass of neurtrals[kg]
    Z_i = [16; 30; 32; 14]; % Ion mass in proton mass units
    m_i = Z_i*m_p;
    q_i = [1;1;1; 1]; % Ion charge

% Gyro frequencies
    q = konst.e;
    w_e  = q.*B./m_e;
    w_p = q.*B./m_p;

    for i = 1:1:length(m_i)
        w_i(:,i)=q_i(i)*w_p/Z_i(i);
    end

    v_i = zeros(length(alt),length(m_i));
    Tr = (Ti+Tn)/2;

% Ion-neutral collision frequency
    for i=1:length(m_i)
        for n = 1:length(m_n)

            if C_in(i,n) == 0
                if i==1 && n==1
                    %OI-O
                    v_in(:,i,n) = (1.7*3.42*10^-17).*Nn(:,n).*((Tr).^0.5).*...
                    (1.08 - 0.139.*log10(Tr) + (4.51*10^-3 ).*(log10(Tr)).^2);
                elseif i==3 && n==3
                    %O2I - O2
                    v_in(Tr>800,i,n) = (2.4*10^-19).*Nn((Tr>800),n).*(Tr(Tr>800)).*(10.4-0.76.*log10((Tr(Tr>800)))).^2;
                    v_in(Tr<=800,i,n) = (4.08*10^-16).*Nn((Tr<=800),n);
                elseif i==5 && n==1
                    %H+-O
                    v_in(Tr>300,i,n) = (6.61*10^-17).*Nn((Tr>300),n).*(Tr(Tr>300).^0.5).*(1-0.047.*log10((Tr(Tr>300)))).^2;
                end
            else
                v_in(:,i,n) = C_in(i,n).*Nn(:,n);
            end

            v_i(:,i) = v_i(:,i) + v_in(:,i,n);

        end
    end

% Electron-neutral Collision frequency
    v_eO  = (8.9*10^-17).*Nn(:,1).*(1+(5.7*10^-4).*Te).*Te.^0.5;
    v_eO2 = (1.82*10^-16).*Nn(:,3).*(1+(3.6*10^-2).*Te.^0.5).*Te.^0.5;
    v_eN2 = (2.33*10^-17).*Nn(:,2).*(1+(1.21*10^-4).*Te).*Te;
    v_en = v_eO + v_eO2 + v_eN2;

% Electron density assignment
    if mode==0
        ne = electronDensity(:);
    else
        ne = iriData.Ne(:);
    end
%     ne(ne<=0)=nan;
%     ne = interp_nans(ne);

% Electron-ion collision frequency
%   v_ei = ((ne.*(10^-6)).*Te.^-1.5).*(34 + 4.18.*log(((Te.^3)./(ne.*(10^-6))))); % Michael C Kelley
    v_ei = ((ne.*(10^-6)).*Te.^-1.5).*(59 + 4.18.*log10(((Te.^3).*(ne)).^3)); %Michel C Kelley?
%     v_ei_1=  54.5*(ne.*(10^-6).*Te.^-1.5).*(C(:,1).*q_i(1).^2 + C(:,2).*q_i(2).^2 + C(:,3).*q_i(3).^2);

% Electron collision frequency
    v_e   = v_en+v_ei;

% Mobility
    for i=1:length(m_i)
        k_i(:,i) = w_i(:,i)./v_i(:,i);
    end
    k_e = w_e./v_e;

% Conductivity
    sigma_P = zeros(length(alt1),1);
    sigma_P1 = zeros(length(alt1),1);
    sigma_H = zeros(length(alt1),1);
    sigma_H1 = zeros(length(alt1),1);
    for i=1:1:length(m_i)
        sigma_P =  sigma_P + C(:,i).*k_i(:,i)./(1+(k_i(:,i).^2));
        sigma_H1 =  sigma_H1 + C(:,i).*(k_i(:,i).^2)./(1+(k_i(:,i).^2));
    end
    sigma_P = (q*ne./B).*(sigma_P + (k_e)./(1+k_e.^2));
    sigma_H = -1*(q*ne./B).*(sigma_H1 - (k_e.^2)./(1+k_e.^2)); % The negative sign is important

% Storing output

    if any(strcmp(outputs,'input')) || any(strcmp(outputs,'all'))
        input.Bmag=B; input.description.Bmag = {'IGRF Field |B|','[nT]'};
        input.Tn=Tn; input.description.Tn = {'MSIS Neutral Temperature','[K]'};
        input.Ti=Ti; input.description.Ti = {'IRI Ion Temperature','[K]'};
        input.Te=Te; input.description.Te = {'IRI Electron Temperature','[K]'};

        input.ionConcentration=C;
        input.description.ionConcentration{1} = {'[O+] Concentration','ratio'};
        input.description.ionConcentration{2} = {'[NO+] Concentration','ratio'};
        input.description.ionConcentration{3} = {'[O2+] Concentration','ratio'};
        input.description.ionConcentration{4} = {'[N+] Concentration','ratio'};
        input.description.ionConcentration{5} = {'[H+] Concentration','ratio'};
        input.description.ionConcentration{6} = {'[He+] Concentration','ratio'};

        input.neutralNumberDensity=Nn;
        input.description.neutralNumberDensity{1}={'O','#/m^3'};
        input.description.neutralNumberDensity{2}={'N_2','#/m^3'};
        input.description.neutralNumberDensity{3}={'O_2','#/m^3'};

        input.totalNeutralNumberDensity=msisData.totalNumberDensity;
        input.description.totalNeutralNumberDensity={'Total neutral density','#/m^3'};

        input.electronDensity = ne;
        input.description.electronDensity = {'Ne','#/m^3'};
    end

    if any(strcmp(outputs,'all'))
        data.ionGyroFrequency = w_i; %rad/s
        data.electronGyroFrequency = w_e; %rad/s
        data.electronCollisionFrequency = v_e; %Hz
        data.electronNeutralCollisionFrequency=v_en;
        data.electronIonCollisionFrequency=v_ei; %Hz
        data.ionNeutralCollisionFrequency = v_i; %Hz
        data.electronMobility = k_e;
        data.ionMobility = k_i;
    end

    data.pedersonConductivity = sigma_P; % Pedersen Conductivity S/m
    data.hallConductivity = sigma_H; % Hall Conductivity Conductivity S/m
    data.altitude = alt1; % Altitude [km]

    if setPlotConductivity
    plot_conductivity(data, mode);
    end

end

function plot_conductivity(data, mode)

resize_figure(figure,100,100);

plot(data.pedersonConductivity,data.altitude,'k');
hold on;
plot(data.hallConductivity,data.altitude,'r');
set(gca,'XScale','log');
ylabel('Km');
xlabel('S/m');

if mode == 1
    modeStr = 'IRI';
else
    modeStr = '';
end
legend([modeStr,' \sigma_P'],[modeStr,' \sigma_H']);

end
