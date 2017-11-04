function [ outputArguments ] = get_conductivity_v1( alt, electronDensity, latitude, longitude, time )
%get_conductivity_v1.m returns the hall and pederson conductivity given an 
%input electron density value, along with some other output parameters.
% Date Modified: 27 Sep 2017
%--------------------------------------------------------------------------
    % Initilization
	if nargin<5
	    time        = datenum([2008 03 26 10 00 00]);
	end
	if nargin<4
	    longitude   = -147.5; % Degrees. 
	end
	if nargin<3
	    latitude    = 65; % Degrees.
    end
       
    coordinateSystem = 'geodetic';
    curlDir='C:\Users\Nithin\Documents\GitHub\energy-height-conversion\Tools\External Tools';
    
    % MSIS
    alt1=linspace(min(alt),max(alt),length(alt));
	[D,T] = msis(time, latitude, longitude, alt1);
    Nn(:,1)     = D(:,1); %O
	Nn(:,2)     = D(:,2); %N2
	Nn(:,3)     = D(:,3); %O2
    
    N = sum(Nn,2); % Total density
    
    % IRI2016
    IRI=iri2012(time, latitude, longitude, alt1, true, coordinateSystem,curlDir, [], [], [], [], [], [], [], [], [], [] ,[] ,[], [], [], 'RBV10/TTS03' );
    %iri2012 function modified to use iri2016

    Tn = T(:,1); % Neutral temperature from MSIS
    Ti = IRI(:,4);
    Te = Ti;
    
    % Concentration of ion species
    C(:,1) = IRI(:,6)/100.000;  %OI
    C(:,2) = IRI(:,10)/100.000; %NOI
    C(:,3) = IRI(:,9)/100.000;  %O2I
    C(:,4) = IRI(:,12)/100.000; %N+
    C(:,5) = IRI(:,7)/100.000;  %H+
%   C(:,6) = IRI(:,8)/100.000;  %He+

% Correcting IRI [The discrepancy in the sum of the ion concentrations]
     C(C<0)=0;
     C_Sum = sum(C,2);
     C(:,1) = C(:,1)./C_Sum;
     C(:,2) = C(:,2)./C_Sum;
     C(:,3) = C(:,3)./C_Sum;
     C(:,4) = C(:,4)./C_Sum;
     C(:,5) = C(:,5)./C_Sum;

    
% B field
    [Bx, By, Bz] = igrf(time, latitude, longitude, alt1, coordinateSystem);
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

% Reduced mass
    m_e = 9.11*10^-31; % Electron mass
    m_p = 1.673*10^-27; % Proton mass
    Z_n = [16; 28; 32];
    m_n = Z_n*m_p; % Mass of neurtrals[kg]
    Z_i = [16; 30; 32; 14]; % Ion mass in proton mass units
    m_i = Z_i*m_p;
    q_i = [1;1;1; 1]; % Ion charge
    
    for thisIon=1:length(m_i)
        for thisNeutral=1:length(m_n)
            u_in(thisIon,thisNeutral) = m_i(thisIon).*m_n(thisNeutral)./...
                (m_i(thisIon)+m_n(thisNeutral));
        end
    end
    
    
% Gyro frequencies
    q = 1.6*10^-19;
    k = 1.28*10^-23;
    w_e  = q.*B./m_e;
    w_p = q.*B./m_p;
    
    for i = 1:1:length(m_i)
        w_i(:,i)=q_i(i)*w_p/Z_i(i);
    end;
    
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
%     ne = electronDensity';
    ne = IRI(:,1);
    ne(ne<=0)=10^6;

% Electron-ion collision frequency
%   v_ei = ((ne.*(10^-6)).*Te.^-1.5).*(34 + 4.18.*log(((Te.^3)./(ne.*(10^-6))))); % Michael C Kelley
    v_ei = ((ne.*(10^-6)).*Te.^-1.5).*(59 + 4.18.*log10(((Te.^3).*(ne)).^3));
    v_ei_1=  54.5*(ne.*(10^-6).*Te.^-1.5).*(C(:,1).*q_i(1).^2 + C(:,2).*q_i(2).^2 + C(:,3).*q_i(3).^2);
    
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
    sigma_H = (q*ne./B).*abs((sigma_H1 - (k_e.^2)./(1+k_e.^2)));

% Storing output
    data.B=B; data.Tn=Tn; data.Ti=Ti; data.C=C; data.Nn=Nn; data.N=N;
    data.w_i = w_i; data.w_e = w_e; data.v_e = v_e; data.v_in = v_in;
    data.v_i = v_i;
    data.sigma_P = sigma_P; data.sigma_H = sigma_H; data.alt = alt1; 
    data.C_in=C_in; 
    data.k_e = k_e;
    data.k_i = k_i;
    data.ne = ne;
    data.v_en=v_en;
    data.v_ei=v_ei; data.v_ei_1=v_ei_1;
    data.Te=Te;
%     data.E = E;
%     data.v_in_approx=v_in_approx;
    outputArguments = data;
    % Mass of ion species
    
end

