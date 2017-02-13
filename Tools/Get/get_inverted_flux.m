function [data] = get_inverted_flux( q, time, alt, energyBin, A, guessFlux, DateNumBeg, DateNumEnd )
%% get_inverted_flux.m Estimate energy spectra by inverting q(z) 
% Using altitude profiles of production rates and the maximum entropy inversion 
% method to estimate the energy spectra.
%--------------------------------------------------------------------------
% Input:
%--------
% q         : Production rates [m^-3 s^-1]   [nHxnT]
% alt       : Altitude [km]                  [nHx1]
% time      : Time [matlab units]            [1xnT]
% energyBin : Energy bin values [eV]         [nEx1]
% DateNumBeg: Initial time 
% DateNumEnd: Final time 
%--------------------------------------------------------------------------
% Output:
%--------
% data.flux          : Differential number flux [m-2 s^-1 eV^-1] (Isotropic flux)
% data.energyFlux    : Differential energy flux [eV m-2 sr-1 s-1 eV-1]
% data.chi2          : Measures the deviation between the measured q and inverted q (q-A*flux)/qSigma
% data.qInverted     : Production rate derived from the estimated flux - A*flux [m-3 s-1]
% data.maxIter       : The maximum number of iterations before convergence

% data.energyBin     : The energy values of each energy spectral bin [eV]
% data.time          : Time vector [s]
% data.alt           : Altitude [km]

% data.A             : Energy deposition matrix [nH x nE] [eV m^-1]
% data.qInput        : The production rate which was input to this function [m^-3 s^-1]

%----------------------------------------------------------------------------
% Modified: 25th Sep 2016 
% Created : 25th Sep 2016
% Author  : Nithin Sivadas
% Ref     : 
% 1. Sergienko and Ivanov 1993 "A new approach to calculate the excitation of atmospheric gases by auroral electron impact"
% 2. Semeter and Kamalabadi 2005 
% 3. Rees 1989
% 4. Wedlund et al "Electron Energy Spectra and Auroral Arcs" JGR 2013
%----------------------------------------------------------------------------


    %% Enetering default values 
    if nargin<8
        DateNumEnd=max(time);
    end

    if nargin<7
        DateNumBeg=min(time);
    end;

    if nargin<6
        kb          = 1.38*10^-23;
        eV          = 1.602E-19;
        T_w         = 10*5000*eV/kb;
%         guessFlux = get_kappa_j([10000,T_w,160, 10^12], energyBin); %[m^-2 s^-1]
        guessFlux = (10^13)*get_kappa_E(energyBin,10000,T_w,160); %[m^-2 s^-1]
    end;
    
    if nargin<5
            % Generating the A matrix
            A =get_energy_dep_matrix(alt,energyBin,...
                65,-147.5,datenum([2008 03 26 10 00 00])); %[m^-1 eV] 
    end;

%     [timeMaxNo] = find_time(time, DateNumEnd);
%     [timeMinNo] = find_time(time, DateNumBeg);

    %% Preparing the data for inversion
    % Cropping the time array and density matrix to that prescribed by the user
    [q, time] = crop_time(q, time, DateNumBeg, DateNumEnd);

    % Removing nans
    q = interp_nans(q);

    %  Removing negative values of production rates if any
    totalNegValuesInQ=sum(abs(q(q<0)./q(q<0)));
    q(q<0)=0;


    %% Implementing the Maximum Entropy Method

    % Setting up the inputs for the mem_solve function
    beta    = 15;
    itime    = 1:1:length(time);
    qSigma   =(var(q)).^0.5;
    noIter = 5000;

    % Generating the weighting function
    W=ones(size(A,1),1);
    W = gaussmf((1:1:size(W,1)),[50 100])';
    W=W./sum(W);

    % Maximum entropy inversion
    for thisTime=itime
            qThis=q(:,thisTime);
            [fluxNew(:,thisTime),qNew(:,thisTime),chi2(thisTime),maxIter(thisTime)] = ...
            mem_solve(qThis, A, beta, guessFlux, (qSigma(thisTime)), noIter, W);

            display([num2str((thisTime)/max(itime)*100), ' %',num2str(chi2(thisTime))]);
    end;

    % Generating the output structure

    data.flux = fluxNew; % Differential number flux [m-2 s^-1 eV^-1] (Isotropic flux)
    data.energyFlux = num_to_energy(fluxNew, time, energyBin); % Differential energy flux [eV m-2 sr-1 s-1 eV-1]
    data.chi2 = chi2; % Measures the deviation between the measured q and inverted q (q-A*flux)/qSigma
    data.qInverted = qNew; % Production rate derived from the estimated flux - A*flux [m-3 s-1]
    data.maxIter = maxIter; % The maximum number of iterations before convergence

    data.energyBin = energyBin; % The energy values of each energy spectral bin [eV]
    data.time = time; % Time vector [s]
    data.alt  = alt; % Altitude [km]

    data.A = A; % energy deposition matrix [H x E] [eV m^-1]
    data.qInput = q; % the production rate which was input to this function [m^-3 s^-1]

end

