function [dE] = get_error_in_energyFlux(dq, A, energyBin, energyFlux, time, Gamma)

%% get_error_in_energyFlux.m Calculates the covariance matrix of energy flux estimates
%-------------------------------------------------------------------------
%  Input 
%-------
%  dq         : Error in the form of std. dev. of production rates [m^-3 s^-1]
%  A          : Matrix of production rates per unit number flux [eV m-1]
%  energyBin  : Electron energy bin values [eV]
%  energyFlux : Differential number flux [eV-1 m-2 s-1 ]
%  time       : Time vector in [matlab units]
%  Gamma      : A constant that weighs the importance of two error terms
%               in -S/Gamma &  0.5 [e]'[C_d^-1][e]
%               If Gamma == inf - maximum error (Default)
%               If Gamma == 0   - least possible error
%-------------------------------------------------------------------------
% Output
%--------
% dE          : Std. deviation in energy flux - worst case [eV m-2 s-1 eV-1]
%%
%----------------------------------------------------------------------------
% Modified: 22nd Sep 2016 
% Created : 22nd Sep 2016
% Author  : Nithin Sivadas
% Ref     : D. L. Hysell 2007
%----------------------------------------------------------------------------
%%



    if nargin < 7
        Gamma = -1;
    end;

    numFlux = energy_to_num(energyFlux, time, energyBin); 

    for itime = 1:1:length(time)

        thisFlux = numFlux(:,itime);
        
        Cq = diag(dq(:,i).^2);

        part1 = A'*inv(Cq)*A; 
        if Gamma ~= -1
            part2 = inv(Gamma*diag(thisFlux));            
        else
            part2 = 0; % Simulating Gamma -> inf
        end;
        Cflux_inv = part1 + part2;        
        Cflux = inv(Cflux_inv);
        
        [X,D]=eig(Cflux,'matrix');
           
        dnumFlux=diag(D);        
        denergyFlux(i,:)=num_to_energy(dnumFlux,1,energyBin); 
     
     end;
     
    dE=denergyFlux.^0.5;

    [isThereNAN, totalNAN] = check_nan(dE);
    
end

