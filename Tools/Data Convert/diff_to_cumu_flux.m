function [ cumuEnergyFlux ] = diff_to_cumu_flux( energyFlux, energyBin )
%diff_to_cumu_flux Converts differential energy flux into cumulative energy
% flux
% Input
% energyFlux : differential energy flux [eV m-2 sr-1 s-1 ev-1] [nExnT]
% energyBin  : energy bin values [eV]
%
% Output
% cumuEnergyFlux : cumulative energy flux [eV m-2 sr-1 s-1] [nExnT]
%

%%
%----------------------------------------------------------------------------
% Modified: 26th Sep 2016 
% Created : 26th Sep 2016
% Author  : Nithin Sivadas
% Ref     : 
%----------------------------------------------------------------------------
%%



noEnergyBin = length(energyBin);
noTime = size(energyFlux,2);
cumuEnergyFlux = zeros(noEnergyBin, noTime);

for thisEnergy=1:1:noEnergyBin
    cumuEnergyFlux(thisEnergy,:) =...
        sum(energyFlux(1:thisEnergy,:).*...
        repmat(energyBin(1:thisEnergy),1,noTime),1);
end;

end

