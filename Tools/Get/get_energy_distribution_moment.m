function [ mean, median ] = get_energy_distribution_moment( flux, energyBin) 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Mean

dE = diff(energyBin);
dE = [dE; dE(end)];
P_x = flux./sum(flux.*dE);
    
% Calculating mean of energy distribution
x = energyBin;
mean = sum(P_x.*x.*dE);

% Calculating median (in progress)
for thisE = 1:1:length(energyBin)
    if thisE==1
        C_x(thisE) = P_x(thisE).*dE(thisE);
    else
        C_x(thisE)= P_x(thisE).*dE(thisE)+C_x(thisE-1);
    end;
end;
median = interp1(C_x,energyBin, 0.5);
        
end

