% Plotting lambda to verify with Semeter 2005

E   =[5000, 20000, 30000]; % eV
hi  =[0:0.01:3.0]; % s/R

pad =2;
%   pad     - type of pitch angle distribution 1 for || 2 for
%             isotropic 

figure;

for i=1:1:3
    Y = lambda(hi,E(i),pad);
    hold on;
    plot(hi,Y);
end;
hold off;
