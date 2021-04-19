%%
clear all;

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

% Time array to determine how the pendulum's length changes
timeArray = linspace(0,11,1000); 

% Time range [t_initial t_final] through which the equation of motion of the pendulum will be solved.
tspan = [0 11]; 

% Starting from 1, the minimum length the pendulum attains at the end of
% the timeArray. 
Lmin = 0.5;

% Initial [theta, dtheta_dt] values to be input to the runge-kutta solver
theta_i = [0, 1];

% Acceleration due to gravity;
g = 9.81; 

% Resolution of the output solution
options = odeset('RelTol',1e-7);

% Runge-kutta solver
[t, theta] = ode45(@(t,theta) simplePendulum(t, theta, timeArray, Lmin,  g), tspan, theta_i, options);

% Calculating the length, and dL/dt of the pendulum at the each time point of
% the output solution
[L,dL] = pendulumLength(t,Lmin,timeArray);

% Conjugate momenta
P = L.^2.*theta(:,2); % Conjugate momenta

% Total energy
E = 0.5.*(L.^2.*theta(:,2).^2 - dL.^2 + 2.*g.*L.*(1-cos(theta(:,1)))); 
% Plotting
hFig = figure;
resize_figure(hFig,80,250);

mov(1:length(t)) = struct('cdata', [], 'colormap', []);

for iTime = 1:1:length(t)
   plotPendulum(hFig,t(1:iTime),theta(1:iTime,:),L(1:iTime),P(1:iTime),E(1:iTime))
   mov(iTime) = getframe(gcf);
end

% Movie
% v=VideoWriter('G:\My Drive\Research\Conferences\RICE Seminar 2021\violatingAdiabaticInvariance.avi','Uncompressed AVI');
% open(v)
% writeVideo(v,mov);
% close(v)

%%

function dTheta = simplePendulum(t, theta, tspan1, Lmin, g)
    [l,dl] = pendulumLength(t,Lmin,tspan1);
    dTheta = [theta(2); -2.*dl.*(theta(2)./l) - (g./l).*sin(theta(1))];
end

function [L,dL] = pendulumLength(t,Lmin,tspan1)
L1 = linspace(1,Lmin,length(tspan1));
% L1(round(length(tspan1)/2)+1:end) = linspace(L1(round(length(tspan1)/2)+1),L1(round(length(tspan1)/2)+1),round(length(tspan1)/2));
% L1 = smooth(L1,400);
% L1 = 1-Lmin+Lmin.*cos(tspan1).^2;
% L = interp1(tspan1,linspace(1,Lmin,length(tspan1)),t);
L = interp1(tspan1,L1,t);
dL1 = diff(L1);
dL1(end+1) = dL1(end);
dL = interp1(tspan1,dL1,t);
% dL = interp1(tspan1,diff(linspace(1,Lmin,length(tspan1)+1)),t);
dL = dL./(mean(diff(tspan1)));
end


function plotPendulum(hFig,t,theta,L,P,E)
    
    p = panel(hFig);
    p.pack({{60}},3);
    p.marginright=20;
    % Plotting the pendulum in cartesian coordinates

    p(1,1).select();
    plot(L.*sin(theta(:,1)),-L.*cos(theta(:,1)),'-r');
    line([0,L(end,1).*sin(theta(end,1))],[0,-L(end,1).*cos(theta(end,1))],'Color','k');
    hold on;
    plot(L(end).*sin(theta(end,1)),-L(end).*cos(theta(end,1)),'.b');
    ylim([-1,0]);
    xlim([-1,1]);

    % Plotting phase space
    p(1,2).select();
    plot(P,theta(:,1),'r');
    hold on;
    plot(P(end),theta(end,1),'.b');
    xlabel('P_{\theta}');
    ylabel('\theta');
    ylim([-0.6,+0.6]);
    xlim([-1,+1]);
    % Plotting total energy of the pendulum
    p(1,3).select();

    yyaxis left
    plot(t,E,'k');
    hold on;
    plot(t(end),E(end),'.b');
    ylim([0,1]);
    xlim([0,11]);
    ylabel('E(t): Total Energy');
    % hold on; 
    % plot(t,g.*L.*(1-cos(theta(:,1))));
    % hold on;
    % plot(t,0.5.*(L.^2.*theta(:,2).^2 - dL.^2));

    yyaxis right
    plot(t,(L.^0.5).*E,'r');
    plot(t(end),(L(end).^0.5).*E(end),'.b');
    ylim([0,1]);
    xlim([0,11]);
    ylabel('E(t)l(t)^{1/2}: Phase space area');
    xlabel('Time (t)');

end

% function [L,dL] = pendulumLength(t,Lmin,tspan1)
%     a=-0.05;
%     
%     L = 1 + a.*t;
%     dL = a;
% end