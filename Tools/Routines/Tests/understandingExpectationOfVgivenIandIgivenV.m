%% Finding why E(V|I) is different from E(I|V)
clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

% Case 1: noisy V and I data (gaussian)

SNR = 10; % In linear unitsl and is equal to (rms(Vs)./rms(Ns)).^2
SNRArray = logspace(0.5,-0.5,1000).*SNR;
signals = linspace(0,20,1000); % DC signal voltage for 500 signals
nRepeat = 500;
Vs = repmat(signals,nRepeat,1); %500 instances of each signal
for i = 1:length(signals)
    Ns(:,i) = wgn(1,nRepeat,signals(i).^2,'linear')./sqrt(0.05.*SNRArray(i));
end

for i = 1:length(signals)
    Ns2(:,i) = wgn(1,nRepeat,signals(i).^2,'linear')./sqrt(0.5.*SNRArray(i));
end

yLimRange = [0,6000];
xLimRange = [0,40];

% figure; plot(rms(Vs).^2./rms(Ns).^2);
V = Vs + Ns;
% V = Vs;

% R = linspace(20,1,500); %Constant
R = 1/300;

Isignal = (signals./R) ;
Is = (Vs./R);
INs = (Ns./R);
% I = (V./R);
I = Is + Ns2./R;
% I = Is + INs;

vRange = [0,max(V(:))];
iRange = [0,max(I(:))];

vBinsTemp = linspace(0,max(V(:)),50);
vBins = vBinsTemp(1:end-1) + diff(vBinsTemp);

for i = 1:1:length(vBinsTemp)-1
    EIgV(i) = mean(I(V>=vBinsTemp(i) & V<vBinsTemp(i+1) & V >=0 )); 
end


iBinsTemp = linspace(0,max(I(:)),20);
iBins = iBinsTemp(1:end-1) + diff(iBinsTemp);

for i = 1:1:length(iBinsTemp)-1
    EVgI(i) = mean(V(I>=iBinsTemp(i) & I<iBinsTemp(i+1) & I >=0 )); 
end

hFig = figure; 
resize_figure(hFig,80,250);

    p = panel(hFig);
    p.pack({{60}},3);
    p.marginright=20;

    p(1,1).select();
    scatter(V(:),I(:),10,'filled','MarkerFaceAlpha',0.01);
    r = corrcoef(V(:),I(:));
    text(0.8,0.1,['r_0_0 = ',num2str(r(2,1),3)],'Units','normalized');
    xlabel('E');
    ylabel('J');
    hold on;
    plot(signals,Isignal);
    xlim(xLimRange);
    ylim(yLimRange);
    legend('Signal + Noise','Signal');
    
    p(1,2).select();
    yyaxis left
    plot(vBins,EIgV); 
    ylabel('\langle J|E \rangle');
    xlim(xLimRange);
    ylim(yLimRange);
%     xlim(vRange);
%     ylim(iRange);
    hold on; 

    yyaxis right
    plot(signals,Isignal);
    ylabel('J');
%     ylim(iRange);
    ylim(yLimRange);
    xlabel('E');

    legend('Signal + Noise', 'Signal');

    p(1,3).select(); 
    plot(EVgI,iBins); 
    hold on; 
    plot(signals,Isignal);
    xlim(xLimRange);
    ylim(yLimRange);
%     xlim(vRange);
%     ylim(iRange);
    xlabel('\langle E|J \rangle');
    ylabel('J');
    legend('Signal + Noise', 'Signal','Location','northeast');
