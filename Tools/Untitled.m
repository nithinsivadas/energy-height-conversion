%% Finding why E(V|I) is different from E(I|V)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

% Case 1: noisy V and I data (gaussian)
Vs = linspace(1,10,100); 
V = Vs ;

R = 10; %Constant

Is = (Vs./R);
I = (V./R) + wgn(100,100,0.1);

figure; scatter(V(:),I(:),10,'filled','MarkerFaceAlpha',0.4);
r = corrcoef(V(:),I(:));
text(0.9,0.1,['r_0_0 = ',num2str(r(2,1))],'Units','normalized');

vBinsTemp = linspace(min(V(:)),max(V(:)),10);
vBins = vBinsTemp(1:end-1) + diff(vBinsTemp);

for i = 1:1:length(vBinsTemp)-1
    EIgV(i) = mean(I(V>=vBinsTemp(i) & V<vBinsTemp(i+1) )); 
end


iBinsTemp = linspace(min(I(:)),max(I(:)),10);
iBins = iBinsTemp(1:end-1) + diff(iBinsTemp);

for i = 1:1:length(iBinsTemp)-1
    EVgI(i) = mean(V(I>=iBinsTemp(i) & I<iBinsTemp(i+1) )); 
end

figure; 

yyaxis left
plot(vBins,EIgV); 
ylabel('\langle J|E \rangle');
xlim([min(V(:)),max(V(:))]);
ylim([min(I(:)),max(I(:))]);
hold on; 

yyaxis right
plot(Vs,Is);
ylabel('J');
ylim([min(I(:)),max(I(:))]);
xlabel('E');

legend('Singal + Noise', 'Signal');

figure; 
plot(EVgI,iBins); 
hold on; 
plot(Vs,Is);
xlim([min(V(:)),max(V(:))]);
ylim([min(I(:)),max(I(:))]);
xlabel('\langle E|J \rangle');
ylabel('J');
legend('Singal + Noise', 'Signal');
