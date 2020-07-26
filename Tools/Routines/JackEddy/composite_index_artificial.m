
%% Canonical Correlation Analysis
inputFileStr='G:\My Drive\Nithin\Job Search\JEPostDoc\Data\2018ja025430-sup-0001-data_set_si-s01.txt';
Table = readtable(inputFileStr);
Table.Properties.VariableNames...
    ([5,6,7,8,9,10,11,12,13,14,15,16,17,18])...
    = {'S1','E1','E1_E1pred','Y1_AL','Y2_AU','Y3_PCI',...
    'Y4_Kp','Y5_am','Y6_Dst','Y7_mPe','Y8_mPi','Y9_Pips',...
    'Y10_Fe130','Y11_Ssince'};
time = datetime(Table.year,1,Table.day,Table.hr,0,0,0);
G = [Table.Y1_AL,Table.Y2_AU,Table.Y4_Kp,Table.Y_13__log_vsw_,Table.Y_14__log_nsw_,Table.Y_17__log_Bmag_];
P = [Table.Y7_mPe,Table.Y8_mPi];

[A,B,r,U,V,stats] = canoncorr(G,P);
[G1,G2,P1,P2] = composite_index(G,P,A,B);
%%
% Manual indices
% Get higher temporal resolution 
% omnih5= 'G:\My Drive\Research\Projects\Data\omniManual.h5';
startTimeStr = '01 Jan 9999 08:00';
endTimeStr = '01 Jan 9999 18:00';
omniManual.time = datenum(startTimeStr):1/(24*60):datenum(endTimeStr);
nTime = length(omniManual.time);
growthPhaseTime = 0.5; % Hours
recoveryPhaseTime = 0.5; % 2 hours
amplitude = 0; %nT
onsetTime = 3; % hours
% omniManual.AL = double_gaussian((omniManual.time-omniManual.time(1))*24,...
%     amplitude,onsetTime,growthPhaseTime,recoveryPhaseTime) ;
% omniManual.AU = double_gaussian((omniManual.time-omniManual.time(1))*24,...
%    0,onsetTime,growthPhaseTime+0.5,recoveryPhaseTime) ;

% omniManual.Kp = linspace(2.098,9,nTime);
% omniManual.AU = linspace(0,300,nTime);
omniManual.AL = linspace(-63.5,-1000,nTime);
% omniManual.Vsw = linspace(250,800,nTime);
% omniManual.nsw = linspace(1,20,nTime);
% omniManual.Bmag = linspace(1,20,nTime);

% omniManual.AL = -63.5.*ones(1,nTime);
omniManual.AU = 52.9.*ones(1,nTime);
omniManual.Kp = 2.098.*ones(1,nTime);
omniManual.Vsw = 429.*ones(1,nTime);
omniManual.nsw = 5.3*ones(1,nTime);
omniManual.Bmag = 5.6.*ones(1,nTime);

omniManual.AL(omniManual.AL==99999)=nan;
omniManual.AU(omniManual.AU==99999)=nan;
input.AL = (log10(1+abs(omniManual.AL'))-1.809)./0.5609;
input.AU = (log10(1+abs(omniManual.AU'))-1.732)./0.3962;
input.Kp = (omniManual.Kp'-2.098)./1.391;
input.Vsw = (log10(omniManual.Vsw')-2.633)./0.09717;
input.nsw = (log10(omniManual.nsw')-0.7315)./0.2983;
input.Bmag =(log10(omniManual.Bmag')-0.7527)./0.1931;

index = ~isnan(input.AL) & ~isnan(input.AU) & ~isnan(input.Kp) &  ~isnan(input.Vsw) & ~isnan(input.nsw) & ~isnan(input.Bmag);

G_new = [input.AL(index), input.AU(index), input.Kp(index), input.Vsw(index), input.nsw(index), input.Bmag(index)];

G_i = G_new*A;

time_i = omniManual.time(index);

g1 = polyfit(G1,P1,1);
g2 = polyfit(G2,P2,1);
P_i = [G_i(:,1).*g1(1)+g1(2), G_i(:,2).*g2(1)+g2(2)];
mP = (P_i*inv(B));
mP(:,1) = (10.^((mP(:,1)*0.3985)+1.291))-0.01;  %electron
mP(:,2) = (10.^((mP(:,2)*0.2626)+0.3381))-0.01;  %ions

%% Plotting CCA Example

tIndx1 = find_time(omniManual.time,startTimeStr);
tIndx2 = find_time(omniManual.time,endTimeStr);
GtIndx1 = find_time(omniManual.time,startTimeStr);
GtIndx2 = find_time(omniManual.time,endTimeStr);
trange = tIndx1:tIndx2;
Gtrange = GtIndx1:GtIndx2; 
G1smooth = G_i(Gtrange,1);
G2smooth = G_i(Gtrange,2);

%
lw=1.5;
figureHandle2 = figure;
resize_figure(figureHandle2,100,400);
p=panel();
p.pack({{70}},{{75} {78} {78} {78}});
p.margintop=15;
p.marginleft=18;
p(1,2).pack(3,1);
p(1,2).marginleft=10;
p(1,3).marginleft=15;
p(1,4).marginleft=15;

p(1,3).pack(3,1);
p(1,4).pack(3,1);
p.fontsize=14;
p.fontname='Times';
p.select('all');

% Panel 1
p(1,1).select();
colormap(inferno);
p(1,1).marginright=30;
p(1,2).de.margin = 4;

% histogram2(G_i(:,1),G_i(:,2),'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none'); 
[H, Xedges, Yedges]  = histcounts2(G1,G2,128,'normalization','countdensity');
[X,Y] = meshgrid(Xedges(1:end-1),Yedges(1:end-1));
%[500,300,100,10]
%[0.15,0.1,0.05,0.01]
[C,h]=contour(X,Y,smoothdata(H')/10,[1000,500,100,10],'LineColor',[0.5,0.5,0.5]);
hcl=clabel(C,h);
% colormap(gca,'Gray');

set(gca,'XTick',[-2.5,0,2.5],'YTick',[-2.5,0,2.5]);
xlim([-3,3]);
ylim([-3,3]);
x1h = xlabel('G*_{(1)}');
x1h.Position(2) = x1h.Position(2)+0.1;
y1h=ylabel('G*_{(2)}');
y1h.Position(1) = y1h.Position(1)+0.2;

hold on;
% G1smooth = smoothdata(G_i(Gtrange,1));
% G2smooth = smoothdata(G_i(Gtrange,2));
scatter(G1smooth ,G2smooth,8,time_i(Gtrange),'.');
c11=colorbar_thin('YLabel','Time');
c11.Ticks=[time_i(Gtrange(1)),time_i(Gtrange(end))]; c11.TickLabels={datestr(c11.Ticks(1),'HH:MM'),datestr(c11.Ticks(2),'HH:MM')};
c11.Limits=[time_i(Gtrange(1)),time_i(Gtrange(end))];
caxis(c11.Limits);
colormap(gca,jet);
text(0.5,1.13,'A) Composite-index Space','Units','normalized','FontSize',18,'FontName','Times','horizontalAlignment','center');

n = length(trange);
cdata = [uint8(jet(n)*255) uint8(ones(n,1))].';

% Panel 2
p(1,2,1,1).select();
p1=plot(datetime(omniManual.time(trange),'ConvertFrom','datenum'),omniManual.AL(trange),'r','LineWidth',lw);
drawnow
set(p1.Edge, 'ColorBinding', 'interpolated', 'ColorData', cdata);    
set(gca,'XTickLabel','');
text(0.8,0.3,'AL','Units','normalized','FontSize',18,'FontName','Times');
text(0.015,0.18,'[nT]','Units','normalized','FontSize',12,'FontName','Times');
text(0.5,1.4,'B) Geomagnetic Indices','Units','normalized','FontSize',18,'FontName','Times','horizontalAlignment','center');

p(1,2,2,1).select();
p2=plot(datetime(omniManual.time(trange),'ConvertFrom','datenum'),omniManual.AU(trange),'r','LineWidth',lw);
% drawnow
% set(p2.Edge, 'ColorBinding', 'interpolated', 'ColorData', cdata);
set(gca,'XTickLabel','')
text(0.8,0.3,'AU','Units','normalized','FontSize',18,'FontName','Times');
text(0.015,0.18,'[nT]','Units','normalized','FontSize',12,'FontName','Times');

p(1,3).de.margin = 4;
p(1,2,3,1).select();
p3=plot(datetime(omniManual.time(trange),'ConvertFrom','datenum'),omniManual.Kp(trange),'r','LineWidth',lw);
% drawnow
% set(p3.Edge, 'ColorBinding', 'interpolated', 'ColorData', cdata);
set(gca,'YTick',[4,6,8,10,12],'YTickLabel',{'4','','8','','12'});
text(0.8,0.3,'Kp','Units','normalized','FontSize',18,'FontName','Times');
text(0.015,0.18,'[a.u]','Units','normalized','FontSize',12,'FontName','Times');

% Panel 3
p(1,3,1,1).select();
omniManual.Vsw=interp_nans(omniManual.Vsw')';
p4 = plot(datetime(omniManual.time(trange),'ConvertFrom','datenum'),omniManual.Vsw(trange),'r','LineWidth',lw);
% drawnow
% set(p4.Edge, 'ColorBinding', 'interpolated', 'ColorData', cdata);
set(gca,'XTickLabel','');
text(0.8,0.8,'V_{sw}','Units','normalized','FontSize',18,'FontName','Times');
text(0.015,0.18,'[ms^{-1}]','Units','normalized','FontSize',12,'FontName','Times');
text(0.5,1.4,'C) Solarwind Parameters','Units','normalized','FontSize',18,'FontName','Times','horizontalAlignment','center');

p(1,3,2,1).select();
omniManual.nsw=interp_nans(omniManual.nsw')';
p5=plot(datetime(omniManual.time(trange),'ConvertFrom','datenum'),omniManual.nsw(trange),'r','LineWidth',lw);
% drawnow
% set(p5.Edge, 'ColorBinding', 'interpolated', 'ColorData', cdata);
set(gca,'XTickLabel','');
text(0.8,0.8,'n_{sw}','Units','normalized','FontSize',18,'FontName','Times');
text(0.015,0.18,'[cm^-^3]','Units','normalized','FontSize',12,'FontName','Times');

p(1,3,3,1).select();
omniManual.Bmag=interp_nans(omniManual.Bmag')';
p6=plot(datetime(omniManual.time(trange),'ConvertFrom','datenum'),omniManual.Bmag(trange),'r','LineWidth',lw);
% drawnow
% set(p6.Edge, 'ColorBinding', 'interpolated', 'ColorData', cdata);
text(0.8,0.8,'B_{mag}','Units','normalized','FontSize',18,'FontName','Times');
text(0.015,0.18,'[nT]','Units','normalized','FontSize',12,'FontName','Times');

p(1,4).de.margin = 4;
p(1,4,1,1).select();
p7 = plot(datetime(time_i(Gtrange),'ConvertFrom','datenum'),mP(Gtrange,1),'k','LineWidth',lw);
hold on;
p8 = plot(datetime(time_i(Gtrange),'ConvertFrom','datenum'),mP(Gtrange,2),'r','LineWidth',lw);
set(gca,'XTickLabel','');
text(0.8,0.8,'mPe','Units','normalized','FontSize',18,'FontName','Times','Color','k');
text(0.6,0.8,'mPi','Units','normalized','FontSize',18,'FontName','Times','Color','r');
text(0.015,0.18,'[GW]','Units','normalized','FontSize',12,'FontName','Times');
% ylim([-1.5,+1.5]);

GComplex = G1smooth + 1i*G2smooth;
GAmp = abs(GComplex);
GPh = angle(GComplex);
Gn = length(Gtrange);
Gcdata = [uint8(jet(Gn)*255) uint8(ones(Gn,1))].';

p(1,4,2,1).select();
p9 = plot(datetime(time_i(Gtrange),'ConvertFrom','datenum'),GAmp,'k','LineWidth',lw);
drawnow
set(p9.Edge, 'ColorBinding', 'interpolated', 'ColorData', Gcdata);
set(gca,'XTickLabel','');
text(0.8,0.8,'|G^*_(_1_)+iG^*_(_2_)|','Units','normalized','FontSize',12,'FontName','Times','Color','k');
text(0.015,0.18,'[a.u.]','Units','normalized','FontSize',12,'FontName','Times');

p(1,4,3,1).select();
scatter(datetime(time_i(Gtrange),'ConvertFrom','datenum'),wrapTo2Pi(GPh)./pi,[],time_i(Gtrange),'.');
ylim([0,2]);
text(0.8,0.8,'arg(G^*_(_1_)+iG^*_(_2_))','Units','normalized','FontSize',12,'FontName','Times','Color','k');
text(0.015,0.18,'[\pi]','Units','normalized','FontSize',12,'FontName','Times');
colormap(gca,jet);

%%
figure;
theta.AL = (atan(A(1,2)./A(1,1)));
theta.AU = (atan(A(2,2)./A(2,1)));
theta.Kp = (atan(A(3,2)./A(3,1)));
theta.Vsw = (atan(A(4,2)./A(4,1)));
theta.nsw = (atan(A(5,2)./A(5,1)));
theta.Bmag = (atan(A(6,2)./A(6,1)));
theta.mPi = atan(-g2(1)*B(1,1))./(g1(1)*B(1,2));
theta.mPe = atan(-g2(1)*B(2,1))./(g1(1)*B(2,2));

%AL 
quiver(0,0,real(exp(1i*theta.AL)),imag(exp(1i*theta.AL)),'k');
text(cos(theta.AL),sin(theta.AL),[num2str(wrapTo360(rad2deg(theta.AL)),3),'° AL']);
ylim([-1.5,+1.5]);xlim([-1.5,+1.5]);
hold on;

%AU 
quiver(0,0,real(exp(1i*theta.AU)),imag(exp(1i*theta.AU)),'k');
text(cos(theta.AU),sin(theta.AU),[num2str(wrapTo360(rad2deg(theta.AU)),3),'° AU']);
hold on;

%Kp 
quiver(0,0,real(exp(1i*theta.Kp)),imag(exp(1i*theta.Kp)),'k');
text(cos(theta.Kp),sin(theta.Kp),[num2str(wrapTo360(rad2deg(theta.Kp)),3),'° Kp']);
hold on;

%Vsw 
quiver(0,0,real(exp(1i*theta.Vsw)),imag(exp(1i*theta.Vsw)),'b');
text(cos(theta.Vsw),sin(theta.Vsw)-0.1,[num2str(wrapTo360(rad2deg(theta.Vsw)),3),'° V_s_w']);
hold on;

%nsw 
quiver(0,0,real(exp(1i*theta.nsw)),imag(exp(1i*theta.nsw)),'b');
text(cos(theta.nsw)-0.3,sin(theta.nsw),[num2str(wrapTo360(rad2deg(theta.nsw)),3),'° n_s_w']);
hold on;

%Bmag 
quiver(0,0,real(exp(1i*theta.Bmag)),imag(exp(1i*theta.Bmag)),'b');
text(cos(theta.Bmag),sin(theta.Bmag),[num2str(wrapTo360(rad2deg(theta.Bmag)),3),'° B_m_a_g']);
hold on;

% mPi
quiver(0,0,real(exp(1i*theta.mPi)),imag(exp(1i*theta.mPi)),'r');
text(cos(theta.mPi),sin(theta.mPi),[num2str(wrapTo360(rad2deg(theta.mPi)),3),'° mP_i']);
hold on;

% mPe
quiver(0,0,real(exp(1i*theta.mPe)),imag(exp(1i*theta.mPe)),'r');
text(cos(theta.mPe)-0.1,sin(theta.mPe)-0.05,[num2str(wrapTo360(rad2deg(theta.mPe)),3),'° mP_e']);

xlabel('G^*_(_1_)');
ylabel('G^*_(_2_)');
title('Influence of geomagnetic indices on G1/G2 Trajectory');

%%
nGSize = 64;
[Gx,Gy] = meshgrid(linspace(-3,3,nGSize),linspace(-3,3,nGSize));
Px = Gx.*g1(1)+g1(2);
Py = Gy.*g2(1)+g2(2);
P_j(:,:,1) = Px;
P_j(:,:,2) = Py;
for i = 1:1:nGSize 
    mP_j(:,i,:)=squeeze(P_j(:,i,:))*inv(B);
end
mP_je = (10.^((squeeze(mP_j(:,:,1))*0.3985)+1.291))-0.01;  %electron
mP_ji = (10.^((squeeze(mP_j(:,:,2))*0.2626)+0.3381))-0.01;  %ions

h1=figure; 
resize_figure(h1,200,250);
colormap(inferno);
subplot(2,2,1);
r=pcolor(Gx,Gy,log10(mP_je)); set(r,'EdgeColor','none');
title('Electron Flux [GW]');
xlabel('G^*_(_1_)');
ylabel('G^*_(_2_)');
colorbar;
% hold on;
subplot(2,2,2);
[pex,pey] = gradient(mP_je);
compass(pex,pey);
hold on;
text(max(abs(pex(:))),min((pey(:))),num2str(wrapTo360(rad2deg(mean(angle(pex(:)+1i*pey(:))))),4));
% quiver(Gx,Gy,pex,pey,2);

subplot(2,2,3);
r2=pcolor(Gx,Gy,log10(mP_ji)); set(r2,'EdgeColor','none');
title('Ion Flux [GW]');
colorbar;
xlabel('G^*_(_1_)');
ylabel('G^*_(_2_)');

subplot(2,2,4);
[pix,piy] = gradient(mP_ji);
compass(pix,piy);
hold on;
text(max(abs(pix(:))),max(abs(piy(:))),num2str(wrapTo360(rad2deg(mean(angle(pix(:)+1i*piy(:))))),4));

function [G1,G2,P1,P2] = composite_index(G,P,A,B)
    Temp1 = G*A;
    Temp2 = P*B;
    G1 = Temp1(:,1);
    G2 = Temp1(:,2);
    P1 = Temp2(:,1);
    P2 = Temp2(:,2);
end

function y = normalize(x)
    y = (x-nanmean(x))./nanstd(x);
end

function y=gauss(x,a,m,s)
    y = a;
    y = y.*exp(-0.5.*((x-m)./s).^2);
end

% %% Plotting composite indices
% figureHandle=figure;
% resize_figure(figureHandle,125,300);
% colormap('inferno');
% 
% p = panel();
% p.pack(1,2);
% p.margintop = 4;
% % p.marginleft = 15;
% p.marginbottom = 20;
% p.fontsize=14;
% p.fontname='Times';
% p(1,1).pack(2,2);
% p(1,2).pack('h',{3/4 []});
% p.select('all');
% p(1,1).de.margin=2;
% p(1,1).marginright = 40;
% % p.de.marginright=30;
% 
% p(1,1,1,1).select();
% histogram2(G1,P1,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none'); 
% % c11=colorbar_thin('YLabel','# of Samples');
% xlim([-2.5,2.5]);
% ylim([-2.5,2.5]);
% % xlabel('G1');
% ylabel('P*_{(1)}');
% set(gca,'XTickLabel','','YTick',[-2,0,2]);
% caxis([0,1200]);
% R11=corrcoef(G1,P1);
% legend(['r_1_1=',num2str(R11(1,2),2)]);
% 
% p(1,1,1,2).select();
% histogram2(G2,P1,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none');
% c12=colorbar_thin('YLabel','# of Samples');
% set(c12,'Position',get(c12,'Position')-[0,0.25,0,0],'Units','normalized');
% c12.Ticks=[0,600,1200]; c12.TickLabels={'0','0.6k','1.2k'};
% 
% xlim([-2.5,2.5]);
% ylim([-2.5,2.5]);
% % xlabel('G2');
% % ylabel('P1');
% set(gca,'YTickLabel','','XTickLabel','');
% caxis([0,1200]);
% R12=corrcoef(G2,P1);
% legend(['r_1_2=',num2str(round(R12(1,2),2),2)]);
% 
% p(1,1,2,1).select();
% histogram2(G1,P2,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none'); 
% % c11=colorbar_thin('YLabel','# of Samples');
% xlim([-2.5,2.5]);
% ylim([-2.5,2.5]);
% set(gca,'YTick',[-2,0,2],'XTick',[-2,0,2]);
% xlabel('G*_{(1)}');
% ylabel('P*_{(2)}');
% caxis([0,1200]);
% R21=corrcoef(G1,P2);
% legend(['r_2_1=',num2str(round(R21(1,2),2),2)]);
% 
% 
% p(1,1,2,2).select();
% histogram2(G2,P2,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none'); 
% % c11=colorbar_thin('YLabel','# of Samples');
% xlim([-2.5,2.5]);
% ylim([-2.5,2.5]);
% xlabel('G*_{(2)}');
% % ylabel('P2');
% set(gca,'YTickLabel','','XTick',[-2,0,2]);
% caxis([0,1200]);
% R22=corrcoef(G2,P2);
% legend(['r_2_2=',num2str(R22(1,2),2)]);
% 
% % figureHandle2=figure;
% % resize_figure(figureHandle2,125,150);
% % p = panel();
% 
% % p.select('all');
% p(1,2).de.margin = 2;
% p(1,2,1).select();
% % cx={'AL','AU','PCI','K_p','Dst','Am'};
% % cx={'AL','AU','PCI','K_p','Am','Dst','P_{ips}','F_{e130}','S_{since}','Tilt','V_{sw}','n_{sw}','sin^2(\theta)^(2/3)','B_z','|B|','\theta_{B_n^3}','R_{qint10}'};
% cx = {'AL','AU','K_p','V_{sw}','n_{sw}','|B|'}; %'sin^2(\theta)^{2/3}'
% x = categorical(cx);
% x = reordercats(x,cx);
% y = A(:,1);
% ylim([-1,1]);
% b1=bar(x,A./max(abs(A)));
% b1(1).FaceColor = [0.5,0.5,0.5];
% b1(2).FaceColor = [1,0.5,0.5];
% ylabel('Coefficients');
% l1=legend('G*_1','G*_2','Location','southwest');
% set(gca,'YTick',[-1,-0.5,0,0.5,1]);
% 
% p(1,2,2).select();
% % cx={'AL','AU','PCI','K_p','Dst','Am'};
% dx={'mP_{e}','mP_{i}'};
% x = categorical(dx);
% x = reordercats(x,dx);
% b2=bar(x,B./max(abs(B)),0.8);
% b2(1).FaceColor = [0.6,0.8,0.6];
% b2(2).FaceColor = [0.5,0.6,0.8];
% % ylabel('Coefficients');
% ylim([-1,1]);
% set(gca,'YTickLabel','','YColor','none','YTick',[-1,-0.5,0,0.5,1]);
% hold on;
% % plot(get(gca,'xlim'),[0 0],'k');
% yline(0,'-k');
% l2=legend([b2(1) b2(2)],'P*_1','P*_2','Location','southwestoutside');
% l2.Position(2) = l1.Position(2);
% l2.Position(1) = l2.Position(1) -0.05;
% %export_fig('G:\My Drive\Nithin\Job Search\JEPostDoc\Data\Final\Fig4_corr_and_coeff.pdf','-r600','-pdf','-nocrop');
