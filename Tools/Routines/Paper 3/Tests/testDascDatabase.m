%% Testing dascData old and new
clear all;

h5New = 'C:\Users\nithin\Documents\GitHub\LargeFiles\dascDatabase.h5';
h5Old = 'C:\Users\nithin\Documents\GitHub\LargeFiles\dascDatabase_copy.h5';

year = 2013;

dascNew.time = h5read(h5New,['/DASC/',num2str(year),'/time']);
dascNew.wavelength = h5read(h5New,['/DASC/',num2str(year),'/wavelength']);
dascNew.waveCodeCell = h5read(h5New,['/DASC/',num2str(year),'/wavelengthCode']);
dascNew.waveCode = conv_wavecode(dascNew.waveCodeCell);
dascNew.wavelength = dascNew.waveCode(dascNew.wavelength);
disp(dascNew.waveCodeCell)
A = [dascNew.time(:)'; dascNew.wavelength(:)'];
A0 = (A(:,A(2,:)==0)');
A428 = sort(A(:,A(2,:)==428)');
A558 = sort(A(:,A(2,:)==558)');
A630 = sort(A(:,A(2,:)==630)');
A1111 = sort(A(:,A(2,:)==1111)');
A715 = sort(A(:,A(2,:)==715)');


dascOld.time = h5read(h5Old,['/DASC/',num2str(year),'/time']);
dascOld.wavelength = h5read(h5Old,['/DASC/',num2str(year),'/wavelength']);

B = [dascOld.time; dascOld.wavelength];
B0 = (B(:,B(2,:)==0)');
B428 = sort(B(:,B(2,:)==428)');
B558 = sort(B(:,B(2,:)==558)');
B630 = sort(B(:,B(2,:)==630)');

%%
figure; 
plot(datetime(A0(:,1),'convertfrom','posixtime'),0*ones(size(A0,1),1),'.r'); 
hold on; 
plot(datetime(B0(:,1),'convertfrom','posixtime'),-10*ones(size(B0,1),1),'.k');

plot(datetime(A428(:,1),'convertfrom','posixtime'),428*ones(size(A428,1),1),'.r'); 
plot(datetime(B428(:,1),'convertfrom','posixtime'),418*ones(size(B428,1),1),'.k'); 

plot(datetime(A558(:,1),'convertfrom','posixtime'),558*ones(size(A558,1),1),'.r'); 
plot(datetime(B558(:,1),'convertfrom','posixtime'),548*ones(size(B558,1),1),'.k'); 

plot(datetime(A630(:,1),'convertfrom','posixtime'),630*ones(size(A630,1),1),'.r'); 
plot(datetime(B630(:,1),'convertfrom','posixtime'),620*ones(size(B630,1),1),'.k'); 

plot(datetime(A1111(:,1),'convertfrom','posixtime'),1111*ones(size(A1111,1),1),'.r'); 

plot(datetime(A715(:,1),'convertfrom','posixtime'),715*ones(size(A715,1),1),'.r'); 



ylim([-100,1200]);
legend('cURL','wget');
title('DASC Data Downloaded using different methods');
set(gca,'YTick',[0,428,558,630,715,1111],'yTickLabel',{'0000','4280','5580','6300','L715','1111'});
ylabel('Wavelength [A^{\circ}]');
set(gcf,'color','w');
export_fig(strcat('G:\My Drive\Research\Projects\Paper 3\Figures\Documentation\dascData\',num2str(year),'.png'),'-r300','-png','-nocrop');

function waveCode = conv_wavecode(waveCodeCell)
    n = length(waveCodeCell(1,:));
    
    for i = 1:n
        if strcmp(string(waveCodeCell(1,i)),"L715")
            waveCode(i,1) = 715;
        else
            waveCode(i,1) = str2num(string(waveCodeCell(1,i)));
        end
    end
    
end