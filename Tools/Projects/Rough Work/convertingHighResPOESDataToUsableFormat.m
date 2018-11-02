%% Temp
load('G:\My Drive\Research\Projects\Paper 2\Data\NOAA17\n1720080326_table.mat');
n=size(n1720080326,1);
di = 1./n;
multiWaitbar('Calculating...',0);
for i = 1:1:n
sec = mod(i-1,4);
    if sec==0
    time(i)=(datenum([num2str(n1720080326{i,1}),' 1 ',num2str(n1720080326{i,2}),' ',num2str(n1720080326{i,3}),':',num2str(n1720080326{i,4}),':',num2str(n1720080326{i,5})]));
    lat(i,1) = n1720080326{i,6};
    lon(i,1) = n1720080326{i,7};
    else
    time(i)= (datenum([num2str(n1720080326{i,1}),' 1 ',num2str(n1720080326{i,2}),' ',num2str(n1720080326{i,3}),':',num2str(n1720080326{i,4}),':',num2str(n1720080326{i,5}+2*sec)]));
    lat(i,1)=nan;
    lon(i,1)=nan;
    end

multiWaitbar('Calculating...','Increment',di);
end

poes.time = time';
poes.lat = interp_nans(lat);
poes.lon = interp_nans(lon);
poes.mep0P1 = n1720080326{:,10};
poes.mep0P2 = n1720080326{:,11};
poes.mep0P3 = n1720080326{:,12};
poes.mep0P4 = n1720080326{:,13};
poes.mep0P5 = n1720080326{:,14};
poes.mep0P6 = n1720080326{:,15};
poes.mep0E1 = n1720080326{:,16};
poes.mep0E2 = n1720080326{:,17};
poes.mep0E3 = n1720080326{:,18};
poes.mep90P1 = n1720080326{:,19};
poes.mep90P2 = n1720080326{:,20};
poes.mep90P3 = n1720080326{:,21};
poes.mep90P4 = n1720080326{:,22};
poes.mep90P5 = n1720080326{:,23};
poes.mep90P6 = n1720080326{:,24};
poes.mep90E1 = n1720080326{:,25};
poes.mep90E2 = n1720080326{:,26};
poes.mep90E3 = n1720080326{:,27};