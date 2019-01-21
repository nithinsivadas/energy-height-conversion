function [] = skymap(YYYY,MM,DD,hh,mm,ss,long,lat,UTdiff)

disp('Please wait,... Processing sky map!');

YYYY=2008;
MM=03;
DD=26;
hh=11;
mm=30;
ss=50;
long=-147.5;
lat=65.12;
UTdiff=+0;

maxmag=5.5;
if YYYY<2000 || isnan(YYYY)==1  || isreal(YYYY)==0
    errordlg('You have entered an invalid Year!, Year must be a number from 2000 and above','Bad Input','modal');
    return
end
if MM<1 || MM>12 || isnan(MM)==1 || isreal(MM)==0
    errordlg('You have entered an invalid Month!, Month must be a number between 1 and 12 inclusive','Bad Input','modal');
    return
end
if (DD<1 || DD>31 || isnan(DD)==1 || isreal(DD)==0) && ismember(MM, [1,3,5,7,8,10,12])==1
    errordlg('You have entered an invalid Day of Month!, This Month does not have the number of days you entered','Bad Input','modal');
    return
end
if (DD<1 || DD>30 || isnan(DD)==1 || isreal(DD)==0) && ismember(MM, [4,6,9,11])==1
    errordlg('You have entered an invalid Day of Month!, This Month does not have the number of days you entered','Bad Input','modal');
    return
end
if (DD<1 || DD>28 || isnan(DD)==1 || isreal(DD)==0) && MM==2 && mod(YYYY,4)~=0 
    errordlg('You have entered an invalid Day of Month!, This Month does not have the number of days you entered','Bad Input','modal');
    return
end
if (DD<1 || DD>29 || isnan(DD)==1 || isreal(DD)==0) && MM==2 && mod(YYYY,4)==0 
    errordlg('You have entered an invalid Day of Month!, This Month does not have the number of days you entered','Bad Input','modal');
    return
end
if hh<0 || hh>24 || isnan(hh)==1 || isreal(hh)==0
    errordlg('You have entered an invalid hour!, hour must be a number from 0 to 24','Bad Input','modal');
    return
end
if mm<0 || mm>60 || isnan(mm)==1 || isreal(mm)==0
    errordlg('You have entered an invalid minute!, minute must be a number from 0 to 60','Bad Input','modal');
    return
end
if ss<0 || ss>60 || isnan(ss)==1 || isreal(ss)==0
    errordlg('You have entered an invalid second!, second must be a number from 0 to 60','Bad Input','modal');
    return
end
if long<-180 || long>180 || isnan(long)==1 || isreal(long)==0
    errordlg('You have entered an invalid longitude!, longitude must be a number between -180 and 180 inclusive','Bad Input','modal');
    return
end
if lat<-90 || lat>90 || isnan(lat)==1 || isreal(lat)==0
    errordlg('You have entered an invalid latitude!, latitude must be a number between -90 and 90 inclusive','Bad Input','modal');
    return
end
if isnan(maxmag)==1 || isreal(maxmag)==0
    errordlg('You have entered an invalid star magnitude!, star magnitude must be a real number','Bad Input','modal');
    return
end


%-----------------------------------------------------------------
% GREENWICH MEAN SIDEREAL TIME
%-----------------------------------------------------------------
dfrac=(hh+mm/60+ss/3600)/24;
dwhole =367*YYYY-fix(7*(YYYY+fix((MM+9)/12))/4)+fix(275*MM/9)+DD-730531.5;
d=dwhole+dfrac;  %no. of days elapsed since J2000.0
dd=d-UTdiff/24;  %no of UT days elapsed since UT J2000.0
T=dd/36525;  %fraction of epoch/century time elapsed

% Meeus formula 11.4 for mean sidereal time at zero longitude (Greenwich Mean Sidereal Time)
GMST=280.46061837+(360.98564736629*dd)+(0.000387933*T^2)-(T^3/38710000); 
while GMST>360
    GMST=GMST-360;
end
while GMST<0
    GMST=GMST+360;
end

%------------------------------------------------------------------
% LOCAL SIDEREAL TIME
%------------------------------------------------------------------
LST=GMST+long;
while LST>360
    LST=LST-360;
end
while LST<0
    LST=LST+360;
end
LST=LST/15;  %LST in hours

if maxmag<6
    fid=fopen('sao6.txt');
else
   fid=fopen('saoNAN.txt'); 
end
M=textscan(fid, '%f %f %f %f %f %f %f %f', 'headerlines', 1);
sao=M{1}; mvis=M{3}; RA=(M{5}+(M{7}*(YYYY-2000)))/15; DEC=(M{6}+(M{8}*(YYYY-2000)));

scrnsz=get(0,'ScreenSize');
figure('Position',[scrnsz(1)+scrnsz(3)*0.1 scrnsz(2)+scrnsz(4)*0.1 scrnsz(3)*0.8 scrnsz(4)*0.8]);
colordef black;
for i=1:length(sao)
    
if mvis(i)<maxmag
HA =(LST - RA(i))*15;
El = asind(sind(lat)*sind(DEC(i))+cosd(lat)*cosd(DEC(i))*cosd(HA));  %http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm
if El>0   %star is in horizon
    if El~=90 %star not at zenith
        AZ1 = (-cosd(DEC(i))*sind(HA))/cosd(El); 
        AZ2 = (cosd(lat)*sind(DEC(i))-sind(lat)*cosd(DEC(i))*cosd(HA))/cosd(El);
        if (AZ1>=0) && (AZ2>=0) %1st Quadrant anticlockwise from North
            [xxx,yyy]=linecirc(AZ2/AZ1,0,0,0,90-El);
            scatter(min(xxx),max(yyy),((maxmag-mvis(i))*2),'*','filled','MarkerEdgeColor','white', 'MarkerFaceColor', 'white'); hold on;                       
        elseif (AZ1>=0) && AZ2<0 %2nd Quadrant anticlockwise from North
            [xxx,yyy]=linecirc(AZ2/AZ1,0,0,0,90-El);
            scatter(min(xxx),min(yyy),((maxmag-mvis(i))*2),'*','filled','MarkerEdgeColor','white', 'MarkerFaceColor', 'white'); hold on;
        elseif AZ1<0 && (AZ2>=0) %4th Quadrant anticlockwise from North
            [xxx,yyy]=linecirc(AZ2/AZ1,0,0,0,90-El);
            scatter(max(xxx),max(yyy),((maxmag-mvis(i))*2),'*','filled','MarkerEdgeColor','white', 'MarkerFaceColor', 'white'); hold on;
        elseif AZ1<0 && AZ2<0 %3rd Quadrant anticlockwise from North
            [xxx,yyy]=linecirc(AZ2/AZ1,0,0,0,90-El);
            scatter(max(xxx),min(yyy),((maxmag-mvis(i))*2),'*','filled','MarkerEdgeColor','white', 'MarkerFaceColor', 'white'); hold on;
        end
    elseif El==90 %star at zenith
        scatter(0,0,((maxmag-mvis(i))*2),'*','filled','MarkerEdgeColor','white', 'MarkerFaceColor', 'white'); hold on;
    end
end
end
end
fxn=['x^2+y^2-8100']; ezplot(fxn,[-90,90]);  %draw outer circle
title(['Sky Map for Longitude ' num2str(long) ', Latitude ' num2str(lat) ', ' num2str(YYYY) '-' num2str(MM) '-' num2str(DD) ', ' num2str(hh) ':' num2str(mm) ':' num2str(ss) ' LT']); ylabel(''); xlabel(''); axis equal tight;
fclose all;