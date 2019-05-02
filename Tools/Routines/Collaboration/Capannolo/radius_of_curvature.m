%% To test Luisa's results of radius of curvature at the equatorial plane

clear all;

dataFile = 'G:\My Drive\Research\Projects\Collaborations\Capannolo Luisa\IRBEM_Test\n18_coords_interpolated.txt';
fileID = fopen(dataFile);
A = textscan(fileID, '%26q %f %f %f','HeaderLines',1,'Delimiter','\t');
fclose(fileID);

time = datenum(A{1},'yyyy mm dd hh MM ss.fff');
alt = datenum(A{2}); % in km
lat = datenum(A{3}); % in deg
lon = datenum(A{4}); % in deg

%% Generating maginput

omniH5='G:\My Drive\Research\Projects\Data\omni.h5';
[maginput,magTime] = generate_maginput(omniH5,datestr(time(1)),datestr(time(end)));

%% Calculate L-shell or equatorial parameters
disp('Filtering...');
kext = find_irbem_magFieldModelNo('TS89');
% kext = find_irbem_magFieldModelNo('TS04storm');

maginput = interp_nans(maginput);
maginput=filter_irbem_maginput(kext,maginput);

%%
xGDZ = [alt, lat, convert_longitude(lon)];
%% Calculate Rc
% for i = 1:100
%     maginput(1,1) = i;
    [timeRc, Rc, divBArr, data] = get_radius_of_curvature(time,magTime,maginput,xGDZ,kext); 
%     kp(i) = i;
%     kpRc(i) = Rc;
%     kpPOSIT(:,i) = data.POSIT;
end
%% Plot
h=figure;
p=create_panels(h,'totalPanelNo',2,'demargin',15,'panelHeight',50);

p(1,1).select();
title(find_irbem_magFieldModelStr(kext));
t = datetime(datevec(datestr(timeRc)));
plot(t,Rc);
ylim([0,20]);
ylabel('R_curve [R_E]');

p(1,2).select();
plot(t,divBArr);
ylim([-300,300]);
ylabel('divB');
%% Function 
%         PARMOD = get_parmod(kext,thisMaginput);
function [timeRc, Rc, divBArr, data]=get_radius_of_curvature(time,magTime,maginput,xGDZ,kext)
    options = [0,0,0,0,0];
    multiWaitbar('Processing...',0);
    dt = 1./length(time);
    for iTime = 1:1:length(time)
%     for iTime = 1:1:1
        multiWaitbar('Processing...','Increment',dt);
        thisTime = time(iTime);
        t = datetime(datevec(time(iTime)));
        thisMaginput = interp1(magTime,maginput,thisTime,'nearest');
        
%         [~,~,~,~,POSIT] = onera_desp_lib_trace_field_line(kext,options,...
%             0,thisTime,xGDZ(iTime,1),xGDZ(iTime,2),xGDZ(iTime,3),...
%             thisMaginput,(6371.2+110)/6371.2);
        
        [Bmin, POSIT] = onera_desp_lib_find_magequator(kext, options, 0,...
            thisTime, xGDZ(iTime,1), xGDZ(iTime,2), xGDZ(iTime,3), thisMaginput);
            
        [Bgeo, B, gradBmag, diffB] = onera_desp_lib_get_bderivs(kext,...
            options, 1, thisTime, POSIT(1), POSIT(2), POSIT(3), thisMaginput);
        
        [grad_par,grad_perp,grad_drift,curvature,Rcurv,curv_drift,curlB,divB] =....
            onera_desp_lib_compute_grad_curv_curl(Bgeo,B,gradBmag,diffB); 
        divBArr(iTime) = divB; 
        Rc(iTime) = Rcurv;
        timeRc(iTime) = thisTime;
        
        data.Bmin(iTime) = Bmin;
        data.POSIT(:,iTime) = POSIT;
        data.Bmag(iTime) = B;
        data.curvDrift(:,iTime) = curv_drift;
    end
%         [BX,BY,BZ] = T96(0,PARMODT96,GEOPACK1.PSI,magEq{iTime}(1),magEq{iTime}(2),magEq{iTime}(3));
%         KE(iTime) = ((((minRc(iTime).*C.RE ).*(C.e).*(BZ*10^-9)).^2).*(2^-7).*(C.me).^-1).*(10^-3).*(C.e).^-1; %keV
    
    multiWaitbar('Processing...',1);
end

function [PARMOD,IOPT,magStr] = get_parmod(magFieldNo,maginput)
    
PARMOD = zeros(10,1);

if magFieldNo==4
    kp = round(maginput(1)/10)+1;
    magStr = 'T89';
    PARMOD(1) = kp;
    IOPT = kp;
elseif magFieldNo==7
    magStr = 'T96';
    PARMOD(1) = maginput(5);
    PARMOD(2) = maginput(2);
    PARMOD(3) = maginput(6);
    PARMOD(4) = maginput(7);    
    IOPT = 0;
elseif magFieldNo==9
    magStr = 'T01';
    PARMOD(1) = maginput(5);
    PARMOD(2) = maginput(2);
    PARMOD(3) = maginput(6);
    PARMOD(4) = maginput(7);    
    PARMOD(5) = maginput(8);    
    PARMOD(6) = maginput(9);
    IOPT=0;
else
    error('GEOPACK version does not contain the specified external field.');
end

end