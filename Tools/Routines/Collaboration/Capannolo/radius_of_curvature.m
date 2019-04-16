%% To test Luisa's results of radius of curvature at the equatorial plane

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
smoothFactor = 10;
kext = find_irbem_magFieldModelNo('TS89');
maginput = interp_nans(maginput);
maginput=filter_irbem_maginput(kext,maginput);

%% Coordinate transform
for iTime = 1:1:length(time)
    thisTime = time(iTime);
    xGEO(iTime,:) = onera_desp_lib_coord_trans([alt(iTime),lat(iTime),...
        convert_longitude(lon(iTime),'360to180')],'gdz2geo',thisTime);
end
%% Calculate Rc
[~, minRc] = get_min_scatter_energy(time,magTime,maginput,xGEO,kext); 
%% Function 
%         PARMOD = get_parmod(kext,thisMaginput);
function [minRc]=get_min_scatter_energy(time,magTime,maginput,xGEO,kext)
    global GEOPACK1;
    C = define_universal_constants;
    PARMOD = zeros(10,1);
    PARMODT96 = zeros(10,1);
    for iTime = 1:1:length(time)
        thisTime = time(iTime);
        t = datetime(datevec(time(iTime)));
        thisMaginput = interp1(magTime,maginput,thisTime,'nearest');
        GEOPACK_RECALC(year(t),day(t,'dayofyear'),hour(t),minute(t),second(t));
        PARMODT96 = get_parmod(find_irbem_magFieldModelNo('TS96'),thisMaginput);
%         PARMODT96(1) = thisMaginput(5);
%         PARMODT96(2) = thisMaginput(2);
%         PARMODT96(3) = thisMaginput(6);
%         PARMODT96(4) = thisMaginput(7);

        xGSM = onera_desp_lib_coord_trans(xGEO(iTime,:),'geo2gsm',thisTime);
        [~,~,~,X1,Y1,Z1,~] = GEOPACK_TRACE (xGSM(1),xGSM(2),xGSM(3),...
            1.0,50,(6371.2+110)/6371.2,0,PARMODT96,'T96','GEOPACK_IGRF_GSM');

        [~,~,~,X2,Y2,Z2,~] = GEOPACK_TRACE (xGSM(1),xGSM(2),xGSM(3),...
            -1.0,50,(6371.2+110)/6371.2,0,PARMODT96,'T96','GEOPACK_IGRF_GSM');

        XX{iTime} = [fliplr(X1),X2];
        YY{iTime} = [fliplr(Y1),Y2];
        ZZ{iTime} = [fliplr(Z1),Z2];

        Kc{iTime} = geopack_find_curvature(XX{iTime},YY{iTime},ZZ{iTime});

        [maxKc,indx] = max(Kc{iTime});
        minRc(iTime) = 1./maxKc;
        magEq{iTime} = [XX{iTime}(indx),YY{iTime}(indx),ZZ{iTime}(indx)];
%         [BX,BY,BZ] = T96(0,PARMODT96,GEOPACK1.PSI,magEq{iTime}(1),magEq{iTime}(2),magEq{iTime}(3));
%         KE(iTime) = ((((minRc(iTime).*C.RE ).*(C.e).*(BZ*10^-9)).^2).*(2^-7).*(C.me).^-1).*(10^-3).*(C.e).^-1; %keV
    end
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