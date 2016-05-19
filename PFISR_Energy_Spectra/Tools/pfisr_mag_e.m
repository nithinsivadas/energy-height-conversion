function [ ne, altitude, t ] = pfisr_mag_e( fileNameStr, beam )
%% PFISR_AVG_E.m 
% calculates the average electron density per altitude
% from the electron density measurements along the different beams of
% PFISR system. 

%% Function variables
% Input  : fileNameStr      - The HDF5 file which contains the PFISR data
% Output : total_electron   - The average electron density per altitude
%                             [m^-3]
%          altitude         - Range of the first PFISR beam [km]
%%

% Invoking the GeoData class to extract range and density from the HDF5 file
    pfisrGD = GeoData(@readMadhdf5,fileNameStr,{'range','popl'});

    % Identifying the range of theta and phi within the data
    theta_min = min(pfisrGD.dataloc(:,3));
    theta_max = max(pfisrGD.dataloc(:,3));

    phi_min = min(pfisrGD.dataloc(:,2));
    phi_max = max(pfisrGD.dataloc(:,2));

    % Identifying the number of beams
    temp = pfisrGD.dataloc(1,1);
    j=1;
    while (pfisrGD.dataloc(j+1,1)==temp)
        j=j+1;
    end;    
    nb=j;

    range_data_l = size(pfisrGD.dataloc);
    time_data_l= size(pfisrGD.data.popl);

    % Averaging electron desnity per range slice
    ne=zeros(range_data_l(1)/nb-1,time_data_l(2));
%     total_electron=zeros(range_data_l(1)/nb-1,time_data_l(2));
    altitude=zeros(range_data_l(1)/nb-1,1);
    for i=1:1:(range_data_l(1)/nb - 1)
        ne(i,:)=(10.^pfisrGD.data.popl((i-1)*nb+beam,1:1:time_data_l(2)));    
        r1=pfisrGD.dataloc((i-1)*nb+beam,1);
        altitude(i)=r1*sin(pfisrGD.dataloc((i-1)*nb+beam,end)*pi/180);
    end;
    t=unixtime2matlab(pfisrGD.times(:,1));

end   

