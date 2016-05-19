function [ total_electron, altitude_grid, t ] = pfisr_avg_e_v2(fileNameStr,altitude_grid)
%% PFISR_AVG_E_V2.m 
% calculates the average electron density per altitude
% from the electron density measurements along the different beams of
% PFISR system. 

%% Function variables
% Input  : fileNameStr      [String] - The HDF5 file which contains the PFISR data
%        : altitude_grid    [N x 1]  - A column vector of the altitude for which the 
%                                      for which in mean electron density is required [in km]                              
% Output : total_electron   - The average electron density per altitude
%                             [m^-3]
%          altitude_grid    - The default altitude (of beam 1) or user
%                             specified altitude points
%%
% Invoking the GeoData class to extract range and density from the HDF5 file
    pfisrGD = GeoData(@readMadhdf5,fileNameStr,{'range','popl'});

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
    altitude=zeros(range_data_l(1)/nb-1,1);
    ne=zeros(range_data_l(1)/nb-1,nb,time_data_l(2));
    
    
    inb=1; %Initial beam number
    % Calculating the altitude for each beam
    for beam=1:1:nb
        for i=1:1:(range_data_l(1)/nb - 1)
            r1=pfisrGD.dataloc((i-1)*nb+inb,1);
            ne(i,beam,1:1:time_data_l(2))=10.^pfisrGD.data.popl((i-1)*nb+beam,1:1:time_data_l(2));
            altitude(i,beam)=r1*sin(pfisrGD.dataloc((i-1)*nb+beam,end)*pi/180);
        end;
    end;
    if (nargin<2) 
    % Assigning the final altitude grid
    altitude_grid=altitude(:,1); %Here 1 - is beam 1, along the magnetic field line
    else
        disp('Using user supplied altitude grid');
    end;
    ne_interp=zeros(size(altitude_grid,1),nb,time_data_l(2));
    for beam=1:1:nb
        ne_interp(:,beam,:)=interp1(altitude(:,beam),ne(:,beam,:),altitude_grid,'linear','extrap');
    end;
    
    ne_sum=zeros(size(ne_interp,1),size(ne_interp,3));
    ne_mean=zeros(size(ne_interp,1),size(ne_interp,3));
    for i=1:1:size(ne_interp,1)
        ne_sum(i,:)=squeeze(nansum(ne_interp(i,:,:),2))';
        ne_mean(i,:)=ne_sum(i,:)./(nb*ones(1,time_data_l(2))-squeeze(sum(isnan(ne_interp(i,:,:)),2))');
    end;
    
    total_electron=ne_mean;

    t=unixtime2matlab(pfisrGD.times(:,1));

end   



