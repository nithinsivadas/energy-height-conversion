function peakIonizationAlt=calculate_peak_altitude_of_ionization(energy,time0,lat,lon,alt)
    % Calculated the peak ionization altitude using the Segienko & Ivanov
    % 1993 model of electron precipitation and MSIS00
    % lat, lon, for each alt
    % time0 just 1 element
    % energy just 1 element
    if nargin<5
        alt = (60:0.1:200)';
    end
    if nargin<4
        lat = 64.5*ones(size(alt));
    end
    if nargin<3
        lon = -147.5*ones(size(alt));
    end
    dE = 10;
    if energy >= 10^3 && energy<=1*10^6
        energyBin = (energy-dE:dE:energy+dE)';
    else
        error('Energy out of acceptable range 1keV to 1MeV');
    end
    
    A = get_energy_dep_matrix(alt, energyBin, lat, lon, time0);
   
    [q_max, ialt] = max(A);
    
    peakIonizationAlt = alt(ialt(2));
    
    
end