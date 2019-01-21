%% Test DASC PokerFlat Projection at 110 km

calFileEl = [initialize_root_path,'LargeFiles',filesep,'DASC',filesep,...
        '20080326',filesep,'PKR_Cal_before_2011',filesep,'PKR_20111006_EL_10deg.FITS'];
    
calFileAz = [initialize_root_path,'LargeFiles',filesep,'DASC',filesep,...
        '20080326',filesep,'PKR_Cal_before_2011',filesep,'PKR_20111006_AZ_10deg.FITS'];
    
azOldRes = fitsread(calFileAz);
elOldRes = fitsread(calFileEl);
        
[data.az, data.el, data.lowAzGradientFilter] = calibrate_DASC_pixels(azOldRes,elOldRes);
[data.minElFilter, data.lat, data.lon, data.alt] = DASC_aer_to_geodetic_v2018(data.az, data.el,...
    0, 110);

[dataOld.minElFilter, dataOld.lat, dataOld.lon, dataOld.alt] = DASC_aer_to_geodetic_v2019(data.az, data.el,...
    0, 110);
