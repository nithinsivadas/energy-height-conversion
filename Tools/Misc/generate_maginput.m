function [maginput,time,header,units] = generate_maginput(omniH5FileStr, timeMin, timeMax)
    %% Generate maginput from omni.h5
    omniTime = unixtime2matlab(h5read(omniH5FileStr,'/Time'));
    minTimeIndx = find_time(omniTime,datestr(timeMin));
    maxTimeIndx = find_time(omniTime,datestr(timeMax));
    deltaTimeIndx = maxTimeIndx - minTimeIndx +1;
    timeIndx = minTimeIndx:1:maxTimeIndx;
    header = ["Kp","SYM_H(Dst)","ProtonDensity",...
        "Vsw","Psw","ByGSM","BzGSM","G1","G2","G3","W1","W2","W3","W4","W5","W6","AL"];
    units = ["a.u.","nT","n/cc","km/s","nPa","nT","nT",...
        "a.u.","a.u.","a.u.","a.u.","a.u.","a.u.","a.u.","a.u.","a.u.","a.u."];
    maginput(:,1) = h5read(omniH5FileStr, '/Indices/Kp', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,2) = h5read(omniH5FileStr, '/Indices/SYM_H', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,3) = h5read(omniH5FileStr, '/ProtonDensity', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,4) = h5read(omniH5FileStr, '/Velocity/V', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,5) = h5read(omniH5FileStr, '/FlowPressure', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,6) = h5read(omniH5FileStr, '/BField/ByGSM', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,7) = h5read(omniH5FileStr, '/BField/BzGSM', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,8) = h5read(omniH5FileStr, '/TSY/G1', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,9) = h5read(omniH5FileStr, '/TSY/G2', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,10) = h5read(omniH5FileStr, '/TSY/G3', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,11) = h5read(omniH5FileStr, '/TSY/W1', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,12) = h5read(omniH5FileStr, '/TSY/W2', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,13) = h5read(omniH5FileStr, '/TSY/W3', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,14) = h5read(omniH5FileStr, '/TSY/W4', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,15) = h5read(omniH5FileStr, '/TSY/W5', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,16) = h5read(omniH5FileStr, '/TSY/W6', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,17) = h5read(omniH5FileStr, '/Indices/AL', [1 minTimeIndx], [1 deltaTimeIndx]);
    maginput(:,18:25) = nan(length(timeIndx),8);
    
    time = omniTime(timeIndx);
end