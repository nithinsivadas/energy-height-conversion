%% Linux create omni script
try 
    h5FileStr = "/media/nithin/PFISR_002_006/Nithin/omni.h5";   
    status = create_omni_HDF5(h5FileStr);
catch ME
    ME
end
exit
