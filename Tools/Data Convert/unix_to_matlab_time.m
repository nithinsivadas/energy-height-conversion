function mTime = unix_to_matlab_time(uTime)
%% unix_to_matlab_time.m Converts unix time to matlab time
mTime = uTime/(24*3600) + datenum('jan-01-1970');