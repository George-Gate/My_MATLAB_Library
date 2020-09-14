% generate a time series that is a sampling of white noise
samplingRate = 1000;  % sampling rate, in unit of Hz
Npts = 2^12;   % number of sampled points



dt = 1/samplingRate;
tList = 0:dt:dt*(Npts-1);