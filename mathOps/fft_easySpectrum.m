function [freq, spectrum, enbw] = fft_easySpectrum(x, samplingRate, fft_window)
% [freq, spectrum, enbw] = fft_easySpectrum(x, samplingRate, fft_window)
% Version 2.3 @ 2021-09-15
% By GeorgeGate
    
    %this is a custom function that helps in plotting the two-sided spectrum
    %x is the signal that is to be transformed
    %samplingRate is the sampling rate

    if nargin <3
        fft_window = 'rect';
    end
    
    N=length(x);
    %this part of the code generates that frequency axis
    if mod(N,2)==0
        k=-N/2:N/2-1; % N even
    else
        k=-(N-1)/2:(N-1)/2; % N odd
    end
    T=N/samplingRate;
    freq=k/T; %the frequency axis

    switch upper(fft_window)
        case 'RECT'
            windowVec = ones(N,1); 
            enbw = 1/T;
        case 'HANN'
            % 时域: w(t)=0.5 - 0.5*cos(2*pi*t/T), 其中T是采样时长
            windowVec = hann(N,'periodic')*2;
            enbw = 1.5/T;
        case 'HAMMING'
            % 时域: w(t)=0.54 - 0.46*cos(2*pi*t/T), 其中T是采样时长
            windowVec = hamming(N,'periodic')*50/27;
            enbw = 1987/1458/T;
        otherwise
            error('Unrecognized FFT window "%s".', fft_window);
    end
    windowVec = reshape(windowVec,size(x));
    

    %takes the fft of the signal, and adjusts the amplitude accordingly
    spectrum=1/N*fft(x.*windowVec); % normalize the data so that the height of peak equals Vpk/2, 
                  % ie. for A*sin(wt), you will get a peak whose height is A/2.
    spectrum=fftshift(spectrum); %shifts the fft data so that it is centered
    
    spectrum = reshape(spectrum, size(x));
    freq = reshape(freq, size(x));
end