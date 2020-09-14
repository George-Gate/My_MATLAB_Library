function [freq, X] = fft_easySpectrum(x, samplingRate, fft_window)
% Version 1.1 @ 2020-06-02
% By GeorgeGate
    
    %this is a custom function that helps in plotting the two-sided spectrum
    %x is the signal that is to be transformed
    %samplingRate is the sampling rate

    if nargin <3
        fft_window = 'rect';
    end
    
    N=length(x);

    switch upper(fft_window)
        case 'RECT'
            windowVec = ones(N,1); 
        case 'HANN'
            % 时域: w(t)=0.5 - 0.5*cos(2*pi*t/T), 其中T是采样时长
            windowVec = hann(N,'periodic')*2;
        case 'HAMMING'
            % 时域: w(t)=0.54 - 0.46*cos(2*pi*t/T), 其中T是采样时长
            windowVec = hamming(N,'periodic')*50/27;
        otherwise
            error('Unrecognized FFT window.');
    end
    windowVec = reshape(windowVec,size(x));
    
    
    %this part of the code generates that frequency axis
    if mod(N,2)==0
        k=-N/2:N/2-1; % N even
    else
        k=-(N-1)/2:(N-1)/2; % N odd
    end
    T=N/samplingRate;
    freq=k/T; %the frequency axis

    %takes the fft of the signal, and adjusts the amplitude accordingly
    X=2/N*fft(x.*windowVec); % normalize the data so that the height of peak equals Vpk, 
                  % ie. for A*sin(wt), you will get a peak whose height is A.
    X=fftshift(X); %shifts the fft data so that it is centered

end