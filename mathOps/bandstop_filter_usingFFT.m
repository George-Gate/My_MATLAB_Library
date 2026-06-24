% filtered_sig = bandstop_filter_usingFFT(sig, samplingRate, stopband)
% Version 1.0 @ 2026-06-22
% By GeorgeGate
% 带阻滤波器，使用fft/ifft，在频谱上进行滤波
% sig: 时域信号，需要是等间隔采样
% samplingRate: 时域信号的采样率
% stopband: 阻带频率范围，阻带内的频谱会被强制设为0
%       stopband = [a,b]          --   阻带为a~b
%       stopband = {[a,b];[c,d]}  --   返回两组滤波后的信号，阻带分别为[a,b]和[c,d]
% filtered_sig: 滤波后的时域信号

function filtered_sig = bandstop_filter_usingFFT(sig, samplingRate, stopband)
    
    if ~iscell(stopband) && length(stopband) == 2
        stopband = {stopband};
    end

    spectrum = fft(sig);
    
    %this part of the code generates that frequency axis of spectrum
    N=length(sig);
    if mod(N,2)==0
        k=-N/2:N/2-1; % N even
    else
        k=-(N-1)/2:(N-1)/2; % N odd
    end
    T=N/samplingRate;
    freq=k/T; %the frequency axis
    freq = ifftshift(freq);
    
    % filter sig
    for ii = 1:length(stopband)
        filterIdx = stopband{ii}(1)<=abs(freq) & abs(freq)<=stopband{ii}(2);
        filtered_spectrum = spectrum;
        filtered_spectrum(filterIdx) = 0;
        filtered_sig{ii} = real(ifft(filtered_spectrum));
    end

    if isscalar(filtered_sig)
        filtered_sig = filtered_sig{1};
    end
end