% filtered_sig = bandpass_filter_usingFFT(sig, samplingRate, passband, trans_width)
% Version 1.0 @ 2026-06-22
% By GeorgeGate
% 带通滤波器，使用fft/ifft，在频谱上进行滤波
% sig: 时域信号，需要是等间隔采样
% samplingRate: 时域信号的采样率
% passband: 通带频率范围，通带以外的频谱会被强制设为0
%       passband = [a,b]          --   通带为a~b,  0<=a<b
%       passband = {[a,b];[c,d]}  --   返回两组滤波后的信号，通带分别为[a,b]和[c,d]
% trans_width: 通带边缘过渡区域的相对宽度，默认值为0.05
%            在通带边缘，频谱使用cos函数平滑过度到零，过渡区域宽度为trans_width*(b-a)
% filtered_sig: 滤波后的时域信号


function filtered_sig = bandpass_filter_usingFFT(sig, samplingRate, passband, trans_width)
    
    if nargin <4
        trans_width = 0.05;
    end
    if ~iscell(passband) && length(passband) == 2
        passband = {passband};
    end

    N=length(sig);
    T=N/samplingRate;

    spectrum = fft(sig);
    
    %this part of the code generates that frequency axis of spectrum
    if mod(N,2)==0
        k=-N/2:N/2-1; % N even
    else
        k=-(N-1)/2:(N-1)/2; % N odd
    end
    freq=k/T; %the frequency axis
    freq = ifftshift(freq);
    
    % filter sig
    for ii = 1:length(passband)
        lb = passband{ii}(1);
        ub = passband{ii}(2);
        f_tran = trans_width*(ub-lb);
        mask_fun = @(f) 1*(lb<=f&f<=ub) + (0.5+0.5*cos(pi*(lb-f)/f_tran)).*(lb-f_tran<=f&f<lb) + (0.5+0.5*cos(pi*(ub-f)/f_tran)).*(ub<f&f<=ub+f_tran);
        filtered_spectrum = spectrum.*mask_fun(abs(freq));
        filtered_sig{ii} = real(ifft(filtered_spectrum));
    end

    if isscalar(filtered_sig)
        filtered_sig = filtered_sig{1};
    end
end