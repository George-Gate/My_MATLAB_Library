function [spectrumData] = overlappingFFT_logscale(data, saRate, fft_len, varargin)
%OVERLAPPINGFFT_LOGSCALE overlapping FFT with logscale frequency list
% Author: GeorgeGate @ 2026-09-04
% Version: 3.1
% 计算PSD的平均值，使用overlapping FFT方式进行平均。频域上的取点满足近似的 logspace，可以有效地平滑高频段的频谱，同时保留低频段的频谱不丢失。
% 
% average method: mean(PSD)
% downsampling method: FIR lowpass
% single-sided
% min df: freq_full(2) = (saRate / 2^floor(log2(length(data)/fft_len))) / fft_len
%         will get best min df when  length(data)/fft_len is 2^n

    p = inputParser;
    p.addRequired('data', @(x)validateattributes(x, {'numeric'}, {'vector', 'finite', 'real'}));    % 原始的时域数据
    p.addRequired('saRate', @(x)validateattributes(x, {'numeric'}, {'scalar', 'positive'}));        % 原始数据的采样率
    p.addRequired('fft_len', @(x)validateattributes(x, {'numeric'}, {'scalar', 'integer', '<=', length(data), '>=', 2}));      % 单次FFT使用的数据长度
    p.addParameter('minGapTime', 1/saRate, @(x)validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));                 % 计算 overlapping FFT 时，相邻两个数据样本之间的最短时间间隔（单位：秒）
    p.addParameter('minGapLen', round(fft_len/2), @(x)validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));  % 计算 overlapping FFT 时，相邻两个数据样本之间的最小采样点间隔。将minGapLen设置成fft_len，可以计算非重叠的FFT平均值
                                                                                                                               % 经过测试，通常 gapLen 的默认值就能保证最优的平均效果，更小的 gapLen 没有太大的意义
    p.addParameter('optimize_fft_len', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));    % 是否在输入的 fft_len 附近寻找一个最优的 fft_len, 使得输出的频谱有最佳的分辨率？程序找到的 最优fft_len 总是比输入的 fft_len 小。
    p.addParameter('fft_window_name', 'Hann', @(x)ismember(upper(x), {'HAMMING','HANN','RECT'}));    % 计算FFT时使用的窗函数 
    
    p.addParameter('antiAliasPassbandRatio', 0.80, @(x)validateattributes(x, {'numeric'}, {'scalar', '>', 0, '<', 1}));   % 抗混叠下采样滤波器，passband edge / NEW Nyquist, default 0.80. The remaining 20% is used as transition band.
    p.addParameter('antiAliasStopbandAtten', 200, @(x)validateattributes(x, {'numeric'},  {'scalar', '>=', 20}));         % 抗混叠下采样滤波器，target FIR stopband attenuation in dB, default 100 dB.

    p.parse(data, saRate, fft_len, varargin{:});
    
    minGapTime = p.Results.minGapTime;
    minGapLen = p.Results.minGapLen;
    antiAliasPassbandRatio = p.Results.antiAliasPassbandRatio;
    antiAliasStopbandAtten = p.Results.antiAliasStopbandAtten;

    % error('Not finished');

    if p.Results.optimize_fft_len
        M = length(data);
        N = fft_len;     % target fft_len
        fft_len = floor( M / (2^nextpow2( M/N )) );
        est_min_df = (saRate / 2^floor(log2(length(data)/fft_len))) / fft_len;
    end
    
    data = reshape(data, 1, []);   % transform to row vector
    
    freq_full = [];
    spectrum_full = [];
    spectrumDensity_full = [];
    
    % 使用 fft_len 对 data 做 overlapping FFT
    % 然后对 data 做下采样操作，使采样率减半，再重复上一条操作。
    % 循环直到 data 的长度小于 fft_len
    iterCount = 0;
    iterationInfo = struct();

    currentValidMaxFreq = inf;  % 记录可信的最大频率范围，第一次迭代所有频率都是可信的
                                % 之后的迭代中由于用了抗混叠滤波器，transition band中的频谱是不可信的
                                
    % The normalized anti-alias filter is identical at every /2 stage,
    % therefore design it once. This both avoids repeated work and ensures
    % that the compensation uses exactly the same FIR coefficients as the
    % actual decimator.
    [hAntiAlias, filterInfo] = designAntiAliasFilter( antiAliasPassbandRatio, antiAliasStopbandAtten);

    while length(data) >= fft_len
        iterCount = iterCount + 1;
        fprintf('Iter %d: ', iterCount);
        dataN = length(data);
        deltaN = max(minGapLen, ceil(minGapTime*saRate));   % overlapping的数据间隔
        averageNo = 1 + floor((dataN-fft_len)/deltaN);
        % 计算 overlapping FFT
        counter = 0;
        startN = 1;
        ps_summer = 0;
        strCnt = 0;
        while (startN+fft_len-1 <= dataN)
            [freq, spectrum, enbw] = fft_easySpectrum(data(startN:startN+fft_len-1), saRate, p.Results.fft_window_name);

            % ----- 挑选非负频部分的频谱，并乘2转换成单边谱 ----------
            posMask  = (freq >= 0);
            freq_pos = reshape(freq(posMask), 1, []);
            ps_seg   = 2*reshape(abs(spectrum(posMask)).^2, 1, []);
            % Interior positive-frequency bins are doubled; DC bin must NOT be doubled.
            ps_seg(1) = ps_seg(1) / 2;
            % -------------------------------------------------------
            % sum 
            ps_summer = ps_summer + ps_seg;
            counter = counter + 1;
            startN = startN + deltaN;
            % print progress
            newStr = sprintf('%d/%d', counter, averageNo);
            fprintf([repmat('\b',1,strCnt) newStr]);
            strCnt = length(newStr);
        end
        fprintf('\n');
        ps_summer = ps_summer / counter;  % overlapping平均功率谱
        ls_summer = sqrt(ps_summer);
        freq = freq_pos;
        % 整合数据
        % 记录新增的频率点，并丢弃掉旧数据中重复频段的点
        % === spectrum stitching uses only the trustworthy
        % passband of decimated data. The transition-band / near-Nyquist
        % portion of this iteration is discarded, and the previous
        % higher-sample-rate spectrum is retained there instead.
        if iterCount == 1
            validMaxFreq = currentValidMaxFreq;
            useCurrent = true(size(freq));
            stitchFreq = freq(end);
            useOld = false(size(freq_full));
        else
            validMaxFreq = min(currentValidMaxFreq, freq(end));
            useCurrent = (freq <= validMaxFreq);
            stitchFreq = freq(find(useCurrent, 1, 'last'));
            useOld = (freq_full > stitchFreq);
        end
        freq_used = freq(useCurrent);
        % ================================================================
        % ====== cumulative FIR response compensation ====================
        % ================================================================
        % Iteration k has passed through (k-1) anti-alias filters.  Because
        % the sampling rate halves after each filter, a physical frequency
        % f is seen at a different normalized frequency by each previous
        % FIR.  Therefore the correct response is the PRODUCT of all stages,
        % not merely the response of the most recent filter.
        if iterCount > 1
            Hcum = cumulativeAntiAliasGain(freq_used, saRate, hAntiAlias, iterCount - 1);
            if any(~isfinite(Hcum)) || any(Hcum <= 0)
                error('Invalid cumulative anti-alias filter response in the trusted passband.');
            end
            ls_summer(useCurrent) = ls_summer(useCurrent) ./ Hcum;
        end
        
        freq_full            = [freq_used,                 freq_full(useOld)];
        spectrum_full        = [ls_summer(useCurrent),            spectrum_full(useOld)];
        spectrumDensity_full = [ls_summer(useCurrent)/sqrt(enbw), spectrumDensity_full(useOld)];


        % 记录迭代信息
        iterationInfo(iterCount).iterID = iterCount;
        iterationInfo(iterCount).saRate = saRate;
        iterationInfo(iterCount).data_len = dataN;
        iterationInfo(iterCount).gap_len = deltaN;
        iterationInfo(iterCount).averageNo = averageNo;
        iterationInfo(iterCount).valid_max_freq = validMaxFreq;
        iterationInfo(iterCount).computed_freq_range = [freq(2),      freq(end)];
        iterationInfo(iterCount).freq_range          = [freq_used(2), freq_used(end)];   % 此次迭代生成的频谱 所覆盖的频率范围
        if iterCount > 1   % 更新上一次迭代的频率范围
            iterationInfo(iterCount-1).freq_range(1) = freq_full(length(freq_used)+1);
        end
        % === frequency-response diagnostics ===
        % Diagnostics for detecting whether a visible spectral step can be
        % explained by residual filter gain mismatch.
        if iterCount > 1
            trustedGain = Hcum;
            minCumGain = min(trustedGain);
            maxCumGain = max(trustedGain);
            maxAmpCorrection = max(1 ./ trustedGain);
            minAmpCorrection = min(1 ./ trustedGain);
        else
            minCumGain = NaN;
            maxCumGain = NaN;
            maxAmpCorrection = NaN;
            minAmpCorrection = NaN;
        end
        iterationInfo(iterCount).n_decimation_stages = iterCount - 1;
        iterationInfo(iterCount).cum_gain_min_trusted = minCumGain;
        iterationInfo(iterCount).cum_gain_max_trusted = maxCumGain;
        iterationInfo(iterCount).amp_correction_min_trusted = minAmpCorrection;
        iterationInfo(iterCount).amp_correction_max_trusted = maxAmpCorrection;
        
        % 降低采样率
        % === proper anti-alias filtering BEFORE downsampling. ===
        % Old code:
        %   data = reshape(data(1:floor(dataN/2)*2), 2, []);
        %   data = mean(data, 1);
        % This is only a 2-tap boxcar LPF and gives poor stopband rejection.
        % High-frequency tones can therefore fold into the low-frequency band.
        evenN = floor(dataN/2)*2;
        data = antiAliasDownsampleBy2(data(1:evenN), hAntiAlias);
        saRate = saRate / 2;

        % The next iteration is trusted only inside the FIR passband.
        currentValidMaxFreq = antiAliasPassbandRatio * (saRate/2);
    end
    
    % pack results
    spectrumData.freq = freq_full;
    spectrumData.spectrum = spectrum_full;
    spectrumData.spectrumDensity = spectrumDensity_full;
    spectrumData.fft_len = fft_len;
    spectrumData.inputs = p.Results;
    spectrumData.inputs = rmfield(spectrumData.inputs, 'data');
    spectrumData.iterationInfo = iterationInfo;
    % anti-alias filter
    spectrumData.antiAliasFilter.coefficients = hAntiAlias;
    spectrumData.antiAliasFilter.order = filterInfo.order;
    spectrumData.antiAliasFilter.beta = filterInfo.beta;
    spectrumData.antiAliasFilter.cutoffNormalizedToOldNyquist = filterInfo.cutoff;
    spectrumData.antiAliasFilter.passbandEdgeNormalizedToOldNyquist = filterInfo.fp;
    spectrumData.antiAliasFilter.stopbandEdgeNormalizedToOldNyquist = filterInfo.fsb;
end



% NOTE: fir1() and kaiser() require Signal Processing Toolbox.
function [h, info] = designAntiAliasFilter(passbandRatio, stopbandAtten)
%DESIGNANTIALIASFILTER Design one normalized FIR used at every /2 stage.
%
% Frequencies are normalized to the OLD Nyquist frequency:
%   new Nyquist  = 0.5
%   passband edge = 0.5 * passbandRatio
%   stopband edge = 0.5

    fp = 0.5 * passbandRatio;
    fsb = 0.5;
    transition = fsb - fp;

    deltaOmega = pi * transition;
    filtOrder = ceil((stopbandAtten - 8) / (2.285 * deltaOmega));
    filtOrder = max(filtOrder, 8);
    if mod(filtOrder, 2) ~= 0
        filtOrder = filtOrder + 1;  % even order => integer group delay
    end

    if stopbandAtten > 50
        beta = 0.1102 * (stopbandAtten - 8.7);
    elseif stopbandAtten >= 21
        beta = 0.5842*(stopbandAtten - 21)^0.4 + ...
               0.07886*(stopbandAtten - 21);
    else
        beta = 0;
    end

    cutoff = (fp + fsb) / 2;
    h = fir1(filtOrder, cutoff, 'low', ...
        kaiser(filtOrder+1, beta), 'scale');
    h = reshape(h, 1, []);
    h = h / sum(h);   % exact unity DC gain

    info.order = filtOrder;
    info.beta = beta;
    info.cutoff = cutoff;
    info.fp = fp;
    info.fsb = fsb;
end


function y = antiAliasDownsampleBy2(x, h)
%ANTIALIASDOWNSAMPLEBY2 Apply the exact FIR h and then decimate by 2.

    x = reshape(x, 1, []);
    filtOrder = numel(h) - 1;
    padLen = filtOrder/2;

    if numel(x) > padLen + 1
        leftPad = x((padLen+1):-1:2);
        rightPad = x((end-1):-1:(end-padLen));
        xPad = [leftPad, x, rightPad];
        xFiltPad = conv(xPad, h, 'same');
        xFilt = xFiltPad(padLen+1 : padLen+numel(x));
    else
        xFilt = conv(x, h, 'same');
    end

    y = xFilt(1:2:end);
end


function Hcum = cumulativeAntiAliasGain(freq, currentSaRate, h, nStages)
%CUMULATIVEANTIALIASGAIN Magnitude response of all previous FIR stages.
%
% At the current sampling rate Fs, the most recent FIR saw the physical
% frequency f with digital angular frequency
%       w0 = pi*f/Fs.
% The preceding FIR saw half that normalized frequency, then one quarter,
% etc. Therefore
%       Hcum(f) = product_{m=0}^{nStages-1} |H(pi*f/(Fs*2^m))|.

    freq = reshape(freq, 1, []);
    Hcum = ones(size(freq));

    for m = 0:(nStages-1)
        omega = pi * freq / (currentSaRate * 2^m);
        Hstage = freqz(h, 1, omega);
        Hcum = Hcum .* abs(reshape(Hstage, 1, []));
    end
end
