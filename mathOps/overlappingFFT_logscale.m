function [freq_full, spectrum_full, spectrumDensity_full, fft_len] = overlappingFFT_logscale(data, saRate, fft_len, varargin)
%OVERLAPPINGFFT_LOGSCALE overlapping FFT with logscale frequency list
% Author: GeorgeGate @ 2022-12-07
% Version: 1.0
% 计算PSD的平均值，使用overlapping FFT方式进行平均。频域上的取点满足近似的 logspace，可以有效地平滑高频段的频谱。
% 
% average method: mean(PSD)
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
    
    p.parse(data, saRate, fft_len, varargin{:});
    
    minGapTime = p.Results.minGapTime;
    minGapLen = p.Results.minGapLen;
        
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
    % 然后对 data 做平均操作，使采样率减半，再重复上一条操作。
    % 循环直到 data 的长度小于 fft_len
    iterCount = 0;
    while length(data) >= fft_len
        iterCount = iterCount+1;
        fprintf('Iter %d: ', iterCount);
        dataN = length(data);
        deltaN = max(minGapLen, ceil(minGapTime*saRate));   % overlapping的数据间隔
        % 计算 overlapping FFT
        counter = 0;
        startN = 1;
        ps_summer = 0;
        strCnt = 0;
        while (startN+fft_len-1 <= dataN)
            [freq, spectrum, enbw] = fft_easySpectrum(data(startN:startN+fft_len-1), saRate, p.Results.fft_window_name);
            ps_summer = ps_summer + 2*abs(spectrum(freq>=0)).^2;  % 只保留正频部分，且乘2转换成单边谱
            counter = counter+1;
            startN = startN + deltaN;
            % print progress
            newStr = sprintf('%d/%d', counter, 1+floor((dataN-fft_len)/deltaN));
            fprintf([repmat('\b',1,strCnt) newStr]);
            strCnt = length(newStr);
        end
        fprintf('\n');
        ps_summer = ps_summer/counter;  % overlapping平均功率谱
        ls_summer = sqrt(ps_summer);
        freq = freq(freq>=0);
        % 整合数据
        % 记录新增的频率点，并丢弃掉旧数据中的重复频率点
        id1 = ceil(length(freq)/2)+1;
        freq_full            = [freq,                 freq_full(id1:end)];
        spectrum_full        = [ls_summer,            spectrum_full(id1:end)];
        spectrumDensity_full = [ls_summer/sqrt(enbw), spectrumDensity_full(id1:end)];
        
        % 降低采样率
        data = reshape(data(1:floor(dataN/2)*2), 2, []);
        data = mean(data, 1);
        saRate = saRate / 2;
    end
end

