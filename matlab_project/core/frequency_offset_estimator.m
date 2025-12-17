function [freq_offset_est, phase_inst, diagnostics] = frequency_offset_estimator(rx_signal, fs, config)
% FREQUENCY_OFFSET_ESTIMATOR 基于四次方差分相位法的频偏估计算法
%
% 适用于 QPSK/4QAM 等恒模调制信号的载波频率偏移估计
%
% 算法原理:
%   1. 四次方法消除调制信息: z[n] = r[n]^4
%      对于QPSK信号, 四次方后符号映射到同一点, 仅保留载波相位
%   2. 相位差分估计频偏: Δφ = angle(z[n] * conj(z[n-L]))
%      频偏估计: Δf = Δφ / (4 * 2π * L * Ts)
%   3. 可选的滑动平均滤波以提高估计精度
%
% 输入参数:
%   rx_signal - 接收信号 (复数列向量或行向量)
%   fs        - 符号率 (Hz), 用于将相位差转换为频率
%   config    - 配置结构体 (可选), 包含:
%               .L           : 差分延迟 (默认: 1)
%               .avg_window  : 平均窗口长度 (默认: 64)
%               .mod_order   : 调制阶数 (默认: 4, 用于QPSK/4QAM)
%               .method      : 'differential' | 'fft' (默认: 'differential')
%
% 输出参数:
%   freq_offset_est - 估计的频率偏移 (Hz)
%   phase_inst      - 瞬时相位序列 (rad), 用于后续处理
%   diagnostics     - 诊断信息结构体
%
% 参考:
%   - Mengali & D'Andrea, "Synchronization Techniques for Digital Receivers"
%   - 数字信号处理训练文档: 频偏估计章节
%
% 作者: DSP课程设计项目
% 日期: 2025年12月

    %% 输入验证与默认配置
    if nargin < 3 || isempty(config)
        config = struct();
    end
    
    % 默认参数
    if ~isfield(config, 'L')
        config.L = 1;  % 差分延迟
    end
    if ~isfield(config, 'avg_window')
        config.avg_window = 64;  % 平均窗口
    end
    if ~isfield(config, 'mod_order')
        config.mod_order = 4;  % QPSK/4QAM
    end
    if ~isfield(config, 'method')
        config.method = 'differential';
    end
    
    % 确保为列向量
    rx_signal = rx_signal(:);
    N = length(rx_signal);
    
    if N < config.L + config.avg_window
        error('信号长度不足: 需要至少 %d 个样本', config.L + config.avg_window);
    end
    
    %% 步骤1: M次方消除调制信息
    % 对于M-PSK/M-QAM信号, M次方后调制信息被消除
    M = config.mod_order;
    z = rx_signal .^ M;
    
    %% 步骤2: 频偏估计
    switch lower(config.method)
        case 'differential'
            % 差分相位法
            [freq_offset_est, phase_diff_avg, diagnostics] = ...
                differential_phase_method(z, fs, M, config);
            
        case 'fft'
            % 基于FFT的频偏估计 (适用于大频偏)
            [freq_offset_est, phase_diff_avg, diagnostics] = ...
                fft_based_method(z, fs, M, config);
            
        otherwise
            error('未知的估计方法: %s', config.method);
    end
    
    %% 步骤3: 计算瞬时相位序列
    % 从四次方信号中提取相位, 再除以M恢复原始相位
    phase_inst = angle(z) / M;
    
    %% 填充诊断信息
    diagnostics.config = config;
    diagnostics.signal_length = N;
    diagnostics.mod_order = M;
    diagnostics.method = config.method;
    diagnostics.z_power_mean = mean(abs(z).^2);
    
end

%% ===================== 内部函数 =====================

function [freq_est, phase_diff_avg, diag] = differential_phase_method(z, fs, M, config)
% 差分相位法频偏估计
%
% 原理: Δφ = angle(z[n] * conj(z[n-L]))
%       由于 z = r^M, 相位被放大M倍
%       故 Δf = Δφ / (M * 2π * L * Ts)
%
% 改进: 使用大延迟L和全局累加以提高小频偏估计精度

    N = length(z);
    
    % 对于高符号率系统，使用较大的延迟L以积累足够相位差
    % 默认 L = 1 仅适用于低符号率或大频偏场景
    L = config.L;
    
    % 自适应延迟: 如果符号率很高，增大L以提高估计精度
    % 目标: L * Δf / fs 产生可测量的相位差 (至少 1e-4 rad)
    if L == 1 && fs > 1e9
        % 对于高速系统，使用更大的延迟累积相位
        L = min(1000, floor(N / 10));  % 最多使用1000个符号延迟
    end
    
    % 方法1: 全局累加法 (更适合小频偏)
    % 使用整个序列计算平均相位变化率
    
    % 将序列分成若干块，计算块间相位差
    block_size = min(config.avg_window, floor(N / 4));
    num_blocks = floor(N / block_size);
    
    if num_blocks >= 2
        % 累加每个块的 z 值
        z_blocks = zeros(num_blocks, 1);
        for blk = 1:num_blocks
            idx_start = (blk - 1) * block_size + 1;
            idx_end = blk * block_size;
            z_blocks(blk) = sum(z(idx_start:idx_end));
        end
        
        % 计算相邻块之间的相位差
        phase_diffs = zeros(num_blocks - 1, 1);
        for k = 1:num_blocks - 1
            phase_diffs(k) = angle(z_blocks(k+1) * conj(z_blocks(k)));
        end
        
        % 平均相位差 (每 block_size 个符号的相位变化)
        phase_diff_per_block = mean(phase_diffs);
        
        % 转换为每符号的相位变化
        phase_diff_per_symbol = phase_diff_per_block / block_size;
        
        % 频率偏移估计
        % phase_diff_per_symbol = M * 2π * Δf / fs
        freq_est = phase_diff_per_symbol * fs / (M * 2 * pi);
        
        phase_diff_avg = phase_diff_per_symbol;
        
        % 诊断信息
        diag.method = 'block_accumulation';
        diag.phase_diff = phase_diffs;
        diag.phase_diff_avg = phase_diff_avg;
        diag.phase_diff_std = std(phase_diffs) / block_size;
        diag.block_size = block_size;
        diag.num_blocks = num_blocks;
        diag.z_diff_power = mean(abs(z_blocks).^2);
    else
        % 回退到简单差分法
        z_delayed = [zeros(L, 1); z(1:end-L)];
        z_diff = z .* conj(z_delayed);
        z_diff_valid = z_diff(L+1:end);
        
        % 全局累加 (而非滑动平均)
        z_sum = sum(z_diff_valid);
        phase_diff_avg = angle(z_sum) / L;
        
        freq_est = phase_diff_avg * fs / (M * 2 * pi);
        
        diag.method = 'simple_differential';
        diag.phase_diff = angle(z_diff_valid);
        diag.phase_diff_avg = phase_diff_avg;
        diag.phase_diff_std = std(angle(z_diff_valid));
        diag.z_diff_power = mean(abs(z_diff_valid).^2);
    end
    
end

function [freq_est, phase_diff_avg, diag] = fft_based_method(z, fs, M, config)
% 基于FFT的频偏估计
%
% 原理: z = r^M 的频谱在 M*Δf 处有明显峰值
%       通过FFT找到峰值位置即可估计频偏

    N = length(z);
    nfft = 2^nextpow2(N * 4);  % 4倍过采样以提高频率分辨率
    
    % 应用窗函数减少频谱泄漏
    win = hamming(N);
    z_windowed = z .* win;
    
    % 计算FFT
    Z = fft(z_windowed, nfft);
    Z_mag = abs(Z);
    
    % 找到峰值位置
    [~, peak_idx] = max(Z_mag(1:nfft/2));
    
    % 将索引转换为频率
    freq_resolution = fs / nfft;
    peak_freq = (peak_idx - 1) * freq_resolution;
    
    % 原始频偏 = peak_freq / M
    freq_est = peak_freq / M;
    
    % 如果频偏超过 fs/2M, 需要进行折叠处理
    if freq_est > fs / (2 * M)
        freq_est = freq_est - fs / M;
    end
    
    % 相位差估计 (与差分法统一输出格式)
    phase_diff_avg = freq_est * M * 2 * pi / fs;
    
    % 诊断信息
    diag.nfft = nfft;
    diag.peak_idx = peak_idx;
    diag.peak_freq = peak_freq;
    diag.freq_resolution = freq_resolution;
    diag.Z_mag = Z_mag(1:min(1024, nfft/2));
    diag.phase_diff_avg = phase_diff_avg;
    diag.phase_diff_std = NaN;  % FFT法不直接计算相位序列
    diag.phase_diff = [];
    
end
