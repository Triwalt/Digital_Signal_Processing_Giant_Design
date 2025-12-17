function [phase_est, signal_corrected, diagnostics] = phase_noise_estimator(rx_signal, config)
% PHASE_NOISE_ESTIMATOR 带相位模糊解决的相位噪声估计算法
%
% 基于 Viterbi-Viterbi (V&V) 算法实现相位噪声估计与补偿
% 解决 M 次方法引入的 2π/M 相位模糊问题
%
% 算法原理:
%   1. M次方法消除调制: z[n] = r[n]^M
%   2. 滑动平均滤波平滑相位: φ_avg = angle(sum(z[n-k:n+k]))
%   3. 相位恢复: θ[n] = φ_avg[n] / M
%   4. 相位解模糊: 使用相位unwrap确保连续性
%   5. 可选: 判决引导细化
%
% 输入参数:
%   rx_signal - 接收信号 (复数列向量或行向量), 通常是均衡后的信号
%   config    - 配置结构体 (可选), 包含:
%               .mod_order      : 调制阶数 (默认: 4, QPSK/4QAM)
%               .win_half       : 滑动窗口半宽 (默认: 5, 总长度 = 2*win_half+1)
%               .constellation  : 参考星座点 (默认: QPSK归一化星座)
%               .phase_unwrap   : 是否进行相位展开 (默认: true)
%               .ambiguity_method: 'unwrap' | 'decision_directed' | 'none'
%
% 输出参数:
%   phase_est        - 估计的相位噪声序列 (rad)
%   signal_corrected - 相位校正后的信号
%   diagnostics      - 诊断信息结构体
%
% 参考:
%   - A. J. Viterbi, "Nonlinear Estimation of PSK-Modulated Carrier Phase"
%
% 作者: DSP课程设计项目
% 日期: 2025年12月

    %% 输入验证与默认配置
    if nargin < 2 || isempty(config)
        config = struct();
    end
    
    % 调制阶数 (用于M次方法)
    if ~isfield(config, 'mod_order')
        config.mod_order = 4;
    end
    M = config.mod_order;
    
    % 滑动窗口半宽
    if ~isfield(config, 'win_half')
        config.win_half = 5;
    end
    
    % 参考星座点
    if ~isfield(config, 'constellation')
        % 默认QPSK/4QAM归一化星座
        config.constellation = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);
    end
    
    % 相位展开开关
    if ~isfield(config, 'phase_unwrap')
        config.phase_unwrap = true;
    end
    
    % 相位模糊解决方法
    if ~isfield(config, 'ambiguity_method')
        config.ambiguity_method = 'unwrap';  % 默认使用unwrap方法
    end
    
    % 判决引导块大小
    if ~isfield(config, 'block_size')
        config.block_size = 100;
    end
    
    % 确保为列向量
    rx_signal = rx_signal(:);
    N = length(rx_signal);
    
    %% 步骤1: M次方法消除调制信息
    % 对于 M-PSK/M-QAM 信号: r[n] = a[n] * exp(j*θ[n])
    % M次方后调制信息被消除，仅保留M倍相位噪声
    z = rx_signal .^ M;
    
    %% 步骤2: 计算星座点固有相位偏移
    % 对于归一化QPSK星座 [(1+j)/sqrt(2), ...], 四次方后角度为π
    % 需要补偿这个固有偏移
    constellation_phase_offset = angle(mean(config.constellation .^ M)) / M;
    
    %% 步骤3: 滑动平均滤波
    win_half = config.win_half;
    win_len = 2 * win_half + 1;
    window = ones(win_len, 1) / win_len;
    
    % 使用卷积实现滑动平均
    z_filtered = conv(z, window, 'same');
    
    %% 步骤4: 提取相位并除以M，减去固有偏移
    phase_raw = angle(z_filtered);
    phase_est_raw = phase_raw / M - constellation_phase_offset;
    
    %% 步骤5: 相位模糊解决
    switch lower(config.ambiguity_method)
        case 'unwrap'
            % 使用相位展开解决模糊 - 这是最稳健的方法
            % 先对 M倍相位 进行 unwrap，再除以 M，减去固有偏移
            phase_unwrapped = unwrap(phase_raw);
            phase_est = phase_unwrapped / M - constellation_phase_offset;
            ambiguity_seq = zeros(N, 1);
            
        case 'decision_directed'
            % 判决引导方法
            [phase_est, ambiguity_seq] = resolve_ambiguity_dd(...
                rx_signal, phase_est_raw, config);
            if config.phase_unwrap
                phase_est = unwrap(phase_est);
            end
            
        case 'differential'
            % 差分方法
            [phase_est, ambiguity_seq] = resolve_ambiguity_diff(phase_est_raw, M);
            if config.phase_unwrap
                phase_est = unwrap(phase_est);
            end
            
        case 'none'
            phase_est = phase_est_raw;
            if config.phase_unwrap
                phase_est = unwrap(phase_est);
            end
            ambiguity_seq = zeros(N, 1);
            
        otherwise
            error('未知的相位模糊解决方法: %s', config.ambiguity_method);
    end
    
    %% 步骤5: 信号相位校正
    signal_corrected = rx_signal .* exp(-1j * phase_est);
    
    %% 步骤6: 初始相位模糊消除 (π/2 模糊)
    % V&V算法会产生 k*π/2 的相位模糊，需要通过判决来消除
    % 测试4个可能的模糊值，选择使EVM最小的那个
    ambiguity_candidates = [0, pi/2, pi, 3*pi/2];
    min_evm = inf;
    best_rotation = 0;
    
    for k = 1:length(ambiguity_candidates)
        rotation = ambiguity_candidates(k);
        test_signal = signal_corrected * exp(-1j * rotation);
        [~, test_evm] = calculate_evm_internal(test_signal, config.constellation);
        if test_evm < min_evm
            min_evm = test_evm;
            best_rotation = rotation;
        end
    end
    
    % 应用最佳旋转
    signal_corrected = signal_corrected * exp(-1j * best_rotation);
    phase_est = phase_est + best_rotation;
    
    %% 填充诊断信息
    diagnostics = struct();
    diagnostics.config = config;
    diagnostics.N = N;
    diagnostics.mod_order = M;
    diagnostics.win_len = win_len;
    diagnostics.phase_raw = phase_raw;
    diagnostics.phase_est_raw = phase_est_raw;
    diagnostics.ambiguity_seq = ambiguity_seq;
    diagnostics.z_filtered_power = mean(abs(z_filtered).^2);
    diagnostics.phase_noise_std = std(diff(phase_est));
    
    % 计算校正前后的EVM
    [~, evm_before] = calculate_evm(rx_signal, config.constellation);
    [~, evm_after] = calculate_evm(signal_corrected, config.constellation);
    diagnostics.evm_before = evm_before;
    diagnostics.evm_after = evm_after;
    
end

%% ===================== 内部函数 =====================

function [phase_resolved, ambiguity_seq] = resolve_ambiguity_dd(rx_signal, phase_raw, config)
% 判决引导相位模糊解决

    M = config.mod_order;
    constellation = config.constellation;
    block_size = config.block_size;
    N = length(rx_signal);
    
    ambiguity_candidates = (0:M-1) * 2 * pi / M;
    
    num_blocks = ceil(N / block_size);
    ambiguity_seq = zeros(N, 1);
    phase_resolved = zeros(N, 1);
    
    prev_ambiguity = 0;
    
    for blk = 1:num_blocks
        idx_start = (blk - 1) * block_size + 1;
        idx_end = min(blk * block_size, N);
        blk_indices = idx_start:idx_end;
        
        sig_blk = rx_signal(blk_indices);
        phase_blk = phase_raw(blk_indices);
        
        min_error = inf;
        best_ambiguity = 0;
        
        for k = 1:M
            amb = ambiguity_candidates(k);
            phase_test = phase_blk + amb;
            sig_corrected = sig_blk .* exp(-1j * phase_test);
            
            [~, error_vec] = hard_decision(sig_corrected, constellation);
            total_error = mean(abs(error_vec).^2);
            
            % 连续性惩罚
            amb_diff = min([abs(amb - prev_ambiguity), ...
                           abs(amb - prev_ambiguity + 2*pi), ...
                           abs(amb - prev_ambiguity - 2*pi)]);
            total_error = total_error + 0.1 * amb_diff;
            
            if total_error < min_error
                min_error = total_error;
                best_ambiguity = amb;
            end
        end
        
        phase_resolved(blk_indices) = phase_blk + best_ambiguity;
        ambiguity_seq(blk_indices) = best_ambiguity;
        prev_ambiguity = best_ambiguity;
    end
    
end

function [phase_resolved, ambiguity_seq] = resolve_ambiguity_diff(phase_est_raw, M)
% 差分相位模糊解决
% 使用相位差来跟踪并解决模糊

    N = length(phase_est_raw);
    ambiguity_seq = zeros(N, 1);
    phase_resolved = zeros(N, 1);
    
    ambiguity_step = 2 * pi / M;
    
    % 初始化第一个点
    phase_resolved(1) = phase_est_raw(1);
    
    for n = 2:N
        % 计算预测相位差
        phase_diff = phase_est_raw(n) - phase_est_raw(n-1);
        
        % 将相位差归一化到 [-π/M, π/M] 范围
        while phase_diff > ambiguity_step / 2
            phase_diff = phase_diff - ambiguity_step;
        end
        while phase_diff < -ambiguity_step / 2
            phase_diff = phase_diff + ambiguity_step;
        end
        
        % 使用修正后的相位差更新
        phase_resolved(n) = phase_resolved(n-1) + phase_diff;
        ambiguity_seq(n) = phase_resolved(n) - phase_est_raw(n);
    end
    
end

function [decided_symbols, error_vec] = hard_decision(signal, constellation)
% 硬判决函数

    N = length(signal);
    decided_symbols = zeros(N, 1);
    error_vec = zeros(N, 1);
    
    for n = 1:N
        distances = abs(signal(n) - constellation);
        [~, min_idx] = min(distances);
        decided_symbols(n) = constellation(min_idx);
        error_vec(n) = signal(n) - decided_symbols(n);
    end
    
end

function [decided_symbols, evm_percent] = calculate_evm(signal, constellation)
% 计算EVM

    [decided_symbols, error_vec] = hard_decision(signal, constellation);
    evm_percent = sqrt(mean(abs(error_vec).^2)) / sqrt(mean(abs(decided_symbols).^2)) * 100;
    
end

function [decided_symbols, evm_percent] = calculate_evm_internal(signal, constellation)
% 内部EVM计算函数 - 用于模糊消除

    [decided_symbols, error_vec] = hard_decision(signal, constellation);
    evm_percent = sqrt(mean(abs(error_vec).^2)) / sqrt(mean(abs(decided_symbols).^2)) * 100;
    
end
