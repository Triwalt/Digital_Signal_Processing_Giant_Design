% run_experiments.m: CMA均衡器参数影响分析实验
%
% 功能:
% 1. 步长因子mu的影响
% 2. 滤波器长度M的影响  
% 3. 信噪比SNR性能曲线
%
% 参考: 技术规格书第5章

clear; clc; close all;

fprintf('========================================\n');
fprintf('CMA均衡器参数影响分析实验\n');
fprintf('========================================\n\n');

%% 实验1: 步长因子mu的影响 (技术规格书 5.1节)

fprintf('实验1: 步长因子mu的影响\n');
fprintf('----------------------------------------\n');

% 固定参数
fixed_config = struct();
fixed_config.numSymbols = 10000;
fixed_config.snr_dB = 25;
fixed_config.M = 4;
fixed_config.grayMap = [0 1 3 2];
fixed_config.channelTaps = [1, 0.5, 0.2];
fixed_config.cma.filterLength = 31;
fixed_config.cma.R2 = 1;
fixed_config.fft.nfft = 512;

% 测试不同的步长 (归一化LMS需要较大步长)
mu_values = [0.001, 0.005, 0.01, 0.05];
results_mu = cell(length(mu_values), 1);

% 定义4QAM符号对应的颜色
color_map = [
    0.8, 0.2, 0.2;  % 符号0: 红色
    0.2, 0.8, 0.2;  % 符号1: 绿色
    0.2, 0.2, 0.8;  % 符号2: 蓝色
    0.9, 0.7, 0.1   % 符号3: 黄色
];

figure('Position', [100, 100, 1400, 400]);
for idx = 1:length(mu_values)
    mu = mu_values(idx);
    fprintf('  测试 mu = %.1e ... ', mu);
    
    % 运行仿真
    result = run_single_simulation(fixed_config, mu, fixed_config.cma.filterLength);
    results_mu{idx} = result;
    
    fprintf('BER = %.2e, EVM = %.2f%%\n', result.ber, result.evm);
    
    % 绘制星座图 (按符号着色)
    subplot(1, length(mu_values), idx);
    hold on; grid on; axis equal;
    
    % 计算符号延迟
    tx_symbol_delay = floor((result.best_delay + result.transient*2) / 2);
    num_plot = min(2000, length(result.equalized));
    
    % 按符号绘制不同颜色的点
    for sym_idx = 0:3
        tx_positions = find(result.grayIndices == sym_idx);
        if ~isempty(tx_positions)
            eq_idx = tx_positions - tx_symbol_delay;
            valid_mask = (eq_idx >= 1) & (eq_idx <= length(result.equalized));
            eq_idx = eq_idx(valid_mask);
            eq_idx = eq_idx(eq_idx <= num_plot);
            
            if ~isempty(eq_idx)
                plot(real(result.equalized(eq_idx)), imag(result.equalized(eq_idx)), '.', ...
                     'Color', color_map(sym_idx+1, :), 'MarkerSize', 3);
            end
        end
    end
    
    % 绘制理想星座点
    constellation = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);
    for sym_idx = 0:3
        plot(real(constellation(sym_idx+1)), imag(constellation(sym_idx+1)), ...
             'o', 'MarkerSize', 10, 'LineWidth', 2, ...
             'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_map(sym_idx+1, :));
    end
    
    xlim([-1.5 1.5]); ylim([-1.5 1.5]);
    title(sprintf('\\mu = %.1e\nBER=%.2e', mu, result.ber));
    xlabel('I'); ylabel('Q');
end
sgtitle('实验1: 步长因子mu的影响', 'FontSize', 12, 'FontWeight', 'bold');

fprintf('\n');

%% 实验2: 滤波器长度M的影响 (技术规格书 5.2节)

fprintf('实验2: 滤波器长度M的影响\n');
fprintf('----------------------------------------\n');

% 使用实验1中较优的步长
optimal_mu = 5e-4;

% 测试不同的滤波器长度
M_values = [11, 21, 31, 51];
results_M = cell(length(M_values), 1);

figure('Position', [100, 100, 1400, 400]);
for idx = 1:length(M_values)
    M = M_values(idx);
    fprintf('  测试 M = %d ... ', M);
    
    % 运行仿真
    result = run_single_simulation(fixed_config, optimal_mu, M);
    results_M{idx} = result;
    
    fprintf('BER = %.2e, EVM = %.2f%%\n', result.ber, result.evm);
    
    % 绘制星座图 (按符号着色)
    subplot(1, length(M_values), idx);
    hold on; grid on; axis equal;
    
    tx_symbol_delay = floor((result.best_delay + result.transient*2) / 2);
    num_plot = min(2000, length(result.equalized));
    
    for sym_idx = 0:3
        tx_positions = find(result.grayIndices == sym_idx);
        if ~isempty(tx_positions)
            eq_idx = tx_positions - tx_symbol_delay;
            valid_mask = (eq_idx >= 1) & (eq_idx <= length(result.equalized));
            eq_idx = eq_idx(valid_mask);
            eq_idx = eq_idx(eq_idx <= num_plot);
            
            if ~isempty(eq_idx)
                plot(real(result.equalized(eq_idx)), imag(result.equalized(eq_idx)), '.', ...
                     'Color', color_map(sym_idx+1, :), 'MarkerSize', 3);
            end
        end
    end
    
    constellation = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);
    for sym_idx = 0:3
        plot(real(constellation(sym_idx+1)), imag(constellation(sym_idx+1)), ...
             'o', 'MarkerSize', 10, 'LineWidth', 2, ...
             'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_map(sym_idx+1, :));
    end
    
    xlim([-1.5 1.5]); ylim([-1.5 1.5]);
    title(sprintf('M = %d\nBER=%.2e', M, result.ber));
    xlabel('I'); ylabel('Q');
end
sgtitle('实验2: 滤波器长度M的影响', 'FontSize', 12, 'FontWeight', 'bold');

fprintf('\n');

%% 实验3: 信噪比SNR性能曲线 (技术规格书 5.3节)

fprintf('实验3: 信噪比SNR性能曲线\n');
fprintf('----------------------------------------\n');

% 使用较优参数
optimal_M = 31;

% SNR范围
snr_range = 10:2:30;
ber_with_eq = zeros(size(snr_range));
ber_without_eq = zeros(size(snr_range));

for idx = 1:length(snr_range)
    snr = snr_range(idx);
    fprintf('  测试 SNR = %d dB ... ', snr);
    
    % 有均衡器
    test_config = fixed_config;
    test_config.snr_dB = snr;
    result = run_single_simulation(test_config, optimal_mu, optimal_M);
    ber_with_eq(idx) = result.ber;
    
    % 无均衡器 (直接解调)
    ber_no_eq = simulate_without_equalizer(test_config);
    ber_without_eq(idx) = ber_no_eq;
    
    fprintf('有均衡: %.2e, 无均衡: %.2e\n', ber_with_eq(idx), ber_no_eq);
end

% 绘制BER曲线
figure('Position', [100, 100, 800, 600]);
semilogy(snr_range, ber_with_eq, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
semilogy(snr_range, ber_without_eq, 'r--s', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('信噪比 SNR (dB)', 'FontSize', 12);
ylabel('误码率 BER', 'FontSize', 12);
title('实验3: BER vs SNR 性能对比', 'FontSize', 14, 'FontWeight', 'bold');
legend('有CMA均衡器', '无均衡器', 'Location', 'best', 'FontSize', 11);
ylim([1e-5 1]);

fprintf('\n========================================\n');
fprintf('所有实验完成!\n');
fprintf('========================================\n');


%% 辅助函数

function result = run_single_simulation(config, mu, filter_length)
    % 运行单次仿真并返回结果
    
    % 生成发射信号
    numBits = config.numSymbols * log2(config.M);
    txBits = randi([0 1], 1, numBits);
    
    txSymbolIndices = zeros(1, config.numSymbols);
    for i = 1:config.numSymbols
        bit1 = txBits(2*i-1);
        bit2 = txBits(2*i);
        txSymbolIndices(i) = bit1 * 2 + bit2;
    end
    
    grayIndices = config.grayMap(txSymbolIndices + 1);
    constellation_map = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);
    txSymbols = constellation_map(grayIndices + 1);
    
    % 信道
    rxSymbols_ISI = filter(config.channelTaps, 1, txSymbols);
    rxSymbols = awgn(rxSymbols_ISI, config.snr_dB, 'measured');
    
    % CMA均衡 - 样本级处理 (归一化LMS)
    M = filter_length;
    R2 = config.cma.R2;
    
    w = zeros(M, 1);
    w(floor(M/2) + 1) = 1;
    
    rxSymbols_col = rxSymbols(:);
    N_total = length(rxSymbols_col);
    equalizedSymbols = zeros(N_total, 1);
    input_buffer = [zeros(M-1, 1); rxSymbols_col];
    
    for n = 1:N_total
        x_vec = input_buffer(n:n+M-1);
        
        % 归一化 (防止权重发散)
        x_power = real(x_vec' * x_vec) + 1e-6;
        
        y_out = w' * x_vec;
        equalizedSymbols(n) = y_out;
        
        % CMA权重更新 (归一化LMS)
        error_signal = (R2 - abs(y_out)^2) * y_out;
        w = w + (mu / x_power) * conj(error_signal) * x_vec;
    end
    
    % 相位校正
    transient_samples = 1000;
    equalizedSymbols_steady = equalizedSymbols(transient_samples+1:end);
    
    pilot_length = min(500, length(equalizedSymbols_steady));
    pilot_symbols = equalizedSymbols_steady(1:pilot_length);
    
    % 修正: 对每个符号找最近的星座点
    pilot_reference = zeros(size(pilot_symbols));
    for k = 1:length(pilot_symbols)
        [~, idx] = min(abs(pilot_symbols(k) - constellation_map));
        pilot_reference(k) = constellation_map(idx);
    end
    
    phase_offset = angle(sum(pilot_reference .* conj(pilot_symbols)));
    equalizedSymbols_corrected = equalizedSymbols_steady * exp(1j * phase_offset);
    
    % 解调
    % 星座点映射: 索引0->(1+j), 1->(1-j), 2->(-1-1j), 3->(-1+1j)
    % 这些索引对应格雷码: [0, 1, 2, 3]
    % 需要反向映射到原始比特索引
    rxSymbolIndices = zeros(length(equalizedSymbols_corrected), 1);
    for k = 1:length(equalizedSymbols_corrected)
        [~, rxSymbolIndices(k)] = min(abs(equalizedSymbols_corrected(k) - constellation_map));
    end
    
    % rxSymbolIndices现在是1-4,对应格雷码索引0-3
    rxGrayIndices = rxSymbolIndices - 1; % 格雷码索引 0-3
    
    % 格雷码反向映射: 找到格雷码在grayMap中的位置,该位置就是原始比特索引
    rxSymbolIndices_original = zeros(size(rxGrayIndices));
    for k = 1:length(rxGrayIndices)
        pos = find(config.grayMap == rxGrayIndices(k));
        if ~isempty(pos)
            rxSymbolIndices_original(k) = pos - 1; % 转换为0-3
        end
    end
    
    rxBits = zeros(1, length(rxSymbolIndices_original) * 2);
    for i = 1:length(rxSymbolIndices_original)
        idx = rxSymbolIndices_original(i);
        rxBits(2*i-1) = bitshift(idx, -1);
        rxBits(2*i) = bitand(idx, 1);
    end
    
    % 计算BER - 使用延迟搜索
    max_search_delay = min(2000, length(txBits) - 100);  % 扩大搜索范围
    best_ber = 1.0;
    best_delay = 0;
    
    for delay_bits = 0:2:max_search_delay
        tx_start = delay_bits + 1;
        rx_start = 1;
        test_length = min(length(txBits) - tx_start + 1, length(rxBits) - rx_start + 1);
        
        if test_length < 1000
            continue;
        end
        
        test_length = min(test_length, 3000);  % 使用更多比特测试
        tx_test = txBits(tx_start : tx_start + test_length - 1);
        rx_test = rxBits(rx_start : rx_start + test_length - 1);
        
        test_errors = sum(tx_test ~= rx_test);
        test_ber = test_errors / test_length;
        
        if test_ber < best_ber
            best_ber = test_ber;
            best_delay = delay_bits;
        end
    end
    
    % 使用最佳延迟对齐
    tx_align_start = best_delay + 1;
    rx_align_start = 1;
    compare_length = min(length(txBits) - tx_align_start + 1, length(rxBits) - rx_align_start + 1);
    
    txBits_aligned = txBits(tx_align_start : tx_align_start + compare_length - 1);
    rxBits_aligned = rxBits(rx_align_start : rx_align_start + compare_length - 1);
    
    bit_errors = sum(txBits_aligned ~= rxBits_aligned);
    ber = bit_errors / compare_length;
    
    % 计算EVM
    tx_symbol_delay = floor((best_delay + transient_samples*2) / 2);
    tx_symbol_start = tx_symbol_delay + 1;
    tx_symbol_end = min(tx_symbol_start + length(equalizedSymbols_corrected) - 1, length(txSymbols));
    txSymbols_ref = txSymbols(tx_symbol_start:tx_symbol_end);
    eq_symbols_trim = equalizedSymbols_corrected(1:(tx_symbol_end - tx_symbol_start + 1));
    
    alpha = (eq_symbols_trim' * txSymbols_ref(:)) / (eq_symbols_trim' * eq_symbols_trim);
    equalizedSymbols_aligned = alpha * eq_symbols_trim;
    
    error_vector = equalizedSymbols_aligned - txSymbols_ref(:);
    evm = sqrt(mean(abs(error_vector).^2)) / sqrt(mean(abs(txSymbols_ref).^2)) * 100;
    
    % 返回结果
    result = struct();
    result.ber = ber;
    result.evm = evm;
    result.equalized = equalizedSymbols_corrected(1:min(2000, length(equalizedSymbols_corrected)));
    result.grayIndices = grayIndices;  % 用于颜色标记
    result.best_delay = best_delay;    % 延迟信息
    result.transient = transient_samples;
end


function ber = simulate_without_equalizer(config)
    % 无均衡器的BER仿真
    
    numBits = config.numSymbols * log2(config.M);
    txBits = randi([0 1], 1, numBits);
    
    txSymbolIndices = zeros(1, config.numSymbols);
    for i = 1:config.numSymbols
        bit1 = txBits(2*i-1);
        bit2 = txBits(2*i);
        txSymbolIndices(i) = bit1 * 2 + bit2;
    end
    
    grayIndices = config.grayMap(txSymbolIndices + 1);
    constellation_map = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);
    txSymbols = constellation_map(grayIndices + 1);
    
    rxSymbols_ISI = filter(config.channelTaps, 1, txSymbols);
    rxSymbols = awgn(rxSymbols_ISI, config.snr_dB, 'measured');
    
    % 直接解调
    transient = 100;
    rxSymbols_demod = rxSymbols(transient+1:end);
    
    rxSymbolIndices = zeros(length(rxSymbols_demod), 1);
    for k = 1:length(rxSymbols_demod)
        [~, rxSymbolIndices(k)] = min(abs(rxSymbols_demod(k) - constellation_map));
    end
    rxGrayIndices = rxSymbolIndices - 1;
    
    % 格雷码反向映射
    rxSymbolIndices_original = zeros(size(rxGrayIndices));
    for k = 1:length(rxGrayIndices)
        pos = find(config.grayMap == rxGrayIndices(k));
        if ~isempty(pos)
            rxSymbolIndices_original(k) = pos - 1;
        end
    end
    
    rxBits = zeros(1, length(rxSymbolIndices_original) * 2);
    for i = 1:length(rxSymbolIndices_original)
        idx = rxSymbolIndices_original(i);
        rxBits(2*i-1) = bitshift(idx, -1);
        rxBits(2*i) = bitand(idx, 1);
    end
    
    % BER - 使用延迟搜索
    max_delay = min(500, length(txBits) - 100);
    best_ber = 1.0;
    best_delay = 0;
    
    for delay_bits = 0:2:max_delay
        tx_start = delay_bits + 1;
        test_len = min(length(txBits) - tx_start + 1, length(rxBits));
        if test_len < 500
            continue;
        end
        test_len = min(test_len, 2000);
        tx_test = txBits(tx_start : tx_start + test_len - 1);
        rx_test = rxBits(1:test_len);
        test_err = sum(tx_test ~= rx_test);
        test_ber = test_err / test_len;
        if test_ber < best_ber
            best_ber = test_ber;
            best_delay = delay_bits;
        end
    end
    
    tx_start = best_delay + 1;
    compare_length = min(length(txBits) - tx_start + 1, length(rxBits));
    
    txBits_aligned = txBits(tx_start : tx_start + compare_length - 1);
    rxBits_aligned = rxBits(1:compare_length);
    
    bit_errors = sum(txBits_aligned ~= rxBits_aligned);
    ber = bit_errors / compare_length;
end