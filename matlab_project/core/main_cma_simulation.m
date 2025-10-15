% main_cma_simulation.m: 4QAM通信系统CMA盲均衡仿真主脚本
%
% 功能:
% 完整的4QAM通信链路仿真,包括:
% 1. 发射端: 信源生成、符号映射、4QAM调制
% 2. 信道: ISI信道、AWGN噪声
% 3. 接收端: CMA盲均衡、相位校正、解调
% 4. 性能评估: BER、EVM、星座图
%
% 参考: 技术规格书第4章

clear; clc; close all;

% 将自建函数目录加入搜索路径
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
customFuncDir = fullfile(scriptDir, '..', 'Gemini_generated_simulation');
if exist(customFuncDir, 'dir')
    addpath(customFuncDir);
else
    warning('自建函数目录未找到: %s', customFuncDir);
end

fprintf('========================================\n');
fprintf('4QAM通信系统CMA盲均衡仿真\n');
fprintf('========================================\n\n');

%% 步骤1: 仿真参数配置 (技术规格书 4.1节)

config = struct();

% 基本参数
config.numSymbols = 30000;          % 仿真符号数量
config.snr_dB = 15;                  % 信噪比 (dB)

% 4QAM调制参数
config.M = 4;                        % 调制阶数
config.grayMap = [0 1 3 2];          % 格雷码映射表 (00->0, 01->1, 11->3, 10->2)

% 信道参数
config.channelTaps = [1, 0.5, 0.2, 0.1, 0.5, 0.7]; % ISI信道冲激响应
config.channelNormalizePower = true; % 是否归一化信道能量
config.channelPhaseOffsetDeg = 30;   % 额外的信道相位偏移 (单位: 度)

% CMA均衡器参数
config.cma.filterLength = 31;        % 均衡器长度 (M)
config.cma.stepSize = 0.01;         % 初始步长 (mu)
config.cma.stepDecay = 1;         % 每遍迭代的步长衰减因子
config.cma.numPasses = 10;            % 遍历接收序列的次数
config.cma.ddStartPass = 5;          % 从第几遍开始切换到判决引导; 0表示禁用
config.cma.R2 = 1;                   % 4QAM的收敛半径 R^2 = 1

% 快速卷积参数
config.conv.blockLength = 256;       % 自建快速卷积块长 (Overlap-Add方法)

fprintf('配置参数:\n');
fprintf('  符号数量: %d\n', config.numSymbols);
fprintf('  SNR: %d dB\n', config.snr_dB);
fprintf('  信道抽头: [%.2f, %.2f, %.2f]\n', config.channelTaps);
fprintf('  信道归一化: %s\n', mat2str(config.channelNormalizePower));
fprintf('  均衡器长度: %d\n', config.cma.filterLength);
fprintf('  CMA步长: %.4f\n', config.cma.stepSize);
fprintf('  信道相位偏移: %.2f 度\n', config.channelPhaseOffsetDeg);
fprintf('  Fast OLA卷积块长: %d\n', config.conv.blockLength);
if config.cma.ddStartPass > 0
    fprintf('  判决引导模式: 第%d遍开始启用\n', config.cma.ddStartPass);
else
    fprintf('  判决引导模式: 禁用\n');
end
fprintf('\n');

%% 步骤2: 发射端实现 (技术规格书 4.2节)

fprintf('发射端处理...\n');

% 2.1 生成随机比特
numBits = config.numSymbols * log2(config.M);
txBits = randi([0 1], 1, numBits);
fprintf('  生成 %d 个比特\n', numBits);

% 2.2 比特到符号映射 (向量化)
bitPairs = reshape(txBits, 2, []).';      % 每行对应一个符号的两位比特
txSymbolIndices = (bitPairs(:, 1) * 2 + bitPairs(:, 2)).';

% 2.3 格雷编码映射
grayIndices = config.grayMap(txSymbolIndices + 1); % MATLAB索引从1开始

% 2.4 4QAM调制
% 4QAM星座点: {(1+j), (1-j), (-1+j), (-1-j)} / sqrt(2)
% 映射: 0->(1+j), 1->(1-j), 2->(-1-j), 3->(-1+j)
constellation_map = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);
txSymbols = constellation_map(grayIndices + 1);

fprintf('  调制 %d 个4QAM符号\n', config.numSymbols);
fprintf('  符号能量: %.4f\n', mean(abs(txSymbols).^2));
fprintf('\n');

%% 步骤3: 信道建模 (技术规格书 4.3节)

fprintf('信道传输...\n');

% 3.0 计算实际信道冲激响应
channelImpulseResponse = config.channelTaps(:).';
if config.channelNormalizePower
    channel_rms = sqrt(sum(abs(channelImpulseResponse).^2));
    if channel_rms > 0
        channelImpulseResponse = channelImpulseResponse / channel_rms;
    end
    fprintf('  信道能量归一化系数: %.4f\n', channel_rms);
end

% 3.1 通过ISI信道 (使用自建快速卷积函数)
blockLength = config.conv.blockLength;
txSymbols_row = txSymbols(:).';
rxSymbols_full = fast_conv_add(txSymbols_row, channelImpulseResponse, blockLength);
rxSymbols_ISI = rxSymbols_full(1:length(txSymbols_row));
fprintf('  使用fast_conv_add模拟ISI信道 (块长=%d)\n', blockLength);

% 3.2 添加可配置的信道相位偏移
phase_offset_rad_channel = config.channelPhaseOffsetDeg * pi / 180;
if abs(config.channelPhaseOffsetDeg) > 1e-6
    rxSymbols_ISI = rxSymbols_ISI * exp(1j * phase_offset_rad_channel);
    fprintf('  施加信道相位偏移: %.2f 度\n', config.channelPhaseOffsetDeg);
end

% 3.3 添加AWGN噪声 (使用自建噪声发生器)
rxSymbols = my_awgn(rxSymbols_ISI, config.snr_dB, 'measured');
fprintf('  添加自建AWGN噪声 (SNR = %d dB)\n', config.snr_dB);
fprintf('\n');

%% 步骤4: CMA均衡器 (样本级实现)

fprintf('CMA均衡处理...\n');

% 4.1 初始化
M = config.cma.filterLength;
mu = config.cma.stepSize;
R2 = config.cma.R2;

% 中心抽头初始化
w = zeros(M, 1);
w(floor(M/2) + 1) = 1;

% 准备输入信号
rxSymbols_col = rxSymbols(:);
N_total = length(rxSymbols_col);

% 输入缓冲区(添加前导零)
rx_padded = [zeros(M-1, 1); rxSymbols_col];

% 输出
equalizedSymbols = zeros(N_total, 1);

fprintf('  初始化完成\n');
fprintf('  均衡器长度: %d, 符号数: %d\n', M, N_total);
fprintf('  遍历次数: %d\n', config.cma.numPasses);

% 4.2 样本级CMA处理
fprintf('  开始样本级CMA处理...\n');

for pass = 1:config.cma.numPasses
    mu_pass = mu * (config.cma.stepDecay)^(pass-1);
    fprintf('  Pass %d/%d 开始... (步长=%.4g)\n', pass, config.cma.numPasses, mu_pass);
    use_decision_directed = config.cma.ddStartPass > 0 && pass >= config.cma.ddStartPass;
    for n = 1:N_total
        % 获取输入向量 (最新样本在前)
        idx = n + M - 1;
        x_vec = rx_padded(idx:-1:idx-M+1);
        
        % 计算均衡器输出
        y_out = w' * x_vec;
        equalizedSymbols(n) = y_out;
        
        % CMA权重更新 (Godard算法, 归一化LMS形式)
        % 避免权重发散,使用归一化
        x_power = real(x_vec' * x_vec) + 1e-6;  % 输入功率 + 小常数避免除零
        mu_norm = mu_pass / x_power;
        if use_decision_directed
            [~, dec_idx_current] = min(abs(y_out - constellation_map));
            decision_symbol = constellation_map(dec_idx_current);
            error_signal = conj(decision_symbol - y_out);
        else
            error_signal = (R2 - abs(y_out)^2) * conj(y_out);  % 复数误差信号
        end
        w = w + mu_norm * error_signal * x_vec;
        
        % 权重剪裁以防止发散
        max_coeff = max(abs(w));
        if max_coeff > 10
            w = w * (10 / max_coeff);
        end
        
        % 显示进度
        if mod(n, 5000) == 0 || n == N_total
            fprintf('    Pass %d: 处理进度 %d/%d (%.1f%%), norm(w)=%.4f, |y|=%.4f\n', ...
                    pass, n, N_total, 100*n/N_total, norm(w), abs(y_out));
        end
    end
end

fprintf('  均衡完成\n');
fprintf('\n');

%% 步骤5: 相位校正与解调 (技术规格书 4.5节)

fprintf('相位校正与解调...\n');

% 5.1 去除瞬态样本
transient_samples = 1000;
equalizedSymbols_steady = equalizedSymbols(transient_samples+1:end);
equalizedSymbols_steady = equalizedSymbols_steady(:);

% 5.2 简单相位校正 - 使用前导符号估计相位偏移
pilot_length = min(500, length(equalizedSymbols_steady));
pilot_symbols = equalizedSymbols_steady(1:pilot_length);
pilot_symbols = pilot_symbols(:);

% 硬判决得到最近星座点 (向量化)
[~, pilot_indices] = min(abs(bsxfun(@minus, pilot_symbols, constellation_map)), [], 2);
pilot_reference = constellation_map(pilot_indices);
pilot_reference = pilot_reference(:);

% 估计相位旋转
sum_val = sum(pilot_reference .* conj(pilot_symbols));
phase_offset = angle(sum_val);
fprintf('  估计相位偏移: %.2f 度\n', phase_offset * 180 / pi);

% 应用相位校正
equalizedSymbols_corrected = equalizedSymbols_steady * exp(1j * phase_offset);
equalizedSymbols_corrected = equalizedSymbols_corrected(:);

% 5.3 硬判决解调 (向量化)
[~, rxSymbolIndices] = min(abs(bsxfun(@minus, equalizedSymbols_corrected, constellation_map)), [], 2);

% rxSymbolIndices现在是1-4,对应格雷码索引0-3
rxGrayIndices = rxSymbolIndices - 1; % 格雷码索引 0-3

% 5.4 格雷码反向映射 (预计算逆映射)
grayInverse = zeros(size(config.grayMap));
grayInverse(config.grayMap + 1) = 0:config.M-1;
rxSymbolIndices_original = grayInverse(rxGrayIndices + 1);

% 5.5 符号到比特
rxBitPairs = zeros(length(rxSymbolIndices_original), 2);
rxBitPairs(:, 1) = bitshift(rxSymbolIndices_original, -1);
rxBitPairs(:, 2) = bitand(rxSymbolIndices_original, 1);
rxBits = reshape(rxBitPairs.', 1, []);

fprintf('  解调完成\n');
fprintf('\n');

%% 步骤6: 性能评估 (技术规格书 4.6节)

fprintf('性能评估...\n');

% 6.1 序列对齐 - 使用互相关搜索最佳延迟
fprintf('  搜索最佳延迟对齐...\n');
max_search_delay = min(6000, length(txBits) - 100);  % 搜索范围
best_ber = 1.0;
best_delay = 0;

for delay_bits = 0:2:max_search_delay  % 以符号为单位 (2比特步长)
    tx_start = delay_bits + 1;
    rx_start = 1;
    test_length = min(length(txBits) - tx_start + 1, length(rxBits) - rx_start + 1);
    
    if test_length < 1000  % 至少需要1000比特来统计
        continue;
    end
    
    % 只用前5000比特快速测试
    test_length = min(test_length, 5000);
    
    tx_test = txBits(tx_start : tx_start + test_length - 1);
    rx_test = rxBits(rx_start : rx_start + test_length - 1);
    
    test_errors = sum(tx_test ~= rx_test);
    test_ber = test_errors / test_length;
    
    if test_ber < best_ber
        best_ber = test_ber;
        best_delay = delay_bits;
    end
end

fprintf('    最佳延迟: %d 比特 (%d 符号)\n', best_delay, best_delay/2);
fprintf('    估计BER: %.4e (基于%d比特测试)\n', best_ber, test_length);

% 使用最佳延迟对齐全部序列
tx_align_start = best_delay + 1;
rx_align_start = 1;
compare_length = min(length(txBits) - tx_align_start + 1, length(rxBits) - rx_align_start + 1);

txBits_aligned = txBits(tx_align_start : tx_align_start + compare_length - 1);
rxBits_aligned = rxBits(rx_align_start : rx_align_start + compare_length - 1);

% 6.2 计算BER
bit_errors = sum(txBits_aligned ~= rxBits_aligned);
ber = bit_errors / compare_length;

fprintf('  比较比特数: %d\n', compare_length);
fprintf('  比特错误数: %d\n', bit_errors);
fprintf('  误码率 (BER): %.4e\n', ber);

% 6.3 计算EVM
% 获取对应的发射符号 (延迟以符号计)，并在邻域内搜索最优对齐以最小化EVM
base_symbol_delay = floor(best_delay / 2);
search_radius = 40;  % 可调邻域范围
candidate_offsets = base_symbol_delay + (-search_radius:search_radius);

evm_best = inf;
alpha_best = 1;
tx_symbol_delay = base_symbol_delay;
txSymbols_ref_best = [];
eqSymbols_ref_best = [];

eqSymbols_eval = equalizedSymbols_corrected(:);
total_tx_symbols = length(txSymbols);

for offset = candidate_offsets
    if offset < 0 || offset >= total_tx_symbols
        continue;
    end
    tx_symbol_start = offset + 1;
    remaining_tx = total_tx_symbols - offset;
    segment_length = min(length(eqSymbols_eval), remaining_tx);
    if segment_length < 1000  % 保证统计意义
        continue;
    end
    eq_segment = eqSymbols_eval(1:segment_length);
    tx_segment = txSymbols(tx_symbol_start:tx_symbol_start + segment_length - 1);

    alpha_candidate = (eq_segment' * tx_segment(:)) / (eq_segment' * eq_segment);
    error_vec = alpha_candidate * eq_segment - tx_segment(:);
    evm_candidate = sqrt(mean(abs(error_vec).^2)) / sqrt(mean(abs(tx_segment).^2)) * 100;

    if evm_candidate < evm_best
        evm_best = evm_candidate;
        alpha_best = alpha_candidate;
        tx_symbol_delay = offset;
        txSymbols_ref_best = tx_segment(:);
        eqSymbols_ref_best = eq_segment;
    end
end

if isempty(txSymbols_ref_best)
    warning('EVM计算时未找到有效对齐，使用基础延迟。');
    tx_symbol_delay = base_symbol_delay;
    tx_symbol_start = tx_symbol_delay + 1;
    tx_symbol_end = min(tx_symbol_start + length(eqSymbols_eval) - 1, total_tx_symbols);
    segment_length = tx_symbol_end - tx_symbol_start + 1;
    txSymbols_ref_best = txSymbols(tx_symbol_start:tx_symbol_end);
    eqSymbols_ref_best = eqSymbols_eval(1:segment_length);
    alpha_best = (eqSymbols_ref_best' * txSymbols_ref_best(:)) / (eqSymbols_ref_best' * eqSymbols_ref_best);
    evm_best = sqrt(mean(abs(alpha_best * eqSymbols_ref_best - txSymbols_ref_best(:)).^2)) / ...
              sqrt(mean(abs(txSymbols_ref_best).^2)) * 100;
end

equalizedSymbols_aligned = alpha_best * eqSymbols_ref_best;
evm = evm_best;
fprintf('  EVM最优对齐符号延迟: %d\n', tx_symbol_delay);

fprintf('  EVM: %.2f%%\n', evm);
fprintf('\n');

%% 可视化 (技术规格书 4.6节)

fprintf('生成可视化结果...\n');

% 定义4QAM符号对应的颜色 (为每个符号分配一个颜色)
% 符号映射: 0->(1+j)/√2, 1->(1-j)/√2, 2->(-1-j)/√2, 3->(-1+j)/√2
color_map = [
    0.8, 0.2, 0.2;  % 符号0: 红色 (第一象限)
    0.2, 0.8, 0.2;  % 符号1: 绿色 (第四象限)
    0.2, 0.2, 0.8;  % 符号2: 蓝色 (第三象限)
    0.9, 0.7, 0.1   % 符号3: 黄色 (第二象限)
];

figure('Position', [100, 100, 1600, 450]);

% 选择用于绘图的样本数量
num_plot_samples = min(3000, length(txSymbols));

% 使用EVM搜索得到的最佳符号延迟对齐各阶段星座

% 1. 接收信号星座图 (颜色=发射符号，位置=接收信号实际值)
subplot(1, 3, 1);
hold on; grid on; axis equal;

% 按发射符号着色，但绘制接收信号的实际位置
for sym_idx = 0:3
    % 找到发射了该符号的位置
    tx_positions = find(grayIndices == sym_idx);
    % 只绘制前 num_plot_samples 个
    rx_idx = tx_positions(tx_positions <= num_plot_samples);
    
    if ~isempty(rx_idx)
        % 关键：颜色来自发射符号(sym_idx)，位置来自接收信号(rxSymbols)
        plot(real(rxSymbols(rx_idx)), imag(rxSymbols(rx_idx)), '.', ...
             'Color', color_map(sym_idx+1, :), 'MarkerSize', 4);
    end
end

% 绘制理想星座点
for sym_idx = 0:3
    plot(real(constellation_map(sym_idx+1)), imag(constellation_map(sym_idx+1)), ...
         'o', 'MarkerSize', 12, 'LineWidth', 2.5, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_map(sym_idx+1, :));
end
xlim([-1.8 1.8]); ylim([-1.8 1.8]);
title('接收信号 (ISI+噪声)');
xlabel('同相分量 (I)'); ylabel('正交分量 (Q)');
legend('符号 00', '符号 01', '符号 11', '符号 10', 'Location', 'northeast');

% 2. CMA均衡后星座图 (颜色=发射符号，位置=均衡后实际值)
subplot(1, 3, 2);
hold on; grid on; axis equal;

for sym_idx = 0:3
    tx_positions = find(grayIndices == sym_idx);
    
    % 计算均衡后对应的索引（考虑延迟）
    eq_idx = tx_positions - tx_symbol_delay;
    % 只取有效范围内的索引
    valid_mask = (eq_idx >= 1) & (eq_idx <= length(equalizedSymbols_steady));
    eq_idx = eq_idx(valid_mask);
    eq_idx = eq_idx(eq_idx <= num_plot_samples);
    
    if ~isempty(eq_idx)
        % 关键：颜色来自发射符号(sym_idx)，位置来自均衡后信号(equalizedSymbols_steady)
        plot(real(equalizedSymbols_steady(eq_idx)), imag(equalizedSymbols_steady(eq_idx)), '.', ...
             'Color', color_map(sym_idx+1, :), 'MarkerSize', 4);
    end
end

% 绘制理想星座点
for sym_idx = 0:3
    plot(real(constellation_map(sym_idx+1)), imag(constellation_map(sym_idx+1)), ...
         'o', 'MarkerSize', 12, 'LineWidth', 2.5, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_map(sym_idx+1, :));
end
xlim([-1.8 1.8]); ylim([-1.8 1.8]);
title('CMA均衡后 (未相位校正)');
xlabel('同相分量 (I)'); ylabel('正交分量 (Q)');

% 3. 相位校正后星座图 (颜色=发射符号，位置=校正后实际值)
subplot(1, 3, 3);
hold on; grid on; axis equal;

for sym_idx = 0:3
    tx_positions = find(grayIndices == sym_idx);
    
    % 计算校正后对应的索引（考虑延迟）
    eq_idx = tx_positions - tx_symbol_delay;
    valid_mask = (eq_idx >= 1) & (eq_idx <= length(equalizedSymbols_corrected));
    eq_idx = eq_idx(valid_mask);
    eq_idx = eq_idx(eq_idx <= num_plot_samples);
    
    if ~isempty(eq_idx)
        % 关键：颜色来自发射符号(sym_idx)，位置来自校正后信号(equalizedSymbols_corrected)
        plot(real(equalizedSymbols_corrected(eq_idx)), imag(equalizedSymbols_corrected(eq_idx)), '.', ...
             'Color', color_map(sym_idx+1, :), 'MarkerSize', 4);
    end
end

% 绘制理想星座点
for sym_idx = 0:3
    plot(real(constellation_map(sym_idx+1)), imag(constellation_map(sym_idx+1)), ...
         'o', 'MarkerSize', 12, 'LineWidth', 2.5, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_map(sym_idx+1, :));
end
xlim([-1.8 1.8]); ylim([-1.8 1.8]);
title(sprintf('相位校正后 (BER=%.2e, EVM=%.2f%%)', ber, evm));
xlabel('同相分量 (I)'); ylabel('正交分量 (Q)');

sgtitle('4QAM CMA盲均衡仿真结果 (颜色=发射符号)', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('可视化完成\n');
fprintf('\n========================================\n');
fprintf('仿真完成!\n');
fprintf('========================================\n');
