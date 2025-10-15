% demo_color_constellation.m: 彩色星座图演示
% 
% 功能: 通过颜色标记展示CMA均衡和相位校正的效果
% 不同符号用不同颜色标记，便于观察各符号的聚合情况

clear; clc; close all;

fprintf('========================================\n');
fprintf('彩色星座图演示 - CMA均衡效果可视化\n');
fprintf('========================================\n\n');

%% 配置参数
config = struct();
config.numSymbols = 5000;  % 使用较少符号以便清晰显示
config.snr_dB = 20;        % 中等SNR以显示均衡效果
config.M = 4;
config.grayMap = [0 1 3 2];
config.channelTaps = [1, 0.5, 0.2];
config.cma.filterLength = 31;
config.cma.stepSize = 0.01;
config.cma.R2 = 1;

fprintf('配置:\n');
fprintf('  符号数: %d\n', config.numSymbols);
fprintf('  SNR: %d dB\n', config.snr_dB);
fprintf('  信道: ISI [1, 0.5, 0.2]\n\n');

%% 发射端
numBits = config.numSymbols * 2;
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

fprintf('发射端: 生成 %d 个4QAM符号\n', config.numSymbols);

%% 信道
rxSymbols_ISI = filter(config.channelTaps, 1, txSymbols);
rxSymbols = awgn(rxSymbols_ISI, config.snr_dB, 'measured');
fprintf('信道: ISI + AWGN (SNR=%ddB)\n', config.snr_dB);

%% CMA均衡
M = config.cma.filterLength;
mu = config.cma.stepSize;
R2 = config.cma.R2;

w = zeros(M, 1);
w(floor(M/2) + 1) = 1;

input_buffer = [zeros(M-1, 1); rxSymbols(:)];
N_total = length(rxSymbols);
equalizedSymbols = zeros(N_total, 1);

for n = 1:N_total
    x_vec = input_buffer(n:n+M-1);
    x_power = real(x_vec' * x_vec) + 1e-6;
    y_out = w' * x_vec;
    equalizedSymbols(n) = y_out;
    error_signal = (R2 - abs(y_out)^2) * y_out;
    w = w + (mu / x_power) * conj(error_signal) * x_vec;
end

fprintf('CMA均衡: 完成 (权重范数=%.4f)\n', norm(w));

%% 相位校正
transient_samples = 500;
equalizedSymbols_steady = equalizedSymbols(transient_samples+1:end);

pilot_length = min(200, length(equalizedSymbols_steady));
pilot_symbols = equalizedSymbols_steady(1:pilot_length);

pilot_reference = zeros(size(pilot_symbols));
for k = 1:length(pilot_symbols)
    [~, idx] = min(abs(pilot_symbols(k) - constellation_map));
    pilot_reference(k) = constellation_map(idx);
end

phase_offset = angle(sum(pilot_reference .* conj(pilot_symbols)));
equalizedSymbols_corrected = equalizedSymbols_steady * exp(1j * phase_offset);

fprintf('相位校正: 偏移=%.2f度\n\n', phase_offset*180/pi);

%% 彩色可视化
% 定义颜色: 每个符号一个颜色
color_map = [
    0.85, 0.15, 0.15;  % 符号00 (0): 鲜红色
    0.15, 0.85, 0.15;  % 符号01 (1): 鲜绿色  
    0.15, 0.15, 0.85;  % 符号11 (2): 鲜蓝色
    0.95, 0.75, 0.05   % 符号10 (3): 金黄色
];

symbol_labels = {'00 (I)', '01 (IV)', '11 (III)', '10 (II)'};

% 计算延迟
tx_symbol_delay = floor((transient_samples) / 2) + floor(M/2) + length(config.channelTaps) - 1;

fprintf('生成彩色星座图...\n');
figure('Position', [50, 50, 1800, 500]);

%% 子图1: 接收信号 (颜色=发射符号，位置=接收实际值)
subplot(1, 3, 1);
hold on; grid on; box on;
num_plot = min(2000, config.numSymbols);

% 按发射符号着色，绘制接收信号的实际位置
for sym_idx = 0:3
    tx_positions = find(grayIndices == sym_idx);
    rx_idx = tx_positions(tx_positions <= num_plot);
    if ~isempty(rx_idx)
        % 颜色=发射符号，位置=接收信号实际值
        plot(real(rxSymbols(rx_idx)), imag(rxSymbols(rx_idx)), '.', ...
             'Color', color_map(sym_idx+1, :), 'MarkerSize', 6);
    end
end

% 绘制理想星座点
for sym_idx = 0:3
    plot(real(constellation_map(sym_idx+1)), imag(constellation_map(sym_idx+1)), ...
         'o', 'MarkerSize', 14, 'LineWidth', 3, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_map(sym_idx+1, :));
end

axis equal; xlim([-2.2 2.2]); ylim([-2.2 2.2]);
xlabel('同相分量 (I)', 'FontSize', 11);
ylabel('正交分量 (Q)', 'FontSize', 11);
title('(a) 接收信号 (ISI + 噪声)', 'FontSize', 12, 'FontWeight', 'bold');
legend(symbol_labels, 'Location', 'northeast', 'FontSize', 9);

%% 子图2: CMA均衡后 (颜色=发射符号，位置=均衡后实际值)
subplot(1, 3, 2);
hold on; grid on; box on;

for sym_idx = 0:3
    tx_positions = find(grayIndices == sym_idx);
    eq_idx = tx_positions - tx_symbol_delay;
    valid_mask = (eq_idx >= 1) & (eq_idx <= length(equalizedSymbols_steady));
    eq_idx = eq_idx(valid_mask);
    eq_idx = eq_idx(eq_idx <= num_plot);
    
    if ~isempty(eq_idx)
        % 颜色=发射符号，位置=均衡后实际值
        plot(real(equalizedSymbols_steady(eq_idx)), imag(equalizedSymbols_steady(eq_idx)), '.', ...
             'Color', color_map(sym_idx+1, :), 'MarkerSize', 6);
    end
end

for sym_idx = 0:3
    plot(real(constellation_map(sym_idx+1)), imag(constellation_map(sym_idx+1)), ...
         'o', 'MarkerSize', 14, 'LineWidth', 3, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_map(sym_idx+1, :));
end

axis equal; xlim([-2.2 2.2]); ylim([-2.2 2.2]);
xlabel('同相分量 (I)', 'FontSize', 11);
ylabel('正交分量 (Q)', 'FontSize', 11);
title('(b) CMA均衡后 (相位未校正)', 'FontSize', 12, 'FontWeight', 'bold');

%% 子图3: 相位校正后 (颜色=发射符号，位置=校正后实际值)
subplot(1, 3, 3);
hold on; grid on; box on;

for sym_idx = 0:3
    tx_positions = find(grayIndices == sym_idx);
    eq_idx = tx_positions - tx_symbol_delay;
    valid_mask = (eq_idx >= 1) & (eq_idx <= length(equalizedSymbols_corrected));
    eq_idx = eq_idx(valid_mask);
    eq_idx = eq_idx(eq_idx <= num_plot);
    
    if ~isempty(eq_idx)
        % 颜色=发射符号，位置=校正后实际值
        plot(real(equalizedSymbols_corrected(eq_idx)), imag(equalizedSymbols_corrected(eq_idx)), '.', ...
             'Color', color_map(sym_idx+1, :), 'MarkerSize', 6);
    end
end

for sym_idx = 0:3
    plot(real(constellation_map(sym_idx+1)), imag(constellation_map(sym_idx+1)), ...
         'o', 'MarkerSize', 14, 'LineWidth', 3, ...
         'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_map(sym_idx+1, :));
end

axis equal; xlim([-2.2 2.2]); ylim([-2.2 2.2]);
xlabel('同相分量 (I)', 'FontSize', 11);
ylabel('正交分量 (Q)', 'FontSize', 11);
title('(c) 相位校正后 (最终输出)', 'FontSize', 12, 'FontWeight', 'bold');

sgtitle(sprintf('CMA盲均衡效果演示 (SNR=%ddB, 不同颜色代表不同符号)', config.snr_dB), ...
        'FontSize', 14, 'FontWeight', 'bold');

fprintf('可视化完成!\n');
fprintf('\n观察要点:\n');
fprintf('  - 左图: 接收信号受ISI影响，各颜色点混叠分散\n');
fprintf('  - 中图: CMA均衡消除ISI，各颜色点聚集但存在相位旋转\n');
fprintf('  - 右图: 相位校正后，各颜色点准确对齐到理想星座点\n');
fprintf('\n========================================\n');
