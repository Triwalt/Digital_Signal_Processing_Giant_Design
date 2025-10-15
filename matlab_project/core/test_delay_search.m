% 测试:使用互相关寻找最佳延迟对齐
clear; clc;

fprintf('==== 互相关延迟搜索测试 ====\n\n');

% 简化参数
numSymbols = 2000;
M = 15;
mu = 0.01;
R2 = 1;
snr_dB = 25;

grayMap = [0 1 3 2];
constellation_map = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);

%% 发射
numBits = numSymbols * 2;
txBits = randi([0 1], 1, numBits);
txSymbolIndices = zeros(1, numSymbols);
for i = 1:numSymbols
    bit1 = txBits(2*i-1);
    bit2 = txBits(2*i);
    decimal_idx = bit1 * 2 + bit2;
    txSymbolIndices(i) = decimal_idx;
end
grayIndices = grayMap(txSymbolIndices + 1);
txSymbols = constellation_map(grayIndices + 1);

%% 信道
channelTaps = [1, 0.5, 0.2];
rxSymbols_ISI = filter(channelTaps, 1, txSymbols);
rxSymbols = awgn(rxSymbols_ISI, snr_dB, 'measured');

%% CMA
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

fprintf('CMA完成, norm(w)=%.4f\n', norm(w));

%% 解调 (无相位校正,先测试)
transient = 200;
equalizedSymbols_steady = equalizedSymbols(transient+1:end);

rxSymbolIndices = zeros(length(equalizedSymbols_steady), 1);
for k = 1:length(equalizedSymbols_steady)
    [~, rxSymbolIndices(k)] = min(abs(equalizedSymbols_steady(k) - constellation_map));
end

rxGrayIndices = rxSymbolIndices - 1;
rxSymbolIndices_original = zeros(size(rxGrayIndices));
for k = 1:length(rxGrayIndices)
    pos = find(grayMap == rxGrayIndices(k));
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

%% 使用互相关搜索最佳延迟
fprintf('\n搜索最佳延迟...\n');
max_delay = 500;  % 搜索范围
best_ber = 1.0;
best_delay = 0;

for delay_bits = 0:2:max_delay  % 以符号为单位搜索(2比特步长)
    tx_start = delay_bits + 1;
    rx_start = 1;
    len = min(length(txBits) - tx_start + 1, length(rxBits) - rx_start + 1);
    
    if len < 100
        continue;
    end
    
    tx_seg = txBits(tx_start : tx_start + len - 1);
    rx_seg = rxBits(rx_start : rx_start + len - 1);
    
    errors = sum(tx_seg ~= rx_seg);
    ber = errors / len;
    
    if ber < best_ber
        best_ber = ber;
        best_delay = delay_bits;
    end
end

fprintf('  最佳延迟: %d bits (%d 符号)\n', best_delay, best_delay/2);
fprintf('  最佳BER: %.4e\n\n', best_ber);

% 使用最佳延迟重新计算
tx_start = best_delay + 1;
rx_start = 1;
len = min(length(txBits) - tx_start + 1, length(rxBits) - rx_start + 1);
tx_aligned = txBits(tx_start : tx_start + len - 1);
rx_aligned = rxBits(rx_start : rx_start + len - 1);

fprintf('对齐后:\n');
fprintf('  TX前20bit: '); disp(tx_aligned(1:20));
fprintf('  RX前20bit: '); disp(rx_aligned(1:20));
fprintf('  BER: %.4e\n', sum(tx_aligned ~= rx_aligned) / len);

fprintf('\n==== 完成 ====\n');
