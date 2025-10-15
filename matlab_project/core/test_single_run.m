% 简化测试 - 只测试一个CMA配置
clear; clc;

fprintf('==== 单次CMA测试 ====\n\n');

% 配置
config = struct();
config.numSymbols = 10000;
config.snr_dB = 25;
config.M = 4;
config.grayMap = [0 1 3 2];
config.channelTaps = [1, 0.5, 0.2];
config.cma.filterLength = 31;
config.cma.R2 = 1;

mu = 0.01;
filter_length = 31;

% 发射
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

fprintf('TX: %d symbols, %d bits\n', config.numSymbols, numBits);

% 信道
rxSymbols_ISI = filter(config.channelTaps, 1, txSymbols);
rxSymbols = awgn(rxSymbols_ISI, config.snr_dB, 'measured');

% CMA
M = filter_length;
R2 = config.cma.R2;

w = zeros(M, 1);
w(floor(M/2) + 1) = 1;

rxSymbols_col = rxSymbols(:);
N_total = length(rxSymbols_col);
equalizedSymbols = zeros(N_total, 1);
input_buffer = [zeros(M-1, 1); rxSymbols_col];

fprintf('CMA Processing...\n');
for n = 1:N_total
    x_vec = input_buffer(n:n+M-1);
    x_power = real(x_vec' * x_vec) + 1e-6;
    y_out = w' * x_vec;
    equalizedSymbols(n) = y_out;
    error_signal = (R2 - abs(y_out)^2) * y_out;
    w = w + (mu / x_power) * conj(error_signal) * x_vec;
    
    if mod(n, 2000) == 0
        fprintf('  %d/%d\n', n, N_total);
    end
end

fprintf('CMA Done, norm(w)=%.4f\n', norm(w));

% 相位校正
transient_samples = 1000;
equalizedSymbols_steady = equalizedSymbols(transient_samples+1:end);

pilot_length = min(500, length(equalizedSymbols_steady));
pilot_symbols = equalizedSymbols_steady(1:pilot_length);

pilot_reference = zeros(size(pilot_symbols));
for k = 1:length(pilot_symbols)
    [~, idx] = min(abs(pilot_symbols(k) - constellation_map));
    pilot_reference(k) = constellation_map(idx);
end

phase_offset = angle(sum(pilot_reference .* conj(pilot_symbols)));
fprintf('Phase offset: %.2f deg\n', phase_offset*180/pi);

equalizedSymbols_corrected = equalizedSymbols_steady * exp(1j * phase_offset);

% 解调
rxSymbolIndices = zeros(length(equalizedSymbols_corrected), 1);
for k = 1:length(equalizedSymbols_corrected)
    [~, rxSymbolIndices(k)] = min(abs(equalizedSymbols_corrected(k) - constellation_map));
end

rxGrayIndices = rxSymbolIndices - 1;

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

fprintf('RX: %d bits demodulated\n', length(rxBits));
fprintf('TX bits: %d, RX bits: %d\n', length(txBits), length(rxBits));

% BER计算 - 延迟搜索
fprintf('\nSearching for best delay...\n');
max_search_delay = min(2000, length(txBits) - 100);
best_ber = 1.0;
best_delay = 0;

for delay_bits = 0:2:max_search_delay
    tx_start = delay_bits + 1;
    rx_start = 1;
    test_length = min(length(txBits) - tx_start + 1, length(rxBits) - rx_start + 1);
    
    if test_length < 1000
        continue;
    end
    
    test_length = min(test_length, 3000);
    tx_test = txBits(tx_start : tx_start + test_length - 1);
    rx_test = rxBits(rx_start : rx_start + test_length - 1);
    
    test_errors = sum(tx_test ~= rx_test);
    test_ber = test_errors / test_length;
    
    if test_ber < best_ber
        best_ber = test_ber;
        best_delay = delay_bits;
        fprintf('  delay=%d bits, BER=%.4f\n', delay_bits, test_ber);
    end
end

fprintf('\nBest delay: %d bits (%d symbols)\n', best_delay, best_delay/2);
fprintf('Best BER (from search): %.4e\n\n', best_ber);

% 使用最佳延迟
tx_align_start = best_delay + 1;
rx_align_start = 1;
compare_length = min(length(txBits) - tx_align_start + 1, length(rxBits) - rx_align_start + 1);

txBits_aligned = txBits(tx_align_start : tx_align_start + compare_length - 1);
rxBits_aligned = rxBits(rx_align_start : rx_align_start + compare_length - 1);

bit_errors = sum(txBits_aligned ~= rxBits_aligned);
ber = bit_errors / compare_length;

fprintf('Final BER: %.4e (%d errors in %d bits)\n', ber, bit_errors, compare_length);

fprintf('\n==== 完成 ====\n');
