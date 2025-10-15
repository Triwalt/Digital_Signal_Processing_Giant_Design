% 完整CMA测试 - 从发射到接收
clear; clc;

fprintf('==== 完整CMA链路测试 ====\n\n');

% 参数
numSymbols = 1000;
M = 15;  % 均衡器长度
mu = 0.001;
R2 = 1;
snr_dB = 25;

% 4QAM参数
grayMap = [0 1 3 2];
constellation_map = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);

%% 发射端
fprintf('1. 发射端\n');
numBits = numSymbols * 2;
txBits = randi([0 1], 1, numBits);

% 比特到符号
txSymbolIndices = zeros(1, numSymbols);
for i = 1:numSymbols
    bit1 = txBits(2*i-1);
    bit2 = txBits(2*i);
    decimal_idx = bit1 * 2 + bit2;
    txSymbolIndices(i) = decimal_idx;
end

% 格雷编码
grayIndices = grayMap(txSymbolIndices + 1);
txSymbols = constellation_map(grayIndices + 1);

fprintf('   符号数: %d\n', numSymbols);
fprintf('   符号能量: %.4f\n', mean(abs(txSymbols).^2));

%% 信道
fprintf('\n2. 信道\n');
channelTaps = [1, 0.5, 0.2];
rxSymbols_ISI = filter(channelTaps, 1, txSymbols);
rxSymbols = awgn(rxSymbols_ISI, snr_dB, 'measured');
fprintf('   ISI信道: [1, 0.5, 0.2]\n');
fprintf('   SNR: %d dB\n', snr_dB);

%% CMA均衡
fprintf('\n3. CMA均衡 (归一化LMS)\n');
w = zeros(M, 1);
w(floor(M/2) + 1) = 1;

input_buffer = [zeros(M-1, 1); rxSymbols(:)];
N_total = length(rxSymbols);
equalizedSymbols = zeros(N_total, 1);

for n = 1:N_total
    x_vec = input_buffer(n:n+M-1);
    
    % 归一化
    x_power = x_vec' * x_vec;
    if x_power < 1e-10
        x_power = 1e-10;
    end
    
    y_out = w' * x_vec;
    equalizedSymbols(n) = y_out;
    
    error_signal = (R2 - abs(y_out)^2) * y_out;
    w = w + (mu / x_power) * error_signal * conj(x_vec);
end

fprintf('   权重范数: %.4f\n', norm(w));
fprintf('   输出能量: %.4f\n', mean(abs(equalizedSymbols).^2));

%% 解调
fprintf('\n4. 解调\n');
transient = 100;
equalizedSymbols_steady = equalizedSymbols(transient+1:end);

% 相位校正
pilot_length = min(200, length(equalizedSymbols_steady));
pilot_symbols = equalizedSymbols_steady(1:pilot_length);

pilot_reference = zeros(size(pilot_symbols));
for k = 1:length(pilot_symbols)
    [~, idx] = min(abs(pilot_symbols(k) - constellation_map));
    pilot_reference(k) = constellation_map(idx);
end

phase_offset = angle(sum(pilot_reference .* conj(pilot_symbols)));
fprintf('   相位偏移: %.2f度\n', phase_offset * 180/pi);

equalizedSymbols_corrected = equalizedSymbols_steady * exp(1j * phase_offset);

% 硬判决
rxSymbolIndices = zeros(length(equalizedSymbols_corrected), 1);
for k = 1:length(equalizedSymbols_corrected)
    [~, rxSymbolIndices(k)] = min(abs(equalizedSymbols_corrected(k) - constellation_map));
end

rxGrayIndices = rxSymbolIndices - 1;

% 格雷逆映射
rxSymbolIndices_original = zeros(size(rxGrayIndices));
for k = 1:length(rxGrayIndices)
    pos = find(grayMap == rxGrayIndices(k));
    if ~isempty(pos)
        rxSymbolIndices_original(k) = pos - 1;
    end
end

% 符号到比特
rxBits = zeros(1, length(rxSymbolIndices_original) * 2);
for i = 1:length(rxSymbolIndices_original)
    idx = rxSymbolIndices_original(i);
    rxBits(2*i-1) = bitshift(idx, -1);
    rxBits(2*i) = bitand(idx, 1);
end

%% BER计算
fprintf('\n5. BER计算\n');
channel_delay = length(channelTaps) - 1;
equalizer_delay = floor(M / 2);
total_delay = channel_delay + equalizer_delay + transient;

tx_align_start = total_delay * 2 + 1;
rx_align_start = 1;
compare_length = min(length(txBits) - tx_align_start + 1, length(rxBits) - rx_align_start + 1);

fprintf('   信道延迟: %d\n', channel_delay);
fprintf('   均衡器延迟: %d\n', equalizer_delay);
fprintf('   瞬态长度: %d\n', transient);
fprintf('   总延迟: %d 符号 = %d 比特\n', total_delay, total_delay*2);

txBits_aligned = txBits(tx_align_start : tx_align_start + compare_length - 1);
rxBits_aligned = rxBits(rx_align_start : rx_align_start + compare_length - 1);

bit_errors = sum(txBits_aligned ~= rxBits_aligned);
ber = bit_errors / compare_length;

fprintf('   比较长度: %d bits\n', compare_length);
fprintf('   误码数: %d\n', bit_errors);
fprintf('   BER: %.4e\n\n', ber);

%% 可视化
figure('Position', [100 100 1200 400]);

subplot(1,3,1);
plot(real(rxSymbols), imag(rxSymbols), 'b.', 'MarkerSize', 3);
grid on; axis equal; axis([-2 2 -2 2]);
title('接收信号 (ISI+噪声)');
xlabel('实部'); ylabel('虚部');

subplot(1,3,2);
plot(real(equalizedSymbols_corrected), imag(equalizedSymbols_corrected), 'r.', 'MarkerSize', 3);
hold on;
plot(real(constellation_map), imag(constellation_map), 'ko', 'MarkerSize', 10, 'LineWidth', 2);
grid on; axis equal; axis([-1.5 1.5 -1.5 1.5]);
title(sprintf('均衡后 (BER=%.2e)', ber));
xlabel('实部'); ylabel('虚部');
legend('均衡符号', '理想星座', 'Location', 'best');

subplot(1,3,3);
plot(abs(equalizedSymbols_corrected));
hold on;
yline(sqrt(2)/sqrt(2), 'r--', 'R=1');
grid on;
title('均衡符号幅度');
xlabel('符号索引'); ylabel('幅度');

fprintf('==== 测试完成 ====\n');
