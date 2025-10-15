% 测试CMA均衡器
clear; clc; close all;

fprintf('CMA均衡器测试\n');
fprintf('================\n\n');

% 简单配置
numSymbols = 1000;
snr_dB = 25;
channelTaps = [1, 0.5, 0.2];
filterLength = 31;
stepSize = 0.001;
R2 = 1;
nfft = 512;

% 生成信号
grayMap = [0 1 3 2];
constellation_map = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);

txBits = randi([0 1], 1, numSymbols * 2);
txSymbolIndices = zeros(1, numSymbols);
for i = 1:numSymbols
    txSymbolIndices(i) = txBits(2*i-1) * 2 + txBits(2*i);
end
grayIndices = grayMap(txSymbolIndices + 1);
txSymbols = constellation_map(grayIndices + 1);

% 信道
rxSymbols_ISI = filter(channelTaps, 1, txSymbols);
rxSymbols = awgn(rxSymbols_ISI, snr_dB, 'measured');

fprintf('信号生成完成\n');
fprintf('  发射能量: %.4f\n', mean(abs(txSymbols).^2));
fprintf('  接收能量: %.4f\n', mean(abs(rxSymbols).^2));

% CMA均衡
M = filterLength;
L = nfft - M + 1;

w = zeros(M, 1);
w(floor(M/2) + 1) = 1;

rxSymbols_col = rxSymbols(:);
N_total = length(rxSymbols_col);
equalizedSymbols = zeros(N_total, 1);
input_buffer = [zeros(M-1, 1); rxSymbols_col];

num_blocks = ceil(N_total / L);
fprintf('\nCMA均衡: %d块\n', num_blocks);

for block_idx = 1:num_blocks
    block_start = (block_idx - 1) * L + 1;
    block_end = min(block_idx * L, N_total);
    block_samples = block_end - block_start + 1;
    
    input_start = block_start;
    input_end = block_end + M - 1;
    input_block = input_buffer(input_start:input_end);
    
    if block_samples == L
        output_block = fast_conv_os(input_block, w, nfft);
        equalizedSymbols(block_start:block_end) = output_block(1:L);
    else
        output_block = conv(input_block, w);
        equalizedSymbols(block_start:block_end) = output_block(1:block_samples);
    end
    
    for n = block_start:block_end
        x_tilde = equalizedSymbols(n);
        y_vec = input_buffer(n:n+M-1);
        error_term = R2 - abs(x_tilde)^2;
        w = w + stepSize * error_term * x_tilde * conj(y_vec);
    end
end

fprintf('均衡完成\n');
fprintf('  均衡后能量: %.4f\n', mean(abs(equalizedSymbols).^2));
fprintf('  权重范数: %.4f\n', norm(w));

% 绘图
figure('Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
plot(real(rxSymbols), imag(rxSymbols), '.', 'MarkerSize', 3);
title('接收信号');
grid on; axis equal; xlim([-2 2]); ylim([-2 2]);

subplot(1, 3, 2);
plot(real(equalizedSymbols), imag(equalizedSymbols), '.', 'MarkerSize', 3);
title('均衡后(全部)');
grid on; axis equal; xlim([-2 2]); ylim([-2 2]);

subplot(1, 3, 3);
plot(real(equalizedSymbols(500:end)), imag(equalizedSymbols(500:end)), '.', 'MarkerSize', 3);
title('均衡后(去除瞬态)');
grid on; axis equal; xlim([-2 2]); ylim([-2 2]);
