% 测试fast_conv_os
clear; clc;

M = 31;
nfft = 512;
L = nfft - M + 1; % 482

% 创建测试信号
N_total = 1000;
input_signal = randn(N_total, 1) + 1j * randn(N_total, 1);
input_buffer = [zeros(M-1, 1); input_signal];

% 创建滤波器
w = zeros(M, 1);
w(floor(M/2) + 1) = 1;

fprintf('测试fast_conv_os\n');
fprintf('M = %d, nfft = %d, L = %d\n', M, nfft, L);
fprintf('N_total = %d\n', N_total);

num_blocks = ceil(N_total / L); % 3块
fprintf('num_blocks = %d\n', num_blocks);

for block_idx = 1:num_blocks
    block_start = (block_idx - 1) * L + 1;
    block_end = min(block_idx * L, N_total);
    block_samples = block_end - block_start + 1;
    
    fprintf('\n块 %d: start=%d, end=%d, samples=%d\n', ...
            block_idx, block_start, block_end, block_samples);
    
    input_start = block_start;
    input_end = block_end + M - 1;
    input_block = input_buffer(input_start:input_end);
    
    fprintf('  input_block长度: %d\n', length(input_block));
    
    if block_samples == L
        output_block = fast_conv_os(input_block, w, nfft);
        fprintf('  输出长度: %d, 期望: %d\n', length(output_block), length(input_block) + M - 1);
        fprintf('  使用前%d个样本\n', L);
    else
        output_block = conv(input_block, w);
        fprintf('  conv输出长度: %d\n', length(output_block));
        fprintf('  使用前%d个样本\n', block_samples);
    end
end
