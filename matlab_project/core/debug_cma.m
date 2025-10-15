% Debug CMA
clear; clc;

% 简单测试
M = 11; % 更短的滤波器
mu = 0.001;
R2 = 1;

% 生成信号
numSymbols = 100;
constellation_map = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);
tx_idx = randi([1 4], 1, numSymbols);
txSymbols = constellation_map(tx_idx);

% 信道
channelTaps = [1, 0.3];
rxSymbols = filter(channelTaps, 1, txSymbols);
rxSymbols = awgn(rxSymbols, 20, 'measured');

% CMA
w = zeros(M, 1);
w(floor(M/2) + 1) = 1;

input_buffer = [zeros(M-1, 1); rxSymbols(:)];
N_total = length(rxSymbols);
equalizedSymbols = zeros(N_total, 1);

fprintf('CMA调试\n');
fprintf('初始权重: w(%d) = 1\n', floor(M/2)+1);

for n = 1:N_total
    x_vec = input_buffer(n:n+M-1);
    y_out = w' * x_vec;
    equalizedSymbols(n) = y_out;
    
    error_signal = (R2 - abs(y_out)^2) * y_out;
    w = w + mu * conj(x_vec) * error_signal;
    
    if n <= 15 || (n >= 50 && n <= 55)
        fprintf('n=%d: x_vec(1)=%.2f%+.2fj, |y|=%.3f, |error|=%.3f\n', ...
                n, real(x_vec(1)), imag(x_vec(1)), abs(y_out), abs(error_signal));
    end
    
    if any(~isfinite(w))
        fprintf('权重变为inf/NaN at n=%d!\n', n);
        break;
    end
end

fprintf('\n最终:\n');
fprintf('  权重范数: %.4f\n', norm(w));
fprintf('  输出能量: %.4f\n', mean(abs(equalizedSymbols).^2));
fprintf('  NaN数量: %d\n', sum(~isfinite(equalizedSymbols)));

% 绘图
figure;
subplot(1,2,1);
plot(real(rxSymbols), imag(rxSymbols), '.');
title('接收');
grid on; axis equal;

subplot(1,2,2);
plot(real(equalizedSymbols), imag(equalizedSymbols), '.');
title('均衡后');
grid on; axis equal;
