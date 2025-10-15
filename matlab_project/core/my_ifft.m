function x = my_ifft(X)
% my_ifft: 实现逆快速傅里叶变换
%
% 实现原理: 利用FFT函数高效实现IFFT
% IDFT(X) = (1/N) * conj(DFT(conj(X)))
%
% 输入:
%   X - 频域信号,复数列向量,长度必须是2的幂
%
% 输出:
%   x - 时域信号,复数列向量
%
% 参考: 技术规格书第2.2节

    % 确保输入为列向量
    if size(X, 2) > 1
        X = X(:);
    end
    
    N = length(X);
    
    % 验证N是2的幂
    if N == 0 || mod(log2(N), 1) ~= 0
        error('输入长度必须是2的整数次幂');
    end
    
    % 步骤1: 计算输入的共轭
    X_conj = conj(X);
    
    % 步骤2: 调用my_fft
    fft_result = my_fft(X_conj);
    
    % 步骤3: 取共轭并除以N
    x = conj(fft_result) / N;
end
