function X = my_fft_mix(x)
% my_fft_mix  用于任意长度一维序列的混合基FFT实现
%   X = my_fft_mix(x) 使用递归混合基Cooley-Tukey算法计算x的N点离散傅里叶变换。
%   输入长度N可以是任意正整数。对于合数N，算法将N分解为N = r * m，
%   并应用radix-r和radix-m分解。对于小N或质数N，回退到O(N^2)的直接DFT计算。
%
%   本实现避免调用MATLAB内置的fft函数，主要用于教学参考，而非追求最高性能。
%
%   输入可以是行向量或列向量；输出将保持与输入相同的方向。
%   本函数支持复数值输入。
%
%   示例:
%       n = 0:14;
%       x = cos(2*pi*3*n/15) + 0.5*sin(2*pi*5*n/15);
%       X1 = my_fft_mix(x);
%       X2 = fft(x);
%       max_err = max(abs(X1 - X2));
%       fprintf('与fft的最大差异: %g\n', max_err);
%
% 作者: 由Cascade生成的混合基FFT实现

    % 确保输入是向量
    if ~isvector(x)
        error('my_fft_mix:InputNotVector', '输入x必须是一维向量');
    end

    % 记住原始方向（行向量或列向量）
    wasColumn = iscolumn(x);

    % 内部统一使用行向量处理
    x = x(:).';

    % 调用递归混合基FFT核心函数
    X = local_fft_mix(x);

    % 恢复原始向量方向
    if wasColumn
        X = X.';  % 如果输入是列向量，将输出也转为列向量
    end
end

function X = local_fft_mix(x)
    % 递归辅助函数，假设x是行向量
    N = length(x);

    % Planner 策略：先判断质数/合数，再根据规模选择算法

    % 1) 质数长度：小质数用 codelet，大质数用 Bluestein
    if isprime(N)
        if N <= 11
            % 小质数：调用专门的 codelet（当前实现为通用DFT，将来可替换为C/MEX）
            X = fft_codelet_small_prime(x);
        else
            % 大质数：使用 Bluestein 算法，将DFT转为卷积，再用混合基FFT实现
            X = bluestein_fft(x);
        end
        return;
    end

    % 2) 极小长度（非质数部分）：直接DFT
    if N <= 16
        X = dft_direct(x);
        return;
    end

    % 3) 合数：使用混合基 Cooley-Tukey（可看作一般的 mixed-radix / PFA 框架）
    % 尝试将N分解为r × m，其中r是N的最小非1非N的因子
    r = smallest_factor(N);

    % 安全性检查：如果意外没有找到合适因子，则退回DFT
    if isempty(r) || r == 1 || r == N
        X = dft_direct(x);
        return;
    end

    % 计算另一个因子m = N/r
    m = N / r;  % 现在N = r * m

    % 将输入序列x重新排列为m行r列的矩阵
    % 每一列包含m个点，共r列
    x_mat = reshape(x, m, r);

    % 步骤1：对每一列进行m点FFT
    % 这对应于Cooley-Tukey算法的分解步骤
    X_col = zeros(m, r);  % 预分配结果矩阵
    for j = 1:r
        % 对第j列进行m点FFT（递归调用）
        % 注意：由于MATLAB的reshape按列优先，所以需要转置
        X_col(:, j) = local_fft_mix(x_mat(:, j).').';
    end

    % 步骤2：乘以旋转因子(twiddle factors)
    % 这些因子用于调整子FFT结果的位置
    k = (0:m-1).';                % 行索引，表示每个m点FFT内部的频率分量
    j = 0:(r-1);                  % 列索引，表示r个m点FFT的编号
    % 计算m×r的旋转因子矩阵
    % 每个元素W_N^(k*j) = exp(-j*2π*k*j/N)
    twiddle = exp(-1j * 2 * pi * (k * j) / N);
    % 应用旋转因子（逐元素相乘）
    X_col = X_col .* twiddle;

    % 步骤3：对每一行进行r点FFT，并将结果映射到最终输出
    X = zeros(1, N);  % 预分配最终输出
    for k1 = 0:(m-1)  % 遍历每一行
        % 获取当前行（包含r个点）
        row = X_col(k1+1, :);

        % 对当前行进行r点FFT（递归调用，这样如果r可以继续分解，会自动处理）
        row_fft = local_fft_mix(row);

        % 将结果映射到全局索引
        % 根据Cooley-Tukey算法，全局频率索引k = k1 + m*k2
        % 其中k1 = 0,...,m-1, k2 = 0,...,r-1
        for k2 = 0:(r-1)
            % 计算全局频率索引（注意MATLAB是1-based索引）
            idx = k1 + m * k2 + 1;
            % 存储结果
            X(idx) = row_fft(k2+1);
        end
    end
end

function X = dft_direct(x)
    % 通用的 O(N^2) 直接DFT实现，作为小N或兜底方案
    N = length(x);
    n = 0:N-1;               % 时域索引 [0,1,2,...,N-1]
    k = n.';                  % 频域索引转置为列向量
    W = exp(-1j * 2 * pi / N);  % 基频旋转因子
    F = W .^ (k * n);         % 构造N×N的DFT矩阵
    X = (F * x(:)).';         % 矩阵乘法计算DFT，并转回行向量
end

function X = fft_codelet_small_prime(x)
    % 小质数长度（如3,5,7,11等）的占位 codelet
    % 当前实现直接调用通用DFT，将来可替换为高度优化的C/MEX实现
    X = dft_direct(x);
end

function X = bluestein_fft(x)
    % Bluestein算法：将N点DFT转换为长度M>=2N-1的卷积，再用混合基FFT计算
    N = length(x);
    n = 0:N-1;

    % 构造 a(n) = x(n) * exp(-j*pi*n^2/N)
    x_row = x(:).';
    a = x_row .* exp(-1j * pi * (n.^2) / N);

    % 构造 chirp 序列 b(n) = exp(j*pi*n^2/N)
    b = exp(1j * pi * (n.^2) / N);

    % 线性卷积长度至少为 2N-1，选择不小于该长度的2的幂
    M = 2^nextpow2(2 * N - 1);

    % 零填充 a 和 b
    a_pad = [a, zeros(1, M - N)];
    b_pad = [b, zeros(1, M - 2 * N + 1), b(N:-1:2)];

    % 使用当前混合基FFT实现卷积：C = IFFT( FFT(a) .* FFT(b) )
    A = my_fft_mix(a_pad);
    B = my_fft_mix(b_pad);
    C = A .* B;
    c = my_ifft_via_fft(C);

    % 取前N个点，并乘以 exp(-j*pi*k^2/N) 得到最终DFT结果
    c = c(1:N);
    X = exp(-1j * pi * (n.^2) / N) .* c;
    X = X(:).';
end

function x = my_ifft_via_fft(X)
    % 使用 my_fft_mix 实现 IFFT：ifft(X) = conj(fft(conj(X))) / N
    N = length(X);
    X_row = X(:).';
    x = conj(my_fft_mix(conj(X_row))) / N;
    x = x(:).';
end

function f = smallest_factor(N)
    % 返回N的最小非平凡因子（2 <= f <= N/2）
    % 如果N是质数，返回空数组[]
    
    % 处理N <= 3的边界情况
    if N <= 3
        f = [];  % 2和3是质数，1没有非平凡因子
        return;
    end
    
    % 首先检查2是否是因子（偶数情况）
    if mod(N, 2) == 0
        f = 2;
        return;
    end
    
    % 检查奇数因子，只需要检查到sqrt(N)
    limit = floor(sqrt(N));  % 最大需要检查的因子
    
    % 从3开始，步进2（只检查奇数）
    for d = 3:2:limit
        if mod(N, d) == 0
            f = d;  % 找到最小因子
            return;
        end
    end
    
    % 如果没找到因子，说明N是质数
    f = [];
end
