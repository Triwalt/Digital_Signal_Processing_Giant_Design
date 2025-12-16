function X = my_fft_mix(x)
% =========================================================================
%            混合基快速傅里叶变换 (Mixed-Radix FFT) 实现
% =========================================================================
%
% 【算法概述】
%   本函数实现了任意长度 N 的混合基 FFT，不要求输入长度是 2 的整数次幂。
%   采用类似 FFTW 库的 Planner 策略，根据输入长度自动选择最优算法:
%   - 合数长度: 使用混合基 Cooley-Tukey 分解
%   - 质数长度: 使用 Bluestein (Chirp-Z) 算法
%   - 小长度: 直接 DFT 计算
%
% 【Planner 策略】
%   ┌────────────────────────────────────────────────────────────────┐
%   │  输入长度 N                                                     │
%   │      │                                                         │
%   │      ├── N 是质数 ───┬── N ≤ 11 ──→ Codelet (直接DFT)          │
%   │      │               └── N > 11 ──→ Bluestein O(N logN)        │
%   │      │                                                         │
%   │      └── N 是合数 ───┬── N ≤ 16 ──→ 直接 DFT O(N²)             │
%   │                      └── N > 16 ──→ 混合基分解 (递归)           │
%   └────────────────────────────────────────────────────────────────┘
%
% 【混合基 Cooley-Tukey 算法】
%   对于合数 N = r × m，将 N 点 DFT 分解为:
%   1. r 个 m 点 DFT (列变换)
%   2. 旋转因子相乘: twiddle[k,j] = exp(-j·2π·k·j/N)
%   3. m 个 r 点 DFT (行变换)
%
%   数学原理:
%       X[k₀ + m·k₁] = Σ(n₁) { W_N^(k₀·n₁) ·
%                      [ Σ(n₀) x[m·n₁+n₀] · W_m^(n₀·k₀) ] } · W_r^(n₁·k₁)
%
% 【Bluestein 算法 (处理质数长度)】
%   将 DFT 转换为卷积: X[k] = exp(-jπk²/N) · (a ⊛ b)[k]
%   其中: a[n] = x[n]·exp(-jπn²/N), b[n] = exp(jπn²/N)
%   卷积用 FFT 加速: O(N log N) 复杂度
%
% 【复杂度分析】
%   N = 2^m       : O(N log N)
%   N = r₁·r₂·... : O(N log N)
%   N 为大质数    : O(N log N) (常数因子较大)
%   N ≤ 16        : O(N²)
%
% 【输入参数】
%   x - 输入信号，复数或实数向量，长度 N 可为任意正整数
%       支持行向量或列向量，输出保持相同方向
%
% 【输出参数】
%   X - x 的 N 点 DFT 结果，方向与输入相同
%
% 【使用示例】
%   n = 0:14;  % N = 15 = 3 × 5
%   x = cos(2*pi*3*n/15) + 0.5*sin(2*pi*5*n/15);
%   X1 = my_fft_mix(x);
%   X2 = fft(x);
%   max_err = max(abs(X1 - X2));
%   fprintf('与 fft() 的最大差异: %g\n', max_err);
%
% 【参考文献】
%   [1] Cooley, J.W., Tukey, J.W. (1965). "An algorithm for the machine
%       calculation of complex Fourier series"
%   [2] Bluestein, L.I. (1970). "A linear filtering approach to the
%       computation of discrete Fourier transform"
%   [3] Frigo, M., Johnson, S.G. (2005). "The Design and Implementation
%       of FFTW3"
%
% =========================================================================

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
    X = X.';
end
end

% =========================================================================
% 递归混合基 FFT 核心函数
% =========================================================================
function X = local_fft_mix(x)
% local_fft_mix: 混合基 FFT 的递归实现
%   根据输入长度 N 自动选择算法:
%   - 小质数 (≤11): codelet
%   - 大质数 (>11): Bluestein
%   - 合数 (≤16): 直接 DFT
%   - 合数 (>16): 混合基 Cooley-Tukey

N = length(x);

% -----------------------------------------------------------------
% 情况 1: 质数长度
% -----------------------------------------------------------------
if isprime(N)
    if N <= 11
        % 小质数: codelet (可替换为优化实现)
        X = fft_codelet_small_prime(x);
    else
        % 大质数: Bluestein 算法
        X = bluestein_fft(x);
    end
    return;
end

% -----------------------------------------------------------------
% 情况 2: 极小合数
% -----------------------------------------------------------------
if N <= 16
    X = dft_direct(x);
    return;
end

% -----------------------------------------------------------------
% 情况 3: 合数 N > 16 - 混合基 Cooley-Tukey
% -----------------------------------------------------------------
r = smallest_factor(N);  % 找最小因子

if isempty(r) || r == 1 || r == N
    X = dft_direct(x);
    return;
end

m = N / r;  % N = r × m

% 步骤 1: 重排为 m×r 矩阵，对每列做 m 点 FFT
x_mat = reshape(x, m, r);
X_col = zeros(m, r);
for j = 1:r
    X_col(:, j) = local_fft_mix(x_mat(:, j).').';
end

% 步骤 2: 乘以旋转因子
k = (0:m-1).';
j_idx = 0:(r-1);
twiddle = exp(-1j * 2 * pi * (k * j_idx) / N);
X_col = X_col .* twiddle;

% 步骤 3: 对每行做 r 点 FFT，映射到输出
X = zeros(1, N);
for k1 = 0:(m-1)
    row = X_col(k1+1, :);
    row_fft = local_fft_mix(row);
    for k2 = 0:(r-1)
        idx = k1 + m * k2 + 1;
        X(idx) = row_fft(k2+1);
    end
end
end

% =========================================================================
% 直接 DFT 实现 - O(N²)
% =========================================================================
function X = dft_direct(x)
% dft_direct: 矩阵乘法实现的直接 DFT
%   X[k] = Σ x[n] · W_N^(nk), W_N = exp(-j·2π/N)

N = length(x);
n = 0:N-1;
k = n.';
W = exp(-1j * 2 * pi / N);
F = W .^ (k * n);
X = (F * x(:)).';
end

% =========================================================================
% 小质数 Codelet
% =========================================================================
function X = fft_codelet_small_prime(x)
% fft_codelet_small_prime: 小质数长度的 FFT
%   当前实现: 直接 DFT
%   可替换为: Winograd 或 C/MEX 优化实现

X = dft_direct(x);
end

% =========================================================================
% Bluestein FFT 算法
% =========================================================================
function X = bluestein_fft(x)
% bluestein_fft: Bluestein (Chirp-Z) 算法
%   将 N 点 DFT 转为卷积，用 FFT 加速
%   复杂度: O(N log N)

N = length(x);
n = 0:N-1;

% 调制序列: a[n] = x[n] · exp(-jπn²/N)
x_row = x(:).';
a = x_row .* exp(-1j * pi * (n.^2) / N);

% Chirp 序列: b[n] = exp(jπn²/N)
b = exp(1j * pi * (n.^2) / N);

% 卷积长度: 2 的幂
M = 2^nextpow2(2 * N - 1);

% 零填充
a_pad = [a, zeros(1, M - N)];
b_pad = [b, zeros(1, M - 2 * N + 1), b(N:-1:2)];

% FFT 卷积
A = my_fft_mix(a_pad);
B = my_fft_mix(b_pad);
C = A .* B;
c = my_ifft_via_fft(C);

% 后处理
c = c(1:N);
X = exp(-1j * pi * (n.^2) / N) .* c;
X = X(:).';
end

% =========================================================================
% 利用 FFT 实现 IFFT
% =========================================================================
function x = my_ifft_via_fft(X)
% my_ifft_via_fft: IFFT(X) = conj(FFT(conj(X))) / N

N = length(X);
X_row = X(:).';
x = conj(my_fft_mix(conj(X_row))) / N;
x = x(:).';
end

% =========================================================================
% 寻找最小因子
% =========================================================================
function f = smallest_factor(N)
% smallest_factor: 返回 N 的最小非平凡因子 (2 ≤ f ≤ √N)
%   若 N 是质数，返回空数组 []

if N <= 3
    f = [];
    return;
end

if mod(N, 2) == 0
    f = 2;
    return;
end

limit = floor(sqrt(N));
for d = 3:2:limit
    if mod(N, d) == 0
        f = d;
        return;
    end
end

f = [];
end
