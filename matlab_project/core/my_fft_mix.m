function X = my_fft_mix(x)
% =========================================================================
%            混合基快速傅里叶变换 (Mixed-Radix FFT) 实现
% =========================================================================
%
% 【算法概述】
%   本函数实现任意长度 N 的 FFT，不要求 N 为 2 的整数次幂。
%   设计思路类似 FFTW 的 "Planner"：先判断 N 的结构，再选择合适的算法路径。
%
%   - N 为合数: 优先使用混合基 Cooley-Tukey 分解，将 N 拆成 r×m 的小 FFT 组合
%   - N 为质数且较大: 使用 Bluestein (Chirp-Z) 将 DFT 转成卷积，再用 FFT 加速
%   - N 很小: 直接 DFT (O(N^2))，避免递归/卷积的常数开销
%
%   与本仓库的 `my_fft.m` (基2 DIT) 相比：
%   - `my_fft.m` 追求 "原位(in-place) + 规则蝶形"，适合 2^m 点与硬件映射
%   - `my_fft_mix.m` 追求 "任意点数可用"，以递归与数据重排换取通用性
%
% 【Planner 策略】
%   ┌────────────────────────────────────────────────────────────────┐
%   │  输入长度 N                                                     │
%   │      │                                                         │
%   │      ├── N 是质数 ───┬── N ≤ 11 ──→ Codelet (直接DFT)           │
%   │      │               └── N > 11 ──→ Bluestein O(N logN)        │
%   │      │                                                         │
%   │      └── N 是合数 ───┬── N ≤ 16 ──→ 直接 DFT O(N²)              │
%   │                      └── N > 16 ──→ 混合基分解 (递归)           │
%   └────────────────────────────────────────────────────────────────┘
%
% 【混合基 Cooley-Tukey 算法】
%   对于合数 N = r × m (本实现取 r 为最小因子，使递归更快收敛)，将 N 点 DFT 分解为:
%   1) "列变换"：把输入 reshape 成 m×r，对每一列做 m 点 FFT
%   2) "旋转因子"：对中间结果乘 twiddle[k,j] = exp(-j·2π·k·j/N)
%   3) "行变换"：对每一行做 r 点 FFT，并按输出索引规则写回一维序列
%
%   可视化示例 (N=80=5×16，本实现取最小因子 r=5, m=16；输入用 x0..x79 表示):
%   ---------------------------------------------------------------------
%   设输入为行向量:
%       x = [x0 x1 x2 ... x79]
%
%   (1) 重排: x_mat = reshape(x, m, r) 得到 16×5 矩阵 (MATLAB 按列填充):
%       第 1 列: [x0  x1  ... x15]'
%       第 2 列: [x16 x17 ... x31]'
%       第 3 列: [x32 x33 ... x47]'
%       第 4 列: [x48 x49 ... x63]'
%       第 5 列: [x64 x65 ... x79]'
%     等价索引关系:
%       x_mat(n0+1, n1+1) = x( m*n1 + n0 + 1 )
%       其中 n0=0..15, n1=0..4。
%
%   (2) 列 FFT: 对每列做 16 点 FFT，得到 X_col(16×5)
%
%   (3) 旋转因子: 对 X_col(k+1, j+1) 乘 twiddle(k,j)
%       twiddle(k,j) = exp(-j*2*pi*(k*j)/80)
%       - k = 0..15 (行索引)
%       - j = 0..4  (列索引)
%       乘法按元素进行，twiddle 的尺寸同为 16×5。
%
%   (4) 行 FFT: 对每一行做 5 点 FFT，得到 row_fft(1×5)
%
%   回写输出: 使用 idx = r*k1 + k2 + 1 (r=5)
%       - k1=0..15 对应 "列 FFT 输出索引"
%       - k2=0..4  对应 "行 FFT 输出索引"
%       输出 X 的下标为:
%         k1=0 时写入: X(1), X(2), X(3), X(4), X(5)      即 k=0,1,2,3,4
%         k1=1 时写入: X(6), X(7), X(8), X(9), X(10)     即 k=5,6,7,8,9
%         ...
%       即连续写入，符合 DFT 输出的自然顺序。
%
%   可视化示例 (N=15=3×5，取 r=3, m=5；输入用符号 a..o 表示):
%       x = [a b c d e  f g h i j  k l m n o]
%       reshape(x, m, r) 得到 5×3 矩阵 (MATLAB 按列填充):
%           [ a  f  k
%             b  g  l
%             c  h  m
%             d  i  n
%             e  j  o ]
%       - 第 1 列 [a b c d e] 做 5 点 FFT
%       - 第 2 列 [f g h i j] 做 5 点 FFT
%       - 第 3 列 [k l m n o] 做 5 点 FFT
%       之后对第 k 行(0..m-1)做 3 点 FFT，并用 idx = k + m*k2 + 1 映射回输出
%
%   数学原理:
%       X[k₀ + m·k₁] = Σ(n₁) { W_N^(k₀·n₁) ·
%                      [ Σ(n₀) x[m·n₁+n₀] · W_m^(n₀·k₀) ] } · W_r^(n₁·k₁)
%
% 【Bluestein 算法 (处理质数长度)】
%   对于较大的质数 N，Cooley-Tukey 无法因式分解。
%   Bluestein 将 N 点 DFT 转换为长度 M(取 2 的幂) 的卷积问题:
%
%       X[k] = exp(-jπk^2/N) · (a ⊛ b)[k]
%       a[n] = x[n] · exp(-jπn^2/N)
%       b[n] = exp(+jπn^2/N)
%
%   其中 (a ⊛ b) 可通过 FFT 卷积实现:
%       (a ⊛ b) = IFFT( FFT(a_pad) .* FFT(b_pad) )
%
%   备注:
%   - Bluestein 会引入零填充与额外 FFT/ IFFT，常数因子较大
%   - 因此本实现仅在 "大质数" 场景使用 Bluestein，小质数直接 DFT 更划算
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

% =====================================================================
% 步骤 0: 获取序列长度
% =====================================================================
% 与 `my_fft.m` 的 "要求 N=2^m" 不同，本函数允许任意 N

N = length(x);

% -----------------------------------------------------------------
% 情况 1: 质数长度
% -----------------------------------------------------------------
% 质数 N 无法做 Cooley-Tukey 因式分解。
% 这里做一个 "常数项" 层面的权衡：
%   - N 很小: 直接 DFT 反而更快
%   - N 较大: Bluestein 通过卷积 + FFT 把复杂度降到 O(N log N)
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
% 虽然合数可以分解，但当 N 很小时，递归/重排/旋转因子的开销可能超过收益
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

% =====================================================================
% Cooley-Tukey 混合基 FFT (DIF 风格)
% =====================================================================
% 对于 N = N1 × N2 (这里 N1=m, N2=r)，DIF 分解公式为：
%   X[k1 + N1*k2] = Σ_{n2=0}^{N2-1} W_{N2}^{n2*k2} 
%                   [ W_N^{n2*k1} Σ_{n1=0}^{N1-1} x[N2*n1 + n2] W_{N1}^{n1*k1} ]
%
% 索引关系：
%   输入: n = N2*n1 + n2 = r*n1 + n2  (n1=0..m-1, n2=0..r-1)
%   输出: k = k1 + N1*k2 = k1 + m*k2  (k1=0..m-1, k2=0..r-1)
%
% 步骤：
% 1. 将输入按 x_mat(n1, n2) = x[r*n1 + n2] 排列成 m×r 矩阵
% 2. 对每列 (固定 n2) 做 m 点 FFT：得到 Σ_{n1} x[r*n1+n2] W_m^{n1*k1}
% 3. 乘以 twiddle: W_N^{n2*k1} = exp(-j*2π*n2*k1/N)
% 4. 对每行 (固定 k1) 做 r 点 FFT：得到 Σ_{n2} [...] W_r^{n2*k2}
% 5. 输出 X[k1 + m*k2]
%
% 关键：输入索引是 x[r*n1 + n2]，而 MATLAB reshape 按列填充
%       reshape(x, r, m) 得到 r×m 矩阵，其中第 (n2+1, n1+1) 元素 = x[r*n1 + n2 + 1]
%       转置后得到 m×r 矩阵，第 (n1+1, n2+1) 元素 = x[r*n1 + n2 + 1]

% 步骤 1: 数据重排
% reshape(x, r, m) 产生 r×m 矩阵:
%   列 n1 包含 [x[r*n1], x[r*n1+1], ..., x[r*n1+r-1]]
% 转置得到 m×r 矩阵:
%   行 n1 包含 [x[r*n1], x[r*n1+1], ..., x[r*n1+r-1]]
%   即 x_mat(n1+1, n2+1) = x[r*n1 + n2 + 1]
x_mat = reshape(x, r, m).';  % m×r 矩阵

% 步骤 2: 列 FFT (m 点)
% 对每列 (固定 n2) 做 m 点 FFT
% X_col(k1+1, n2+1) = Σ_{n1=0}^{m-1} x_mat(n1+1, n2+1) * W_m^{n1*k1}
X_col = zeros(m, r);
for n2 = 1:r
    X_col(:, n2) = local_fft_mix(x_mat(:, n2).').';
end

% 步骤 3: 旋转因子 (Twiddle) 相乘
% twiddle(k1, n2) = W_N^{n2*k1} = exp(-j*2π*n2*k1/N)
% 注意：n2 是从 0 开始的数学索引，对应 MATLAB 的列索引 n2+1
k1_idx = (0:m-1).';  % m×1
n2_idx = 0:(r-1);     % 1×r
twiddle = exp(-1j * 2 * pi * (k1_idx * n2_idx) / N);  % m×r
X_col = X_col .* twiddle;

% 步骤 4: 行 FFT (r 点) + 输出索引映射
% 对每行 (固定 k1) 做 r 点 FFT
% 输出 X[k1 + m*k2] = row_fft(k2+1)
X = zeros(1, N);
for k1 = 0:(m-1)
    row = X_col(k1+1, :);
    row_fft = local_fft_mix(row);
    for k2 = 0:(r-1)
        idx = k1 + m * k2 + 1;  % k = k1 + N1*k2 = k1 + m*k2
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

