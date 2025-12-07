# 算法参考手册

本文档详细介绍项目中实现的核心算法原理与实现细节。

---

## 目录

1. [FFT算法](#1-fft算法)
2. [快速卷积算法](#2-快速卷积算法)
3. [CMA盲均衡算法](#3-cma盲均衡算法)
4. [4QAM调制解调](#4-4qam调制解调)

---

## 1. FFT算法

### 1.1 数学基础

离散傅里叶变换(DFT)定义：
```
X[k] = Σ(n=0 to N-1) x[n] · e^(-j·2π·k·n/N)
```

直接计算复杂度为O(N²)，FFT将其降低到O(N log N)。

本项目提供两种FFT实现：
- **`my_fft.m`** (core/): 基2-FFT，要求输入长度为2的幂，效率最高
- **`my_fft_mix.m`** (legacy/): 混合基FFT，支持任意长度输入

### 1.2 Radix-2 DIT算法

`core/my_fft.m` 采用**基2按时间抽取(DIT)迭代蝶形运算**实现：

**算法步骤**：
1. **比特反转重排**: 将输入按比特反转顺序重新排列
2. **迭代蝶形运算**: log₂(N)级蝶形运算
3. **旋转因子乘法**: 每级使用不同的旋转因子

**蝶形运算**：
```
A' = A + W·B
B' = A - W·B
```
其中 W = e^(-j·2π·k/N) 为旋转因子。

### 1.3 实现特点

**文件**: `core/my_fft.m`

```matlab
function X = my_fft(x)
    % 比特反转重排
    X = x(bitrev_idx);
    
    % 迭代蝶形运算
    for stage = 1:num_stages
        for each butterfly
            temp = X(top);
            X(top) = temp + W * X(bot);
            X(bot) = temp - W * X(bot);
        end
    end
end
```

**优化技术**：
- 旋转因子预计算并缓存
- 利用对称性: W_N^(N-k) = conj(W_N^k)
- 利用周期性: W_N^(k+N) = W_N^k
- 向量化蝶形运算

### 1.4 IFFT实现

利用FFT实现IFFT：
```
IFFT(X) = (1/N) · conj(FFT(conj(X)))
```

**文件**: `core/my_ifft.m`

### 1.5 性能指标

| 指标 | 数值 |
|------|------|
| 精度误差 | ~1e-14 (机器精度) |
| 复杂度 | O(N log N) |
| 内存 | O(N) |

### 1.6 混合基FFT (任意点数)

**文件**: `legacy/my_fft_mix.m`, `legacy/my_ifft_mix.m`

对于非2的幂次长度，使用混合基Cooley-Tukey算法。

#### 算法策略 (Planner)

根据输入长度N选择不同策略：

```
N的类型判断
    │
    ├─ 质数N ─┬─ N ≤ 11 → 直接DFT (codelet)
    │         └─ N > 11 → Bluestein算法
    │
    ├─ 合数N ≤ 16 → 直接DFT
    │
    └─ 合数N > 16 → 混合基分解 (N = r × m)
```

#### 混合基分解

对于合数N = r × m（r为N的最小非平凡因子）：

**步骤1**: 将输入重排为m×r矩阵
```matlab
x_mat = reshape(x, m, r);  % m行r列
```

**步骤2**: 对每列进行m点FFT（递归调用）
```matlab
for j = 1:r
    X_col(:, j) = local_fft_mix(x_mat(:, j));
end
```

**步骤3**: 乘以旋转因子
```matlab
twiddle = exp(-1j * 2 * pi * (k * j) / N);  % m×r矩阵
X_col = X_col .* twiddle;
```

**步骤4**: 对每行进行r点FFT，映射到输出
```matlab
for k1 = 0:(m-1)
    row_fft = local_fft_mix(X_col(k1+1, :));
    for k2 = 0:(r-1)
        idx = k1 + m * k2 + 1;  % 全局频率索引
        X(idx) = row_fft(k2+1);
    end
end
```

#### Bluestein算法 (大质数)

将N点DFT转换为卷积问题：

1. **Chirp调制**: `a(n) = x(n) · exp(-jπn²/N)`
2. **构造chirp序列**: `b(n) = exp(jπn²/N)`
3. **卷积**: 补零到2的幂次长度M，用FFT计算 `c = IFFT(FFT(a) · FFT(b))`
4. **解调**: `X(k) = exp(-jπk²/N) · c(k)`

#### 复杂度分析

| 输入长度 | 算法 | 复杂度 |
|----------|------|--------|
| N = 2^k | Radix-2 | O(N log N) |
| N = 合数 | 混合基 | O(N log N) |
| N = 质数 | Bluestein | O(N log N) |
| N ≤ 16 | 直接DFT | O(N²) |

#### 使用示例

```matlab
% 任意长度FFT
x = randn(1, 100);  % 100点（非2的幂）
X = my_fft_mix(x);

% 验证
X_ref = fft(x);
max_error = max(abs(X - X_ref));  % ~1e-14
```

---

## 2. 快速卷积算法

### 2.1 原理

利用卷积定理：时域卷积等于频域乘积
```
y = x * h  ⟺  Y = X · H
```

### 2.2 重叠保留法 (Overlap-Save)

**文件**: `core/fast_conv_os.m`

**算法步骤**：
1. 在输入前端补(M-1)个零
2. 提取长度为N的重叠块（相邻块重叠M-1个样本）
3. FFT变换到频域
4. 频域相乘
5. IFFT变换回时域
6. 丢弃前(M-1)个混叠样本，保留后L=N-M+1个有效样本

**参数关系**：
- N: FFT点数（必须是2的幂）
- M: 滤波器长度
- L = N - M + 1: 每块有效输出样本数

```matlab
function y = fast_conv_os(x_long, h, nfft)
    L = nfft - M + 1;           % 有效数据块长度
    H = my_fft(h_padded);       % 滤波器频域表示
    
    for each block
        X_block = my_fft(current_block);
        Y_block = X_block .* H;
        y_block = my_ifft(Y_block);
        valid_samples = y_block(M:end);  % 丢弃混叠部分
    end
end
```

### 2.3 重叠相加法 (Overlap-Add)

**文件**: `legacy/fast_conv_add.m`

**算法步骤**：
1. 将输入分成不重叠的L点块
2. 每块补零到N点
3. FFT → 频域乘积 → IFFT
4. 相邻块输出在重叠区域相加

### 2.4 方法对比

| 特性 | 重叠保留法 | 重叠相加法 |
|------|------------|------------|
| 输入块 | 重叠 | 不重叠 |
| 输出处理 | 丢弃混叠 | 重叠相加 |
| 适用场景 | 实时/流处理 | 批处理 |
| 内存使用 | 较低 | 较高 |

---

## 3. CMA盲均衡算法

### 3.1 算法原理

恒模算法(CMA)是一种**盲均衡**技术，利用信号的恒模特性消除ISI，无需训练序列。

**代价函数**：
```
J = E[(|y(n)|² - R²)²]
```
其中：
- y(n): 均衡器输出
- R²: 恒模半径（对于归一化4QAM，R² = 1）

**权重更新**：
```
w(n+1) = w(n) + μ · e(n) · x*(n)
e(n) = y(n) · (R² - |y(n)|²)
```

### 3.2 实现架构

本项目采用**混合处理架构**：

```
输入信号 → 输入缓冲区
              ↓
        ┌─────────────┐
        │  块级滤波    │ ← 快速卷积 O(N log N)
        │ (fast_conv) │
        └─────────────┘
              ↓
        均衡输出
              ↓
        ┌─────────────┐
        │ 逐样本更新   │ ← CMA权重更新
        │ (权重更新)   │
        └─────────────┘
              ↓
        更新后的权重
```

**优势**：
- 滤波操作利用FFT加速：O(N²) → O(N log N)
- 权重更新保持样本级精度
- 速度提升10-50倍

### 3.3 关键参数

| 参数 | 符号 | 典型值 | 影响 |
|------|------|--------|------|
| 步长 | μ | 1e-4 ~ 1e-2 | 大→快但不稳，小→慢但稳 |
| 滤波器长度 | M | 信道长度×2-3 | 太小无法均衡 |
| 恒模半径 | R² | 1 (4QAM) | 取决于调制方式 |
| 遍历次数 | - | 1-10 | 多遍提高收敛质量 |

### 3.4 相位模糊问题

CMA对相位不敏感，均衡后存在相位旋转。解决方法：

**导频辅助相位估计**：
```matlab
% 硬判决得到参考符号
pilot_reference = nearest_constellation_point(pilot_symbols);

% 估计相位偏移
phase_offset = angle(sum(pilot_reference .* conj(pilot_symbols)));

% 校正
corrected = equalized * exp(1j * phase_offset);
```

### 3.5 归一化LMS变体

为提高稳定性，使用归一化步长：
```matlab
x_power = x_vec' * x_vec + epsilon;
mu_norm = mu / x_power;
w = w + mu_norm * error * conj(x_vec);
```

### 3.6 判决引导模式

收敛后可切换到判决引导(DD)模式提高性能：
```matlab
if pass >= dd_start_pass
    decision = nearest_constellation_point(y_out);
    error = conj(decision - y_out);
else
    error = (R2 - abs(y_out)^2) * conj(y_out);
end
```

### 3.7 2x2 MIMO CMA (双极化处理)

**文件**: `homework1/CMA_homework.m`, `homework1/CMA_reference.m`

用于处理双极化光纤通信信号（如`TD_TRdata.mat`），需要2x2 MIMO结构联合均衡。

**信号模型**：
```
[输出X]   [hxx  hxy] [输入X]
[输出Y] = [hyx  hyy] [输入Y]
```

**均衡器结构**：
```
输出X = conv(输入X, xx) + conv(输入Y, xy)
输出Y = conv(输入X, yx) + conv(输入Y, yy)
```

**频域实现** (`CMA_homework.m`)：
```matlab
% 抽头FFT
XX_fft = my_fft(pad_xx);
XY_fft = my_fft(pad_xy);
YX_fft = my_fft(pad_yx);
YY_fft = my_fft(pad_yy);

% 输入FFT
X_fft = my_fft(pad_input_x);
Y_fft = my_fft(pad_input_y);

% 频域卷积
conv_xx = my_ifft(X_fft .* XX_fft);
conv_xy = my_ifft(Y_fft .* XY_fft);
conv_yx = my_ifft(X_fft .* YX_fft);
conv_yy = my_ifft(Y_fft .* YY_fft);

% 组合输出
block_out_x = conv_xx + conv_xy;
block_out_y = conv_yx + conv_yy;
```

**权重更新**：
```matlab
errx = Rx - abs(block_out_x).^2;
erry = Ry - abs(block_out_y).^2;

xx_accu = xx_accu + step * errx .* block_out_x * conj(input_x);
xy_accu = xy_accu + step * errx .* block_out_x * conj(input_y);
yx_accu = yx_accu + step * erry .* block_out_y * conj(input_x);
yy_accu = yy_accu + step * erry .* block_out_y * conj(input_y);
```

**参数配置**：
| 参数 | 值 | 说明 |
|------|-----|------|
| `seg_len` | 32 | 每块符号数 |
| `tap_len` | 33 | FIR滤波器抽头数 |
| `par_num` | 4 | 并行处理块数 |
| `step` | 2^-8 | 步长 |
| `Rx`, `Ry` | 2 | 恒模半径 |

### 3.8 V&V相位恢复

CMA均衡后存在相位漂移，使用Viterbi-Viterbi算法恢复：

```matlab
% 四次方消除调制
z = y.^4;

% 滑动平均滤波
win = ones(2*win_half+1, 1) / (2*win_half+1);
f = conv(z, win, 'same');

% 估计相位
theta = angle(f) / 4;

% 相位校正
y_corrected = y .* exp(-1j * theta);
```

---

## 4. 4QAM调制解调

### 4.1 星座映射

4QAM星座点（归一化功率为1）：
```
符号0: (1+j)/√2   → 比特 00
符号1: (1-j)/√2   → 比特 01
符号2: (-1-j)/√2  → 比特 11
符号3: (-1+j)/√2  → 比特 10
```

### 4.2 格雷编码

使用格雷码映射减少误码扩散：
```matlab
grayMap = [0 1 3 2];  % 00→0, 01→1, 11→3, 10→2
constellation_map = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);
```

### 4.3 调制流程

```matlab
% 比特到符号索引
txSymbolIndices = bit1 * 2 + bit2;

% 格雷编码
grayIndices = grayMap(txSymbolIndices + 1);

% 映射到星座点
txSymbols = constellation_map(grayIndices + 1);
```

### 4.4 解调流程

```matlab
% 硬判决：找最近星座点
[~, rxSymbolIndices] = min(abs(received - constellation_map));

% 格雷解码
rxGrayIndices = rxSymbolIndices - 1;
rxOriginalIndices = grayInverse(rxGrayIndices + 1);

% 符号到比特
bit1 = bitshift(rxOriginalIndices, -1);
bit2 = bitand(rxOriginalIndices, 1);
```

---

## 5. 性能指标计算

### 5.1 误码率 (BER)

```matlab
BER = sum(txBits ~= rxBits) / length(txBits);
```

### 5.2 误差矢量幅度 (EVM)

```matlab
error_vector = equalized - decided_symbols;
EVM = sqrt(mean(abs(error_vector).^2)) / sqrt(mean(abs(decided_symbols).^2)) * 100;
```

### 5.3 延迟对齐

由于信道和均衡器引入延迟，需要对齐发射和接收序列：
```matlab
% 搜索最佳延迟
for delay = 0:max_delay
    test_ber = calculate_ber(txBits(delay+1:end), rxBits);
    if test_ber < best_ber
        best_delay = delay;
    end
end
```

---

## 参考文献

1. Cooley, J.W. and Tukey, J.W., "An Algorithm for the Machine Calculation of Complex Fourier Series"
2. Godard, D.N., "Self-Recovering Equalization and Carrier Tracking in Two-Dimensional Data Communication Systems"
3. Oppenheim, A.V. and Schafer, R.W., "Discrete-Time Signal Processing"
4. Haykin, S., "Adaptive Filter Theory"

---

*最后更新: 2025年12月*
