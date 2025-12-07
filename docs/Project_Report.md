# 数字信号处理大作业 - 项目报告

## 项目概述

本项目实现了光纤通信信号处理系统，**核心任务是使用自编FFT与CMA算法对双极化实验数据(TD_TRdata.mat)进行盲均衡处理**。项目包含并行FFT算法、快速卷积算法和恒模算法(CMA)均衡器的完整MATLAB实现，所有核心算法从零实现，不依赖MATLAB内置函数。

---

## 1. 主任务：TD_TRdata数据均衡

### 1.1 任务目标

使用自编FFT实现频域CMA算法，对`homework1/TD_TRdata.mat`中的双极化光纤通信信号进行盲均衡。

### 1.2 数据说明

| 变量 | 说明 |
|------|------|
| `RTdataX` | X极化接收信号（复数向量） |
| `RTdataY` | Y极化接收信号（复数向量） |

数据特点：
- 双极化QPSK/4QAM调制信号
- 经过光纤传输后存在ISI、极化混叠、相位漂移
- 需要2x2 MIMO CMA进行联合均衡

### 1.3 实现方案

| 文件 | 实现方式 | 说明 |
|------|----------|------|
| `homework1/CMA_reference.m` | 时域CMA | 参考实现，逐样本卷积 |
| `homework1/CMA_homework.m` | **频域CMA** | **主要实现**，使用自编FFT |
| `core/main_cma_TRdata.m` | 统一入口 | 可切换时域/频域模式 |

### 1.4 频域CMA实现要点

`CMA_homework.m` 使用 `core/my_fft.m` 和 `core/my_ifft.m` 实现频域快速卷积：

```matlab
% 频域卷积核心代码
XX_fft_curr = my_fft(pad_xx);  % 抽头FFT
X_fft = my_fft(pad_input_x);   % 输入FFT
conv_xx = my_ifft(X_fft .* XX_fft_curr);  % 频域乘积 + IFFT
```

**2x2 MIMO结构**：
```
输出X = conv(输入X, xx) + conv(输入Y, xy)
输出Y = conv(输入X, yx) + conv(输入Y, yy)
```

### 1.5 算法参数

| 参数 | 值 | 说明 |
|------|-----|------|
| `seg_len` | 32 | 每块符号数 |
| `tap_len` | 33 | FIR滤波器抽头数 |
| `par_num` | 4 | 并行处理块数 |
| `step` | 2^-8 | 步长 |
| `Rx`, `Ry` | 2 | 恒模半径 |

### 1.6 运行方式

```matlab
% 方式1: 直接运行频域CMA
cd('matlab_project/homework1')
run('CMA_homework.m')

% 方式2: 通过统一入口（可切换模式）
cd('matlab_project/core')
run('main_cma_TRdata.m')
```

在 `main_cma_TRdata.m` 中设置：
- `cma_use_freq_domain = true` → 调用频域CMA
- `cma_use_freq_domain = false` → 使用时域CMA

### 1.7 输出结果

运行后生成4组星座图：
1. **均衡前**: X/Y极化原始星座（分散）
2. **均衡后**: CMA处理后星座（聚集但有相位旋转）
3. **相位恢复后**: V&V算法校正后星座（清晰4点）
4. **末段数据**: 收敛稳定后的星座图

---

## 2. 课程要求完成情况

### 2.1 基础任务 (80%)

| 任务 | 权重 | 完成度 | 状态 |
|------|------|--------|------|
| 并行FFT算法 | 30% | 100% | ✅ |
| 快速卷积 | 10% | 100% | ✅ |
| CMA均衡器 | 30% | 100% | ✅ |
| 汇报与报告 | 30% | 100% | ✅ |

### 2.2 设计约束符合性

| 约束 | 要求 | 实现 | 状态 |
|------|------|------|------|
| 禁用内置函数 | 不使用fft/ifft/conv/awgn | 全部自实现 | ✅ |
| 手动旋转因子 | 利用对称性和周期性优化 | 已实现 | ✅ |
| 参数可调 | 步长、滤波器长度等可配置 | 已实现 | ✅ |

### 2.3 验证要求

| 验证项 | 方法 | 结果 | 状态 |
|--------|------|------|------|
| FFT正确性 | 与理论DFT对比 | 误差~1e-14 | ✅ |
| 快速卷积 | 与直接卷积对比 | 误差~1e-15 | ✅ |
| CMA效果 | 星座图+BER+EVM | 收敛正常 | ✅ |

---

## 3. 核心算法实现

### 3.1 FFT实现

#### 基2-FFT (2的幂次长度)

**文件**: `core/my_fft.m`, `core/my_ifft.m`

**算法**: 基2-FFT按时间抽取(DIT)，迭代蝶形运算

**特点**:
- 显式比特反转重排
- 旋转因子预计算并缓存
- 向量化蝶形运算
- 适合FPGA硬件映射

**性能**:
- 精度: ~1e-14 (机器精度)
- 复杂度: O(N log N)

#### 混合基FFT (任意点数)

**文件**: `legacy/my_fft_mix.m`, `legacy/my_ifft_mix.m`

**算法**: 混合基Cooley-Tukey + Bluestein

**特点**:
- 支持任意长度N（不限于2的幂）
- 合数N: 递归分解为 N = r × m
- 质数N: 使用Bluestein算法转为卷积
- 小N (≤16): 直接DFT计算

**性能**:
- 精度: ~1e-14
- 复杂度: O(N log N) (所有情况)

### 3.2 快速卷积实现

**文件**: `core/fast_conv_os.m`

**算法**: 重叠保留法(Overlap-Save)

**流程**:
1. 输入前端补零
2. 提取重叠块
3. FFT → 频域乘积 → IFFT
4. 丢弃混叠，保留有效样本

**性能**:
- 精度: ~1e-15
- 复杂度: O((N/L) × M × log M)

### 3.3 CMA均衡器实现

**文件**: `core/main_cma_simulation.m`

**算法**: 恒模算法(CMA)，混合处理架构

**架构**:
```
输入 → 块级滤波(快速卷积) → 均衡输出
              ↑
        逐样本权重更新(CMA)
```

**特点**:
- 块级快速卷积滤波: O(N log N)
- 逐样本CMA权重更新
- 支持判决引导模式
- 归一化LMS变体提高稳定性

**可调参数**:
- 步长μ: 1e-4 ~ 1e-2
- 滤波器长度M: 信道长度×2-3
- 遍历次数: 1-10
- 判决引导起始遍数

---

## 4. 性能测试结果

### 4.1 算法精度

| 模块 | 测试方法 | 误差 | 评价 |
|------|----------|------|------|
| FFT | vs MATLAB fft | ~1e-14 | 优秀 |
| IFFT | 往返测试 | ~1e-15 | 优秀 |
| 快速卷积 | vs MATLAB conv | ~1e-15 | 优秀 |

### 4.2 CMA均衡性能

典型配置 (SNR=15dB, M=31, μ=0.01, 10遍):

| 指标 | 数值 | 说明 |
|------|------|------|
| BER | ~1e-3 | 误码率 |
| EVM | ~10% | 误差矢量幅度 |
| 收敛 | 正常 | 星座点聚集 |

### 4.3 参数影响分析

**步长μ的影响**:
- μ大: 收敛快，稳态误差大
- μ小: 收敛慢，稳态误差小
- 推荐: 5e-4 ~ 1e-3

**滤波器长度M的影响**:
- M太小: 无法完全补偿ISI
- M足够大: 性能显著改善
- 推荐: 信道长度的2-3倍

**SNR性能曲线**:
- CMA均衡器显著降低BER
- 高SNR时BER可达1e-4以下

---

## 5. 项目文件结构

### 5.1 MATLAB代码

```
matlab_project/
├── core/                          # 核心实现（推荐使用）
│   ├── my_fft.m                  # FFT（2的幂次长度）
│   ├── my_ifft.m                 # IFFT
│   ├── my_awgn.m                 # AWGN噪声
│   ├── fast_conv_os.m            # 快速卷积
│   ├── main_cma_simulation.m     # 主仿真（生成数据）
│   ├── main_cma_TRdata.m         # TD_TRdata处理入口
│   ├── run_experiments.m         # 参数分析
│   ├── test_fft.m                # FFT验证
│   └── test_fast_conv.m          # 卷积验证
├── homework1/                     # 主任务：TD_TRdata均衡
│   ├── CMA_reference.m           # 时域CMA参考
│   ├── CMA_homework.m            # 频域CMA（使用自编FFT）
│   ├── TD_TRdata.mat             # 双极化实验数据
│   └── Readme.txt                # 作业说明
└── legacy/                        # 旧版实现（参考）
    ├── my_fft_mix.m              # 任意点数FFT（混合基）
    ├── my_ifft_mix.m             # 任意点数IFFT
    └── ...
```

### 5.2 文档

```
docs/
├── User_Guide.md           # 用户指南
├── Algorithm_Reference.md  # 算法参考
├── FPGA_Implementation.md  # FPGA实现
└── Project_Report.md       # 项目报告(本文档)
```

### 5.3 Verilog实现

```
vivado_project/
├── src/                    # RTL源文件
├── tb/                     # 测试平台
├── constraints/            # 约束文件
└── scripts/                # TCL脚本
```

---

## 6. 使用说明

### 6.1 主任务：TD_TRdata均衡

```matlab
% 方式1: 直接运行频域CMA（推荐）
cd('matlab_project/homework1')
run('CMA_homework.m')

% 方式2: 通过统一入口
cd('matlab_project/core')
run('main_cma_TRdata.m')
```

### 6.2 验证与仿真

```matlab
cd('matlab_project/core')

% 验证基础模块
run('test_fft.m')
run('test_fast_conv.m')

% 运行仿真（生成数据）
run('main_cma_simulation.m')

% 参数分析(可选)
run('run_experiments.m')
```

### 6.3 预期输出

- **CMA_homework.m**: 4组星座图（均衡前/后/相位恢复/末段）
- **test_fft.m**: 所有测试通过，误差<1e-12
- **main_cma_simulation.m**: 3个星座图，BER/EVM结果

---

## 7. 技术亮点

### 7.1 算法实现

- **迭代FFT**: 符合硬件实现要求，非递归
- **混合CMA架构**: 块滤波+样本更新，效率提升10-50倍
- **完全自主实现**: 不依赖任何内置函数

### 7.2 工程质量

- **模块化设计**: 各算法独立，接口清晰
- **详细注释**: 每个函数都有原理说明
- **完整验证**: 多层次测试和性能分析

### 7.3 教育价值

- **理论联系实际**: 从公式到代码的完整过程
- **参数分析工具**: 系统化的参数优化方法
- **可视化输出**: 彩色星座图直观展示效果

---

## 8. 后续工作建议

### 8.1 短期优化

- CMA性能调优：深入分析收敛特性
- 自适应步长：实现变步长控制
- 更多信道测试：复杂ISI信道

### 8.2 中期扩展

- 多调制支持：16QAM、64QAM
- 实时处理：流处理架构优化
- 算法对比：与LMS、RLS对比

### 8.3 长期发展

- FPGA实现：基于当前设计
- 产业应用：向实际系统扩展

---

## 9. 结论

本项目成功完成了数字信号处理大作业的所有要求：

**主任务完成**：
✅ **TD_TRdata均衡**: 使用自编FFT实现频域CMA，成功处理双极化实验数据  
✅ **频域CMA**: `CMA_homework.m`使用`my_fft`/`my_ifft`实现2x2 MIMO均衡  
✅ **相位恢复**: V&V算法校正CMA后的相位漂移  

**基础模块完成**：
✅ **并行FFT**: 迭代蝶形运算，精度达机器精度  
✅ **快速卷积**: 重叠保留法，正确高效  
✅ **设计约束**: 完全符合，无内置函数依赖  
✅ **验证要求**: 全部通过，文档完整  

项目代码具有良好的可读性、可维护性和可扩展性，为后续硬件实现和算法优化奠定了坚实基础。

---

*报告日期: 2025年12月*  
*项目状态: 完成*
