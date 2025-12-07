# 用户指南

本指南帮助您快速上手CMA盲均衡器MATLAB仿真系统。

---

## 目录

1. [环境要求](#1-环境要求)
2. [快速开始](#2-快速开始)
3. [运行仿真](#3-运行仿真)
4. [理解输出结果](#4-理解输出结果)
5. [自定义参数](#5-自定义参数)
6. [彩色星座图说明](#6-彩色星座图说明)
7. [常见问题](#7-常见问题)

---

## 1. 环境要求

### MATLAB环境
- **版本**: MATLAB R2018b 或更高
- **工具箱**: 无需额外工具箱（所有函数从零实现）
- **内存**: 推荐至少4GB RAM

### 项目结构
```
matlab_project/
├── core/                    # 核心实现（推荐使用）
│   ├── my_fft.m            # FFT实现（2的幂次长度）
│   ├── my_ifft.m           # IFFT实现
│   ├── my_awgn.m           # AWGN噪声
│   ├── fast_conv_os.m      # 快速卷积
│   ├── main_cma_simulation.m   # 主仿真（生成数据）
│   ├── main_cma_TRdata.m   # TD_TRdata处理入口
│   ├── run_experiments.m   # 参数分析实验
│   ├── test_fft.m          # FFT验证
│   └── test_fast_conv.m    # 卷积验证
├── homework1/              # 主任务：TD_TRdata均衡
│   ├── CMA_reference.m     # 时域CMA参考
│   ├── CMA_homework.m      # 频域CMA（使用自编FFT）
│   ├── TD_TRdata.mat       # 双极化实验数据
│   └── Readme.txt          # 作业说明
└── legacy/                 # 旧版实现（仅供参考）
```

---

## 2. 快速开始

### 步骤1: 切换到工作目录
```matlab
cd('matlab_project/core')
```

### 步骤2: 验证FFT/IFFT
```matlab
run('test_fft.m')
```

预期输出：
```
========================================
FFT和IFFT验证测试
========================================
--- 测试长度 N = 64 ---
  测试1: 单位冲激... 通过 (误差: 1.23e-14)
  ...
========================================
所有测试通过! ✓
========================================
```

### 步骤3: 验证快速卷积
```matlab
run('test_fast_conv.m')
```

### 步骤4: 运行主仿真
```matlab
run('main_cma_simulation.m')
```

---

## 3. 运行仿真

### 3.1 主仿真脚本

`main_cma_simulation.m` 执行完整的4QAM CMA盲均衡仿真：

1. 生成4QAM符号（格雷编码）
2. 通过ISI信道
3. 添加AWGN噪声
4. CMA盲均衡
5. 相位校正
6. 解调与性能评估

**运行命令**:
```matlab
run('main_cma_simulation.m')
```

**输出**:
- 控制台：配置参数、处理进度、BER/EVM结果
- 图形：3个星座图对比（接收/均衡/校正）

### 3.2 参数影响分析

`run_experiments.m` 执行三组实验：

| 实验 | 内容 | 输出 |
|------|------|------|
| 实验1 | 步长μ的影响 | 4个星座图对比 |
| 实验2 | 滤波器长度M的影响 | 4个星座图对比 |
| 实验3 | SNR性能曲线 | BER vs SNR曲线 |

**运行命令**:
```matlab
run('run_experiments.m')
```

**运行时间**: 约3-5分钟

### 3.3 TD_TRdata实验数据处理（主任务）

本项目的核心任务是使用自编FFT实现频域CMA，对`TD_TRdata.mat`进行盲均衡。

#### 方式1: 直接运行频域CMA（推荐）
```matlab
cd('matlab_project/homework1')
run('CMA_homework.m')
```

#### 方式2: 通过统一入口
```matlab
cd('matlab_project/core')
run('main_cma_TRdata.m')
```

在`main_cma_TRdata.m`中可切换模式：
- `cma_use_freq_domain = true` → 频域CMA（使用自编FFT）
- `cma_use_freq_domain = false` → 时域CMA（参考实现）

#### 输出说明
运行后生成4组星座图：
1. **均衡前**: X/Y极化原始星座（分散，无法区分符号）
2. **均衡后**: CMA处理后星座（聚集但有相位旋转）
3. **相位估计**: V&V算法估计的相位漂移曲线
4. **末段数据**: 收敛稳定后的清晰星座图

---

## 4. 理解输出结果

### 4.1 性能指标

#### BER (Bit Error Rate) - 误码率
- **定义**: 错误比特数 / 总比特数
- **目标值**: < 1e-3（可接受）, < 1e-4（良好）

#### EVM (Error Vector Magnitude) - 误差矢量幅度
- **定义**: 实际符号与理想符号的均方根误差（百分比）
- **目标值**: < 10%（可接受）, < 5%（良好）

### 4.2 星座图解读

| 阶段 | 图位置 | 预期现象 |
|------|--------|----------|
| 接收信号 | 左图 | 星座点严重分散，无法区分符号 |
| CMA均衡后 | 中图 | 星座点聚集，但存在相位旋转 |
| 相位校正后 | 右图 | 清晰的4个点云，对齐到理想位置 |

---

## 5. 自定义参数

### 5.1 主要参数

编辑 `main_cma_simulation.m` 中的 `config` 结构体：

```matlab
config = struct();

% 基本参数
config.numSymbols = 300000;     % 符号数量
config.snr_dB = 15;             % 信噪比 (dB)

% 信道参数
config.channelTaps = [1, 0.5, 0.2, 0.1, 0.5, 0.7];  % ISI信道

% CMA均衡器参数
config.cma.filterLength = 31;   % 均衡器长度
config.cma.stepSize = 0.01;     % 步长
config.cma.numPasses = 10;      % 遍历次数
config.cma.R2 = 1;              % 恒模半径

% 快速卷积参数
config.conv.blockLength = 256;  % 块长度
```

### 5.2 参数选择指南

| 参数 | 推荐范围 | 影响 |
|------|----------|------|
| `numSymbols` | 10000-300000 | 越大越准确，但越慢 |
| `snr_dB` | 15-30 | 信噪比，影响噪声水平 |
| `filterLength` | 信道长度×2-3 | 太小无法均衡，太大增加计算量 |
| `stepSize` | 1e-4 ~ 1e-2 | 大则快但不稳，小则慢但稳 |
| `numPasses` | 1-10 | 多遍迭代提高收敛质量 |

### 5.3 测试不同信道

```matlab
% 轻微ISI
config.channelTaps = [1, 0.2];

% 中等ISI
config.channelTaps = [1, 0.5, 0.2];

% 严重ISI
config.channelTaps = [1, 0.8, 0.5, 0.3];

% 复数信道
config.channelTaps = [1, 0.3*exp(1j*pi/6), 0.1*exp(-1j*pi/4)];
```

---

## 6. 彩色星座图说明

### 6.1 颜色映射

不同4QAM符号用不同颜色标记，便于观察均衡效果：

| 比特 | 符号 | 星座点 | 颜色 | 象限 |
|------|------|--------|------|------|
| 00 | 0 | (1+j)/√2 | 🔴 红色 | I |
| 01 | 1 | (1-j)/√2 | 🟢 绿色 | IV |
| 11 | 2 | (-1-j)/√2 | 🔵 蓝色 | III |
| 10 | 3 | (-1+j)/√2 | 🟡 黄色 | II |

### 6.2 观察要点

- **接收信号**: 各颜色点混叠分散 → ISI影响严重
- **CMA均衡后**: 各颜色点聚集但有旋转 → 均衡有效，存在相位模糊
- **相位校正后**: 各颜色点准确对齐 → 可正确解调

---

## 7. 常见问题

### Q1: "输入长度必须是2的整数次幂"错误
**解决**: 检查FFT块长度是否为2的幂（如256, 512, 1024）

### Q2: BER很高或星座图混乱
**解决**:
```matlab
config.cma.stepSize = 1e-4;      % 减小步长
config.cma.filterLength = 51;    % 增加滤波器长度
config.snr_dB = 30;              % 提高SNR
```

### Q3: 仿真运行很慢
**解决**:
```matlab
config.numSymbols = 5000;        % 减少符号数
config.conv.blockLength = 512;   % 增大块长度
```

### Q4: 星座图有旋转未校正
**解决**: 增加瞬态样本数或检查相位校正代码

### Q5: 如何保存结果
```matlab
saveas(gcf, 'results.png');
save('data.mat', 'ber', 'evm', 'equalizedSymbols_corrected');
```

---

## 附录：文件功能速查

| 文件 | 功能 | 运行方式 |
|------|------|----------|
| `homework1/CMA_homework.m` | **频域CMA（主任务）** | `run('CMA_homework.m')` |
| `core/main_cma_TRdata.m` | TD_TRdata处理入口 | `run('main_cma_TRdata.m')` |
| `core/test_fft.m` | FFT/IFFT验证 | `run('test_fft.m')` |
| `core/test_fast_conv.m` | 快速卷积验证 | `run('test_fast_conv.m')` |
| `core/main_cma_simulation.m` | 仿真数据生成 | `run('main_cma_simulation.m')` |
| `core/run_experiments.m` | 参数分析 | `run('run_experiments.m')` |

---

*最后更新: 2025年12月*
