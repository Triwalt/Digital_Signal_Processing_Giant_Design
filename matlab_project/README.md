# MATLAB项目 - CMA盲均衡算法实现

## 项目结构

```
matlab_project/
├── core/                          # 核心实现(符合技术规格书)
│   ├── my_fft.m                  # 自定义FFT实现(迭代蝶形运算)
│   ├── my_ifft.m                 # 自定义IFFT实现
│   ├── fast_conv_os.m            # 快速卷积(重叠保留法)
│   ├── test_fft.m                # FFT/IFFT验证脚本
│   ├── test_fast_conv.m          # 快速卷积验证脚本
│   ├── main_cma_simulation.m     # 主仿真脚本
│   └── run_experiments.m         # 参数影响分析实验
└── Gemini_generated_simulation/   # 旧版实现(已废弃)
```

## 快速开始

### 1. 验证基础模块

```matlab
% 在MATLAB中切换到core目录
cd('matlab_project/core')

% 验证FFT/IFFT实现
run('test_fft.m')

% 验证快速卷积实现
run('test_fast_conv.m')
```

### 2. 运行主仿真

```matlab
% 运行完整的4QAM CMA盲均衡仿真
run('main_cma_simulation.m')
```

### 3. 参数影响分析

```matlab
% 运行参数影响分析实验
run('run_experiments.m')
```

## 核心模块说明

### my_fft.m
- **算法**: 基2-FFT,按时间抽取(DIT),迭代实现
- **特点**: 
  - 比特反转重排
  - 旋转因子预计算并利用对称性
  - 原位蝶形运算
- **输入**: 列向量,长度必须是2的幂
- **输出**: DFT结果

### my_ifft.m
- **算法**: 利用FFT实现IFFT: IDFT(X) = (1/N) * conj(DFT(conj(X)))
- **输入**: 频域信号(列向量,长度为2的幂)
- **输出**: 时域信号

### fast_conv_os.m
- **算法**: 重叠保留法(Overlap-Save)
- **特点**:
  - 分块处理长信号
  - 频域快速卷积
  - 去除混叠部分
- **输入**: 
  - x_long: 长输入信号
  - h: 滤波器冲激响应
  - nfft: FFT点数(2的幂)
- **输出**: 线性卷积结果

### main_cma_simulation.m
- **功能**: 完整的4QAM通信系统仿真
- **流程**:
  1. 发射端: 比特生成 → 格雷编码 → 4QAM调制
  2. 信道: ISI信道 → AWGN噪声
  3. 接收端: CMA盲均衡(混合架构) → 相位校正 → 解调
  4. 性能评估: BER、EVM、星座图
- **核心算法**: CMA混合处理架构
  - 逐样本权重更新
  - 块滤波(使用快速卷积)

### run_experiments.m
- **功能**: 参数影响分析
- **实验**:
  1. 步长μ的影响
  2. 滤波器长度M的影响
  3. SNR性能曲线(有/无均衡器对比)

## 技术规格书符合性

本实现严格遵循 `CMA均衡器MATLAB代码生成技术规格书.md` 的要求:

- ✅ **第2章**: FFT/IFFT从零实现,迭代蝶形运算,旋转因子优化
- ✅ **第3章**: 快速卷积(重叠保留法)
- ✅ **第4章**: 完整的4QAM通信链路仿真
  - 中心抽头初始化
  - 混合处理架构(逐样本更新+块滤波)
  - 相位校正
  - 性能评估
- ✅ **第5章**: 参数影响分析实验(μ, M, SNR)

## 关键参数说明

### CMA均衡器参数
- `filterLength`: 均衡器长度,建议为信道长度的2-3倍
- `stepSize (μ)`: 步长因子,影响收敛速度和稳态误差
  - 大步长: 收敛快,但稳态误差大
  - 小步长: 收敛慢,但稳态误差小
  - 推荐范围: 1e-4 ~ 1e-3
- `R2`: 恒模半径,对于归一化4QAM,R² = 1

### 仿真参数
- `numSymbols`: 符号数量,建议≥10000以观察收敛
- `snr_dB`: 信噪比(dB)
- `channelTaps`: ISI信道冲激响应

## 性能指标

- **BER (Bit Error Rate)**: 误码率
- **EVM (Error Vector Magnitude)**: 误差矢量幅度,反映星座点聚集度

## 注意事项

1. 所有信号向量应为**列向量**
2. FFT长度必须是**2的整数次幂**
3. 均衡器长度应**≥信道长度**,建议为信道长度的2-3倍
4. CMA存在**相位模糊性**,需要后续相位校正
5. 瞬态样本应被丢弃后再进行性能评估

## 作者

根据技术规格书要求实现

## 参考文献

见 `CMA均衡器MATLAB代码生成技术规格书.md`
