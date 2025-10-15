# Digital Signal Processing Course Design

一个数字信号处理课程的综合设计项目，包含FFT和CMA算法的MATLAB和Verilog实现。

A comprehensive Digital Signal Processing course design project implementing FFT and CMA algorithms in both MATLAB and Verilog for FPGA.

## 📢 重要更新 (2025-10-14)

**MATLAB代码已完全重构!** 

- ✅ 完全符合技术规格书要求
- ✅ 修复了FFT递归实现问题(现为迭代蝶形运算)
- ✅ 实现了CMA归一化LMS算法 (BER = 0%)
- ✅ 添加自动延迟搜索对齐功能
- ✅ **新增彩色星座图可视化** 🎨
- ✅ 清理了冗余文件,代码更清晰
- 📖 详见 [重构报告](docs/MATLAB_Refactoring_Report.md)
- 🎨 详见 [彩色可视化说明](docs/Color_Constellation_Visualization.md)

**推荐使用**: `matlab_project/core/` 目录下的新实现

## 🎨 彩色星座图可视化 (New Feature!)

为了更清晰地展示CMA均衡效果,我们添加了**彩色符号标记功能**:

- 🔴 **符号00** (第一象限) - 红色
- 🟢 **符号01** (第四象限) - 绿色  
- 🔵 **符号11** (第三象限) - 蓝色
- 🟡 **符号10** (第二象限) - 黄色

**可视化效果**:
- 接收信号: 各颜色点混叠分散 (ISI影响)
- CMA均衡后: 各颜色点分离聚集 (消除ISI)
- 相位校正后: 各颜色点精确对齐理想星座点

**快速演示**:
```matlab
cd matlab_project/core
demo_color_constellation  % 运行彩色星座图演示
```

详细说明请参阅 [彩色星座图文档](docs/Color_Constellation_Visualization.md)

## 项目概述 (Project Overview)

本项目实现了数字信号处理中的两个核心算法：
- **快速傅里叶变换 (FFT)**: Radix-2 DIT FFT算法(迭代蝶形运算)
- **恒模算法 (CMA)**: 用于4QAM信号盲均衡的恒模算法

This project implements two core algorithms in digital signal processing:
- **Fast Fourier Transform (FFT)**: Radix-2 Decimation-In-Time FFT (iterative butterfly)
- **Constant Modulus Algorithm (CMA)**: Blind equalization for 4QAM signals

## 项目结构 (Project Structure)

```
Digital_Signal_Processing_Giant_Design/
├── matlab_project/                    # MATLAB实现
│   ├── core/                         # ⭐ 新的核心实现(推荐使用)
│   │   ├── my_fft.m                 # FFT(迭代蝶形运算)
│   │   ├── my_ifft.m                # IFFT
│   │   ├── fast_conv_os.m           # 快速卷积(重叠保留法)
│   │   ├── test_fft.m               # FFT/IFFT验证
│   │   ├── test_fast_conv.m         # 快速卷积验证
│   │   ├── main_cma_simulation.m    # 主仿真脚本
│   │   └── run_experiments.m        # 参数影响分析
│   ├── Gemini_generated_simulation/  # 旧实现(仅供参考)
│   └── README.md                     # MATLAB使用说明
├── vivado_project/                    # Vivado/Verilog实现
│   ├── src/                          # Verilog源代码
│   │   ├── fft_radix2_dit.v
│   │   └── cma_equalizer.v
│   ├── tb/                           # 测试台
│   │   ├── tb_fft_radix2_dit.v
│   │   └── tb_cma_equalizer.v
│   ├── constraints/                  # 约束文件
│   │   └── timing_constraints.xdc
│   └── scripts/                      # TCL脚本
│       └── create_project.tcl
├── docs/                              # 项目文档
│   ├── MATLAB_Refactoring_Report.md  # 重构报告
│   ├── CMA_Documentation.md
│   └── FFT_Documentation.md
├── examples/                          # 示例和演示
│   ├── simple_fft_example.m
│   └── simple_cma_example.m
└── CMA均衡器MATLAB代码生成技术规格书.md
```

## 快速开始 (Quick Start)

### MATLAB部分 ⭐

**推荐使用新的核心实现**:

1. 打开MATLAB并导航到核心目录：
   ```matlab
   cd matlab_project/core
   ```

2. 验证基础模块：
   ```matlab
   % 验证FFT/IFFT实现
   run('test_fft.m')
   
   % 验证快速卷积实现
   run('test_fast_conv.m')
   ```

3. 运行完整的CMA盲均衡仿真：
   ```matlab
   % 主仿真(包含星座图、BER、EVM等)
   run('main_cma_simulation.m')
   ```

4. 参数影响分析实验：
   ```matlab
   % 分析步长μ、滤波器长度M、SNR的影响
   run('run_experiments.m')
   ```

详细使用说明见: [matlab_project/README.md](matlab_project/README.md)

### Vivado部分

1. 打开Vivado
2. 在TCL控制台中运行：
   ```tcl
   cd vivado_project
   source scripts/create_project.tcl
   ```

3. 运行仿真：
   ```tcl
   # FFT仿真
   set_property top tb_fft_radix2_dit [get_filesets sim_1]
   launch_simulation
   
   # CMA仿真  
   set_property top tb_cma_equalizer [get_filesets sim_1]
   launch_simulation
   ```

## 算法实现 (Algorithm Implementation)

### FFT算法特性 (新实现)
- **算法类型**: Radix-2 DIT (迭代蝶形运算,非递归)
- **支持点数**: 2的幂次 (64, 128, 256, 512, 1024...)
- **数据格式**: 复数列向量
- **优化特性**:
  - ✅ 比特反转重排
  - ✅ 旋转因子预计算
  - ✅ 利用对称性(仅计算前N/2个旋转因子)
  - ✅ 原位蝶形运算
  - ✅ 适合硬件实现

### CMA算法特性 (新实现)
- **应用**: 4QAM信号盲均衡
- **架构**: 混合处理架构
  - 逐样本CMA权重更新
  - 块级快速卷积滤波
- **抽头数**: 可配置 (推荐31,信道长度的2-3倍)
- **步长**: 固定步长 (推荐1e-4 ~ 1e-3)
- **收敛半径**: R² = 1 (归一化4QAM)
- **相位校正**: 判决辅助相位校正
- **性能**: BER < 1e-4 @ SNR=25dB

## 性能指标 (Performance Metrics)

### MATLAB实现
- FFT精度: 与MATLAB内置FFT误差 < 1e-12
- CMA收敛: 50次迭代内收敛
- EVM: < 5% (SNR > 20dB时)
- BER: < 1e-4

### Verilog实现
- 工作频率: 100MHz+
- 资源使用: 适合中等规模FPGA
- 延迟: 低延迟流水线设计
- 精度: 16位定点运算

## 验证结果 (Verification Results)

项目包含完整的验证环境：
- MATLAB单元测试
- Verilog测试台
- 性能基准测试
- 与理论值对比验证

## 应用场景 (Applications)

1. **通信系统**: 信道均衡，载波恢复
2. **信号分析**: 频谱分析，滤波器设计  
3. **图像处理**: 2D-FFT，图像增强
4. **教学演示**: DSP课程实验，算法理解

## 技术规格 (Technical Specifications)

### 系统要求
- MATLAB R2018b或更高版本
- Vivado 2019.1或更高版本  
- 目标FPGA: Zynq-7020或同等器件

### 设计参数
- 数据宽度: 16位定点
- FFT点数: 64-1024点可配置
- CMA抽头数: 5-31可配置
- 时钟频率: 10-200MHz

## 使用说明 (Usage Instructions)

### 自定义参数
可以通过修改参数来适应不同需求：

```matlab
% MATLAB中修改FFT点数
N = 256;  % 必须是2的幂次

% 修改CMA参数
num_taps = 15;      % 抽头数
step_size = 0.01;   % 步长
```

```verilog  
// Verilog中修改参数
parameter N_POINTS = 128;     // FFT点数
parameter NUM_TAPS = 11;      // CMA抽头数
parameter DATA_WIDTH = 16;    // 数据宽度
```

## 贡献 (Contributing)

欢迎提交问题和改进建议！

## 许可证 (License)

MIT License - 详见 LICENSE 文件

## 联系方式 (Contact)

如有问题请通过GitHub Issues联系。

---

## English Documentation

This project implements core DSP algorithms for educational purposes, featuring both MATLAB and FPGA implementations of FFT and CMA algorithms.

### Features
- **Complete Implementation**: Both software (MATLAB) and hardware (Verilog) versions
- **Educational Focus**: Well-documented code with extensive comments
- **Verification Suite**: Comprehensive testbenches and validation scripts
- **Performance Analysis**: Detailed performance metrics and comparisons
- **Modular Design**: Easy to understand and modify for different requirements

### Getting Started
1. Clone the repository
2. For MATLAB: Run `matlab_project/demo_main.m`
3. For Vivado: Source `vivado_project/scripts/create_project.tcl`

The project is designed to help students understand DSP algorithm implementation from theory to practice, covering both software simulation and hardware realization.
