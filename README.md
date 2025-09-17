# Digital Signal Processing Course Design

一个数字信号处理课程的综合设计项目，包含FFT和CMA算法的MATLAB和Verilog实现。

A comprehensive Digital Signal Processing course design project implementing FFT and CMA algorithms in both MATLAB and Verilog for FPGA.

## 项目概述 (Project Overview)

本项目实现了数字信号处理中的两个核心算法：
- **快速傅里叶变换 (FFT)**: Radix-2 抽取时间FFT算法
- **常模算法 (CMA)**: 用于盲均衡的常模算法

This project implements two core algorithms in digital signal processing:
- **Fast Fourier Transform (FFT)**: Radix-2 Decimation-In-Time FFT algorithm  
- **Constant Modulus Algorithm (CMA)**: Blind equalization algorithm

## 项目结构 (Project Structure)

```
Digital_Signal_Processing_Giant_Design/
├── matlab_project/           # MATLAB实现
│   ├── src/                 # 源代码
│   │   ├── fft_implementation.m
│   │   └── cma_algorithm.m
│   ├── test/                # 测试脚本
│   │   ├── test_fft.m
│   │   └── test_cma.m
│   ├── docs/                # 文档
│   └── demo_main.m          # 主演示脚本
├── vivado_project/          # Vivado/Verilog实现
│   ├── src/                 # Verilog源代码
│   │   ├── fft_radix2_dit.v
│   │   └── cma_equalizer.v
│   ├── tb/                  # 测试台
│   │   ├── tb_fft_radix2_dit.v
│   │   └── tb_cma_equalizer.v
│   ├── constraints/         # 约束文件
│   │   └── timing_constraints.xdc
│   └── scripts/             # TCL脚本
│       └── create_project.tcl
├── docs/                    # 项目文档
└── examples/                # 示例和演示
```

## 快速开始 (Quick Start)

### MATLAB部分

1. 打开MATLAB并导航到项目目录
2. 运行主演示脚本：
   ```matlab
   matlab_project/demo_main.m
   ```

3. 运行单独的测试：
   ```matlab
   % FFT测试
   cd matlab_project/test
   test_fft
   
   % CMA测试
   test_cma
   ```

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

### FFT算法特性
- **算法类型**: Radix-2 抽取时间 (DIT)
- **支持点数**: 2的幂次 (64, 128, 256, 512, 1024...)
- **数据格式**: 复数 (实部+虚部)
- **优化**: 位反序优化，蝶形运算优化

### CMA算法特性  
- **应用**: QPSK信号盲均衡
- **抽头数**: 可配置 (默认11)
- **步长**: 自适应步长控制
- **收敛**: 快速收敛，低稳态误差

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
