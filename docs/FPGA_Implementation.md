# FPGA实现指南

本文档描述FFT和CMA均衡器的FPGA/Verilog实现策略与验证计划。

---

## 目录

1. [设计概述](#1-设计概述)
2. [FFT硬件架构](#2-fft硬件架构)
3. [CMA均衡器架构](#3-cma均衡器架构)
4. [验证计划](#4-验证计划)
5. [资源评估](#5-资源评估)

---

## 1. 设计概述

### 1.1 目标与约束

- **FFT点数**: 32点与128点，参数化设计
- **实现要求**: 禁止使用厂商FFT IP核，基于Verilog自主实现
- **目标器件**: xc7z010clg400-1 (Zynq-7010)
- **性能目标**: 100 MHz时钟下持续流式输出

### 1.2 项目文件结构

```
vivado_project/
├── src/                    # RTL源文件
│   ├── fft_radix2_dit.v   # FFT模块
│   └── cma_equalizer.v    # CMA均衡器
├── tb/                     # 测试平台
│   ├── tb_fft_radix2_dit.v
│   └── tb_cma_equalizer.v
├── constraints/            # 约束文件
├── scripts/                # TCL脚本
│   ├── create_project.tcl
│   ├── run_fft_sim.tcl
│   └── report_fft_util.tcl
└── reports/                # 报告输出
```

---

## 2. FFT硬件架构

### 2.1 架构选型

采用 **Radix-2² SDF**（单路径延迟反馈）流水线结构：

- 输入输出均为一拍一个复数样点，保证持续吞吐
- 各stage结构相同，便于参数化复用
- 相比传统Radix-2 DIT，更易做寄存器平衡和DSP映射

### 2.2 数据格式

| 位置 | 格式 | 说明 |
|------|------|------|
| 外部接口 | Q2.13 (16位) | 1符号 + 2整数 + 13小数 |
| 旋转因子 | Q1.14 (16位) | 高精度系数 |
| 乘法结果 | 32位 | 16×16，截断前暂存 |
| 内部累加 | 扩展位宽 | 防止溢出 |

### 2.3 模块划分

```
fft_top
├── butterfly_stage[0..log2(N)-1]  # 蝶形运算级
│   ├── butterfly_unit             # 蝶形运算器
│   ├── delay_line                 # 延迟线
│   ├── twiddle_rom                # 旋转因子ROM
│   └── scaling_ctrl               # 缩放控制
├── bit_reversal                   # 比特反转
└── io_adapter                     # 接口适配
```

### 2.4 旋转因子管理

**生成方式**:
- 使用MATLAB/Python脚本预生成
- 输出文件: `twiddle_real_XXX.mem`, `twiddle_imag_XXX.mem`
- 格式: 十六进制补码，便于`$readmemh`读取

**地址编码**:
```
addr = {stage_id, butterfly_idx[stage_id:0]}
```

### 2.5 流水线策略

- 每级蝶形后加寄存器
- DSP48 + LUT路径在一个时钟周期内完成
- Stage数量: 32点→5级，128点→7级
- 延迟: ≈ stage数 + 延迟线总长度

---

## 3. CMA均衡器架构

### 3.1 模块参数

```verilog
module cma_equalizer #(
    parameter DATA_WIDTH = 16,      // 数据位宽
    parameter NUM_TAPS = 11,        // 均衡器抽头数
    parameter FRAC_BITS = 12,       // 小数位数
    parameter STEP_SIZE = 16'h0010  // 步长(定点)
)(
    input clk, rst_n, enable,
    input signed [DATA_WIDTH-1:0] data_in_real, data_in_imag,
    input data_valid,
    output signed [DATA_WIDTH-1:0] data_out_real, data_out_imag,
    output data_out_valid,
    output [31:0] error_magnitude
);
```

### 3.2 架构组成

1. **输入缓冲**: 移位寄存器实现抽头延迟线
2. **权重存储**: 寄存器组存储自适应滤波器权重
3. **复数乘法器**: 计算滤波器输出
4. **误差计算**: CMA误差 e = y·(R² - |y|²)
5. **权重更新**: 梯度下降更新 w = w + μ·e·x*

### 3.3 定点考虑

- 输入/输出: Q(DATA_WIDTH-FRAC_BITS).FRAC_BITS
- 权重: 与数据相同格式
- 内部计算: 扩展精度防止溢出
- 饱和算术: 处理溢出情况

---

## 4. 验证计划

### 4.1 MATLAB参考数据生成

**脚本**: `matlab_project/scripts/generate_fft_vectors.m`

**生成内容**:
- 多组激励: 单频正弦、双频组合、白噪声、冲激、随机QAM
- 定点量化: Q2.13格式
- 黄金结果: MATLAB FFT输出

**运行方式**:
```matlab
cd('matlab_project/scripts');
summary = generate_fft_vectors('PointSizes', [32 128], ...
                               'WordLength', 16, ...
                               'FracLength', 13);
```

**输出文件**:
```
matlab_project/data/fft_vectors/<N>/
├── input_real_<N>.mem
├── input_imag_<N>.mem
├── output_real_<N>.mem
├── output_imag_<N>.mem
└── metadata.json
```

### 4.2 RTL仿真

**Testbench结构**:
```
tb_fft32_file.sv
├── mem_loader          # 读取.mem文件
├── fft_stimulus        # 生成激励
├── fft_scoreboard      # 比较结果
└── logger              # 记录日志
```

**验收标准**:
- 所有激励场景误差 ≤ 1 LSB
- 输出延迟与预期一致
- 仿真退出状态为PASS

**运行命令** (Windows):
```cmd
vivado -mode batch -source vivado_project\scripts\run_fft_sim.tcl -tclargs -n 32
```

### 4.3 综合与资源分析

**运行命令**:
```tcl
source scripts/report_fft_util.tcl -n 32
source scripts/report_fft_util.tcl -n 128
```

**输出报告**:
- `reports/fft32_utilization.rpt`
- `reports/fft128_utilization.rpt`

**验证点**:
- LUT、FF、DSP、BRAM数量
- Worst Negative Slack ≥ 0 (100 MHz)

### 4.4 性能评估

**数值统计**:
- 均方误差 (MSE)
- 信噪比 (SNR)
- 谱峰失真 (SFDR)

**长序列测试**:
- 随机激励1k帧
- 验证持续输出不丢拍
- 统计throughput和latency

---

## 5. 资源评估

### 5.1 FFT资源估算

| 资源 | 32点 | 128点 | 说明 |
|------|------|-------|------|
| LUT | ~1500 | ~3000 | 逻辑单元 |
| FF | ~1200 | ~2500 | 寄存器 |
| DSP48 | ~6 | ~14 | 乘法器 |
| BRAM | ~1 | ~2 | 旋转因子ROM |

### 5.2 CMA资源估算 (11抽头)

| 资源 | 数量 | 说明 |
|------|------|------|
| LUT | ~2000 | 逻辑单元 |
| FF | ~1500 | 寄存器 |
| DSP48 | ~8 | 复数乘法 |
| BRAM | ~0 | 权重用寄存器 |

### 5.3 时序目标

- 目标频率: 100 MHz
- 时钟周期: 10 ns
- 关键路径: 复数乘法 + 累加

---

## 6. 实施路径

### 6.1 里程碑

| 阶段 | 交付物 | 验证方式 |
|------|--------|----------|
| 设计策略 | 本文档 | 审阅 |
| MATLAB脚本 | 激励/黄金数据 | 自检 |
| 32点仿真 | fft32_sim.log | 误差PASS |
| 128点仿真 | fft128_sim.log | 误差PASS |
| 资源报告 | utilization.rpt | 综合PASS |
| 性能分析 | 指标表格 | 达标 |

### 6.2 上板准备（可选）

1. **顶层封装**: AXI-Stream/AXI-Lite适配层
2. **硬件验证**: PS端生成激励，DMA推入
3. **结果捕获**: 逻辑分析仪或BRAM捕获

---

## 附录：Vivado常用命令

### 创建项目
```tcl
create_project dsp_design ./build -part xc7z010clg400-1
add_files -norecurse {src/fft_radix2_dit.v src/cma_equalizer.v}
add_files -fileset sim_1 -norecurse {tb/tb_fft_radix2_dit.v}
```

### 运行仿真
```tcl
set_property top tb_fft_radix2_dit [get_filesets sim_1]
launch_simulation
run all
```

### 综合与实现
```tcl
launch_runs synth_1
wait_on_run synth_1
launch_runs impl_1 -to_step write_bitstream
wait_on_run impl_1
```

---

*最后更新: 2025年12月*
