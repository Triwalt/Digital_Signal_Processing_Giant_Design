# DSP 项目的 Claude Code 配置

## 项目概述

这是一个数字信号处理（DSP）课程设计项目，在 MATLAB 与 Verilog（用于 FPGA 部署）中实现 FFT（快速傅里叶变换）与 CMA（恒模算法）。

## Claude Code 配置项

### 自动批准设置

以下工具被配置为自动批准以最大化效率：

- **文件操作**：Read、Write、Edit、MultiEdit
- **搜索操作**：Glob、Grep
- **系统操作**：Bash 命令
- **任务管理**：Task、TodoWrite
- **网络操作**：WebFetch、WebSearch
- **IDE 集成**：MCP commands、NotebookEdit
- **规划**：ExitPlanMode

### 效率优化

- `confirmCommands: false` - 执行操作无需确认提示
- `autoRunSuggestibleTasks: true` - 自动运行可建议的任务
- `autoRunBashCommands: true` - 自动运行 bash 命令
- `autoRunPythonCode: true` - 在笔记本中自动运行 Python 代码
- `mode: auto` - 自动模式选择
- `verbose: false` - 降低输出冗长度

## DSP 项目特定指南

### 算法实现规范

#### FFT 实现指南

1. **基2 DIT 结构**：遵循既定的蝶形计算模式
2. **位反转优化**：使用现有的 `bit_reverse` 函数进行输入重排
3. **复数处理**：MATLAB 使用 `1j` 作为虚数单位，Verilog 使用分离的实部/虚部
4. **性能目标**：
   - MATLAB：相对内置 FFT 的误差 < 1e-12
   - Verilog：100MHz+ 工作频率，16 位定点精度

#### CMA 算法指南

1. **盲均衡**：聚焦于 QPSK 信号恢复
2. **LMS 风格更新**：遵循随机梯度下降方法
3. **中心抽头初始化**：中心抽头初始化为 1，其余为 0
4. **收敛标准**：目标在 50 次迭代内收敛

### 代码结构模式

#### MATLAB 代码风格

```matlab
% Function naming: lowercase_with_underscores
function [output1, output2] = function_name(input1, input2)
    % Extensive comments explaining algorithm theory
    % Input validation
    % Core implementation
    % Output validation
end
```

#### Verilog 代码风格

```verilog
// Module naming: lowercase_with_underscores
module module_name #(
    parameter PARAM1 = default_value,
    parameter PARAM2 = default_value
)(
    input wire clk,
    input wire rst_n,
    // Additional ports
);
    // Local parameters
    // Signal declarations
    // Implementation logic
endmodule
```

### 测试与验证规范

#### MATLAB 测试

- **功能测试**：正弦信号、随机噪声
- **性能测试**：不同规模下的计时基准
- **精度测试**：与 MATLAB 内置函数对比
- **边界情况**：最小规模、非法输入

#### Verilog 测试

- **测试平台结构**：为每个模块单独编写 testbench
- **信号生成**：包含噪声的真实输入模式
- **性能验证**：时序与资源利用
- **FPGA 验证**：满足目标器件的约束

### 常见开发任务

#### 添加新的 DSP 算法

1. 在相应目录创建实现文件：
   - MATLAB：`matlab_project/src/`
   - Verilog：`vivado_project/src/`
2. 遵循既定命名规范
3. 添加详尽注释解释算法原理
4. 创建相应测试文件
5. 更新主演示脚本

#### 运行测试

```bash
# MATLAB tests
cd matlab_project/test
matlab -batch "test_fft"
matlab -batch "test_cma"

# Verilog simulation (in Vivado TCL console)
cd vivado_project
source scripts/create_project.tcl
```

#### 性能优化

1. **MATLAB**：向量化、预分配、选择高效算法
2. **Verilog**：流水线优化、资源共享、时序收敛

### 项目特定命令

#### MATLAB 开发

```matlab
% Run main demonstration
matlab_project/demo_main

% Individual algorithm tests
cd matlab_project/test
test_fft
test_cma
```

#### Verilog 开发

```tcl
# In Vivado TCL console
cd vivado_project
source scripts/create_project.tcl

# Run simulations
set_property top tb_fft_radix2_dit [get_filesets sim_1]
launch_simulation

set_property top tb_cma_equalizer [get_filesets sim_1]
launch_simulation
```

### 错误处理与调试

#### 常见问题

1. **FFT 尺寸校验**：大小必须为 2 的幂
2. **CMA 收敛**：检查步长与抽头数
3. **Verilog 时序**：满足建立/保持时间要求
4. **MATLAB 精度**：浮点与定点差异

#### 调试策略

1. **增量式测试**：逐个组件进行测试
2. **可视化**：使用绘图进行信号分析
3. **对比**：与参考实现进行比较
4. **日志**：添加诊断输出以便排障

### 文档规范

#### 代码注释

- 在函数/模块层面解释算法原理
- 对实现步骤进行逐步注释
- 参数说明与使用示例
- 性能特征与限制

#### 文件头模板

```matlab
% Function: algorithm_name
% Description: Brief description of the algorithm
% Author: [Name]
% Date: [Date]
% Parameters:
%   input1 - Description of input1
% Returns:
%   output1 - Description of output1
% Algorithm: Reference to theoretical basis
```

```verilog
// Module: module_name
// Description: Brief description of the module
// Author: [Name]
// Date: [Date]
// Parameters:
//   PARAM1 - Description of parameter1
// Ports:
//   clk - Clock input
//   rst_n - Active-low reset
```

### 性能指标

#### MATLAB 性能

- **FFT 精度**：相对内置 FFT 的误差 < 1e-12
- **CMA 收敛**：< 50 次迭代
- **EVM**：< 5%（SNR > 20dB）
- **BER**：< 1e-4

#### Verilog 性能

- **工作频率**：100MHz+
- **资源占用**：适配 Zynq-7020
- **时延**：优化的流水线设计
- **精度**：16 位定点

### 集成指南

#### 跨平台一致性

1. **算法等价性**：MATLAB 与 Verilog 实现应产生等价结果
2. **参数映射**：参数名称与取值保持一致
3. **测试向量一致**：两个平台使用相同测试模式
4. **性能相关性**：结果应在统计意义上等价

#### 版本控制

1. **提交规范**：清晰、具描述性的提交信息
2. **分支策略**：新算法使用特性分支
3. **提交前测试**：运行完整测试套件
4. **文档更新**：保持 README 与 docs 同步更新

此配置可在确保准确性、性能与教学价值的同时，帮助利用 Claude Code 高效开展 DSP 算法开发。