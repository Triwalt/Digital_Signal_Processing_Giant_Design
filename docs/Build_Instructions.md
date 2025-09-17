# Build and Test Instructions

## MATLAB Environment Setup

### Prerequisites
- MATLAB R2018b or later
- Signal Processing Toolbox
- Communications System Toolbox (recommended)

### Running MATLAB Code

1. **Navigate to project directory:**
   ```matlab
   cd 'path/to/Digital_Signal_Processing_Giant_Design'
   ```

2. **Run main demonstration:**
   ```matlab
   matlab_project/demo_main
   ```

3. **Run individual tests:**
   ```matlab
   % FFT tests
   cd matlab_project/test
   test_fft
   
   % CMA tests  
   test_cma
   ```

4. **Run simple examples:**
   ```matlab
   % Simple FFT example
   examples/simple_fft_example
   
   % Simple CMA example
   examples/simple_cma_example
   ```

## Vivado Environment Setup

### Prerequisites
- Vivado 2019.1 or later
- Target FPGA: Zynq-7020 or equivalent

### Creating Vivado Project

1. **Open Vivado and navigate to project directory:**
   ```tcl
   cd vivado_project
   ```

2. **Source the project creation script:**
   ```tcl
   source scripts/create_project.tcl
   ```

3. **Alternatively, create project manually:**
   ```tcl
   create_project dsp_giant_design ./vivado_project_build -part xc7z020clg484-1
   add_files -norecurse {src/fft_radix2_dit.v src/cma_equalizer.v}
   add_files -fileset sim_1 -norecurse {tb/tb_fft_radix2_dit.v tb/tb_cma_equalizer.v}
   ```

### Running Simulations

1. **FFT Simulation:**
   ```tcl
   set_property top tb_fft_radix2_dit [get_filesets sim_1]
   launch_simulation
   run all
   ```

2. **CMA Simulation:**
   ```tcl
   set_property top tb_cma_equalizer [get_filesets sim_1]
   launch_simulation  
   run all
   ```

### Synthesis and Implementation

1. **Run Synthesis:**
   ```tcl
   reset_run synth_1
   launch_runs synth_1
   wait_on_run synth_1
   ```

2. **Run Implementation:**
   ```tcl
   launch_runs impl_1 -to_step write_bitstream
   wait_on_run impl_1
   ```

## Testing Procedures

### MATLAB Testing

1. **Verify FFT Implementation:**
   ```matlab
   % Test against MATLAB's built-in FFT
   x = randn(1, 64);
   X_custom = fft_implementation(x);
   X_matlab = fft(x);
   error = max(abs(X_custom - X_matlab));
   fprintf('Max error: %.2e\n', error);
   ```

2. **Verify CMA Implementation:**
   ```matlab
   % Test convergence with known channel
   symbols = (2*randi([0,1], 100, 1) - 1) + 1j*(2*randi([0,1], 100, 1) - 1);
   channel = [1, 0.5];
   received = conv(symbols, channel, 'same');
   [eq_out, weights, errors] = cma_algorithm(received, 0.01, 7, 20);
   final_error = errors(end);
   fprintf('Final CMA error: %.2e\n', final_error);
   ```

### Verilog Testing

1. **FFT Testbench:**
   - Generates sinusoidal input
   - Verifies frequency domain output
   - Checks for proper peak detection
   - Validates timing and control signals

2. **CMA Testbench:**
   - Generates QPSK symbols
   - Applies channel distortion
   - Monitors convergence
   - Analyzes output constellation

### Expected Results

#### MATLAB Results:
- **FFT Accuracy**: Error < 1e-12 compared to MATLAB FFT
- **CMA Convergence**: Error reduces by 2-3 orders of magnitude
- **Processing Speed**: Custom FFT ~10x slower than MATLAB (acceptable for educational purposes)

#### Verilog Results:
- **FFT Simulation**: Correct frequency peaks at expected bins
- **CMA Simulation**: Error magnitude decreases over time
- **Synthesis**: Should meet timing at 100MHz for moderate FPGA

## Troubleshooting

### Common MATLAB Issues:

1. **Path Problems:**
   ```matlab
   addpath('matlab_project/src');
   ```

2. **Memory Issues with Large FFTs:**
   - Reduce FFT size for testing
   - Clear workspace regularly

3. **Plotting Issues:**
   ```matlab
   close all; % Clear existing figures
   ```

### Common Vivado Issues:

1. **Simulation Not Starting:**
   - Check testbench is set as top module
   - Verify all source files are added

2. **Synthesis Errors:**
   - Check for missing signals in sensitivity lists
   - Verify parameter consistency

3. **Timing Violations:**
   - Reduce clock frequency
   - Add pipeline stages

### Performance Validation:

1. **FFT Validation:**
   ```matlab
   % Compare with known transform pairs
   x = [1, 0, 0, 0]; % Impulse
   X = fft_implementation(x);
   % Should give [1, 1, 1, 1]
   ```

2. **CMA Validation:**
   ```matlab
   % Test with no channel distortion
   symbols = [1+1j, -1+1j, -1-1j, 1-1j]; % QPSK
   [eq_out, ~, ~] = cma_algorithm(symbols, 0.01, 3, 10);
   % Output should match input after convergence
   ```

## Build Automation

### MATLAB Batch Testing:
```matlab
% Create test suite
results = struct();
results.fft_test = run_fft_tests();
results.cma_test = run_cma_tests();
save('test_results.mat', 'results');
```

### Vivado Batch Processing:
```tcl
# Automated build script
source scripts/create_project.tcl
launch_runs synth_1
wait_on_run synth_1
launch_runs impl_1 -to_step write_bitstream
wait_on_run impl_1
```

This completes the build and test instructions for the DSP course design project.