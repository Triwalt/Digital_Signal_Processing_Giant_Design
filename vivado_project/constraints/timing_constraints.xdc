# Timing constraints for DSP modules
# Create clock constraint
create_clock -period 10.000 -name sys_clk [get_ports clk]

# Set input delay constraints
set_input_delay -clock sys_clk -max 2.000 [get_ports {data_in_real data_in_imag data_valid start enable}]
set_input_delay -clock sys_clk -min 0.500 [get_ports {data_in_real data_in_imag data_valid start enable}]

# Set output delay constraints  
set_output_delay -clock sys_clk -max 2.000 [get_ports {data_out_real data_out_imag data_out_valid fft_done error_magnitude}]
set_output_delay -clock sys_clk -min 0.500 [get_ports {data_out_real data_out_imag data_out_valid fft_done error_magnitude}]

# Set reset constraints
set_false_path -from [get_ports rst_n]

# Set clock groups (if multiple clocks exist)
# set_clock_groups -asynchronous -group [get_clocks sys_clk]