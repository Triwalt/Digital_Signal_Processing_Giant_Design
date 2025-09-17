# Vivado TCL Script for DSP Project Setup
# This script creates and configures the Vivado project for the DSP course design

# Set project variables
set project_name "dsp_giant_design"
set project_dir "./vivado_project_build"
set part_name "xc7z020clg484-1"  # Zynq-7020 FPGA (common for educational boards)

# Create project
create_project $project_name $project_dir -part $part_name -force

# Add source files
add_files -norecurse {
    ./src/fft_radix2_dit.v
    ./src/cma_equalizer.v
}

# Add testbench files
add_files -fileset sim_1 -norecurse {
    ./tb/tb_fft_radix2_dit.v
    ./tb/tb_cma_equalizer.v
}

# Add constraints (if any)
# add_files -fileset constrs_1 -norecurse ./constraints/timing_constraints.xdc

# Set top module for synthesis
set_property top fft_radix2_dit [current_fileset]

# Set top module for simulation
set_property top tb_fft_radix2_dit [get_filesets sim_1]

# Project settings
set_property target_language Verilog [current_project]
set_property simulator_language Verilog [current_project]

# Set synthesis strategy
set_property strategy {Vivado Synthesis Defaults} [get_runs synth_1]

# Set implementation strategy  
set_property strategy {Vivado Implementation Defaults} [get_runs impl_1]

puts "Project created successfully!"
puts "Source files added: FFT and CMA modules"
puts "Testbenches added for both modules"
puts "Ready for simulation and synthesis"

# Optional: Run synthesis
# reset_run synth_1
# launch_runs synth_1
# wait_on_run synth_1

# Optional: Run simulation
# reset_run sim_1
# launch_simulation

puts "Project setup complete. Use the following commands:"
puts "  - To run FFT simulation: set_property top tb_fft_radix2_dit [get_filesets sim_1]; launch_simulation"
puts "  - To run CMA simulation: set_property top tb_cma_equalizer [get_filesets sim_1]; launch_simulation"
puts "  - To run synthesis: launch_runs synth_1"