`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: DSP Course Design
// Engineer: Student
// 
// Create Date: 2025
// Design Name: CMA Equalizer
// Module Name: cma_equalizer
// Project Name: Digital Signal Processing Giant Design
// Target Devices: FPGA
// Tool Versions: Vivado
// Description: 
//   Constant Modulus Algorithm (CMA) equalizer implementation in Verilog
//   Performs blind equalization for communication signals
//
//////////////////////////////////////////////////////////////////////////////////

module cma_equalizer #(
    parameter DATA_WIDTH = 16,       // Data width for real and imaginary parts
    parameter NUM_TAPS = 11,         // Number of equalizer taps
    parameter FRAC_BITS = 12,        // Number of fractional bits for fixed-point
    parameter STEP_SIZE = 1024       // Step size (mu) in fixed-point format
)(
    input clk,
    input rst_n,
    input enable,
    input signed [DATA_WIDTH-1:0] data_in_real,
    input signed [DATA_WIDTH-1:0] data_in_imag,
    input data_valid,
    
    output reg signed [DATA_WIDTH-1:0] data_out_real,
    output reg signed [DATA_WIDTH-1:0] data_out_imag,
    output reg data_out_valid,
    output reg [31:0] error_magnitude
);

    // Internal registers and wires
    reg signed [DATA_WIDTH-1:0] input_buffer_real [0:NUM_TAPS-1];
    reg signed [DATA_WIDTH-1:0] input_buffer_imag [0:NUM_TAPS-1];
    reg signed [DATA_WIDTH-1:0] weights_real [0:NUM_TAPS-1];
    reg signed [DATA_WIDTH-1:0] weights_imag [0:NUM_TAPS-1];
    
    wire signed [2*DATA_WIDTH-1:0] mult_result_real [0:NUM_TAPS-1];
    wire signed [2*DATA_WIDTH-1:0] mult_result_imag [0:NUM_TAPS-1];
    
    reg signed [DATA_WIDTH+$clog2(NUM_TAPS)-1:0] output_real_sum;
    reg signed [DATA_WIDTH+$clog2(NUM_TAPS)-1:0] output_imag_sum;
    
    reg signed [DATA_WIDTH-1:0] equalizer_out_real;
    reg signed [DATA_WIDTH-1:0] equalizer_out_imag;
    
    reg signed [2*DATA_WIDTH-1:0] magnitude_squared;
    reg signed [DATA_WIDTH-1:0] error_real;
    reg signed [DATA_WIDTH-1:0] error_imag;
    
    // Constants
    localparam R2 = (1 << FRAC_BITS);  // R^2 = 1 in fixed-point
    
    integer i;
    
    // Complex multiplication: (a+jb) * (c+jd) = (ac-bd) + j(ad+bc)
    genvar tap_idx;
    generate
        for (tap_idx = 0; tap_idx < NUM_TAPS; tap_idx = tap_idx + 1) begin : gen_multipliers
            assign mult_result_real[tap_idx] = 
                input_buffer_real[tap_idx] * weights_real[tap_idx] - 
                input_buffer_imag[tap_idx] * weights_imag[tap_idx];
            assign mult_result_imag[tap_idx] = 
                input_buffer_real[tap_idx] * weights_imag[tap_idx] + 
                input_buffer_imag[tap_idx] * weights_real[tap_idx];
        end
    endgenerate
    
    // Initialize weights (center tap to 1, others to 0)
    initial begin
        for (i = 0; i < NUM_TAPS; i = i + 1) begin
            if (i == NUM_TAPS/2) begin
                weights_real[i] = (1 << FRAC_BITS);  // 1.0 in fixed-point
                weights_imag[i] = 0;
            end else begin
                weights_real[i] = 0;
                weights_imag[i] = 0;
            end
        end
    end
    
    // Main processing pipeline
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            // Reset all registers
            for (i = 0; i < NUM_TAPS; i = i + 1) begin
                input_buffer_real[i] <= 0;
                input_buffer_imag[i] <= 0;
            end
            output_real_sum <= 0;
            output_imag_sum <= 0;
            equalizer_out_real <= 0;
            equalizer_out_imag <= 0;
            data_out_real <= 0;
            data_out_imag <= 0;
            data_out_valid <= 0;
            magnitude_squared <= 0;
            error_real <= 0;
            error_imag <= 0;
            error_magnitude <= 0;
        end else if (enable && data_valid) begin
            
            // Shift input buffer
            for (i = NUM_TAPS-1; i > 0; i = i - 1) begin
                input_buffer_real[i] <= input_buffer_real[i-1];
                input_buffer_imag[i] <= input_buffer_imag[i-1];
            end
            input_buffer_real[0] <= data_in_real;
            input_buffer_imag[0] <= data_in_imag;
            
            // Compute equalizer output (sum of tap outputs)
            output_real_sum = 0;
            output_imag_sum = 0;
            for (i = 0; i < NUM_TAPS; i = i + 1) begin
                output_real_sum = output_real_sum + (mult_result_real[i] >>> FRAC_BITS);
                output_imag_sum = output_imag_sum + (mult_result_imag[i] >>> FRAC_BITS);
            end
            
            // Saturate and assign outputs
            if (output_real_sum > ((1 << (DATA_WIDTH-1)) - 1))
                equalizer_out_real <= (1 << (DATA_WIDTH-1)) - 1;
            else if (output_real_sum < -(1 << (DATA_WIDTH-1)))
                equalizer_out_real <= -(1 << (DATA_WIDTH-1));
            else
                equalizer_out_real <= output_real_sum[DATA_WIDTH-1:0];
                
            if (output_imag_sum > ((1 << (DATA_WIDTH-1)) - 1))
                equalizer_out_imag <= (1 << (DATA_WIDTH-1)) - 1;
            else if (output_imag_sum < -(1 << (DATA_WIDTH-1)))
                equalizer_out_imag <= -(1 << (DATA_WIDTH-1));
            else
                equalizer_out_imag <= output_imag_sum[DATA_WIDTH-1:0];
            
            // Compute magnitude squared: |y|^2
            magnitude_squared <= (equalizer_out_real * equalizer_out_real) + 
                               (equalizer_out_imag * equalizer_out_imag);
            
            // Compute CMA error: y*(R^2 - |y|^2)
            error_real <= (equalizer_out_real * (R2 - (magnitude_squared >>> FRAC_BITS))) >>> FRAC_BITS;
            error_imag <= (equalizer_out_imag * (R2 - (magnitude_squared >>> FRAC_BITS))) >>> FRAC_BITS;
            
            // Update weights: w = w + mu * conj(error) * x
            for (i = 0; i < NUM_TAPS; i = i + 1) begin
                weights_real[i] <= weights_real[i] + 
                    ((STEP_SIZE * (error_real * input_buffer_real[i] + error_imag * input_buffer_imag[i])) >>> (2*FRAC_BITS));
                weights_imag[i] <= weights_imag[i] + 
                    ((STEP_SIZE * (error_imag * input_buffer_real[i] - error_real * input_buffer_imag[i])) >>> (2*FRAC_BITS));
            end
            
            // Output
            data_out_real <= equalizer_out_real;
            data_out_imag <= equalizer_out_imag;
            data_out_valid <= 1;
            
            // Error magnitude for monitoring
            error_magnitude <= (error_real * error_real) + (error_imag * error_imag);
            
        end else begin
            data_out_valid <= 0;
        end
    end

endmodule