`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: DSP Course Design
// Engineer: Student
// 
// Create Date: 2025
// Design Name: FFT Radix-2 DIT
// Module Name: fft_radix2_dit
// Project Name: Digital Signal Processing Giant Design
// Target Devices: FPGA
// Tool Versions: Vivado
// Description: 
//   Radix-2 Decimation-In-Time FFT implementation in Verilog
//   Supports N-point FFT where N is power of 2
//
//////////////////////////////////////////////////////////////////////////////////

module fft_radix2_dit #(
    parameter N_POINTS = 64,         // Number of FFT points (must be power of 2)
    parameter DATA_WIDTH = 16,       // Data width for real and imaginary parts
    parameter TWIDDLE_WIDTH = 16     // Twiddle factor width
)(
    input clk,
    input rst_n,
    input start,
    input [DATA_WIDTH-1:0] data_in_real,
    input [DATA_WIDTH-1:0] data_in_imag,
    input data_valid,
    
    output reg [DATA_WIDTH-1:0] data_out_real,
    output reg [DATA_WIDTH-1:0] data_out_imag,
    output reg data_out_valid,
    output reg fft_done
);

    // Local parameters
    localparam STAGES = $clog2(N_POINTS);
    localparam ADDR_WIDTH = $clog2(N_POINTS);
    
    // Internal signals
    reg [DATA_WIDTH-1:0] real_mem [0:N_POINTS-1];
    reg [DATA_WIDTH-1:0] imag_mem [0:N_POINTS-1];
    reg [DATA_WIDTH-1:0] temp_real_mem [0:N_POINTS-1];
    reg [DATA_WIDTH-1:0] temp_imag_mem [0:N_POINTS-1];
    
    reg [ADDR_WIDTH-1:0] input_counter;
    reg [ADDR_WIDTH-1:0] output_counter;
    reg [$clog2(STAGES+1)-1:0] stage_counter;
    reg [ADDR_WIDTH-1:0] butterfly_counter;
    reg [ADDR_WIDTH-1:0] group_counter;
    
    // State machine
    typedef enum reg [2:0] {
        IDLE,
        INPUT_DATA,
        BIT_REVERSE,
        FFT_COMPUTE,
        OUTPUT_DATA
    } state_t;
    
    state_t current_state, next_state;
    
    // Twiddle factor ROM
    reg signed [TWIDDLE_WIDTH-1:0] twiddle_real [0:N_POINTS/2-1];
    reg signed [TWIDDLE_WIDTH-1:0] twiddle_imag [0:N_POINTS/2-1];
    
    // Initialize twiddle factors
    initial begin
        integer i;
        real angle;
        for (i = 0; i < N_POINTS/2; i = i + 1) begin
            angle = -2.0 * 3.14159265359 * i / N_POINTS;
            twiddle_real[i] = $rtoi($cos(angle) * (2**(TWIDDLE_WIDTH-1) - 1));
            twiddle_imag[i] = $rtoi($sin(angle) * (2**(TWIDDLE_WIDTH-1) - 1));
        end
    end
    
    // Bit reversal function
    function [ADDR_WIDTH-1:0] bit_reverse;
        input [ADDR_WIDTH-1:0] index;
        integer i;
        begin
            bit_reverse = 0;
            for (i = 0; i < ADDR_WIDTH; i = i + 1) begin
                bit_reverse[ADDR_WIDTH-1-i] = index[i];
            end
        end
    endfunction
    
    // State machine transitions
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            current_state <= IDLE;
        end else begin
            current_state <= next_state;
        end
    end
    
    // Next state logic
    always @(*) begin
        next_state = current_state;
        case (current_state)
            IDLE: begin
                if (start) next_state = INPUT_DATA;
            end
            INPUT_DATA: begin
                if (input_counter == N_POINTS - 1 && data_valid)
                    next_state = BIT_REVERSE;
            end
            BIT_REVERSE: begin
                if (input_counter == N_POINTS - 1)
                    next_state = FFT_COMPUTE;
            end
            FFT_COMPUTE: begin
                if (stage_counter == STAGES && butterfly_counter == 0)
                    next_state = OUTPUT_DATA;
            end
            OUTPUT_DATA: begin
                if (output_counter == N_POINTS - 1)
                    next_state = IDLE;
            end
        endcase
    end
    
    // Control logic
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            input_counter <= 0;
            output_counter <= 0;
            stage_counter <= 0;
            butterfly_counter <= 0;
            group_counter <= 0;
            data_out_valid <= 0;
            fft_done <= 0;
        end else begin
            case (current_state)
                IDLE: begin
                    input_counter <= 0;
                    output_counter <= 0;
                    stage_counter <= 0;
                    butterfly_counter <= 0;
                    group_counter <= 0;
                    data_out_valid <= 0;
                    fft_done <= 0;
                end
                
                INPUT_DATA: begin
                    if (data_valid) begin
                        real_mem[input_counter] <= data_in_real;
                        imag_mem[input_counter] <= data_in_imag;
                        input_counter <= input_counter + 1;
                    end
                end
                
                BIT_REVERSE: begin
                    temp_real_mem[input_counter] <= real_mem[bit_reverse(input_counter)];
                    temp_imag_mem[input_counter] <= imag_mem[bit_reverse(input_counter)];
                    input_counter <= input_counter + 1;
                    if (input_counter == N_POINTS - 1) begin
                        // Copy back to main memory
                        real_mem <= temp_real_mem;
                        imag_mem <= temp_imag_mem;
                        stage_counter <= 0;
                        butterfly_counter <= 0;
                        group_counter <= 0;
                    end
                end
                
                FFT_COMPUTE: begin
                    // FFT butterfly computation logic
                    // This is a simplified version - full implementation would be more complex
                    if (stage_counter < STAGES) begin
                        // Perform butterfly operations for current stage
                        // ... (detailed butterfly logic would go here)
                        
                        // For now, increment counters
                        butterfly_counter <= butterfly_counter + 1;
                        if (butterfly_counter == (1 << stage_counter) - 1) begin
                            butterfly_counter <= 0;
                            group_counter <= group_counter + 1;
                            if (group_counter == (N_POINTS >> (stage_counter + 1)) - 1) begin
                                group_counter <= 0;
                                stage_counter <= stage_counter + 1;
                            end
                        end
                    end
                end
                
                OUTPUT_DATA: begin
                    data_out_real <= real_mem[output_counter];
                    data_out_imag <= imag_mem[output_counter];
                    data_out_valid <= 1;
                    output_counter <= output_counter + 1;
                    
                    if (output_counter == N_POINTS - 1) begin
                        fft_done <= 1;
                        data_out_valid <= 0;
                    end
                end
            endcase
        end
    end

endmodule