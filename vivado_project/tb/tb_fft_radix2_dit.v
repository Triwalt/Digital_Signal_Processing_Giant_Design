`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: DSP Course Design
// Engineer: Student
// 
// Create Date: 2025
// Design Name: FFT Testbench
// Module Name: tb_fft_radix2_dit
// Project Name: Digital Signal Processing Giant Design
// Description: 
//   Testbench for FFT Radix-2 DIT module
//
//////////////////////////////////////////////////////////////////////////////////

module tb_fft_radix2_dit;

    // Parameters
    parameter N_POINTS = 64;
    parameter DATA_WIDTH = 16;
    parameter TWIDDLE_WIDTH = 16;
    parameter CLK_PERIOD = 10; // 100MHz clock
    
    // Testbench signals
    reg clk;
    reg rst_n;
    reg start;
    reg [DATA_WIDTH-1:0] data_in_real;
    reg [DATA_WIDTH-1:0] data_in_imag;
    reg data_valid;
    
    wire [DATA_WIDTH-1:0] data_out_real;
    wire [DATA_WIDTH-1:0] data_out_imag;
    wire data_out_valid;
    wire fft_done;
    
    // Test data storage
    reg signed [DATA_WIDTH-1:0] test_data_real [0:N_POINTS-1];
    reg signed [DATA_WIDTH-1:0] test_data_imag [0:N_POINTS-1];
    reg signed [DATA_WIDTH-1:0] output_data_real [0:N_POINTS-1];
    reg signed [DATA_WIDTH-1:0] output_data_imag [0:N_POINTS-1];
    
    integer i, j;
    integer input_count, output_count;
    
    // Instantiate Unit Under Test (UUT)
    fft_radix2_dit #(
        .N_POINTS(N_POINTS),
        .DATA_WIDTH(DATA_WIDTH),
        .TWIDDLE_WIDTH(TWIDDLE_WIDTH)
    ) uut (
        .clk(clk),
        .rst_n(rst_n),
        .start(start),
        .data_in_real(data_in_real),
        .data_in_imag(data_in_imag),
        .data_valid(data_valid),
        .data_out_real(data_out_real),
        .data_out_imag(data_out_imag),
        .data_out_valid(data_out_valid),
        .fft_done(fft_done)
    );
    
    // Clock generation
    initial begin
        clk = 0;
        forever #(CLK_PERIOD/2) clk = ~clk;
    end
    
    // Test stimulus
    initial begin
        // Initialize signals
        rst_n = 0;
        start = 0;
        data_in_real = 0;
        data_in_imag = 0;
        data_valid = 0;
        input_count = 0;
        output_count = 0;
        
        // Generate test data: simple sinusoid
        for (i = 0; i < N_POINTS; i = i + 1) begin
            // Generate sinusoid: cos(2*pi*4*i/N) + j*sin(2*pi*4*i/N)
            test_data_real[i] = $rtoi($cos(2.0 * 3.14159265359 * 4 * i / N_POINTS) * (2**(DATA_WIDTH-2) - 1));
            test_data_imag[i] = $rtoi($sin(2.0 * 3.14159265359 * 4 * i / N_POINTS) * (2**(DATA_WIDTH-2) - 1));
        end
        
        // Reset
        #(10*CLK_PERIOD);
        rst_n = 1;
        #(5*CLK_PERIOD);
        
        $display("Starting FFT Test...");
        $display("Input data (first 8 samples):");
        for (i = 0; i < 8; i = i + 1) begin
            $display("Sample %d: Real=%d, Imag=%d", i, test_data_real[i], test_data_imag[i]);
        end
        
        // Start FFT
        start = 1;
        #CLK_PERIOD;
        start = 0;
        
        // Feed input data
        #(2*CLK_PERIOD);
        for (i = 0; i < N_POINTS; i = i + 1) begin
            @(posedge clk);
            data_in_real = test_data_real[i];
            data_in_imag = test_data_imag[i];
            data_valid = 1;
            #CLK_PERIOD;
            data_valid = 0;
            #CLK_PERIOD;
        end
        
        // Wait for FFT completion
        wait(fft_done);
        $display("FFT computation completed!");
        
        // Collect output data
        output_count = 0;
        while (output_count < N_POINTS) begin
            @(posedge clk);
            if (data_out_valid) begin
                output_data_real[output_count] = data_out_real;
                output_data_imag[output_count] = data_out_imag;
                $display("Output %d: Real=%d, Imag=%d", output_count, data_out_real, data_out_imag);
                output_count = output_count + 1;
            end
        end
        
        $display("FFT Test completed!");
        $display("Output data (first 8 samples):");
        for (i = 0; i < 8; i = i + 1) begin
            $display("Sample %d: Real=%d, Imag=%d", i, output_data_real[i], output_data_imag[i]);
        end
        
        // Find peak frequency bin
        integer max_magnitude, max_bin;
        integer magnitude;
        max_magnitude = 0;
        max_bin = 0;
        
        for (i = 0; i < N_POINTS/2; i = i + 1) begin
            magnitude = (output_data_real[i] * output_data_real[i]) + 
                       (output_data_imag[i] * output_data_imag[i]);
            if (magnitude > max_magnitude) begin
                max_magnitude = magnitude;
                max_bin = i;
            end
        end
        
        $display("Peak frequency bin: %d (expected: 4)", max_bin);
        
        // Test 2: Impulse response
        $display("\nTesting impulse response...");
        
        // Reset for new test
        rst_n = 0;
        #(5*CLK_PERIOD);
        rst_n = 1;
        #(5*CLK_PERIOD);
        
        // Generate impulse
        for (i = 0; i < N_POINTS; i = i + 1) begin
            if (i == 0) begin
                test_data_real[i] = (1 << (DATA_WIDTH-2));
                test_data_imag[i] = 0;
            end else begin
                test_data_real[i] = 0;
                test_data_imag[i] = 0;
            end
        end
        
        // Start FFT
        start = 1;
        #CLK_PERIOD;
        start = 0;
        
        // Feed input data
        #(2*CLK_PERIOD);
        for (i = 0; i < N_POINTS; i = i + 1) begin
            @(posedge clk);
            data_in_real = test_data_real[i];
            data_in_imag = test_data_imag[i];
            data_valid = 1;
            #CLK_PERIOD;
            data_valid = 0;
            #CLK_PERIOD;
        end
        
        // Wait for completion
        wait(fft_done);
        $display("Impulse FFT completed!");
        
        // Collect output
        output_count = 0;
        while (output_count < N_POINTS) begin
            @(posedge clk);
            if (data_out_valid) begin
                output_data_real[output_count] = data_out_real;
                output_data_imag[output_count] = data_out_imag;
                output_count = output_count + 1;
            end
        end
        
        $display("Impulse response (first 8 samples):");
        for (i = 0; i < 8; i = i + 1) begin
            $display("Sample %d: Real=%d, Imag=%d", i, output_data_real[i], output_data_imag[i]);
        end
        
        $display("\nAll FFT tests completed!");
        $finish;
    end
    
    // Monitor for debugging
    initial begin
        $monitor("Time=%t, State=%s, Input_count=%d, Output_count=%d", 
                 $time, uut.current_state.name(), input_count, output_count);
    end
    
    // Timeout
    initial begin
        #(1000000); // 1ms timeout
        $display("ERROR: Test timeout!");
        $finish;
    end

endmodule