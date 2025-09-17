`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: DSP Course Design
// Engineer: Student
// 
// Create Date: 2025
// Design Name: CMA Testbench
// Module Name: tb_cma_equalizer
// Project Name: Digital Signal Processing Giant Design
// Description: 
//   Testbench for CMA Equalizer module
//
//////////////////////////////////////////////////////////////////////////////////

module tb_cma_equalizer;

    // Parameters
    parameter DATA_WIDTH = 16;
    parameter NUM_TAPS = 11;
    parameter FRAC_BITS = 12;
    parameter STEP_SIZE = 41;  // 0.01 in fixed-point (41/4096 ≈ 0.01)
    parameter CLK_PERIOD = 10;
    parameter NUM_SYMBOLS = 200;
    
    // Testbench signals
    reg clk;
    reg rst_n;
    reg enable;
    reg signed [DATA_WIDTH-1:0] data_in_real;
    reg signed [DATA_WIDTH-1:0] data_in_imag;
    reg data_valid;
    
    wire signed [DATA_WIDTH-1:0] data_out_real;
    wire signed [DATA_WIDTH-1:0] data_out_imag;
    wire data_out_valid;
    wire [31:0] error_magnitude;
    
    // Test data storage
    reg signed [DATA_WIDTH-1:0] qpsk_symbols_real [0:NUM_SYMBOLS-1];
    reg signed [DATA_WIDTH-1:0] qpsk_symbols_imag [0:NUM_SYMBOLS-1];
    reg signed [DATA_WIDTH-1:0] channel_output_real [0:NUM_SYMBOLS-1];
    reg signed [DATA_WIDTH-1:0] channel_output_imag [0:NUM_SYMBOLS-1];
    reg signed [DATA_WIDTH-1:0] equalizer_output_real [0:NUM_SYMBOLS-1];
    reg signed [DATA_WIDTH-1:0] equalizer_output_imag [0:NUM_SYMBOLS-1];
    
    // Channel model coefficients (simple multipath)
    reg signed [DATA_WIDTH-1:0] channel_h0_real, channel_h0_imag;
    reg signed [DATA_WIDTH-1:0] channel_h1_real, channel_h1_imag;
    reg signed [DATA_WIDTH-1:0] channel_h2_real, channel_h2_imag;
    
    integer i, j;
    integer symbol_count, output_count;
    integer seed = 123456;
    
    // Instantiate Unit Under Test (UUT)
    cma_equalizer #(
        .DATA_WIDTH(DATA_WIDTH),
        .NUM_TAPS(NUM_TAPS),
        .FRAC_BITS(FRAC_BITS),
        .STEP_SIZE(STEP_SIZE)
    ) uut (
        .clk(clk),
        .rst_n(rst_n),
        .enable(enable),
        .data_in_real(data_in_real),
        .data_in_imag(data_in_imag),
        .data_valid(data_valid),
        .data_out_real(data_out_real),
        .data_out_imag(data_out_imag),
        .data_out_valid(data_out_valid),
        .error_magnitude(error_magnitude)
    );
    
    // Clock generation
    initial begin
        clk = 0;
        forever #(CLK_PERIOD/2) clk = ~clk;
    end
    
    // Random number generator function (simplified LFSR)
    function integer random_bit;
        input integer seed_val;
        begin
            seed = (seed * 1103515245 + 12345) & 32'h7fffffff;
            random_bit = (seed >> 16) & 1;
        end
    endfunction
    
    // QPSK symbol generation
    task generate_qpsk_symbols;
        begin
            for (i = 0; i < NUM_SYMBOLS; i = i + 1) begin
                // Generate random QPSK symbols {+1, -1} + j{+1, -1}
                if (random_bit(seed))
                    qpsk_symbols_real[i] = (1 << (FRAC_BITS-2));  // +0.25 in fixed-point
                else
                    qpsk_symbols_real[i] = -(1 << (FRAC_BITS-2)); // -0.25 in fixed-point
                    
                if (random_bit(seed))
                    qpsk_symbols_imag[i] = (1 << (FRAC_BITS-2));  // +0.25 in fixed-point
                else
                    qpsk_symbols_imag[i] = -(1 << (FRAC_BITS-2)); // -0.25 in fixed-point
            end
        end
    endtask
    
    // Channel simulation (simple multipath with 3 taps)
    task simulate_channel;
        reg signed [2*DATA_WIDTH-1:0] temp_real, temp_imag;
        begin
            // Channel coefficients: h = [0.8, 0.4, 0.2]
            channel_h0_real = $rtoi(0.8 * (1 << FRAC_BITS));   // 0.8
            channel_h0_imag = 0;
            channel_h1_real = $rtoi(0.4 * (1 << FRAC_BITS));   // 0.4
            channel_h1_imag = 0;
            channel_h2_real = $rtoi(0.2 * (1 << FRAC_BITS));   // 0.2
            channel_h2_imag = 0;
            
            // Convolve symbols with channel
            for (i = 0; i < NUM_SYMBOLS; i = i + 1) begin
                temp_real = 0;
                temp_imag = 0;
                
                // h0 * x[i]
                temp_real = temp_real + (qpsk_symbols_real[i] * channel_h0_real - qpsk_symbols_imag[i] * channel_h0_imag);
                temp_imag = temp_imag + (qpsk_symbols_real[i] * channel_h0_imag + qpsk_symbols_imag[i] * channel_h0_real);
                
                // h1 * x[i-1]
                if (i >= 1) begin
                    temp_real = temp_real + (qpsk_symbols_real[i-1] * channel_h1_real - qpsk_symbols_imag[i-1] * channel_h1_imag);
                    temp_imag = temp_imag + (qpsk_symbols_real[i-1] * channel_h1_imag + qpsk_symbols_imag[i-1] * channel_h1_real);
                end
                
                // h2 * x[i-2]
                if (i >= 2) begin
                    temp_real = temp_real + (qpsk_symbols_real[i-2] * channel_h2_real - qpsk_symbols_imag[i-2] * channel_h2_imag);
                    temp_imag = temp_imag + (qpsk_symbols_real[i-2] * channel_h2_imag + qpsk_symbols_imag[i-2] * channel_h2_real);
                end
                
                // Scale down and store
                channel_output_real[i] = temp_real >>> FRAC_BITS;
                channel_output_imag[i] = temp_imag >>> FRAC_BITS;
            end
        end
    endtask
    
    // Test stimulus
    initial begin
        // Initialize signals
        rst_n = 0;
        enable = 0;
        data_in_real = 0;
        data_in_imag = 0;
        data_valid = 0;
        symbol_count = 0;
        output_count = 0;
        
        $display("Starting CMA Equalizer Test...");
        
        // Generate test data
        generate_qpsk_symbols();
        $display("Generated %d QPSK symbols", NUM_SYMBOLS);
        
        // Simulate channel
        simulate_channel();
        $display("Applied channel distortion");
        
        // Display first few symbols
        $display("Original symbols (first 8):");
        for (i = 0; i < 8; i = i + 1) begin
            $display("Symbol %d: Real=%d, Imag=%d", i, qpsk_symbols_real[i], qpsk_symbols_imag[i]);
        end
        
        $display("Channel output (first 8):");
        for (i = 0; i < 8; i = i + 1) begin
            $display("Sample %d: Real=%d, Imag=%d", i, channel_output_real[i], channel_output_imag[i]);
        end
        
        // Reset and start
        #(10*CLK_PERIOD);
        rst_n = 1;
        #(5*CLK_PERIOD);
        enable = 1;
        
        // Feed symbols to equalizer
        for (i = 0; i < NUM_SYMBOLS; i = i + 1) begin
            @(posedge clk);
            data_in_real = channel_output_real[i];
            data_in_imag = channel_output_imag[i];
            data_valid = 1;
            
            // Collect output if valid
            if (data_out_valid) begin
                equalizer_output_real[output_count] = data_out_real;
                equalizer_output_imag[output_count] = data_out_imag;
                output_count = output_count + 1;
                
                // Print some outputs for monitoring
                if (output_count <= 10 || output_count % 20 == 0) begin
                    $display("Output %d: Real=%d, Imag=%d, Error_mag=%d", 
                             output_count-1, data_out_real, data_out_imag, error_magnitude);
                end
            end
            
            #CLK_PERIOD;
            data_valid = 0;
            #(2*CLK_PERIOD); // Some spacing between symbols
        end
        
        // Continue collecting any remaining outputs
        for (i = 0; i < 50; i = i + 1) begin
            @(posedge clk);
            if (data_out_valid) begin
                equalizer_output_real[output_count] = data_out_real;
                equalizer_output_imag[output_count] = data_out_imag;
                output_count = output_count + 1;
                $display("Final Output %d: Real=%d, Imag=%d, Error_mag=%d", 
                         output_count-1, data_out_real, data_out_imag, error_magnitude);
            end
        end
        
        $display("CMA Equalization completed!");
        $display("Total outputs collected: %d", output_count);
        
        // Analyze convergence by looking at last few outputs
        if (output_count >= 10) begin
            $display("Final equalized symbols (last 10):");
            for (i = output_count-10; i < output_count; i = i + 1) begin
                $display("Symbol %d: Real=%d, Imag=%d", i, 
                         equalizer_output_real[i], equalizer_output_imag[i]);
            end
        end
        
        // Simple constellation analysis
        integer pos_real_count, neg_real_count, pos_imag_count, neg_imag_count;
        pos_real_count = 0; neg_real_count = 0;
        pos_imag_count = 0; neg_imag_count = 0;
        
        // Analyze last half of outputs (assuming convergence)
        for (i = output_count/2; i < output_count; i = i + 1) begin
            if (equalizer_output_real[i] > 0) pos_real_count = pos_real_count + 1;
            else neg_real_count = neg_real_count + 1;
            
            if (equalizer_output_imag[i] > 0) pos_imag_count = pos_imag_count + 1;
            else neg_imag_count = neg_imag_count + 1;
        end
        
        $display("Constellation analysis (last half of outputs):");
        $display("Positive real: %d, Negative real: %d", pos_real_count, neg_real_count);
        $display("Positive imag: %d, Negative imag: %d", pos_imag_count, neg_imag_count);
        
        $display("\nCMA Equalizer test completed successfully!");
        $finish;
    end
    
    // Error monitoring
    always @(posedge clk) begin
        if (data_out_valid && enable) begin
            // Log error magnitude every 10 symbols
            if (output_count % 10 == 0) begin
                $display("Time %t: Symbol %d, Error magnitude = %d", $time, output_count, error_magnitude);
            end
        end
    end
    
    // Timeout
    initial begin
        #(5000000); // 5ms timeout
        $display("ERROR: Test timeout!");
        $finish;
    end

endmodule