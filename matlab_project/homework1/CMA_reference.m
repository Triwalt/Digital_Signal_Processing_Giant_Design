%% 说明: 以下代码基于时域CMA，完成数据的解调
st = dbstack;
if numel(st) == 1
    clc;
    clear;
    close all;
end

%% Load Datas 
load TD_TRdata.mat;
SigX_in = RTdataX;
SigY_in = RTdataY;
figure();
subplot(121);plot(RTdataX,'.r');title('X-pol before CMA');xlim([-2 2]);ylim([-2 2]); grid on;
subplot(122);plot(RTdataY,'.b');title('Y-pol before CMA');xlim([-2 2]);ylim([-2 2]); grid on;

%% Parameters
seg_len     = 32;                               % block corresponding data length
tap_len     = 33;                               % FIR filter tap length
clk_dly     = 1;                                % Tap update delay in clock cycles (simulate FPGA application scenarios)
par_num     = 4;                                % Number of blocks processed in parallel per clock cycle
step        = 2^-8;                             % Step size for tap updates
Rx          = 2;                                % Convergence radius for X polarization
Ry          = 2;                                % Convergence radius for Y polarization
seg_num     = floor(length(SigX_in)/seg_len);   % Total number of blocks in the data
clk_num     = ceil(seg_num/par_num);            % Total number of clock cycles to run

%% Prepare Input Signals for Time-Domain Convolution
total_len = seg_len * seg_num;
SigX_in = SigX_in(1:total_len);
SigY_in = SigY_in(1:total_len);

SigX_in_padded = [zeros(tap_len-1, 1); SigX_in(:)];
SigY_in_padded = [zeros(tap_len-1, 1); SigY_in(:)];

%% Taps Initialization
xx = zeros(tap_len, clk_num + clk_dly);
yy = zeros(tap_len, clk_num + clk_dly);
xy = zeros(tap_len, clk_num + clk_dly);
yx = zeros(tap_len, clk_num + clk_dly);

center_tap = (tap_len + 1) / 2;
xx(center_tap, :) = 1;
yy(center_tap, :) = 1;

%% Time Domain CMA
SigX_out = zeros(total_len, 1);
SigY_out = zeros(total_len, 1);
errx_vec = zeros(total_len, 1);
erry_vec = zeros(total_len, 1);

xx_accu = zeros(tap_len, 1);
yy_accu = zeros(tap_len, 1);
xy_accu = zeros(tap_len, 1);
yx_accu = zeros(tap_len, 1);

for k = 1:seg_num
    clk_pos = ceil(k / par_num);
    
    xx_curr = xx(:, clk_pos);
    yy_curr = yy(:, clk_pos);
    xy_curr = xy(:, clk_pos);
    yx_curr = yx(:, clk_pos);

    for i = 1:seg_len
        current_idx = (k - 1) * seg_len + i;
        
        padded_idx = current_idx + tap_len - 1;
        input_x = SigX_in_padded(padded_idx:-1:padded_idx - tap_len + 1);
        input_y = SigY_in_padded(padded_idx:-1:padded_idx - tap_len + 1);
        
        SigX_out(current_idx) = xx_curr.' * input_x + xy_curr.' * input_y;
        SigY_out(current_idx) = yx_curr.' * input_x + yy_curr.' * input_y;
        
        errx_vec(current_idx) = Rx - abs(SigX_out(current_idx))^2;
        erry_vec(current_idx) = Ry - abs(SigY_out(current_idx))^2;
        
        xx_accu = xx_accu + step * errx_vec(current_idx) * SigX_out(current_idx) * conj(input_x);
        yy_accu = yy_accu + step * erry_vec(current_idx) * SigY_out(current_idx) * conj(input_y);
        xy_accu = xy_accu + step * errx_vec(current_idx) * SigX_out(current_idx) * conj(input_y);
        yx_accu = yx_accu + step * erry_vec(current_idx) * SigY_out(current_idx) * conj(input_x);
    end
    
    if mod(k, par_num) == 0
        update_clk_pos = clk_pos + clk_dly;
        
        % Update taps with the accumulated adjustments, averaged over the parallel blocks
        xx(:, update_clk_pos) = xx(:, update_clk_pos - 1) + xx_accu / par_num;
        yy(:, update_clk_pos) = yy(:, update_clk_pos - 1) + yy_accu / par_num;
        xy(:, update_clk_pos) = xy(:, update_clk_pos - 1) + xy_accu / par_num;
        yx(:, update_clk_pos) = yx(:, update_clk_pos - 1) + yx_accu / par_num;
    
        % Reset accumulators for the next set of blocks
        xx_accu = zeros(tap_len, 1);
        yy_accu = zeros(tap_len, 1);
        xy_accu = zeros(tap_len, 1);
        yx_accu = zeros(tap_len, 1);
    end
end

CMAdataX = SigX_out.';
CMAdataY = SigY_out.';

final_xx = xx(:, clk_num);
final_yy = yy(:, clk_num);
final_xy = xy(:, clk_num);
final_yx = yx(:, clk_num);
Tap_Matx = [[final_xx; final_xy] [final_yx; final_yy]];


errx_reshaped = reshape(errx_vec, seg_len, seg_num);
erry_reshaped = reshape(erry_vec, seg_len, seg_num);

%% Figure and Results
figure();
subplot(121);plot(CMAdataX,'.r');title('Xpol after CMA');xlim([-2 2]);ylim([-2 2]);
subplot(122);plot(CMAdataY,'.b');title('Ypol after CMA');xlim([-2 2]);ylim([-2 2]);

CMAdataX=CMAdataX(end-131072:end);
CMAdataY=CMAdataY(end-131072:end);

figure();
subplot(121);plot(CMAdataX,'.r');title('Last part of Xpol after CMA');xlim([-2 2]);ylim([-2 2]);
subplot(122);plot(CMAdataY,'.b');title('Last part of Ypol after CMA');xlim([-2 2]);ylim([-2 2]);