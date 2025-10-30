%% 说明: 使用自行编写的FFT，将CMA_reference.m中的核心算法在频域下实现
clc;
clear;
close all;

%% 自建FFT/IFFT路径
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
coreDir = fullfile(scriptDir, '..', 'core');
if exist(coreDir, 'dir')
    addpath(coreDir);
else
    warning('核心函数目录不存在: %s', coreDir);
end

%% Load Datas from Timing Recovery
load TD_TRdata.mat;
SigX_in = RTdataX(:);
SigY_in = RTdataY(:);

figure();
subplot(121);plot(RTdataX,'.r');title('Xpol before CMA');xlim([-2 2]);ylim([-2 2]);
subplot(122);plot(RTdataY,'.b');title('Ypol before CMA');xlim([-2 2]);ylim([-2 2]);

%% CMA参数设置 (与时域版本保持一致)
seg_len     = 32;
tap_len     = 33;
clk_dly     = 1;
par_num     = 4;
step        = 2^-8;
Rx          = 2;
Ry          = 2;

seg_num     = floor(length(SigX_in) / seg_len);
clk_num     = ceil(seg_num / par_num);

%% 数据预处理
total_len = seg_len * seg_num;
SigX_in = SigX_in(1:total_len);
SigY_in = SigY_in(1:total_len);

SigX_in_padded = [zeros(tap_len - 1, 1); SigX_in];
SigY_in_padded = [zeros(tap_len - 1, 1); SigY_in];

%% 初始化均衡器抽头 (双极化)
xx = zeros(tap_len, clk_num + clk_dly);
yy = zeros(tap_len, clk_num + clk_dly);
xy = zeros(tap_len, clk_num + clk_dly);
yx = zeros(tap_len, clk_num + clk_dly);

center_tap = (tap_len + 1) / 2;
xx(center_tap, :) = 1;
yy(center_tap, :) = 1;

%% 输出缓存
SigX_out = zeros(total_len, 1);
SigY_out = zeros(total_len, 1);
errx_vec = zeros(total_len, 1);
erry_vec = zeros(total_len, 1);

%% 抽头增量缓存
xx_accu = zeros(tap_len, 1);
yy_accu = zeros(tap_len, 1);
xy_accu = zeros(tap_len, 1);
yx_accu = zeros(tap_len, 1);

%% 频域卷积参数准备
block_input_len = seg_len + tap_len - 1;
conv_len = block_input_len + tap_len - 1;
nfft_conv = 2^nextpow2(conv_len);

pad_input_x = zeros(nfft_conv, 1);
pad_input_y = zeros(nfft_conv, 1);
pad_xx = zeros(nfft_conv, 1);
pad_xy = zeros(nfft_conv, 1);
pad_yx = zeros(nfft_conv, 1);
pad_yy = zeros(nfft_conv, 1);

XX_fft_curr = [];
XY_fft_curr = [];
YX_fft_curr = [];
YY_fft_curr = [];
prev_clk_pos = -1;

%% 第一步&第二步: 输入与抽头转换至频域, 频域乘积实现快速卷积
for k = 1:seg_num
    clk_pos = ceil(k / par_num);

    xx_curr = xx(:, clk_pos);
    yy_curr = yy(:, clk_pos);
    xy_curr = xy(:, clk_pos);
    yx_curr = yx(:, clk_pos);

    block_start = (k - 1) * seg_len + 1;
    block_end   = block_start + seg_len - 1;

    % 提取包含前导样本的分段数据
    block_input_x = SigX_in_padded(block_start:block_end + tap_len - 1);
    block_input_y = SigY_in_padded(block_start:block_end + tap_len - 1);

    pad_input_x(:) = 0;
    pad_input_y(:) = 0;
    pad_input_x(1:block_input_len) = block_input_x;
    pad_input_y(1:block_input_len) = block_input_y;

    if clk_pos ~= prev_clk_pos
        pad_xx(:) = 0;
        pad_xy(:) = 0;
        pad_yx(:) = 0;
        pad_yy(:) = 0;

        pad_xx(1:tap_len) = xx_curr;
        pad_xy(1:tap_len) = xy_curr;
        pad_yx(1:tap_len) = yx_curr;
        pad_yy(1:tap_len) = yy_curr;

        XX_fft_curr = my_fft(pad_xx);
        XY_fft_curr = my_fft(pad_xy);
        YX_fft_curr = my_fft(pad_yx);
        YY_fft_curr = my_fft(pad_yy);
        prev_clk_pos = clk_pos;
    end

    X_fft = my_fft(pad_input_x);
    Y_fft = my_fft(pad_input_y);

    conv_xx = my_ifft(X_fft .* XX_fft_curr);
    conv_xy = my_ifft(Y_fft .* XY_fft_curr);
    conv_yx = my_ifft(X_fft .* YX_fft_curr);
    conv_yy = my_ifft(Y_fft .* YY_fft_curr);

    % 截取有效输出 (丢弃前tap_len-1个因重叠产生的样本)
    block_out_x = conv_xx(tap_len:tap_len + seg_len - 1) + conv_xy(tap_len:tap_len + seg_len - 1);
    block_out_y = conv_yx(tap_len:tap_len + seg_len - 1) + conv_yy(tap_len:tap_len + seg_len - 1);

    SigX_out(block_start:block_end) = block_out_x;
    SigY_out(block_start:block_end) = block_out_y;

    %% 第三步: 计算误差并在频域卷积的基础上更新抽头
    errx_block = Rx - abs(block_out_x).^2;
    erry_block = Ry - abs(block_out_y).^2;

    errx_vec(block_start:block_end) = errx_block;
    erry_vec(block_start:block_end) = erry_block;

    grad_x = step * (errx_block .* block_out_x);
    grad_y = step * (erry_block .* block_out_y);

    % 生成输入矩阵用于梯度累积
    x_mat = toeplitz(block_input_x(tap_len:end), block_input_x(tap_len:-1:1));
    y_mat = toeplitz(block_input_y(tap_len:end), block_input_y(tap_len:-1:1));

    xx_accu = xx_accu + x_mat' * grad_x;
    xy_accu = xy_accu + y_mat' * grad_x;
    yx_accu = yx_accu + x_mat' * grad_y;
    yy_accu = yy_accu + y_mat' * grad_y;

    if mod(k, par_num) == 0
        update_clk_pos = clk_pos + clk_dly;

        xx(:, update_clk_pos) = xx(:, update_clk_pos - 1) + xx_accu / par_num;
        yy(:, update_clk_pos) = yy(:, update_clk_pos - 1) + yy_accu / par_num;
        xy(:, update_clk_pos) = xy(:, update_clk_pos - 1) + xy_accu / par_num;
        yx(:, update_clk_pos) = yx(:, update_clk_pos - 1) + yx_accu / par_num;

        xx_accu(:) = 0;
        yy_accu(:) = 0;
        xy_accu(:) = 0;
        yx_accu(:) = 0;
    end
end

%% 输出处理与可视化
CMAdataX = SigX_out.';
CMAdataY = SigY_out.';

final_xx = xx(:, clk_num);
final_yy = yy(:, clk_num);
final_xy = xy(:, clk_num);
final_yx = yx(:, clk_num);
Tap_Matx = [[final_xx; final_xy] [final_yx; final_yy]];

errx_reshaped = reshape(errx_vec, seg_len, seg_num);
erry_reshaped = reshape(erry_vec, seg_len, seg_num);

figure();
subplot(121);plot(CMAdataX,'.r');title('Xpol after CMA');xlim([-2 2]);ylim([-2 2]);
subplot(122);plot(CMAdataY,'.b');title('Ypol after CMA');xlim([-2 2]);ylim([-2 2]);

CMAdataX = CMAdataX(end-131072:end);
CMAdataY = CMAdataY(end-131072:end);

figure();
subplot(121);plot(CMAdataX,'.r');title('Last part of Xpol after CMA');xlim([-2 2]);ylim([-2 2]);
subplot(122);plot(CMAdataY,'.b');title('Last part of Ypol after CMA');xlim([-2 2]);ylim([-2 2]);