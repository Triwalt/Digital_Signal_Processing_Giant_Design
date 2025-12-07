%% main_cma_TRdata.m: TD_TRdata.mat 双极化数据均衡处理入口
%
% 功能说明:
%   对 homework1/TD_TRdata.mat 中的双极化信号 (RTdataX, RTdataY) 进行
%   2x2 MIMO CMA 盲均衡处理，支持时域和频域两种实现方式。
%
% 配置开关:
%   cma_use_freq_domain  - true:  调用 homework1/CMA_homework.m (频域CMA，使用自编FFT)
%                        - false: 使用本文件内的时域CMA实现 (与CMA_reference.m一致)
%   enable_vv_phase_lock - true:  启用V&V相位恢复算法
%                        - false: 跳过相位恢复
%   enable_plot_output   - true:  输出相位漂移估计图
%                        - false: 不输出相位漂移估计图
%
% 输出:
%   - 均衡前星座图 (X/Y极化)
%   - 均衡后星座图 (含相位恢复)
%   - 相位漂移估计曲线
%   - 末段稳定数据星座图
%
% 相关文件:
%   - homework1/CMA_homework.m   : 频域CMA实现 (使用 core/my_fft.m)
%   - homework1/CMA_reference.m  : 时域CMA参考实现
%   - homework1/TD_TRdata.mat    : 双极化实验数据

clc;
clear;
close all;

%% ==================== 配置开关 ====================
cma_use_freq_domain  = false;   % true: 频域CMA (自编FFT), false: 时域CMA
enable_vv_phase_lock = true;    % true: 启用V&V相位恢复
enable_plot_output   = true;    % true: 启用图形输出
%% =================================================
if cma_use_freq_domain
    scriptDir = fileparts(mfilename('fullpath'));
    if isempty(scriptDir)
        scriptDir = pwd;
    end
    projectRoot = fileparts(scriptDir);
    hwDir = fullfile(projectRoot, 'homework1');
    if ~exist(hwDir, 'dir')
        error('homework1 目录不存在: %s', hwDir);
    end
    addpath(hwDir);
    % 配置开关传递给CMA_homework.m (变量在workspace中预设)
    run(fullfile(hwDir, 'CMA_homework.m'));
    return;
end

%% 路径与数据加载
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

projectRoot = fileparts(scriptDir);                 % e.g. .../matlab_project
dataDir     = fullfile(projectRoot, 'homework1');

dataFile = fullfile(dataDir, 'TD_TRdata.mat');
if ~exist(dataFile, 'file')
    error('找不到数据文件: %s', dataFile);
end

S = load(dataFile);
if ~isfield(S, 'RTdataX') || ~isfield(S, 'RTdataY')
    error('TD_TRdata.mat 中必须包含变量 RTdataX 和 RTdataY');
end

SigX_in = S.RTdataX(:);    % 转为列向量
SigY_in = S.RTdataY(:);

fprintf('Loaded TD_TRdata.mat from: %s\n', dataFile);
fprintf('  Length of RTdataX: %d samples\n', length(SigX_in));
fprintf('  Length of RTdataY: %d samples\n', length(SigY_in));

%% 均衡前星座图
figure();
subplot(1, 2, 1);
plot(SigX_in, '.r');
title('X-pol before CMA'); grid on; axis equal;
xlim([-2 2]); ylim([-2 2]);

subplot(1, 2, 2);
plot(SigY_in, '.b');
title('Y-pol before CMA'); grid on; axis equal;
xlim([-2 2]); ylim([-2 2]);

%% CMA 参数设置 (与 CMA_reference.m 保持一致)
seg_len     = 32;      % 每块处理的符号数
tap_len     = 33;      % FIR 滤波器抽头长度
clk_dly     = 1;       % 抽头更新时钟延迟
par_num     = 4;       % 每个时钟并行处理的块数
step        = 2^-8;    % 步长
Rx          = 2;       % X 极化收敛半径
Ry          = 2;       % Y 极化收敛半径

seg_num     = floor(length(SigX_in) / seg_len);  % 有效块数
clk_num     = ceil(seg_num / par_num);           % 时钟数

fprintf('\nCMA configuration:\n');
fprintf('  seg_len = %d, tap_len = %d\n', seg_len, tap_len);
fprintf('  par_num = %d, clk_dly = %d\n', par_num, clk_dly);
fprintf('  step = %.5f, Rx = %.1f, Ry = %.1f\n', step, Rx, Ry);
fprintf('  seg_num = %d, clk_num = %d\n', seg_num, clk_num);

%% 输入信号预处理
% 截断到整块长度
total_len = seg_len * seg_num;
SigX_in = SigX_in(1:total_len);
SigY_in = SigY_in(1:total_len);

% 在前面补 tap_len-1 个零, 便于构造输入向量
a = tap_len - 1;
SigX_in_padded = [zeros(a, 1); SigX_in(:)];
SigY_in_padded = [zeros(a, 1); SigY_in(:)];

%% 抽头矩阵初始化 (2x2 MIMO)
xx = zeros(tap_len, clk_num + clk_dly);   % X <- X
yy = zeros(tap_len, clk_num + clk_dly);   % Y <- Y
xy = zeros(tap_len, clk_num + clk_dly);   % X <- Y
yx = zeros(tap_len, clk_num + clk_dly);   % Y <- X

center_tap = (tap_len + 1) / 2;
xx(center_tap, :) = 1;
yy(center_tap, :) = 1;

%% 时域 CMA 主循环
SigX_out = zeros(total_len, 1);
SigY_out = zeros(total_len, 1);
errx_vec = zeros(total_len, 1);
erry_vec = zeros(total_len, 1);

xx_accu = zeros(tap_len, 1);
yy_accu = zeros(tap_len, 1);
xy_accu = zeros(tap_len, 1);
yx_accu = zeros(tap_len, 1);

fprintf('\nStart time-domain 2x2 CMA ...\n');

for k = 1:seg_num
    clk_pos = ceil(k / par_num);

    % 当前时钟下的 4 组抽头
    xx_curr = xx(:, clk_pos);
    yy_curr = yy(:, clk_pos);
    xy_curr = xy(:, clk_pos);
    yx_curr = yx(:, clk_pos);

    for i = 1:seg_len
        current_idx = (k - 1) * seg_len + i;

        % 对应到补零后的索引
        padded_idx = current_idx + a;
        input_x = SigX_in_padded(padded_idx:-1:padded_idx - tap_len + 1);
        input_y = SigY_in_padded(padded_idx:-1:padded_idx - tap_len + 1);

        % 线性组合得到双极化输出
        SigX_out(current_idx) = xx_curr.' * input_x + xy_curr.' * input_y;
        SigY_out(current_idx) = yx_curr.' * input_x + yy_curr.' * input_y;

        % CMA 误差
        errx_vec(current_idx) = Rx - abs(SigX_out(current_idx))^2;
        erry_vec(current_idx) = Ry - abs(SigY_out(current_idx))^2;

        % 梯度累积 (与 CMA_reference 形式一致)
        xx_accu = xx_accu + step * errx_vec(current_idx) * SigX_out(current_idx) * conj(input_x);
        yy_accu = yy_accu + step * erry_vec(current_idx) * SigY_out(current_idx) * conj(input_y);
        xy_accu = xy_accu + step * errx_vec(current_idx) * SigX_out(current_idx) * conj(input_y);
        yx_accu = yx_accu + step * erry_vec(current_idx) * SigY_out(current_idx) * conj(input_x);
    end

    % 每处理 par_num 个块, 更新一次下一时钟位置的抽头
    if mod(k, par_num) == 0
        update_clk_pos = clk_pos + clk_dly;

        xx(:, update_clk_pos) = xx(:, update_clk_pos - 1) + xx_accu / par_num;
        yy(:, update_clk_pos) = yy(:, update_clk_pos - 1) + yy_accu / par_num;
        xy(:, update_clk_pos) = xy(:, update_clk_pos - 1) + xy_accu / par_num;
        yx(:, update_clk_pos) = yx(:, update_clk_pos - 1) + yx_accu / par_num;

        % 清零累积量
        xx_accu(:) = 0;
        yy_accu(:) = 0;
        xy_accu(:) = 0;
        yx_accu(:) = 0;
    end

    if mod(k, 100) == 0 || k == seg_num
        fprintf('  Processed segment %d / %d (%.1f%%)\n', ...
            k, seg_num, 100 * k / seg_num);
    end
end

fprintf('CMA finished.\n');

%% 输出整理
CMAdataX = SigX_out.';   % 行向量, 便于绘图
CMAdataY = SigY_out.';

yX = CMAdataX(:);
yY = CMAdataY(:);
N = length(yX);
zX = yX.^4;
zY = yY.^4;
win_half = 3;
win = ones(2*win_half+1, 1) / (2*win_half+1);
fX = conv(zX, win, 'same');
fY = conv(zY, win, 'same');
thetaX = angle(fX) / 4;
thetaY = angle(fY) / 4;
if enable_vv_phase_lock
    yX_corr = yX .* exp(-1j * thetaX);
    yY_corr = yY .* exp(-1j * thetaY);
    CMAdataX = yX_corr.';
    CMAdataY = yY_corr.';
end

final_xx = xx(:, clk_num);
final_yy = yy(:, clk_num);
final_xy = xy(:, clk_num);
final_yx = yx(:, clk_num);
Tap_Matx = [[final_xx; final_xy] [final_yx; final_yy]]; %#ok<NASGU>

errx_reshaped = reshape(errx_vec, seg_len, seg_num); %#ok<NASGU>
erry_reshaped = reshape(erry_vec, seg_len, seg_num); %#ok<NASGU>

%% 均衡后星座图
figure();
subplot(1, 2, 1);
plot(CMAdataX, '.r');
title('X-pol after CMA'); grid on; axis equal;
xlim([-2 2]); ylim([-2 2]);

subplot(1, 2, 2);
plot(CMAdataY, '.b');
title('Y-pol after CMA'); grid on; axis equal;
xlim([-2 2]); ylim([-2 2]);

%% 相位漂移估计图 (由 enable_plot_output 控制)
if enable_plot_output
    figure();
    subplot(2, 1, 1);
    plot(-499:0, unwrap(thetaX(end-499:end)));  % Show last 500 points
    xlabel('Sample Index (relative to end)');
    title('Estimated phase drift of X-pol (last 500 samples)'); grid on;

    subplot(2, 1, 2);
    plot(-499:0, unwrap(thetaY(end-499:end)));  % Show last 500 points
    xlabel('Sample Index (relative to end)');
    title('Estimated phase drift of Y-pol (last 500 samples)'); grid on;
end

%% 末尾稳定部分星座图
num_show = min(131072, length(CMAdataX) - 1);
if num_show > 0
    CMAdataX_tail = CMAdataX(end-num_show:end);
    CMAdataY_tail = CMAdataY(end-num_show:end);

    figure();
    subplot(1, 2, 1);
    plot(CMAdataX_tail, '.r');
    title('Last part of X-pol after CMA'); grid on; axis equal;
    xlim([-2 2]); ylim([-2 2]);

    subplot(1, 2, 2);
    plot(CMAdataY_tail, '.b');
    title('Last part of Y-pol after CMA'); grid on; axis equal;
    xlim([-2 2]); ylim([-2 2]);
end

fprintf('\nProcessing of TD_TRdata with time-domain 2x2 CMA is complete.\n');
