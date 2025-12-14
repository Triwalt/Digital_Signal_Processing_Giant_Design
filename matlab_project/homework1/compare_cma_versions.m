%% compare_cma_versions.m
% 比较时域CMA (Reference) 和频域CMA (Homework) 的输出结果
% 特别关注尾部误差

clc;
clear;
close all;

%% 1. 运行 Reference 版本 (时域)
fprintf('正在运行 Reference 版本 (时域)...\n');
run('CMA_reference.m');

% 保存 Reference 结果
ref_CMAdataX = CMAdataX;
ref_CMAdataY = CMAdataY;
ref_SigX_out = SigX_out;
ref_SigY_out = SigY_out;

% 清除可能冲突的变量，保留 Reference 结果
% 注意：CMA_homework.m 会重新加载数据，所以我们可以清除大部分变量
% 但要保留 ref_* 变量
keep_vars = {'ref_CMAdataX', 'ref_CMAdataY', 'ref_SigX_out', 'ref_SigY_out'};
all_vars = who;
vars_to_clear = setdiff(all_vars, keep_vars);
clear(vars_to_clear{:});

%% 2. 运行 Homework 版本 (频域)
fprintf('正在运行 Homework 版本 (频域)...\n');

% 设置配置开关
enable_vv_phase_lock = false; % 关闭相位校正以进行公平比较
enable_plot_output = false;   % 关闭绘图以加快速度

run('CMA_homework.m');

% 获取 Homework 结果
hw_CMAdataX = CMAdataX;
hw_CMAdataY = CMAdataY;
hw_SigX_out = SigX_out;
hw_SigY_out = SigY_out;

%% 3. 比较结果
fprintf('正在比较结果...\n');

% 确保长度一致
% 注意: CMA_reference.m 可能会截断 CMAdataX，所以我们使用 SigX_out (全长)
len = min(length(ref_SigX_out), length(hw_SigX_out));
ref_data_X = ref_SigX_out(1:len);
ref_data_Y = ref_SigY_out(1:len);
hw_data_X = hw_SigX_out(1:len);
hw_data_Y = hw_SigY_out(1:len);

% 计算误差
err_X = abs(ref_data_X - hw_data_X);
err_Y = abs(ref_data_Y - hw_data_Y);

% 整体误差统计
max_err_X = max(err_X);
max_err_Y = max(err_Y);
mse_X = mean(err_X.^2);
mse_Y = mean(err_Y.^2);

fprintf('整体误差统计:\n');
fprintf('Max Error X: %e\n', max_err_X);
fprintf('Max Error Y: %e\n', max_err_Y);
fprintf('MSE X: %e\n', mse_X);
fprintf('MSE Y: %e\n', mse_Y);

%% 4. 尾部误差比较
% 比较最后 1000 个点 (或者更多，取决于数据长度)
tail_len = 1000;
if len < tail_len
    tail_len = len;
end

tail_err_X = err_X(end-tail_len+1:end);
tail_err_Y = err_Y(end-tail_len+1:end);

max_tail_err_X = max(tail_err_X);
max_tail_err_Y = max(tail_err_Y);
mse_tail_X = mean(tail_err_X.^2);
mse_tail_Y = mean(tail_err_Y.^2);

fprintf('\n尾部 (%d 点) 误差统计:\n', tail_len);
fprintf('Max Tail Error X: %e\n', max_tail_err_X);
fprintf('Max Tail Error Y: %e\n', max_tail_err_Y);
fprintf('MSE Tail X: %e\n', mse_tail_X);
fprintf('MSE Tail Y: %e\n', mse_tail_Y);

%% 5. 绘图
figure('Name', 'CMA Comparison Error');
subplot(2,1,1);
plot(err_X);
title('Error Magnitude X (Ref vs Homework)');
xlabel('Sample Index');
ylabel('Error');
grid on;

subplot(2,1,2);
plot(err_Y);
title('Error Magnitude Y (Ref vs Homework)');
xlabel('Sample Index');
ylabel('Error');
grid on;

% 绘制尾部对比
figure('Name', 'CMA Tail Comparison');
subplot(2,1,1);
plot(real(ref_data_X(end-100:end)), 'b.-'); hold on;
plot(real(hw_data_X(end-100:end)), 'r--');
title('Tail Real Part X (Last 100 samples)');
legend('Reference', 'Homework');
grid on;

subplot(2,1,2);
plot(imag(ref_data_X(end-100:end)), 'b.-'); hold on;
plot(imag(hw_data_X(end-100:end)), 'r--');
title('Tail Imag Part X (Last 100 samples)');
legend('Reference', 'Homework');
grid on;

fprintf('\n比较完成。\n');
