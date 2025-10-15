% test_fast_conv.m: 快速卷积函数的验证脚本
%
% 功能:
% 验证fast_conv_os与MATLAB内置conv的一致性
%
% 参考: 技术规格书第3.2节

clear; clc;
fprintf('========================================\n');
fprintf('快速卷积验证测试\n');
fprintf('========================================\n\n');

% 测试配置
tolerance = 1e-10;
all_passed = true;

% 测试1: 长信号与短滤波器
fprintf('测试1: 长信号与短滤波器\n');
rng(42);
x_long = randn(10000, 1);
h = randn(31, 1);

fprintf('  计算基准结果 (conv)... ');
tic;
y_ref = conv(x_long, h);
t_conv = toc;
fprintf('完成 (%.3f秒)\n', t_conv);

fprintf('  计算测试结果 (fast_conv_os)... ');
tic;
y_test = fast_conv_os(x_long, h, 512);
t_fast = toc;
fprintf('完成 (%.3f秒)\n', t_fast);

fprintf('  加速比: %.2fx\n', t_conv / t_fast);

% 验证长度
if length(y_test) ~= length(y_ref)
    fprintf('  长度不匹配: 期望%d, 得到%d\n', length(y_ref), length(y_test));
    all_passed = false;
else
    fprintf('  长度匹配: %d ✓\n', length(y_ref));
end

% 验证数值精度
max_error = max(abs(y_ref - y_test));
fprintf('  最大误差: %.2e\n', max_error);

if max_error < tolerance
    fprintf('  结果验证: 通过 ✓\n');
else
    fprintf('  结果验证: 失败 ✗\n');
    all_passed = false;
end

fprintf('\n');

% 测试2: 不同FFT大小
fprintf('测试2: 不同FFT大小的影响\n');
x_test = randn(5000, 1);
h_test = randn(25, 1);
y_ref = conv(x_test, h_test);

nfft_values = [128, 256, 512, 1024];
fprintf('  FFT大小\t时间(秒)\t最大误差\n');
fprintf('  ----------------------------------------\n');

for nfft = nfft_values
    tic;
    y_test = fast_conv_os(x_test, h_test, nfft);
    t_elapsed = toc;
    max_error = max(abs(y_ref - y_test));
    fprintf('  %4d\t\t%.4f\t\t%.2e\n', nfft, t_elapsed, max_error);
end

fprintf('\n');

% 测试3: 边界情况
fprintf('测试3: 边界情况\n');

% 3a: 滤波器长度为1
fprintf('  3a: 滤波器长度为1... ');
x_test = randn(100, 1);
h_test = randn(1, 1);
y_ref = conv(x_test, h_test);
y_test = fast_conv_os(x_test, h_test, 64);
max_error = max(abs(y_ref - y_test));
if max_error < tolerance
    fprintf('通过 (误差: %.2e)\n', max_error);
else
    fprintf('失败 (误差: %.2e)\n', max_error);
    all_passed = false;
end

% 3b: 输入长度小于FFT大小
fprintf('  3b: 输入长度小于FFT大小... ');
x_test = randn(50, 1);
h_test = randn(10, 1);
y_ref = conv(x_test, h_test);
y_test = fast_conv_os(x_test, h_test, 128);
max_error = max(abs(y_ref - y_test));
if max_error < tolerance
    fprintf('通过 (误差: %.2e)\n', max_error);
else
    fprintf('失败 (误差: %.2e)\n', max_error);
    all_passed = false;
end

fprintf('\n');

% 总结
fprintf('========================================\n');
if all_passed
    fprintf('所有测试通过! ✓\n');
else
    fprintf('部分测试失败! ✗\n');
end
fprintf('========================================\n');
