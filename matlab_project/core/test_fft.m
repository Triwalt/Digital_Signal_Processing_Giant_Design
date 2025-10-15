% test_fft.m: FFT和IFFT函数的验证脚本
%
% 功能:
% 1. 验证my_fft与MATLAB内置fft的一致性
% 2. 验证my_ifft的正确性
% 3. 测试多种信号类型和长度
%
% 参考: 技术规格书第2.3节

clear; clc;
fprintf('========================================\n');
fprintf('FFT和IFFT验证测试\n');
fprintf('========================================\n\n');

% 测试配置
test_lengths = [64, 256, 1024];
tolerance = 1e-12;
all_passed = true;

for N = test_lengths
    fprintf('--- 测试长度 N = %d ---\n', N);
    
    % 测试1: 单位冲激
    fprintf('  测试1: 单位冲激... ');
    x_impulse = zeros(N, 1);
    x_impulse(1) = 1;
    
    Y_custom = my_fft(x_impulse);
    Y_builtin = fft(x_impulse);
    max_error = max(abs(Y_custom - Y_builtin));
    
    if max_error < tolerance
        fprintf('通过 (误差: %.2e)\n', max_error);
    else
        fprintf('失败 (误差: %.2e)\n', max_error);
        all_passed = false;
    end
    
    % 测试2: 直流信号
    fprintf('  测试2: 直流信号... ');
    x_dc = ones(N, 1);
    
    Y_custom = my_fft(x_dc);
    Y_builtin = fft(x_dc);
    max_error = max(abs(Y_custom - Y_builtin));
    
    if max_error < tolerance
        fprintf('通过 (误差: %.2e)\n', max_error);
    else
        fprintf('失败 (误差: %.2e)\n', max_error);
        all_passed = false;
    end
    
    % 测试3: 单频正弦波
    fprintf('  测试3: 单频正弦波... ');
    k = 5; % 频率bin
    n = (0:N-1)';
    x_sine = exp(1j * 2 * pi * k * n / N);
    
    Y_custom = my_fft(x_sine);
    Y_builtin = fft(x_sine);
    max_error = max(abs(Y_custom - Y_builtin));
    
    if max_error < tolerance
        fprintf('通过 (误差: %.2e)\n', max_error);
    else
        fprintf('失败 (误差: %.2e)\n', max_error);
        all_passed = false;
    end
    
    % 测试4: 高斯白噪声
    fprintf('  测试4: 高斯白噪声... ');
    rng(42); % 固定随机种子
    x_noise = randn(N, 1) + 1j * randn(N, 1);
    
    Y_custom = my_fft(x_noise);
    Y_builtin = fft(x_noise);
    max_error = max(abs(Y_custom - Y_builtin));
    
    if max_error < tolerance
        fprintf('通过 (误差: %.2e)\n', max_error);
    else
        fprintf('失败 (误差: %.2e)\n', max_error);
        all_passed = false;
    end
    
    % 测试5: IFFT验证(往返变换)
    fprintf('  测试5: IFFT往返变换... ');
    x_test = randn(N, 1) + 1j * randn(N, 1);
    
    Y_fft = my_fft(x_test);
    x_reconstructed = my_ifft(Y_fft);
    max_error = max(abs(x_reconstructed - x_test));
    
    if max_error < tolerance
        fprintf('通过 (误差: %.2e)\n', max_error);
    else
        fprintf('失败 (误差: %.2e)\n', max_error);
        all_passed = false;
    end
    
    fprintf('\n');
end

% 总结
fprintf('========================================\n');
if all_passed
    fprintf('所有测试通过! ✓\n');
else
    fprintf('部分测试失败! ✗\n');
end
fprintf('========================================\n');
