function y = my_awgn(x, snr_db, power_mode)
% my_awgn: 向信号添加高斯白噪声以达到指定SNR
%
% 输入:
%   x         - 输入信号
%   snr_db    - 信噪比 (dB)
%   power_mode - 'measured' 测量信号功率, 或数值指定功率
%
% 输出:
%   y         - 添加噪声后的信号
%
% 算法:
% 1. 根据power_mode计算信号功率
% 2. 根据SNR要求确定噪声功率
% 3. 生成复高斯白噪声
% 4. 缩放噪声以达到目标SNR
% 5. 将噪声添加到信号

    % 确保输入为行向量以保持一致性
    if size(x, 1) > 1
        x = x.';
        transpose_output = true;
    else
        transpose_output = false;
    end
    
    % 计算信号功率
    if strcmp(power_mode, 'measured')
        % 测量实际信号功率
        signal_power = mean(abs(x).^2);
    else
        % 使用指定的信号功率 (如未指定则假设归一化为1)
        signal_power = 1;
    end
    
    % 将SNR从dB转换为线性值
    snr_linear = 10^(snr_db / 10);
    
    % 计算所需噪声功率
    noise_power = signal_power / snr_linear;
    
    % 生成复高斯白噪声
    % 对于复信号: 噪声功率在实部和虚部之间平均分配
    if ~isreal(x)
        % 复噪声: 每个分量的方差 = noise_power/2
        noise_real = sqrt(noise_power/2) * randn(size(x));
        noise_imag = sqrt(noise_power/2) * randn(size(x));
        noise = noise_real + 1j * noise_imag;
    else
        % 实噪声: 方差 = noise_power
        noise = sqrt(noise_power) * randn(size(x));
    end
    
    % 将噪声添加到信号
    y = x + noise;
    
    % 如需要则恢复原始方向
    if transpose_output
        y = y.';
    end
end
