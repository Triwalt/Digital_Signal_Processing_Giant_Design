function y = fast_conv_os(x_long, h, nfft)
% fast_conv_os: 使用重叠保留法(Overlap-Save)实现快速卷积
%
% 算法原理:
% 1. 将长输入信号分成重叠的块
% 2. 每块与滤波器在频域相乘
% 3. 丢弃每块的前M-1个样本(混叠部分)
% 4. 保留后L=nfft-M+1个有效样本
%
% 输入:
%   x_long - 长输入信号序列(列向量)
%   h      - FIR滤波器核(列向量)
%   nfft   - FFT点数,必须是2的整数次幂
%
% 输出:
%   y - x_long与h线性卷积的结果
%
% 参考: 技术规格书第3.1节

    % 确保列向量
    if size(x_long, 2) > 1
        x_long = x_long(:);
    end
    if size(h, 2) > 1
        h = h(:);
    end
    
    M = length(h);
    Lx = length(x_long);
    
    % 参数验证
    if mod(log2(nfft), 1) ~= 0
        error('nfft必须是2的整数次幂');
    end
    
    if nfft < M
        error('nfft必须大于等于滤波器长度M');
    end
    
    % 计算有效数据块长度
    L = nfft - M + 1;
    
    % 滤波器频域表示
    h_padded = [h; zeros(nfft - M, 1)];
    H = my_fft(h_padded);
    
    % 在输入前端补M-1个零
    x_padded = [zeros(M-1, 1); x_long];
    Lx_padded = length(x_padded);
    
    % 计算需要的块数
    num_blocks = ceil(Lx_padded / L);
    
    % 初始化输出
    y = [];
    
    % 处理每个块
    for i = 1:num_blocks
        % 提取当前块(长度nfft,相邻块重叠M-1个样本)
        block_start = (i-1) * L + 1;
        block_end = min(block_start + nfft - 1, Lx_padded);
        
        % 提取并补零到nfft长度
        current_block = x_padded(block_start:block_end);
        if length(current_block) < nfft
            current_block = [current_block; zeros(nfft - length(current_block), 1)];
        end
        
        % FFT
        X_block = my_fft(current_block);
        
        % 频域相乘
        Y_block = X_block .* H;
        
        % IFFT
        y_block = my_ifft(Y_block);
        
        % 保留有效部分(丢弃前M-1个样本)
        valid_samples = y_block(M:end);
        
        % 拼接到输出
        y = [y; valid_samples];
    end
    
    % 调整输出长度为标准卷积长度
    final_length = Lx + M - 1;
    y = y(1:final_length);
end
