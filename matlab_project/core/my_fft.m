function X = my_fft(x)
% my_fft: 实现基于迭代蝶形运算的基2-FFT算法 (按时间抽取 DIT)
%
% 算法特点:
% 1. 迭代实现(非递归),适合硬件实现
% 2. 比特反转重排输入
% 3. 旋转因子预计算并利用对称性优化
% 4. 原位蝶形运算
%
% 输入:
%   x - 复数列向量,长度N必须是2的整数次幂
%
% 输出:
%   X - x的N点DFT结果,复数列向量
%
% 参考: 技术规格书第2.1节

    % 确保输入为列向量
    if size(x, 2) > 1
        x = x(:);
    end
    
    N = length(x);
    
    % 验证N是2的幂
    if N == 0 || mod(log2(N), 1) ~= 0
        error('输入长度必须是2的整数次幂');
    end
    
    % 计算级数
    num_stages = log2(N);
    
    % 步骤1: 比特反转重排
    X = bit_reverse(x, N);
    
    % 步骤2: 预计算旋转因子 (利用对称性,仅计算前N/2个)
    % W_N^k = exp(-j*2*pi*k/N), k = 0, 1, ..., N/2-1
    twiddle_factors = exp(-1j * 2 * pi * (0:N/2-1) / N);
    
    % 步骤3: 迭代蝶形运算 (log2(N)级)
    for stage = 1:num_stages
        % 当前级的参数
        butterfly_span = 2^stage;          % 蝶形运算跨度
        num_butterflies = N / butterfly_span; % 每组蝶形数量
        half_span = butterfly_span / 2;    % 半跨度
        
        % 对每一组进行蝶形运算
        for group = 0:(num_butterflies - 1)
            % 组的起始索引
            group_start = group * butterfly_span + 1;
            
            % 对组内每个蝶形单元进行计算
            for butterfly = 0:(half_span - 1)
                % 计算蝶形的上下节点索引
                idx_top = group_start + butterfly;
                idx_bot = idx_top + half_span;
                
                % 计算旋转因子索引 (利用对称性和周期性)
                % 旋转因子为 W_N^(butterfly * N / butterfly_span)
                twiddle_idx = mod(butterfly * N / butterfly_span, N);
                
                % 利用对称性获取旋转因子
                if twiddle_idx < N/2
                    W = twiddle_factors(twiddle_idx + 1);
                else
                    % W_N^(k+N/2) = -W_N^k
                    W = -twiddle_factors(twiddle_idx - N/2 + 1);
                end
                
                % 蝶形运算
                temp = X(idx_top);
                X(idx_top) = temp + W * X(idx_bot);
                X(idx_bot) = temp - W * X(idx_bot);
            end
        end
    end
end


function y = bit_reverse(x, N)
% bit_reverse: 对输入向量按照比特反转顺序重新排列
%
% 输入:
%   x - 输入向量
%   N - 向量长度 (必须是2的幂)
%
% 输出:
%   y - 比特反转后的向量

    num_bits = log2(N);
    y = zeros(size(x));
    
    for i = 0:(N-1)
        % 计算i的比特反转索引
        reversed_idx = 0;
        temp = i;
        
        for bit = 0:(num_bits-1)
            reversed_idx = bitshift(reversed_idx, 1);
            reversed_idx = reversed_idx + mod(temp, 2);
            temp = bitshift(temp, -1);
        end
        
        % 重新排列 (MATLAB索引从1开始)
        y(reversed_idx + 1) = x(i + 1);
    end
end
