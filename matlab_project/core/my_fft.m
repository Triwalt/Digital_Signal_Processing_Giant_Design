function X = my_fft(x)
% my_fft: 实现基于迭代蝶形运算的基2-FFT算法 (按时间抽取 DIT)
%
% 优化要点:
% 1. 迭代实现, 适合硬件映射
% 2. 缓存比特反转索引与旋转因子, 避免重复计算
% 3. 利用向量化蝶形运算, 减少MATLAB循环开销
%
% 输入:
%   x - 复数列向量, 长度N必须是2的整数次幂
%
% 输出:
%   X - x的N点DFT结果, 复数列向量

    % 确保输入为列向量
    if size(x, 2) > 1
        x = x(:);
    end

    N = length(x);

    % 验证N是2的幂
    if N == 0 || mod(log2(N), 1) ~= 0
        error('输入长度必须是2的整数次幂');
    end

    % 获取/更新缓存
    cache = get_fft_cache(N);

    % 步骤1: 比特反转重排
    X = x(cache.bitrev_idx);

    % 步骤2: 迭代蝶形运算
    for stage = 1:cache.num_stages
        span = cache.stage_span(stage);
        half_span = span / 2;
        W = cache.stage_twiddles{stage};

        X = reshape(X, span, []);

        top_val = X(1:half_span, :);
        bot_val = X(half_span + 1:end, :);

        twiddle_mul = W .* bot_val;
        X(1:half_span, :) = top_val + twiddle_mul;
        X(half_span + 1:end, :) = top_val - twiddle_mul;

        X = X(:);
    end
end


function cache = get_fft_cache(N)
% get_fft_cache: 返回给定N的FFT缓存(比特反转索引和旋转因子)

    persistent cached_N cached_data

    if isempty(cached_N)
        cached_N = [];
        cached_data = {};
    end

    idx = find(cached_N == N, 1);
    if isempty(idx)
        num_stages = log2(N);

        bitrev_idx = compute_bitrev_indices(N);
        twiddle_full = exp(-1j * 2 * pi * (0:(N/2 - 1)) / N).';

        stage_span = zeros(num_stages, 1);
        stage_twiddles = cell(num_stages, 1);

        for stage = 1:num_stages
            span = 2^stage;
            half_span = span / 2;
            stride = N / span;
            idx_vec = (0:half_span-1) * stride + 1;
            stage_span(stage) = span;
            stage_twiddles{stage} = twiddle_full(idx_vec);
        end

        cache = struct( ...
            'N', N, ...
            'num_stages', num_stages, ...
            'bitrev_idx', bitrev_idx, ...
            'stage_span', stage_span, ...
            'stage_twiddles', {stage_twiddles} ...
        );

        cached_N(end + 1) = N;
        cached_data{end + 1} = cache;
    else
        cache = cached_data{idx};
    end
end


function idx = compute_bitrev_indices(N)
% compute_bitrev_indices: 生成长度为N的比特反转索引

    num_bits = log2(N);
    values = (0:N-1).';
    reversed = zeros(N, 1);

    for bit = 1:num_bits
        reversed = bitshift(reversed, 1) + bitand(values, 1);
        values = bitshift(values, -1);
    end

    idx = reversed + 1;
end
