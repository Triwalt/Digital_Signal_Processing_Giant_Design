function x = my_ifft_mix(X)
% my_ifft_mix  与 my_fft_mix 对应的任意长度一维逆FFT实现
%   x = my_ifft_mix(X) 计算频域序列 X 的逆离散傅里叶变换，得到时域序列 x。
%   输入长度 N 可以是任意正整数，内部调用 my_fft_mix，并利用
%   ifft 的共轭关系：ifft(X) = conj(fft(conj(X))) / N。
%
%   本实现与 my_fft_mix 使用相同的 Planner / 算法组合：
%   - 合数：混合基 Cooley-Tukey (mixed-radix)
%   - 小质数：占位 codelet（当前为直接DFT，将来可替换为C/MEX）
%   - 大质数：Bluestein 算法
%
%   输入可以是行向量或列向量；输出将保持与输入相同的方向。
%   本函数支持复数值输入。

    % 确保输入是向量
    if ~isvector(X)
        error('my_ifft_mix:InputNotVector', '输入X必须是一维向量');
    end

    % 记住原始方向（行向量或列向量）
    wasColumn = iscolumn(X);

    % 内部统一使用行向量处理
    X = X(:).';

    % 利用共轭关系: ifft(X) = conj(fft(conj(X))) / N
    N = length(X);
    X_conj = conj(X);
    x_tmp = my_fft_mix(X_conj);  % 这里复用 my_fft_mix 的 Planner
    x_row = conj(x_tmp) / N;

    % 恢复原始方向
    if wasColumn
        x = x_row.';
    else
        x = x_row;
    end
end
