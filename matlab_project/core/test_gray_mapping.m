% 测试格雷码映射
clear; clc;

grayMap = [0 1 3 2];
constellation_map = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);

fprintf('格雷码映射测试\n');
fprintf('================\n\n');

% 测试所有4种比特组合
for b1 = 0:1
    for b2 = 0:1
        % 发射端
        bits = [b1, b2];
        decimal_idx = b1 * 2 + b2; % 0-3
        gray_idx = grayMap(decimal_idx + 1); % 格雷码
        tx_symbol = constellation_map(gray_idx + 1);
        
        fprintf('比特 [%d,%d]:\n', b1, b2);
        fprintf('  十进制索引: %d\n', decimal_idx);
        fprintf('  格雷码: %d\n', gray_idx);
        fprintf('  星座点索引(数组): %d\n', gray_idx + 1);
        fprintf('  星座点: (%.2f%+.2fj)\n', real(tx_symbol), imag(tx_symbol));
        
        % 接收端(硬判决)
        [~, rx_constellation_idx] = min(abs(tx_symbol - constellation_map));
        rx_gray_idx = rx_constellation_idx - 1;
        
        fprintf('  接收星座点索引: %d\n', rx_constellation_idx);
        fprintf('  接收格雷码: %d\n', rx_gray_idx);
        
        % 方法1: ismember (旧)
        [~, pos1] = ismember(rx_gray_idx, grayMap);
        rx_decimal_idx1 = pos1 - 1;
        rx_b1_1 = bitshift(rx_decimal_idx1, -1);
        rx_b2_1 = bitand(rx_decimal_idx1, 1);
        
        fprintf('  方法1(ismember): 十进制索引=%d, 比特=[%d,%d]', ...
                rx_decimal_idx1, rx_b1_1, rx_b2_1);
        if rx_b1_1 == b1 && rx_b2_1 == b2
            fprintf(' ✓\n');
        else
            fprintf(' ✗\n');
        end
        
        % 方法2: find (新)
        pos2 = find(grayMap == rx_gray_idx);
        rx_decimal_idx2 = pos2 - 1;
        rx_b1_2 = bitshift(rx_decimal_idx2, -1);
        rx_b2_2 = bitand(rx_decimal_idx2, 1);
        
        fprintf('  方法2(find): 十进制索引=%d, 比特=[%d,%d]', ...
                rx_decimal_idx2, rx_b1_2, rx_b2_2);
        if rx_b1_2 == b1 && rx_b2_2 == b2
            fprintf(' ✓\n');
        else
            fprintf(' ✗\n');
        end
        
        fprintf('\n');
    end
end
