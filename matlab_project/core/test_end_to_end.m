% 简化的端到端测试
clear; clc;

fprintf('简化端到端测试\n');
fprintf('================\n\n');

% 配置
numSymbols = 100;
grayMap = [0 1 3 2];
constellation_map = [1+1j, 1-1j, -1-1j, -1+1j] / sqrt(2);

% 发射端
txBits = randi([0 1], 1, numSymbols * 2);

txSymbolIndices = zeros(1, numSymbols);
for i = 1:numSymbols
    bit1 = txBits(2*i-1);
    bit2 = txBits(2*i);
    txSymbolIndices(i) = bit1 * 2 + bit2;
end

grayIndices = grayMap(txSymbolIndices + 1);
txSymbols = constellation_map(grayIndices + 1);

fprintf('发射: %d 符号\n', numSymbols);

% 接收端(无噪声,直接硬判决)
rxSymbolIndices = zeros(length(txSymbols), 1);
for k = 1:length(txSymbols)
    [~, rxSymbolIndices(k)] = min(abs(txSymbols(k) - constellation_map));
end
rxGrayIndices = rxSymbolIndices - 1;

fprintf('接收格雷码示例: [%d %d %d %d %d]\n', rxGrayIndices(1:5));
fprintf('发射格雷码示例: [%d %d %d %d %d]\n', grayIndices(1:5));

% 方法1: ismember
[~, rxSymbolIndices_original1] = ismember(rxGrayIndices, grayMap);
rxSymbolIndices_original1 = rxSymbolIndices_original1 - 1;

rxBits1 = zeros(1, length(rxSymbolIndices_original1) * 2);
for i = 1:length(rxSymbolIndices_original1)
    idx = rxSymbolIndices_original1(i);
    rxBits1(2*i-1) = bitshift(idx, -1);
    rxBits1(2*i) = bitand(idx, 1);
end

bit_errors1 = sum(txBits ~= rxBits1);
ber1 = bit_errors1 / length(txBits);
fprintf('\n方法1 (ismember): BER = %.4f (%d / %d)\n', ber1, bit_errors1, length(txBits));

% 方法2: find
rxSymbolIndices_original2 = zeros(size(rxGrayIndices));
for k = 1:length(rxGrayIndices)
    pos = find(grayMap == rxGrayIndices(k));
    if ~isempty(pos)
        rxSymbolIndices_original2(k) = pos - 1;
    end
end

rxBits2 = zeros(1, length(rxSymbolIndices_original2) * 2);
for i = 1:length(rxSymbolIndices_original2)
    idx = rxSymbolIndices_original2(i);
    rxBits2(2*i-1) = bitshift(idx, -1);
    rxBits2(2*i) = bitand(idx, 1);
end

bit_errors2 = sum(txBits ~= rxBits2);
ber2 = bit_errors2 / length(txBits);
fprintf('方法2 (find): BER = %.4f (%d / %d)\n', ber2, bit_errors2, length(txBits));

% 检查原始索引
fprintf('\n发射索引示例: [%d %d %d %d %d]\n', txSymbolIndices(1:5));
fprintf('接收索引1示例: [%d %d %d %d %d]\n', rxSymbolIndices_original1(1:5)');
fprintf('接收索引2示例: [%d %d %d %d %d]\n', rxSymbolIndices_original2(1:5)');
