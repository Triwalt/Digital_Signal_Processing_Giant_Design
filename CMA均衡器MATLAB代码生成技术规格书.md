

# **CMA均衡器MATLAB代码生成技术规格书**

**文档目标：** 本文档为AI编程助手提供一份完整的技术规格说明，旨在引导其逐步、精确地生成一个基于MATLAB的4QAM通信系统仿真程序。该程序的核心是恒模盲均衡算法（CMA），并且所有关键计算模块（如FFT、快速卷积）均需从零开始实现。请严格遵循本文档定义的模块、函数签名、算法步骤和验证协议。

## **第1章：理论基础与关键参数定义**

本章提供项目所需的理论背景，并明确定义将在代码中使用的关键常量和参数。

### **1.1 理论概述**

1. **信道失真与符号间干扰 (ISI)**：数字信号在物理信道（如光纤）中传输时，会因色散等效应产生失真，导致脉冲展宽，形成符号间干扰（ISI），从而增加误码率（BER）1。  
2. **信道均衡**：为对抗ISI，接收端使用均衡器（一种可调滤波器）来补偿信道失真 3。  
3. **盲均衡与CMA**：盲均衡技术无需训练序列，仅利用接收信号的统计特性进行自适应均衡，提高了传输效率 5。恒模算法（CMA）是其中最经典的一种，适用于具有恒定包络特性的调制信号（如4QAM）6。  
4. CMA核心思想：CMA通过最小化其输出信号模值的波动来恢复信号。其代价函数定义为：

   其中 x\~\[n\]=wH\[n\]y\[n\] 是均衡器输出，w\[n\] 是滤波器系数向量。  
5. CMA权重更新法则：采用随机梯度下降法，得到逐样本的权重更新公式：

   注意：为与作业要求一致，实际编程时可采用等价形式：

6. **相位模糊性**：CMA代价函数仅与模值有关，与相位无关。因此，算法收敛后，输出星座图可能存在一个固有的相位旋转（对于4QAM，可能为 ）8。这必须在后续通过专门的载波恢复或相位校正环节来解决。

### **1.2 实现关键参数定义**

* **调制方式**：4QAM。  
* 收敛半径 R2​：对于归一化能量的4QAM星座点（模值为1），R2​ 的计算公式为 R2​=E\[∣ak​∣4\]/E\[∣ak​∣2\]。由于 ∣ak​∣2=1，因此：

## **第2章：模块一：自定义FFT与IFFT函数的实现**

**任务**：创建两个MATLAB函数，my\_fft.m 和 my\_ifft.m，用于实现快速傅里叶变换及其逆变换。严禁调用MATLAB内置的fft和ifft函数。

### **2.1 函数 my\_fft.m**

* **函数签名**：function X \= my\_fft(x)  
* **输入**：x \- 复数列向量，其长度  必须是2的整数次幂。  
* **输出**：X \- x的N点DFT结果，复数列向量。  
* **实现步骤**：  
  1. **输入验证**：检查输入x是否为列向量，其长度N是否为2的幂。  
  2. **比特反转置换**：创建一个辅助函数 bit\_reverse(x, N)，对输入向量x的元素按其索引的比特反转形式进行重新排序 10。  
  3. **旋转因子预计算（优化）**：  
     * 根据作业要求，必须手动计算并优化旋转因子 12。  
     * 利用对称性 ，仅需预计算并存储前  个旋转因子：twiddle\_factors \= exp(-1j \* 2 \* pi \* (0:N/2-1) / N);  
  4. **迭代蝶形运算**：  
     * 采用迭代方式实现按时间抽取（DIT）的基2-FFT算法。  
     * 循环 log2​N 级。在每一级中，根据蝶形运算公式更新数据：

     * 在运算中，通过索引映射从预计算的twiddle\_factors表中高效获取所需的旋转因子。

### **2.2 函数 my\_ifft.m**

* **函数签名**：function x \= my\_ifft(X)  
* 实现：可以通过复用my\_fft函数高效实现。根据IDFT的定义：

  因此，my\_ifft的实现逻辑为：  
  1. 获取输入X的长度N。  
  2. 计算X的共轭conj(X)。  
  3. 调用my\_fft(conj(X))。  
  4. 将结果取共轭，并除以N。

### **2.3 验证协议 (test\_fft.m)**

**任务**：创建一个名为test\_fft.m的脚本，用于严格验证my\_fft和my\_ifft的正确性。

1. **生成测试向量**：创建多个不同长度（如64, 256, 1024）的测试信号，包括单位冲激、直流信号、单频正弦波和高斯白噪声。  
2. **FFT验证**：  
   * 对每个测试向量x\_test，分别计算 Y\_custom \= my\_fft(x\_test) 和 Y\_builtin \= fft(x\_test)。  
   * 计算最大绝对误差 max\_error \= max(abs(Y\_custom \- Y\_builtin))。  
   * 使用assert(max\_error \< 1e-12)来断言自定义函数与内置函数结果的一致性。  
3. **IFFT验证**：  
   * 对Y\_builtin，计算 x\_reconstructed \= my\_ifft(Y\_builtin)。  
   * 计算最大绝对误差 max\_error \= max(abs(x\_reconstructed \- x\_test))。  
   * 使用assert(max\_error \< 1e-12)来断言逆变换的正确性。

## **第3章：模块二：基于快速卷积的滤波实现**

**任务**：使用已实现的my\_fft和my\_ifft函数，创建一个通过\*\*重叠保留法（Overlap-Save）\*\*实现长序列滤波的函数。

### **3.1 函数 fast\_conv\_os.m**

* **函数签名**：function y \= fast\_conv\_os(x\_long, h, nfft)  
* **输入**：  
  * x\_long：长输入信号序列（列向量）。  
  * h：短的FIR滤波器核（列向量）。  
  * nfft：FFT点数，必须是2的整数次幂。  
* **输出**：y \- x\_long与h线性卷积的结果。  
* **实现步骤**：  
  1. **参数计算与验证**：  
     * 获取滤波器长度 M \= length(h) 和输入信号长度 Lx \= length(x\_long)。  
     * 验证 nfft 是否为2的幂，且 nfft \>= M。  
     * 计算有效数据块长度 L \= nfft \- M \+ 1。  
  2. **滤波器频域表示**：  
     * 将h补零至nfft长度。  
     * 计算一次滤波器的FFT：H \= my\_fft(h\_padded)。  
  3. **分块与卷积**：  
     * 初始化一个输入缓冲区，并对x\_long前端补M-1个零。  
     * 循环处理x\_long：  
       a. 从x\_long中取出nfft长度的重叠数据块（相邻块重叠M-1个样本）14。  
       b. 计算该块的FFT：X\_block \= my\_fft(current\_block)。  
       c. 在频域相乘：Y\_block \= X\_block.\* H。  
       d. 反变换回时域：y\_block \= my\_ifft(Y\_block)。  
       e. 保留有效部分：根据重叠保留法，丢弃y\_block的前M-1个无效样本，保留后L个有效样本 14。  
       f. 将有效样本拼接到最终的输出向量y中。  
  4. **长度调整**：确保最终输出y的长度与conv(x\_long, h)的结果长度一致（即Lx \+ M \- 1）。

### **3.2 验证协议 (test\_fast\_conv.m)**

**任务**：创建一个名为test\_fast\_conv.m的脚本，验证fast\_conv\_os的正确性。

1. **生成测试数据**：创建长随机输入x\_long（如长度10000）和短随机滤波器h（如长度31）。  
2. **计算基准**：使用MATLAB内置函数计算精确结果：y\_ref \= conv(x\_long, h)。  
3. **计算待测结果**：调用y\_test \= fast\_conv\_os(x\_long, h, 512)（其中512是一个合适的nfft值）。  
4. **比较与断言**：  
   * 截取y\_test使其与y\_ref长度一致。  
   * 计算最大绝对误差 max\_error \= max(abs(y\_ref \- y\_test))。  
   * 使用assert(max\_error \< 1e-12)进行断言。

## **第4章：模块三：主仿真脚本的构建**

**任务**：创建主脚本main\_cma\_simulation.m，集成所有模块，搭建并运行一个完整的4QAM通信链路仿真。

### **4.1 步骤1：仿真参数配置**

在脚本开头定义一个config结构体，集中管理所有参数：

Matlab

config.numSymbols \= 200000;      % 仿真符号数量  
config.snr\_dB \= 25;               % 信噪比 (dB)  
% 4QAM调制参数  
config.M \= 4;  
config.grayMap \= \[0 1 3 2\];       % 格雷码映射表  
% 信道参数  
config.channelTaps \= \[1, 0.5, 0.2\]; % 示例ISI信道  
% CMA均衡器参数  
config.cma.filterLength \= 31;     % 均衡器长度 (M)  
config.cma.stepSize \= 0.001;      % 步长 (mu)  
config.cma.R2 \= 1;                % 4QAM的收敛半径  
% 快速卷积参数  
config.fft.nfft \= 512;            % FFT点数 (必须是2的幂)

### **4.2 步骤2：发射端实现**

1. **信源**：生成config.numSymbols \* log2(config.M)个随机比特。  
2. **符号映射**：将比特流两位一组，转换为十进制索引。  
3. **格雷编码**：使用config.grayMap将索引映射为格雷码顺序。  
4. **4QAM调制**：使用qammod函数（或手动实现）将格雷码索引映射为复数星座点。确保输出符号能量归一化。

### **4.3 步骤3：信道建模**

1. **ISI信道**：使用filter(config.channelTaps, 1, txSymbols)使发射信号通过信道。  
2. **AWGN**：使用awgn函数向信号中添加高斯白噪声，信噪比由config.snr\_dB指定。

### **4.4 步骤4：CMA均衡器实现（混合处理架构）**

这是项目的核心。必须采用一种混合架构，解耦逐样本的权重更新和基于块的滤波操作。

1. **初始化**：  
   * 初始化均衡器权重向量w，采用**中心抽头**方式：w \= zeros(config.cma.filterLength, 1); w(floor(config.cma.filterLength/2)+1) \= 1; 17。  
   * 初始化输入和输出缓冲区。  
2. **主处理循环（逐样本）**：  
   * 循环遍历每一个接收到的符号样本n。  
   * **获取均衡器输出**：从**输出缓冲区**中取出当前样本x\_tilde(n)。  
   * **CMA权重更新**：  
     * 获取对应的输入信号向量y\_vec(n)（即均衡器FIR滤波器的输入抽头）。  
     * 根据CMA更新公式计算误差并更新权重w：  
       error\_term \= config.cma.R2 \- abs(x\_tilde(n))^2;  
       w \= w \+ config.cma.stepSize \* error\_term \* x\_tilde(n) \* conj(y\_vec(n));  
   * **触发块滤波**：  
     * 每当接收到L \= config.fft.nfft \- config.cma.filterLength \+ 1个**新**样本时，触发一次快速卷积。  
     * 调用fast\_conv\_os函数，使用**当前最新**的权重向量w和输入缓冲区中的数据块进行滤波。  
     * 将滤波结果的有效部分存入**输出缓冲区**。

### **4.5 步骤5：相位校正与解调**

1. **相位校正**：在均衡器收敛后（例如，忽略前一部分输出），实现一个简单的判决辅助（Decision-Directed）相位校正环路，以消除CMA的相位模糊。  
   * 对均衡器输出x\_tilde(n)做硬判决，得到最近的理想星座点a\_hat(n)。  
   * 计算相位误差：phase\_error \= angle(a\_hat(n) \* conj(x\_tilde(n)))。  
   * 通过一个简单的一阶低通滤波器平滑相位误差，得到相位校正值。  
   * 对x\_tilde(n)进行反向旋转校正。  
2. **解调**：对相位校正后的信号进行硬判决，得到解调符号索引。  
3. **格雷码反向映射**：将解调索引通过格雷码逆映射转换回原始索引。  
4. **符号到比特**：将索引转换为比特流。

### **4.6 步骤6：性能评估**

1. **对齐序列**：由于信道和均衡器引入延迟，必须在计算BER前，将解调出的比特流与原始发射比特流进行对齐。  
2. **计算BER**：比较对齐后的序列，计算误码率。  
3. **可视化**：  
   * 使用scatterplot绘制“均衡前”、“均衡后（CMA输出）”和“相位校正后”的星座图。  
   * 在MATLAB命令窗口打印最终的BER结果。

## **第5章：模块四：参数影响分析实验**

**任务**：创建一个脚本run\_experiments.m，用于系统性地分析关键参数对均衡效果的影响，并生成对比图表。

### **5.1 实验1：步长因子  的影响**

1. **设置**：固定滤波器长度M（如31）和SNR（如25dB）。  
2. **执行**：在一个循环中，使用不同的步长值mu\_values \= \[1e-4, 1e-3, 5e-3\]运行main\_cma\_simulation.m。  
3. **输出**：  
   * 为每个$\\mu$值绘制收敛后的星座图。  
   * 绘制代价函数（或MSE）随符号数变化的收敛曲线图，对比不同$\\mu$值下的收敛速度和稳态误差。  
   * **分析**：观察并记录大步长收敛快但稳态误差大，小步长收敛慢但稳态误差小的现象 5。

### **5.2 实验2：滤波器长度 M 的影响**

1. **设置**：固定步长$\\mu$（使用实验1中找到的较优值）和SNR。  
2. **执行**：在一个循环中，使用不同的滤波器长度M\_values \= （应与信道长度length(config.channelTaps)进行比较）运行仿真。  
3. **输出**：  
   * 为每个M值绘制收敛后的星座图和最终的BER。  
   * **分析**：观察当M过短时，星座图无法完全“开眼”，BER性能差。当M足够长时（如信道长度的2-3倍），性能显著改善 2。

### **5.3 实验3：信噪比（SNR）性能曲线**

1. **设置**：固定步长$\\mu$和滤波器长度M。  
2. **执行**：在一个循环中，遍历一个SNR范围snr\_range\_dB \= 10:2:30，对每个SNR值运行仿真并记录BER。同时，也计算无均衡器情况下的BER作为基准。  
3. **输出**：  
   * 绘制两条BER vs. SNR曲线（有均衡器 vs. 无均衡器）在同一张图上，使用semilogy函数。  
   * **分析**：定量展示CMA均衡器在存在ISI的情况下对系统性能的显著提升。

#### **引用的著作**

1. Blind equalization using the constant modulus criterion \- School of Electrical and Computer Engineering, 访问时间为 十月 14, 2025， [https://people.ece.cornell.edu/johnson/pubLinxEE/IEEE-Proc-10-1998](https://people.ece.cornell.edu/johnson/pubLinxEE/IEEE-Proc-10-1998)  
2. adaptive filter algorithms for channel equalization \- DiVA portal, 访问时间为 十月 14, 2025， [https://www.diva-portal.org/smash/get/diva2:1311153/FULLTEXT01.pdf](https://www.diva-portal.org/smash/get/diva2:1311153/FULLTEXT01.pdf)  
3. 大作业要求250917.pdf  
4. Adaptive Equalization Algorithms:An Overview \- The Science and Information (SAI) Organization, 访问时间为 十月 14, 2025， [https://thesai.org/Downloads/Volume2No3/Paper%2011-%20Adaptive%20Equalization%20Algorithms%20An%20Overview.pdf](https://thesai.org/Downloads/Volume2No3/Paper%2011-%20Adaptive%20Equalization%20Algorithms%20An%20Overview.pdf)  
5. OPTIMAL STEP-SIZE CONSTANT MODULUS ALGORITHM \- Laboratoire I3S, 访问时间为 十月 14, 2025， [https://webusers.i3s.unice.fr/\~zarzoso/biblio/rep04cma.pdf](https://webusers.i3s.unice.fr/~zarzoso/biblio/rep04cma.pdf)  
6. An Introduction to Constant Modulus Algorithm (CMA) \- Wireless Pi, 访问时间为 十月 14, 2025， [https://wirelesspi.com/an-introduction-to-constant-modulus-algorithm-cma/](https://wirelesspi.com/an-introduction-to-constant-modulus-algorithm-cma/)  
7. 4\. THE CONSTANT MODULUS ALGORITHM, 访问时间为 十月 14, 2025， [https://sps.ewi.tudelft.nl/Education/courses/et4147/sheets/cma\_leus.pdf](https://sps.ewi.tudelft.nl/Education/courses/et4147/sheets/cma_leus.pdf)  
8. CMA on Modem Data Signal, 访问时间为 十月 14, 2025， [http://bard.ece.cornell.edu/downloads/tutorials/cmamdm/cmamdm.html](http://bard.ece.cornell.edu/downloads/tutorials/cmamdm/cmamdm.html)  
9. Adaptive Equalizers \- MATLAB & Simulink \- MathWorks, 访问时间为 十月 14, 2025， [https://es.mathworks.com/help/comm/ug/adaptive-equalizers.html](https://es.mathworks.com/help/comm/ug/adaptive-equalizers.html)  
10. Convert Fast Fourier Transform (FFT) to Fixed Point \- MATLAB & Simulink \- MathWorks, 访问时间为 十月 14, 2025， [https://www.mathworks.com/help/fixedpoint/ug/convert-fast-fourier-transform-fft-to-fixed-point.html](https://www.mathworks.com/help/fixedpoint/ug/convert-fast-fourier-transform-fft-to-fixed-point.html)  
11. FFT \- Fast Fourier transform (FFT) of input \- Simulink \- MathWorks, 访问时间为 十月 14, 2025， [https://www.mathworks.com/help/dsp/ref/fft.html](https://www.mathworks.com/help/dsp/ref/fft.html)  
12. fft (MATLAB Function Reference), 访问时间为 十月 14, 2025， [https://math.jhu.edu/\~shiffman/370/help/techdoc/ref/fft.html](https://math.jhu.edu/~shiffman/370/help/techdoc/ref/fft.html)  
13. fft \- Fast Fourier transform \- MATLAB \- MathWorks, 访问时间为 十月 14, 2025， [https://www.mathworks.com/help/matlab/ref/fft.html](https://www.mathworks.com/help/matlab/ref/fft.html)  
14. Overlap-Add/Save \- MATLAB & Simulink \- MathWorks, 访问时间为 十月 14, 2025， [https://www.mathworks.com/help/dsp/ug/overlap-add-save.html](https://www.mathworks.com/help/dsp/ug/overlap-add-save.html)  
15. Overlap Save Method \- File Exchange \- MATLAB Central \- MathWorks, 访问时间为 十月 14, 2025， [https://www.mathworks.com/matlabcentral/fileexchange/41337-overlap-save-method](https://www.mathworks.com/matlabcentral/fileexchange/41337-overlap-save-method)  
16. Overlap–save method \- Wikipedia, 访问时间为 十月 14, 2025， [https://en.wikipedia.org/wiki/Overlap%E2%80%93save\_method](https://en.wikipedia.org/wiki/Overlap%E2%80%93save_method)  
17. Concurrent Modified Constant Modulus Algorithm and Decision Directed Scheme With Barzilai-Borwein Method \- PMC \- PubMed Central, 访问时间为 十月 14, 2025， [https://pmc.ncbi.nlm.nih.gov/articles/PMC8222806/](https://pmc.ncbi.nlm.nih.gov/articles/PMC8222806/)  
18. Adaptive Equalizers \- MATLAB & Simulink \- MathWorks, 访问时间为 十月 14, 2025， [https://www.mathworks.com/help/comm/ug/adaptive-equalizers.html](https://www.mathworks.com/help/comm/ug/adaptive-equalizers.html)  
19. LMS Algorithm Step Size Adjustment for Fast Convergence \- Biblioteka Nauki, 访问时间为 十月 14, 2025， [https://bibliotekanauki.pl/articles/177252.pdf](https://bibliotekanauki.pl/articles/177252.pdf)  
20. SV4: Equalization and Adaptive Filters, 访问时间为 十月 14, 2025， [https://people.ee.ethz.ch/\~isistaff/labs/SV4/SV4.pdf](https://people.ee.ethz.ch/~isistaff/labs/SV4/SV4.pdf)