homework1 - 双极化CMA盲均衡作业
=====================================

文件说明:
- TD_TRdata.mat     : 双极化实验数据 (RTdataX, RTdataY)
- CMA_reference.m   : 时域2x2 CMA参考实现
- CMA_homework.m    : 频域CMA实现（使用自编FFT）

运行方式:
1. 直接运行: run('CMA_homework.m')
2. 或通过入口: run('../core/main_cma_TRdata.m')

CMA_homework.m 实现要点:
- 使用 core/my_fft.m 和 core/my_ifft.m 实现频域快速卷积
- 2x2 MIMO结构处理双极化信号
- 包含V&V相位恢复算法

参数配置:
- seg_len = 32      : 每块符号数
- tap_len = 33      : FIR滤波器抽头数
- step = 2^-8       : 步长
- Rx = Ry = 2       : 恒模半径