# 单一规范代表与 mixed 有序区域 IBP

2026-09-17：3e-3g / 5e-3g 当前输出只保留一个规范代表，仍同时保留 raw 与 IBP 两种等价表示。

- mixed C 补上每条 pure Fourier 传播子的 (-1)，使 I^(2*n+3)*(-1)^n=I^3。
- 原 s 输出携带自身的 f_polynomials；alpha_measure 明确三参数单纯形。
- representative.canonical_topology 与实际 C/J 坐标一致，并保留旋转、腿/alpha置换和尺度关系。
- --tgv-ibp 现在支持 mixed；使用完整端点排序域，保留 TGV-pure 接触面。
- 边界共同因子由 corrected C_raw/(4*(1-2*eps)) 得到；非退化图保留 recoil。
- 先限制 W、S 与排序域，再应用局部有序域超立方映射和 Jacobian。
- 每条有序链的首端、末端、相邻接触面的通量分别是 -c_first、c_last、c_left-c_right。
- raw 文件不含 Euler 补偿。下游准备时 Gamma(1-eps)*Gamma(1-2eps) 对应一次 exp(-3*EulerGamma*eps)。
- status=complete 只表示 IBP 分解已生成；evaluation_status=not_evaluated，formal=not_run。
- type package 保留 J（兼容 raw cube）、J_raw_s 和 J_ibp；旧的 TGV/mixed 包必须重新生成。
- CLI 若有任一图生成失败，非零退出且不汇总可能残留的旧文件。

5e-3g 验证范围：44图共196项，含106 cusp、77 contact、13 recoil。31个退化图无recoil。
01、02、13在共同收敛域做过raw与完整IBP和的数值比较；这是数值证据，不是形式证明。
生成从 1/2<Re(eps)<1、z>0 的共同收敛条带建立等式，之后按亚纯延拓处理 eps=0。

Stage 3 尚未自动消费多项和。必须逐项求解后按 C_ratio_to_raw 汇总，不能丢掉 contact 或重复乘 coefficient。
本次未计算44图的I洛朗系数，未修改独立的2ne成果。

命令：
    montopar --source ../qgraf-pipeline/outputs/3loop/orbit_definitions.json --type 5e-3g --tgv-ibp --package-initial
