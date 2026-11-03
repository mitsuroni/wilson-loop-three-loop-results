# 6e-37：[1,2] / [2,1] 族解析解

日期：2026-09-17。当前 prepared 0.3 根；完整 ε 表达与 ε⁻⁶ 至 ε⁰ 已计算检查，formal=not_run。

## 输入与原始次序

- [原输入](../../prepared/6e.json)：6e-37，包 SHA-256：9699532d1ee3a4027c2eb93c961a5b467a9bceda229c0204c13775cfa2451406。
- a1=1/2，sym=2，κ=-1/16；C 已补偿一次 Euler 因子。
- 原外到内顺序：s11 → s31 → s22 → s42 → s21 → s41。
- 参数：z=z1>0；完整 I 先取实 0<ε<1，固定 z 在 ε=0 展开。

| 变量 | 下界 | 上界 |
|---|---|---|
| s11 | 0 | 1 |
| s31 | 0 | 1 |
| s22 | 0 | 1 |
| s42 | 0 | 1 |
| s21 | 0 | s22 |
| s41 | 0 | s42 |

原被积函数：

$$
\left(- z - 1\right)^{3} \left(s_{11} s_{31} + z \left(1 - s_{11}\right) \left(1 - s_{31}\right)\right)^{\varepsilon - 1} \left(s_{21} s_{41} z + \left(1 - s_{21}\right) \left(1 - s_{41}\right)\right)^{\varepsilon - 1} \left(s_{22} s_{42} z + \left(1 - s_{22}\right) \left(1 - s_{42}\right)\right)^{\varepsilon - 1}
$$

## 分组、母核与完整结果

13 通道只有独立的 s11、s31；24 通道有两条有序链 s21<s22、s41<s42。六个端点互异。
逐个传播子重建确定本图配对为 id，排列 [0, 1]。
所有核的角点奇性在 Re(ε)>−1 可积，放宽到立方体后也绝对可积；据此才能分组。
反射 s11→1−s11、s31→1−s31 将 tilde D 变为 D。H0、Qid、Qswap 均沿用已有定义。

$$
J=H_0 Q_{\mathrm{id}} z^{\frac{3 \varepsilon}{2}}
$$

$$
I=- \frac{G^{3} H_0 Q_{\mathrm{id}}}{16 \cosh{\left(\frac{3 \ell \varepsilon}{2} \right)}}
$$

[公共方法](../../families/mixed-opposite/method.md)给出通用 J→I、crossing、可积界、Meijer G 和已积完的收敛级数。
[不可变结果 r001](../../results/6e-37/r001.json)保存代入 Meijer G 后的完整表达与级数定义。

## Laurent 系数

| ε 幂 | 解析系数 |
|---:|---|
| -6 | $0$ |
| -5 | $0$ |
| -4 | $0$ |
| -3 | $0$ |
| -2 | $0$ |
| -1 | $0$ |
| 0 | $- \frac{M_{0,1}^{3}}{96}$ |

$$I=-\frac{M_{0,1}(z)^3}{96}+O(\varepsilon)=\frac{(\ln^2z+\pi^2)^3}{768}+O(\varepsilon).$$

零极点来自 ε=0 邻域的全纯性；本批不需要新的基函数。两个图的有限项相等，但完整 ε 函数不同。
M01=τ00=−(ln²z+π²)/2；完整 Hlog 坐标保存在结果 JSON。

## 检查与边界

[实际计算检查](../../../../workflow/checks/mixed-opposite/runs/20260917T021108676259Z/checks.json)；[复算说明](../../../../workflow/checks/mixed-opposite/README.md)。
原顺序整数 ε 积分、终止双级数、非整数六维 Sobol、独立系数求积及 O(ε) 余项缩放均已检查。
已有母核输出按明确版本与哈希复用，没有重跑相同 HyperInt 积分。数值一致性不是形式证明。
未修改 Kernel/Lean，未执行本批形式全链、复 z/i0 延拓或一般正幂系数交付。
