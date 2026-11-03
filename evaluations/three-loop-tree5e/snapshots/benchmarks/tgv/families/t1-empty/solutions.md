# T1[] 的 I：完整 ε 表示与 Laurent 系数

七个三圈积分01、02、03、04、17、18、19均已得到当前唯一I规范下的ε⁻⁶…ε⁰，余项O(ε)，固定z>0。formal=not_run。

[统一推导与完整ε定义](method.md) · [计算及核验记录](../../../../workflow/checks/tgv-t1-empty/README.md)

01–04与z无关；17–19仅通过 R=cosh(ε log(z)/2)/cosh(3ε log(z)/2) 依赖z。系数只含π、ζ(3)、ζ(5)及log(z)的多项式，不新增M/τ函数。

全ε恒等式（图号均指5e-3g）：

$$I_{02}=2I_{03},\qquad I_{04}=I_{01}-I_{03},\qquad I_{18}=I_{17}-I_{19}.$$

这些关系来自完整有序域的积分，不只是某一阶系数拟合。下列四个独立解表连同上述线性关系给出全部49个系数；各图JSON另逐项展开存储。

| 图 | 最高极点 | 完整ε形式 |
|---|---:|---|
| [5e-3g-01](../../5e-3g/results/5e-3g-01/r001.json) | ε⁻6 | Gamma |
| [5e-3g-02](../../5e-3g/results/5e-3g-02/r001.json) | ε⁻6 | Gamma＋3F2 |
| [5e-3g-03](../../5e-3g/results/5e-3g-03/r001.json) | ε⁻6 | Gamma＋3F2 |
| [5e-3g-04](../../5e-3g/results/5e-3g-04/r001.json) | ε⁻6 | Gamma＋3F2 |
| [5e-3g-17](../../5e-3g/results/5e-3g-17/r001.json) | ε⁻4 | Gamma |
| [5e-3g-18](../../5e-3g/results/5e-3g-18/r001.json) | ε⁻4 | Gamma＋已逐项积完的显式收敛和 |
| [5e-3g-19](../../5e-3g/results/5e-3g-19/r001.json) | ε⁻3 | Gamma＋已逐项积完的显式收敛和 |

记 $\ell=\log z$，$I_g=\sum_{k=-6}^0c_{g,k}\varepsilon^k+O(\varepsilon)$。

## 5e-3g-01

| k | cₖ |
|---:|---|
| -6 | $\frac{1}{576}$ |
| -5 | $- \frac{1}{288}$ |
| -4 | $\frac{-48 + 11 \pi^{2}}{6912}$ |
| -3 | $- \frac{- 10 \zeta\left(3\right) + \pi^{2} + 16}{1152}$ |
| -2 | $\frac{-512 - 32 \pi^{2} - 64 \zeta\left(3\right) + 13 \pi^{4}}{18432}$ |
| -1 | $- \frac{- 780 \pi^{2} \zeta\left(3\right) - 3024 \zeta\left(5\right) + 960 \zeta\left(3\right) + 19 \pi^{4} + 480 \pi^{2} + 7680}{138240}$ |
| 0 | $\frac{-2580480 - 161280 \pi^{2} - 6384 \pi^{4} - 322560 \zeta\left(3\right) - 20160 \pi^{2} \zeta\left(3\right) - 48384 \zeta\left(5\right) + 342720 \zeta^{2}\left(3\right) + 5893 \pi^{6}}{23224320}$ |

## 5e-3g-03

| k | cₖ |
|---:|---|
| -6 | $\frac{1}{1152}$ |
| -5 | $0$ |
| -4 | $- \frac{5 \pi^{2}}{13824}$ |
| -3 | $\frac{7 \zeta\left(3\right)}{384}$ |
| -2 | $- \frac{127 \pi^{4}}{184320}$ |
| -1 | $- \frac{- 6012 \zeta\left(5\right) + 175 \pi^{2} \zeta\left(3\right)}{23040}$ |
| 0 | $- \frac{- 440640 \zeta^{2}\left(3\right) + 3581 \pi^{6}}{6635520}$ |

## 5e-3g-17

| k | cₖ |
|---:|---|
| -6 | $0$ |
| -5 | $0$ |
| -4 | $\frac{\pi^{2}}{384}$ |
| -3 | $- \frac{\pi^{2} + 9 \zeta\left(3\right)}{192}$ |
| -2 | $- \frac{\pi^{2}}{384}\ell^{2} + \frac{- 240 \pi^{2} + 2160 \zeta\left(3\right) + 77 \pi^{4}}{23040}$ |
| -1 | $\frac{\pi^{2} + 9 \zeta\left(3\right)}{192}\ell^{2} + - \frac{- 2160 \zeta\left(3\right) + 165 \pi^{2} \zeta\left(3\right) + 240 \pi^{2} + 37 \pi^{4} + 3780 \zeta\left(5\right)}{11520}$ |
| 0 | $\frac{11 \pi^{2}}{4608}\ell^{4} + - \frac{- 240 \pi^{2} + 2160 \zeta\left(3\right) + 77 \pi^{4}}{23040}\ell^{2} + \frac{- 24864 \pi^{4} - 1451520 \zeta^{2}\left(3\right) - 161280 \pi^{2} - 50400 \pi^{2} \zeta\left(3\right) + 1451520 \zeta\left(3\right) + 2540160 \zeta\left(5\right) + 9419 \pi^{6}}{3870720}$ |

## 5e-3g-19

| k | cₖ |
|---:|---|
| -6 | $0$ |
| -5 | $0$ |
| -4 | $0$ |
| -3 | $- \frac{- 9 \zeta\left(3\right) + \pi^{2}}{192}$ |
| -2 | $- \frac{- 2160 \zeta\left(3\right) + 7 \pi^{4} + 240 \pi^{2}}{23040}$ |
| -1 | $\frac{- 9 \zeta\left(3\right) + \pi^{2}}{192}\ell^{2} + - \frac{- 3780 \zeta\left(5\right) - 2160 \zeta\left(3\right) - 15 \pi^{2} \zeta\left(3\right) + 240 \pi^{2} + 37 \pi^{4}}{11520}$ |
| 0 | $\frac{- 2160 \zeta\left(3\right) + 7 \pi^{4} + 240 \pi^{2}}{23040}\ell^{2} + - \frac{- 1270080 \zeta\left(5\right) - 725760 \zeta\left(3\right) - 362880 \zeta^{2}\left(3\right) + 25200 \pi^{2} \zeta\left(3\right) + 737 \pi^{6} + 80640 \pi^{2} + 12432 \pi^{4}}{1935360}$ |

## 核验与边界

新鲜HyperInt积分得到F的6阶及Z的4阶矩；发散检查开启。F的反向积分次序与裸TGV Gamma控制均通过；七图原始核、换元Jacobian、唯一C/I连接逐项核对。

21个高精度I截断余项样本符合O(ε)，另有6个有限矩余项检查、7组原始七维J的随机化QMC比较以及完整Z级数的正尾界检查。QMC最大相对差约1.05%，经验误差比较通过，不是高精度或严格数值包围。

18/19的完整ε表达保留显式收敛级数，未声称有限Gamma闭式；未证明物理复z/i0延拓或Lean整链。原分类清单与输入阶段工作区是历史快照，当前求解状态以本解表和结果r001为准。
