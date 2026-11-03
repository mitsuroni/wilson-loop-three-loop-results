# T2[] 的 I：完整 ε 表示与 Laurent 系数

七个三圈积分05、06、07、08、20、21、31均已得到完整ε解析表达与ε⁻⁶…ε⁰共49项（40项非零），余项O(ε)，固定z>0。formal=not_run。

[原始J的反射推导及全部符号定义](method.md) · [计算核验](../../../../workflow/checks/tgv-t2-empty/README.md)

精确反射保持z、原始J及完整C/I规范。下表的相等关系在完整ε下成立；图号均指5e-3g。

| T2 | 相等的T1 | 最高极点 | 完整ε表达 |
|---|---|---:|---|
| [05](../../5e-3g/results/5e-3g-05/r001.json) | 04 | ε⁻6 | $\lambda Gh_1/(576e^5)$ |
| [06](../../5e-3g/results/5e-3g-06/r001.json) | 03 | ε⁻6 | $\lambda G(h_0-h_1)/(576e^5)$ |
| [07](../../5e-3g/results/5e-3g-07/r001.json) | 02 | ε⁻6 | $2I_{06}$ |
| [08](../../5e-3g/results/5e-3g-08/r001.json) | 01 | ε⁻6 | $GU/(576e^6)$ |
| [20](../../5e-3g/results/5e-3g-20/r001.json) | 19 | ε⁻3 | $\lambda GRQ/(64e^4)$ |
| [21](../../5e-3g/results/5e-3g-21/r001.json) | 18 | ε⁻4 | $I_{31}-I_{20}$ |
| [31](../../5e-3g/results/5e-3g-31/r001.json) | 17 | ε⁻4 | $GRU(1-V)/(128e^6)$ |

这里 $e=\varepsilon$，$G=e^{-3\gamma_Ee}\Gamma(1-e)\Gamma(1-2e)$，$\lambda=1-2e$，$R=\cosh(e\log z/2)/\cosh(3e\log z/2)$；$U,V,h_0,h_1,Q$的Beta/3F2/收敛和定义见[推导第2节](method.md#2-自包含的完整ε表达)。

05–08与z无关；20/21/31只通过R依赖z。20/21的完整ε表达保留已逐项积完的显式收敛级数，其余为Gamma或Gamma＋3F2。Laurent系数不需要新M/τ函数。

$$I_{07}=2I_{06},\qquad I_{08}=I_{05}+I_{06},\qquad I_{31}=I_{20}+I_{21}.$$

令 $\ell=\log z$，$I_g=\sum_{k=-6}^{0}c_{g,k}\varepsilon^k+O(\varepsilon)$。以下逐图列出全部系数。

## 5e-3g-05

| k | cₖ |
|---:|---|
| -6 | $\frac{1}{1152}$ |
| -5 | $- \frac{1}{288}$ |
| -4 | $\frac{-32 + 9 \pi^{2}}{4608}$ |
| -3 | $- \frac{\pi^{2} + 11 \zeta\left(3\right) + 16}{1152}$ |
| -2 | $\frac{-5120 - 320 \pi^{2} - 640 \zeta\left(3\right) + 257 \pi^{4}}{184320}$ |
| -1 | $- \frac{- 1830 \pi^{2} \zeta\left(3\right) + 960 \zeta\left(3\right) + 19 \pi^{4} + 480 \pi^{2} + 7680 + 33048 \zeta\left(5\right)}{138240}$ |
| 0 | $\frac{-5160960 - 2399040 \zeta^{2}\left(3\right) - 322560 \pi^{2} - 12768 \pi^{4} - 645120 \zeta\left(3\right) - 40320 \pi^{2} \zeta\left(3\right) - 96768 \zeta\left(5\right) + 36853 \pi^{6}}{46448640}$ |

## 5e-3g-06

| k | cₖ |
|---:|---|
| -6 | $\frac{1}{1152}$ |
| -5 | $0$ |
| -4 | $- \frac{5 \pi^{2}}{13824}$ |
| -3 | $\frac{7 \zeta\left(3\right)}{384}$ |
| -2 | $- \frac{127 \pi^{4}}{184320}$ |
| -1 | $- \frac{- 6012 \zeta\left(5\right) + 175 \pi^{2} \zeta\left(3\right)}{23040}$ |
| 0 | $- \frac{- 440640 \zeta^{2}\left(3\right) + 3581 \pi^{6}}{6635520}$ |

## 5e-3g-07

| k | cₖ |
|---:|---|
| -6 | $\frac{1}{576}$ |
| -5 | $0$ |
| -4 | $- \frac{5 \pi^{2}}{6912}$ |
| -3 | $\frac{7 \zeta\left(3\right)}{192}$ |
| -2 | $- \frac{127 \pi^{4}}{92160}$ |
| -1 | $- \frac{- 6012 \zeta\left(5\right) + 175 \pi^{2} \zeta\left(3\right)}{11520}$ |
| 0 | $- \frac{- 440640 \zeta^{2}\left(3\right) + 3581 \pi^{6}}{3317760}$ |

## 5e-3g-08

| k | cₖ |
|---:|---|
| -6 | $\frac{1}{576}$ |
| -5 | $- \frac{1}{288}$ |
| -4 | $\frac{-48 + 11 \pi^{2}}{6912}$ |
| -3 | $- \frac{- 10 \zeta\left(3\right) + \pi^{2} + 16}{1152}$ |
| -2 | $\frac{-512 - 32 \pi^{2} - 64 \zeta\left(3\right) + 13 \pi^{4}}{18432}$ |
| -1 | $- \frac{- 780 \pi^{2} \zeta\left(3\right) - 3024 \zeta\left(5\right) + 960 \zeta\left(3\right) + 19 \pi^{4} + 480 \pi^{2} + 7680}{138240}$ |
| 0 | $\frac{-2580480 - 161280 \pi^{2} - 6384 \pi^{4} - 322560 \zeta\left(3\right) - 20160 \pi^{2} \zeta\left(3\right) - 48384 \zeta\left(5\right) + 342720 \zeta^{2}\left(3\right) + 5893 \pi^{6}}{23224320}$ |

## 5e-3g-20

| k | cₖ |
|---:|---|
| -6 | $0$ |
| -5 | $0$ |
| -4 | $0$ |
| -3 | $- \frac{- 9 \zeta\left(3\right) + \pi^{2}}{192}$ |
| -2 | $- \frac{- 2160 \zeta\left(3\right) + 7 \pi^{4} + 240 \pi^{2}}{23040}$ |
| -1 | $\frac{- 9 \zeta\left(3\right) + \pi^{2}}{192}\ell^{2} + - \frac{- 3780 \zeta\left(5\right) - 2160 \zeta\left(3\right) - 15 \pi^{2} \zeta\left(3\right) + 240 \pi^{2} + 37 \pi^{4}}{11520}$ |
| 0 | $\frac{- 2160 \zeta\left(3\right) + 7 \pi^{4} + 240 \pi^{2}}{23040}\ell^{2} + - \frac{- 1270080 \zeta\left(5\right) - 725760 \zeta\left(3\right) - 362880 \zeta^{2}\left(3\right) + 25200 \pi^{2} \zeta\left(3\right) + 737 \pi^{6} + 80640 \pi^{2} + 12432 \pi^{4}}{1935360}$ |

## 5e-3g-21

| k | cₖ |
|---:|---|
| -6 | $0$ |
| -5 | $0$ |
| -4 | $\frac{\pi^{2}}{384}$ |
| -3 | $- \frac{3 \zeta\left(3\right)}{32}$ |
| -2 | $- \frac{\pi^{2}}{384}\ell^{2} + \frac{7 \pi^{4}}{1920}$ |
| -1 | $\frac{3 \zeta\left(3\right)}{32}\ell^{2} + - \frac{\pi^{2} \zeta\left(3\right) + 42 \zeta\left(5\right)}{64}$ |
| 0 | $\frac{11 \pi^{2}}{4608}\ell^{4} + - \frac{7 \pi^{4}}{1920}\ell^{2} + \frac{- 725760 \zeta^{2}\left(3\right) + 3631 \pi^{6}}{1290240}$ |

## 5e-3g-31

| k | cₖ |
|---:|---|
| -6 | $0$ |
| -5 | $0$ |
| -4 | $\frac{\pi^{2}}{384}$ |
| -3 | $- \frac{\pi^{2} + 9 \zeta\left(3\right)}{192}$ |
| -2 | $- \frac{\pi^{2}}{384}\ell^{2} + \frac{- 240 \pi^{2} + 2160 \zeta\left(3\right) + 77 \pi^{4}}{23040}$ |
| -1 | $\frac{\pi^{2} + 9 \zeta\left(3\right)}{192}\ell^{2} + - \frac{- 2160 \zeta\left(3\right) + 165 \pi^{2} \zeta\left(3\right) + 240 \pi^{2} + 37 \pi^{4} + 3780 \zeta\left(5\right)}{11520}$ |
| 0 | $\frac{11 \pi^{2}}{4608}\ell^{4} + - \frac{- 240 \pi^{2} + 2160 \zeta\left(3\right) + 77 \pi^{4}}{23040}\ell^{2} + \frac{- 24864 \pi^{4} - 1451520 \zeta^{2}\left(3\right) - 161280 \pi^{2} - 50400 \pi^{2} \zeta\left(3\right) + 1451520 \zeta\left(3\right) + 2540160 \zeta\left(5\right) + 9419 \pi^{6}}{3870720}$ |

## 核验范围

七对原核、完整端点排序、alpha单纯形、Jacobian、颜色与C/I连接精确通过；四种错误映射被拒绝。旧T1的49系数及124项源/快照哈希已复验。此次没有重跑HyperInt；完整函数等式将已有余项证据直接传到T2。

新增14组T2原始七维J的Sobol检查：两个(ε,z)点，每图每点8次scramble、每次2^17点。最大相对差约2.949%，全部通过预设经验误差判据；不是严格数值包围。

未新增Lean/MCP整链收据或证明物理复z/i0延拓。原分类及输入阶段工作区保留历史快照，当前成果以本页和各图r001为准。
