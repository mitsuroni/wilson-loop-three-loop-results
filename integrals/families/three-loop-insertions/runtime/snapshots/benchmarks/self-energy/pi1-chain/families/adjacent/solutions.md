# CH_A：2e-Pi1+Pi1-01 完整I与有限项（现行r002）

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1-chain/families/adjacent/solutions.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。


以下I已采用物理1PI插入约定，组装系数为+1。G、U和各母函数保持原定义；J及有序域不变。系数表覆盖ε⁻⁶至ε⁰，未显示的位置为零，余项O(ε)。

## 2e-Pi1+Pi1-01

[机器结果](../../results/2e-Pi1+Pi1-01/r002.json)；[现行prepared](../../prepared/2e-Pi1+Pi1-r002.json)。

$$
I=\frac{\left(- \varepsilon - 2\right) e^{- 3 \gamma_E \varepsilon} \Gamma\left(1 - 3 \varepsilon\right) \Gamma^{4}\left(1 - \varepsilon\right) \Gamma^{2}\left(\varepsilon + 1\right)}{1152 \varepsilon^{4} \left(1 - 2 \varepsilon\right)^{2} \left(2 \varepsilon + 1\right) \Gamma^{2}\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right)}
$$

| ε幂次 | I 的系数 |
|---:|---|
| -4 | $- \frac{1}{576}$ |
| -3 | $- \frac{5}{1152}$ |
| -2 | $- \frac{1}{64} - \frac{\pi^{2}}{2304}$ |
| -1 | $- \frac{5}{144} - \frac{7 \zeta\left(3\right)}{576} - \frac{5 \pi^{2}}{4608}$ |
| 0 | $- \frac{7}{72} - \frac{\pi^{2}}{256} - \frac{35 \zeta\left(3\right)}{1152} - \frac{67 \pi^{4}}{276480}$ |

### 颜色 CA^2*CF

| ε幂次 | 完整颜色系数 |
|---:|---|
| -4 | $- \frac{25}{1296}$ |
| -3 | $- \frac{395}{7776}$ |
| -2 | $- \frac{707}{3888} - \frac{25 \pi^{2}}{5184}$ |
| -1 | $- \frac{29023}{69984} - \frac{175 \zeta\left(3\right)}{1296} - \frac{395 \pi^{2}}{31104}$ |
| 0 | $- \frac{120923}{104976} - \frac{707 \pi^{2}}{15552} - \frac{2765 \zeta\left(3\right)}{7776} - \frac{335 \pi^{4}}{124416}$ |

### 颜色 CA*CF*TF*nf

| ε幂次 | 完整颜色系数 |
|---:|---|
| -4 | $\frac{5}{162}$ |
| -3 | $\frac{67}{972}$ |
| -2 | $\frac{5 \pi^{2}}{648} + \frac{61}{243}$ |
| -1 | $\frac{67 \pi^{2}}{3888} + \frac{35 \zeta\left(3\right)}{162} + \frac{4577}{8748}$ |
| 0 | $\frac{67 \pi^{4}}{15552} + \frac{469 \zeta\left(3\right)}{972} + \frac{61 \pi^{2}}{972} + \frac{19597}{13122}$ |

### 颜色 CF*(TF*nf)^2

| ε幂次 | 完整颜色系数 |
|---:|---|
| -4 | $- \frac{1}{81}$ |
| -3 | $- \frac{11}{486}$ |
| -2 | $- \frac{7}{81} - \frac{\pi^{2}}{324}$ |
| -1 | $- \frac{703}{4374} - \frac{7 \zeta\left(3\right)}{81} - \frac{11 \pi^{2}}{1944}$ |
| 0 | $- \frac{3179}{6561} - \frac{7 \pi^{2}}{324} - \frac{77 \zeta\left(3\right)}{486} - \frac{67 \pi^{4}}{38880}$ |

由旧定义到现行定义的线性变换、颜色和系数连接已核对；未重新积分，不继承新形式证明。formal=not_run。
