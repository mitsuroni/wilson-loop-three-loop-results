# 2e-Pi2-01：相邻3/4，三个颜色的完整I至有限项（现行r002）

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi2/families/adjacent/solutions.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。


以下I已采用物理1PI插入约定，组装系数为+1。G、U和各母函数保持原定义；J及有序域不变。系数表覆盖ε⁻⁶至ε⁰，未显示的位置为零，余项O(ε)。

## 2e-Pi2-01

[机器结果](../../results/2e-Pi2-01/r002.json)；[现行prepared](../../prepared/2e-Pi2-r002.json)。

$$
I=- \frac{\left(- \varepsilon - 2\right) e^{- 3 \gamma_E \varepsilon} \Gamma\left(1 - 3 \varepsilon\right) \Gamma^{4}\left(1 - \varepsilon\right) \Gamma^{2}\left(\varepsilon + 1\right)}{1152 \varepsilon^{4} \left(1 - 2 \varepsilon\right)^{2} \left(2 \varepsilon + 1\right) \Gamma^{2}\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right)}
$$

| ε幂次 | I 的系数 |
|---:|---|
| -4 | $\frac{1}{576}$ |
| -3 | $\frac{5}{1152}$ |
| -2 | $\frac{\pi^{2}}{2304} + \frac{1}{64}$ |
| -1 | $\frac{5 \pi^{2}}{4608} + \frac{7 \zeta\left(3\right)}{576} + \frac{5}{144}$ |
| 0 | $\frac{67 \pi^{4}}{276480} + \frac{35 \zeta\left(3\right)}{1152} + \frac{\pi^{2}}{256} + \frac{7}{72}$ |

### 颜色 CA^2*CF

| ε幂次 | 完整颜色系数 |
|---:|---|
| -4 | $- \frac{25}{1728}$ |
| -3 | $- \frac{179}{5184}$ |
| -2 | $- \frac{2941}{15552} - \frac{25 \pi^{2}}{6912} + \frac{\zeta\left(3\right)}{144}$ |
| -1 | $- \frac{52453}{93312} - \frac{179 \pi^{2}}{20736} + \frac{\pi^{4}}{8640} + \frac{61 \zeta\left(3\right)}{1728}$ |
| 0 | $- \frac{1187401}{559872} - \frac{2941 \pi^{2}}{62208} + \frac{\pi^{2} \zeta\left(3\right)}{576} + \frac{71 \pi^{4}}{276480} + \frac{7 \zeta\left(5\right)}{144} + \frac{605 \zeta\left(3\right)}{5184}$ |

### 颜色 CA*CF*TF*nf

| ε幂次 | 完整颜色系数 |
|---:|---|
| -4 | $\frac{5}{432}$ |
| -3 | $\frac{7}{324}$ |
| -2 | $\frac{5 \pi^{2}}{1728} + \frac{\zeta\left(3\right)}{18} + \frac{199}{1944}$ |
| -1 | $\frac{7 \pi^{2}}{1296} + \frac{7 \zeta\left(3\right)}{144} + \frac{\pi^{4}}{1080} + \frac{4751}{23328}$ |
| 0 | $\frac{223 \pi^{4}}{207360} + \frac{\pi^{2} \zeta\left(3\right)}{72} + \frac{199 \pi^{2}}{7776} + \frac{7 \zeta\left(5\right)}{18} + \frac{73 \zeta\left(3\right)}{216} + \frac{97655}{139968}$ |

### 颜色 CF^2*TF*nf

| ε幂次 | 完整颜色系数 |
|---:|---|
| -3 | $\frac{1}{72}$ |
| -2 | $\frac{23}{216} - \frac{\zeta\left(3\right)}{9}$ |
| -1 | $- \frac{5 \zeta\left(3\right)}{27} - \frac{\pi^{4}}{540} + \frac{\pi^{2}}{288} + \frac{331}{648}$ |
| 0 | $- \frac{557 \zeta\left(3\right)}{648} - \frac{7 \zeta\left(5\right)}{9} - \frac{\pi^{2} \zeta\left(3\right)}{36} - \frac{\pi^{4}}{324} + \frac{23 \pi^{2}}{864} + \frac{7993}{3888}$ |

由旧定义到现行定义的线性变换、颜色和系数连接已核对；未重新积分，不继承新形式证明。formal=not_run。
