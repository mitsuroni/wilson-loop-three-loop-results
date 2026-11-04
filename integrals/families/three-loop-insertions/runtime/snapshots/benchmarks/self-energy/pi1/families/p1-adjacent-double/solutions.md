# P1：完整 I 与有限项（现行r002）

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1/families/p1-adjacent-double/solutions.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。


以下I已采用物理1PI插入约定，组装系数为+1。G、U和各母函数保持原定义；J及有序域不变。系数表覆盖ε⁻⁶至ε⁰，未显示的位置为零，余项O(ε)。

## 4e-Pi1-01

[机器结果](../../results/4e-Pi1-01/r002.json)；[现行prepared](../../prepared/4e-Pi1-r002.json)。

$$
I=- \frac{e^{- 3 \gamma_E \varepsilon} \Gamma^{3}\left(1 - \varepsilon\right)}{576 \varepsilon^{5} \left(1 - 2 \varepsilon\right) \left(\varepsilon + 1\right)}
$$

| ε幂次 | I 的系数 |
|---:|---|
| -5 | $- \frac{1}{576}$ |
| -4 | $- \frac{1}{576}$ |
| -3 | $- \frac{1}{192} - \frac{\pi^{2}}{2304}$ |
| -2 | $- \frac{5}{576} - \frac{\pi^{2}}{2304} - \frac{\zeta\left(3\right)}{576}$ |
| -1 | $- \frac{11}{576} - \frac{\pi^{2}}{768} - \frac{19 \pi^{4}}{276480} - \frac{\zeta\left(3\right)}{576}$ |
| 0 | $- \frac{7}{192} - \frac{5 \pi^{2}}{2304} - \frac{19 \pi^{4}}{276480} - \frac{\zeta\left(3\right)}{192} - \frac{\pi^{2} \zeta\left(3\right)}{2304} - \frac{\zeta\left(5\right)}{960}$ |

## 4e-Pi1-02

[机器结果](../../results/4e-Pi1-02/r002.json)；[现行prepared](../../prepared/4e-Pi1-r002.json)。

$$
I=- \frac{e^{- 3 \gamma_E \varepsilon} \Gamma^{3}\left(1 - \varepsilon\right)}{576 \varepsilon^{5} \left(1 - 2 \varepsilon\right) \left(\varepsilon + 1\right)}
$$

| ε幂次 | I 的系数 |
|---:|---|
| -5 | $- \frac{1}{576}$ |
| -4 | $- \frac{1}{576}$ |
| -3 | $- \frac{1}{192} - \frac{\pi^{2}}{2304}$ |
| -2 | $- \frac{5}{576} - \frac{\pi^{2}}{2304} - \frac{\zeta\left(3\right)}{576}$ |
| -1 | $- \frac{11}{576} - \frac{\pi^{2}}{768} - \frac{19 \pi^{4}}{276480} - \frac{\zeta\left(3\right)}{576}$ |
| 0 | $- \frac{7}{192} - \frac{5 \pi^{2}}{2304} - \frac{19 \pi^{4}}{276480} - \frac{\zeta\left(3\right)}{192} - \frac{\pi^{2} \zeta\left(3\right)}{2304} - \frac{\zeta\left(5\right)}{960}$ |

由旧定义到现行定义的线性变换、颜色和系数连接已核对；未重新积分，不继承新形式证明。formal=not_run。
