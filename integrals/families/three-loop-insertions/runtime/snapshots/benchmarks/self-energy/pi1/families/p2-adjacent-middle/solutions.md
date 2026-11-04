# P2：完整 I 与有限项（现行r002）

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1/families/p2-adjacent-middle/solutions.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。


以下I已采用物理1PI插入约定，组装系数为+1。G、U和各母函数保持原定义；J及有序域不变。系数表覆盖ε⁻⁶至ε⁰，未显示的位置为零，余项O(ε)。

## 4e-Pi1-05

[机器结果](../../results/4e-Pi1-05/r002.json)；[现行prepared](../../prepared/4e-Pi1-r002.json)。

$$
I=- \frac{\left(- \frac{\Gamma\left(\varepsilon + 1\right) \Gamma\left(2 \varepsilon + 1\right)}{\Gamma\left(3 \varepsilon + 1\right)} + 1\right) e^{- 3 \gamma_E \varepsilon} \cosh{\left(\frac{\varepsilon \log{\left(z \right)}}{2} \right)} \Gamma^{3}\left(1 - \varepsilon\right)}{128 \varepsilon^{5} \left(1 - 2 \varepsilon\right) \left(\varepsilon + 1\right) \cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)}}
$$

| ε幂次 | I 的系数 |
|---:|---|
| -3 | $- \frac{\pi^{2}}{384}$ |
| -2 | $- \frac{\pi^{2}}{384} + \frac{3 \zeta\left(3\right)}{64}$ |
| -1 | $\frac{\pi^{2} \log{\left(z \right)}^{2}}{384} - \frac{37 \pi^{4}}{23040} - \frac{\pi^{2}}{128} + \frac{3 \zeta\left(3\right)}{64}$ |
| 0 | $- \frac{3 \log{\left(z \right)}^{2} \zeta\left(3\right)}{64} + \frac{\pi^{2} \log{\left(z \right)}^{2}}{384} - \frac{37 \pi^{4}}{23040} - \frac{5 \pi^{2}}{384} - \frac{5 \pi^{2} \zeta\left(3\right)}{768} + \frac{9 \zeta\left(3\right)}{64} + \frac{21 \zeta\left(5\right)}{64}$ |

## 4e-Pi1-06

[机器结果](../../results/4e-Pi1-06/r002.json)；[现行prepared](../../prepared/4e-Pi1-r002.json)。

$$
I=- \frac{\left(- \frac{\Gamma\left(\varepsilon + 1\right) \Gamma\left(2 \varepsilon + 1\right)}{\Gamma\left(3 \varepsilon + 1\right)} + 1\right) e^{- 3 \gamma_E \varepsilon} \cosh{\left(\frac{\varepsilon \log{\left(z \right)}}{2} \right)} \Gamma^{3}\left(1 - \varepsilon\right)}{128 \varepsilon^{5} \left(1 - 2 \varepsilon\right) \left(\varepsilon + 1\right) \cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)}}
$$

| ε幂次 | I 的系数 |
|---:|---|
| -3 | $- \frac{\pi^{2}}{384}$ |
| -2 | $- \frac{\pi^{2}}{384} + \frac{3 \zeta\left(3\right)}{64}$ |
| -1 | $\frac{\pi^{2} \log{\left(z \right)}^{2}}{384} - \frac{37 \pi^{4}}{23040} - \frac{\pi^{2}}{128} + \frac{3 \zeta\left(3\right)}{64}$ |
| 0 | $- \frac{3 \log{\left(z \right)}^{2} \zeta\left(3\right)}{64} + \frac{\pi^{2} \log{\left(z \right)}^{2}}{384} - \frac{37 \pi^{4}}{23040} - \frac{5 \pi^{2}}{384} - \frac{5 \pi^{2} \zeta\left(3\right)}{768} + \frac{9 \zeta\left(3\right)}{64} + \frac{21 \zeta\left(5\right)}{64}$ |

由旧定义到现行定义的线性变换、颜色和系数连接已核对；未重新积分，不继承新形式证明。formal=not_run。
