# P4：04/09的完整I与有限项（现行r002）

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1/families/p4-opposite-bare/solutions.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。


以下I已采用物理1PI插入约定，组装系数为+1。G、U和各母函数保持原定义；J及有序域不变。系数表覆盖ε⁻⁶至ε⁰，未显示的位置为零，余项O(ε)。

## 4e-Pi1-04

[机器结果](../../results/4e-Pi1-04/r002.json)；[现行prepared](../../prepared/4e-Pi1-r002.json)。

$$
I=- \frac{\left(\frac{\left(\frac{z^{\frac{\varepsilon}{2}} \left(-1 - \frac{1}{z}\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, \varepsilon, - 2 \varepsilon &  \\0, \varepsilon, \varepsilon &  \end{matrix} \middle| {\frac{1}{z}} \right)}}{\Gamma\left(1 - \varepsilon\right) \Gamma\left(\varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)} - \frac{z^{- \frac{\varepsilon}{2}} \left(- z - 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, \varepsilon, - 2 \varepsilon &  \\0, \varepsilon, \varepsilon &  \end{matrix} \middle| {z} \right)}}{\Gamma\left(1 - \varepsilon\right) \Gamma\left(\varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)}\right) \sinh{\left(\varepsilon \log{\left(z \right)} \right)}}{2} + \left(- \frac{z^{\frac{\varepsilon}{2}} \left(-1 - \frac{1}{z}\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, \varepsilon, - 2 \varepsilon &  \\0, \varepsilon, \varepsilon &  \end{matrix} \middle| {\frac{1}{z}} \right)}}{2 \Gamma\left(1 - \varepsilon\right) \Gamma\left(\varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)} + \frac{z^{- \frac{\varepsilon}{2}} \left(- z - 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, 0, \varepsilon &  \\0, \varepsilon, \varepsilon &  \end{matrix} \middle| {z} \right)}}{\Gamma\left(1 - \varepsilon\right) \Gamma^{2}\left(\varepsilon + 1\right)} - \frac{z^{- \frac{\varepsilon}{2}} \left(- z - 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, \varepsilon, - 2 \varepsilon &  \\0, \varepsilon, \varepsilon &  \end{matrix} \middle| {z} \right)}}{2 \Gamma\left(1 - \varepsilon\right) \Gamma\left(\varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)}\right) \cosh{\left(\varepsilon \log{\left(z \right)} \right)}\right) e^{- 3 \gamma_E \varepsilon} \Gamma^{3}\left(1 - \varepsilon\right)}{128 \varepsilon^{3} \left(1 - 2 \varepsilon\right) \left(\varepsilon + 1\right) \cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)}}
$$

| ε幂次 | I 的系数 |
|---:|---|
| -2 | $- \frac{tau0q1}{128} + \frac{tau2q1}{128}$ |
| -1 | $\frac{eta2q1}{128} - \frac{tau0q1}{128} - \frac{tau0q2}{128} + \frac{tau2q1}{128} + \frac{tau2q2}{128}$ |
| 0 | $\frac{5 L^{2} tau0q1}{1024} - \frac{5 L^{2} tau2q1}{1024} + \frac{eta2q1}{128} + \frac{eta2q2}{128} - \frac{3 tau0q1}{128} - \frac{\pi^{2} tau0q1}{512} - \frac{tau0q2}{128} - \frac{tau0q3}{128} + \frac{\pi^{2} tau2q1}{512} + \frac{3 tau2q1}{128} + \frac{tau2q2}{128} + \frac{tau2q3}{128}$ |

## 4e-Pi1-09

[机器结果](../../results/4e-Pi1-09/r002.json)；[现行prepared](../../prepared/4e-Pi1-r002.json)。

$$
I=- \frac{\left(\frac{\left(\frac{z^{\frac{\varepsilon}{2}} \left(-1 - \frac{1}{z}\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, \varepsilon, - 2 \varepsilon &  \\0, \varepsilon, \varepsilon &  \end{matrix} \middle| {\frac{1}{z}} \right)}}{\Gamma\left(1 - \varepsilon\right) \Gamma\left(\varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)} - \frac{z^{- \frac{\varepsilon}{2}} \left(- z - 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, \varepsilon, - 2 \varepsilon &  \\0, \varepsilon, \varepsilon &  \end{matrix} \middle| {z} \right)}}{\Gamma\left(1 - \varepsilon\right) \Gamma\left(\varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)}\right) \sinh{\left(\varepsilon \log{\left(z \right)} \right)}}{2} + \left(- \frac{z^{\frac{\varepsilon}{2}} \left(-1 - \frac{1}{z}\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, \varepsilon, - 2 \varepsilon &  \\0, \varepsilon, \varepsilon &  \end{matrix} \middle| {\frac{1}{z}} \right)}}{2 \Gamma\left(1 - \varepsilon\right) \Gamma\left(\varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)} + \frac{z^{- \frac{\varepsilon}{2}} \left(- z - 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, 0, \varepsilon &  \\0, \varepsilon, \varepsilon &  \end{matrix} \middle| {z} \right)}}{\Gamma\left(1 - \varepsilon\right) \Gamma^{2}\left(\varepsilon + 1\right)} - \frac{z^{- \frac{\varepsilon}{2}} \left(- z - 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, \varepsilon, - 2 \varepsilon &  \\0, \varepsilon, \varepsilon &  \end{matrix} \middle| {z} \right)}}{2 \Gamma\left(1 - \varepsilon\right) \Gamma\left(\varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)}\right) \cosh{\left(\varepsilon \log{\left(z \right)} \right)}\right) e^{- 3 \gamma_E \varepsilon} \Gamma^{3}\left(1 - \varepsilon\right)}{128 \varepsilon^{3} \left(1 - 2 \varepsilon\right) \left(\varepsilon + 1\right) \cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)}}
$$

| ε幂次 | I 的系数 |
|---:|---|
| -2 | $- \frac{tau0q1}{128} + \frac{tau2q1}{128}$ |
| -1 | $\frac{eta2q1}{128} - \frac{tau0q1}{128} - \frac{tau0q2}{128} + \frac{tau2q1}{128} + \frac{tau2q2}{128}$ |
| 0 | $\frac{5 L^{2} tau0q1}{1024} - \frac{5 L^{2} tau2q1}{1024} + \frac{eta2q1}{128} + \frac{eta2q2}{128} - \frac{3 tau0q1}{128} - \frac{\pi^{2} tau0q1}{512} - \frac{tau0q2}{128} - \frac{tau0q3}{128} + \frac{\pi^{2} tau2q1}{512} + \frac{3 tau2q1}{128} + \frac{tau2q2}{128} + \frac{tau2q3}{128}$ |

由旧定义到现行定义的线性变换、颜色和系数连接已核对；未重新积分，不继承新形式证明。formal=not_run。
