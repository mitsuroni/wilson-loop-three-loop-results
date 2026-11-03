# P6：4e-Pi1-11 的完整 I 与有限项（现行r002）

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1/families/p6-opposite-distinct/solutions.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。


以下I已采用物理1PI插入约定，组装系数为+1。G、U和各母函数保持原定义；J及有序域不变。系数表覆盖ε⁻⁶至ε⁰，未显示的位置为零，余项O(ε)。

## 4e-Pi1-11

[机器结果](../../results/4e-Pi1-11/r002.json)；[现行prepared](../../prepared/4e-Pi1-r002.json)。

$$
I=- \frac{z^{- \frac{\varepsilon}{2}} \left(- z - 1\right) \left(\frac{z^{- \varepsilon} \left(- z - 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, 0, 2 \varepsilon &  \\0, 2 \varepsilon, 2 \varepsilon &  \end{matrix} \middle| {z} \right)}}{\Gamma\left(1 - 2 \varepsilon\right) \Gamma^{2}\left(2 \varepsilon + 1\right)} + \frac{\cosh{\left(\varepsilon \log{\left(z \right)} \right)}}{2 \varepsilon \left(\varepsilon + 1\right)}\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, 0, \varepsilon &  \\0, \varepsilon, \varepsilon &  \end{matrix} \middle| {z} \right)} e^{- 3 \gamma_E \varepsilon} \Gamma^{2}\left(1 - \varepsilon\right)}{64 \varepsilon \left(1 - 2 \varepsilon\right) \cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)} \Gamma^{2}\left(\varepsilon + 1\right)}
$$

| ε幂次 | I 的系数 |
|---:|---|
| -2 | $- \frac{T_{0}}{128}$ |
| -1 | $- \frac{T_{0}^{2}}{64} - \frac{T_{0}}{128} - \frac{T_{1}}{128}$ |
| 0 | $\frac{5 L^{2} T_{0}}{1024} - \frac{T_{0}^{2}}{32} - \frac{3 T_{0} T_{1}}{64} - \frac{3 T_{0}}{128} - \frac{\pi^{2} T_{0}}{512} - \frac{T_{1}}{128} - \frac{T_{2}}{128}$ |

由旧定义到现行定义的线性变换、颜色和系数连接已核对；未重新积分，不继承新形式证明。formal=not_run。
