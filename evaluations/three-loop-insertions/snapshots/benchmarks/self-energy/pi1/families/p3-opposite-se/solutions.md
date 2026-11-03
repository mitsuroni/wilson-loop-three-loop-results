# P3：03/10的完整I与有限项（现行r002）

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1/families/p3-opposite-se/solutions.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。


以下I已采用物理1PI插入约定，组装系数为+1。G、U和各母函数保持原定义；J及有序域不变。系数表覆盖ε⁻⁶至ε⁰，未显示的位置为零，余项O(ε)。

## 4e-Pi1-03

[机器结果](../../results/4e-Pi1-03/r002.json)；[现行prepared](../../prepared/4e-Pi1-r002.json)。

$$
I=- \frac{\left(\frac{\left(\frac{z^{\varepsilon} \left(-1 - \frac{1}{z}\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, - \varepsilon, 2 \varepsilon &  \\0, 2 \varepsilon, 2 \varepsilon &  \end{matrix} \middle| {\frac{1}{z}} \right)}}{\Gamma\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)} - \frac{z^{- \varepsilon} \left(- z - 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, - \varepsilon, 2 \varepsilon &  \\0, 2 \varepsilon, 2 \varepsilon &  \end{matrix} \middle| {z} \right)}}{\Gamma\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)}\right) \sinh{\left(\frac{\varepsilon \log{\left(z \right)}}{2} \right)}}{2} + \left(- \frac{z^{\varepsilon} \left(-1 - \frac{1}{z}\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, - \varepsilon, 2 \varepsilon &  \\0, 2 \varepsilon, 2 \varepsilon &  \end{matrix} \middle| {\frac{1}{z}} \right)}}{2 \Gamma\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)} + \frac{z^{- \varepsilon} \left(- z - 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, 0, 2 \varepsilon &  \\0, 2 \varepsilon, 2 \varepsilon &  \end{matrix} \middle| {z} \right)}}{\Gamma\left(1 - 2 \varepsilon\right) \Gamma^{2}\left(2 \varepsilon + 1\right)} - \frac{z^{- \varepsilon} \left(- z - 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, - \varepsilon, 2 \varepsilon &  \\0, 2 \varepsilon, 2 \varepsilon &  \end{matrix} \middle| {z} \right)}}{2 \Gamma\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)}\right) \cosh{\left(\frac{\varepsilon \log{\left(z \right)}}{2} \right)}\right) e^{- 3 \gamma_E \varepsilon} \Gamma^{3}\left(1 - \varepsilon\right)}{32 \varepsilon^{3} \left(1 - 2 \varepsilon\right) \cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)}} - \frac{\left(\frac{\left(- \frac{3 \Gamma\left(\varepsilon + 1\right) \Gamma\left(2 \varepsilon + 1\right)}{\Gamma\left(3 \varepsilon + 1\right)} + 3\right) \cosh{\left(\frac{\varepsilon \log{\left(z \right)}}{2} \right)}}{\cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)}} + 1\right) e^{- 3 \gamma_E \varepsilon} \Gamma^{3}\left(1 - \varepsilon\right)}{384 \varepsilon^{4} \left(1 - 2 \varepsilon\right) \left(\varepsilon + 1\right)}
$$

| ε幂次 | I 的系数 |
|---:|---|
| -4 | $- \frac{1}{384}$ |
| -3 | $- \frac{1}{384}$ |
| -2 | $\frac{M_{11}}{64} - \frac{5 \pi^{2}}{1536} - \frac{1}{128}$ |
| -1 | $- \frac{L^{2} M_{01}}{256} + \frac{M_{03}}{256} + \frac{M_{11}}{32} + \frac{5 M_{21}}{256} - \frac{5 \pi^{2}}{1536} - \frac{5}{384} + \frac{17 \zeta\left(3\right)}{384}$ |
| 0 | $- \frac{L^{2} M_{01}}{128} - \frac{11 L^{2} M_{11}}{512} + \frac{\pi^{2} L^{2}}{384} + \frac{M_{03}}{128} + \frac{\pi^{2} M_{11}}{256} + \frac{M_{11}}{16} + \frac{13 M_{13}}{1536} + \frac{5 M_{21}}{128} + \frac{19 M_{31}}{1536} - \frac{7 \pi^{4}}{4096} - \frac{5 \pi^{2}}{512} - \frac{11}{384} + \frac{17 \zeta\left(3\right)}{384}$ |

## 4e-Pi1-10

[机器结果](../../results/4e-Pi1-10/r002.json)；[现行prepared](../../prepared/4e-Pi1-r002.json)。

$$
I=- \frac{\left(\frac{\left(\frac{z^{\varepsilon} \left(-1 - \frac{1}{z}\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, - \varepsilon, 2 \varepsilon &  \\0, 2 \varepsilon, 2 \varepsilon &  \end{matrix} \middle| {\frac{1}{z}} \right)}}{\Gamma\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)} - \frac{z^{- \varepsilon} \left(- z - 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, - \varepsilon, 2 \varepsilon &  \\0, 2 \varepsilon, 2 \varepsilon &  \end{matrix} \middle| {z} \right)}}{\Gamma\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)}\right) \sinh{\left(\frac{\varepsilon \log{\left(z \right)}}{2} \right)}}{2} + \left(- \frac{z^{\varepsilon} \left(-1 - \frac{1}{z}\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, - \varepsilon, 2 \varepsilon &  \\0, 2 \varepsilon, 2 \varepsilon &  \end{matrix} \middle| {\frac{1}{z}} \right)}}{2 \Gamma\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)} + \frac{z^{- \varepsilon} \left(- z - 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, 0, 2 \varepsilon &  \\0, 2 \varepsilon, 2 \varepsilon &  \end{matrix} \middle| {z} \right)}}{\Gamma\left(1 - 2 \varepsilon\right) \Gamma^{2}\left(2 \varepsilon + 1\right)} - \frac{z^{- \varepsilon} \left(- z - 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, - \varepsilon, 2 \varepsilon &  \\0, 2 \varepsilon, 2 \varepsilon &  \end{matrix} \middle| {z} \right)}}{2 \Gamma\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right) \Gamma\left(3 \varepsilon + 1\right)}\right) \cosh{\left(\frac{\varepsilon \log{\left(z \right)}}{2} \right)}\right) e^{- 3 \gamma_E \varepsilon} \Gamma^{3}\left(1 - \varepsilon\right)}{32 \varepsilon^{3} \left(1 - 2 \varepsilon\right) \cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)}} - \frac{\left(\frac{\left(- \frac{3 \Gamma\left(\varepsilon + 1\right) \Gamma\left(2 \varepsilon + 1\right)}{\Gamma\left(3 \varepsilon + 1\right)} + 3\right) \cosh{\left(\frac{\varepsilon \log{\left(z \right)}}{2} \right)}}{\cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)}} + 1\right) e^{- 3 \gamma_E \varepsilon} \Gamma^{3}\left(1 - \varepsilon\right)}{384 \varepsilon^{4} \left(1 - 2 \varepsilon\right) \left(\varepsilon + 1\right)}
$$

| ε幂次 | I 的系数 |
|---:|---|
| -4 | $- \frac{1}{384}$ |
| -3 | $- \frac{1}{384}$ |
| -2 | $\frac{M_{11}}{64} - \frac{5 \pi^{2}}{1536} - \frac{1}{128}$ |
| -1 | $- \frac{L^{2} M_{01}}{256} + \frac{M_{03}}{256} + \frac{M_{11}}{32} + \frac{5 M_{21}}{256} - \frac{5 \pi^{2}}{1536} - \frac{5}{384} + \frac{17 \zeta\left(3\right)}{384}$ |
| 0 | $- \frac{L^{2} M_{01}}{128} - \frac{11 L^{2} M_{11}}{512} + \frac{\pi^{2} L^{2}}{384} + \frac{M_{03}}{128} + \frac{\pi^{2} M_{11}}{256} + \frac{M_{11}}{16} + \frac{13 M_{13}}{1536} + \frac{5 M_{21}}{128} + \frac{19 M_{31}}{1536} - \frac{7 \pi^{4}}{4096} - \frac{5 \pi^{2}}{512} - \frac{11}{384} + \frac{17 \zeta\left(3\right)}{384}$ |

由旧定义到现行定义的线性变换、颜色和系数连接已核对；未重新积分，不继承新形式证明。formal=not_run。
