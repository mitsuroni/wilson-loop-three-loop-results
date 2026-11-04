# CH_O：2e-Pi1+Pi1-02 完整I与有限项（现行r002）

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1-chain/families/opposite/solutions.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。


以下I已采用物理1PI插入约定，组装系数为+1。G、U和各母函数保持原定义；J及有序域不变。系数表覆盖ε⁻⁶至ε⁰，未显示的位置为零，余项O(ε)。

## 2e-Pi1+Pi1-02

[机器结果](../../results/2e-Pi1+Pi1-02/r002.json)；[现行prepared](../../prepared/2e-Pi1+Pi1-r002.json)。

$$
I=\frac{z^{- \frac{3 \varepsilon}{2}} \left(z + 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, 0, 3 \varepsilon &  \\0, 3 \varepsilon, 3 \varepsilon &  \end{matrix} \middle| {z} \right)} e^{- 3 \gamma_E \varepsilon} \Gamma^{4}\left(1 - \varepsilon\right) \Gamma^{2}\left(\varepsilon + 1\right)}{128 \varepsilon^{2} \left(1 - 2 \varepsilon\right)^{2} \cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)} \Gamma^{2}\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right) \Gamma^{2}\left(3 \varepsilon + 1\right)} - \frac{e^{- 3 \gamma_E \varepsilon} \Gamma\left(1 - 3 \varepsilon\right) \Gamma^{4}\left(1 - \varepsilon\right) \Gamma^{2}\left(\varepsilon + 1\right)}{384 \varepsilon^{3} \left(1 - 2 \varepsilon\right)^{2} \left(2 \varepsilon + 1\right) \Gamma^{2}\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right)}
$$

| ε幂次 | I 的系数 |
|---:|---|
| -3 | $- \frac{1}{384}$ |
| -2 | $- \frac{M_{01}}{128} - \frac{1}{192}$ |
| -1 | $- \frac{M_{01}}{32} - \frac{3 M_{11}}{256} - \frac{1}{48} - \frac{\pi^{2}}{1536}$ |
| 0 | $\frac{9 L^{2} M_{01}}{1024} - \frac{3 M_{01}}{32} - \frac{\pi^{2} M_{01}}{512} - \frac{3 M_{03}}{1024} - \frac{3 M_{11}}{64} - \frac{9 M_{21}}{1024} - \frac{1}{24} - \frac{7 \zeta\left(3\right)}{384} - \frac{\pi^{2}}{768}$ |

### 颜色 CA^2*CF

| ε幂次 | 完整颜色系数 |
|---:|---|
| -3 | $- \frac{25}{864}$ |
| -2 | $- \frac{25 M_{01}}{288} - \frac{5}{81}$ |
| -1 | $- \frac{155 M_{01}}{432} - \frac{25 M_{11}}{192} - \frac{209}{864} - \frac{25 \pi^{2}}{3456}$ |
| 0 | $\frac{25 L^{2} M_{01}}{256} - \frac{947 M_{01}}{864} - \frac{25 \pi^{2} M_{01}}{1152} - \frac{25 M_{03}}{768} - \frac{155 M_{11}}{288} - \frac{25 M_{21}}{256} - \frac{5845}{11664} - \frac{175 \zeta\left(3\right)}{864} - \frac{5 \pi^{2}}{324}$ |

### 颜色 CA*CF*TF*nf

| ε幂次 | 完整颜色系数 |
|---:|---|
| -3 | $\frac{5}{108}$ |
| -2 | $\frac{5 M_{01}}{36} + \frac{13}{162}$ |
| -1 | $\frac{14 M_{01}}{27} + \frac{5 M_{11}}{24} + \frac{5 \pi^{2}}{432} + \frac{109}{324}$ |
| 0 | $- \frac{5 L^{2} M_{01}}{32} + \frac{5 \pi^{2} M_{01}}{144} + \frac{161 M_{01}}{108} + \frac{5 M_{03}}{96} + \frac{7 M_{11}}{9} + \frac{5 M_{21}}{32} + \frac{13 \pi^{2}}{648} + \frac{35 \zeta\left(3\right)}{108} + \frac{899}{1458}$ |

### 颜色 CF*(TF*nf)^2

| ε幂次 | 完整颜色系数 |
|---:|---|
| -3 | $- \frac{1}{54}$ |
| -2 | $- \frac{M_{01}}{18} - \frac{2}{81}$ |
| -1 | $- \frac{5 M_{01}}{27} - \frac{M_{11}}{12} - \frac{19}{162} - \frac{\pi^{2}}{216}$ |
| 0 | $\frac{L^{2} M_{01}}{16} - \frac{M_{01}}{2} - \frac{\pi^{2} M_{01}}{72} - \frac{M_{03}}{48} - \frac{5 M_{11}}{18} - \frac{M_{21}}{16} - \frac{133}{729} - \frac{7 \zeta\left(3\right)}{54} - \frac{\pi^{2}}{162}$ |

由旧定义到现行定义的线性变换、颜色和系数连接已核对；未重新积分，不继承新形式证明。formal=not_run。
