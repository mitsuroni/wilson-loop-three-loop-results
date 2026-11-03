# 2e-Pi2-02：对边2/4，三个颜色的完整I至有限项（现行r002）

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi2/families/opposite/solutions.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。


以下I已采用物理1PI插入约定，组装系数为+1。G、U和各母函数保持原定义；J及有序域不变。系数表覆盖ε⁻⁶至ε⁰，未显示的位置为零，余项O(ε)。

## 2e-Pi2-02

[机器结果](../../results/2e-Pi2-02/r002.json)；[现行prepared](../../prepared/2e-Pi2-r002.json)。

$$
I=- \frac{z^{- \frac{3 \varepsilon}{2}} \left(z + 1\right) {G_{3, 3}^{3, 3}\left(\begin{matrix} 0, 0, 3 \varepsilon &  \\0, 3 \varepsilon, 3 \varepsilon &  \end{matrix} \middle| {z} \right)} e^{- 3 \gamma_E \varepsilon} \Gamma^{4}\left(1 - \varepsilon\right) \Gamma^{2}\left(\varepsilon + 1\right)}{128 \varepsilon^{2} \left(1 - 2 \varepsilon\right)^{2} \cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)} \Gamma^{2}\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right) \Gamma^{2}\left(3 \varepsilon + 1\right)} + \frac{e^{- 3 \gamma_E \varepsilon} \Gamma\left(1 - 3 \varepsilon\right) \Gamma^{4}\left(1 - \varepsilon\right) \Gamma^{2}\left(\varepsilon + 1\right)}{384 \varepsilon^{3} \left(1 - 2 \varepsilon\right)^{2} \left(2 \varepsilon + 1\right) \Gamma^{2}\left(1 - 2 \varepsilon\right) \Gamma\left(2 \varepsilon + 1\right)}
$$

| ε幂次 | I 的系数 |
|---:|---|
| -3 | $\frac{1}{384}$ |
| -2 | $\frac{M_{01}}{128} + \frac{1}{192}$ |
| -1 | $\frac{M_{01}}{32} + \frac{3 M_{11}}{256} + \frac{\pi^{2}}{1536} + \frac{1}{48}$ |
| 0 | $- \frac{9 L^{2} M_{01}}{1024} + \frac{\pi^{2} M_{01}}{512} + \frac{3 M_{01}}{32} + \frac{3 M_{03}}{1024} + \frac{3 M_{11}}{64} + \frac{9 M_{21}}{1024} + \frac{\pi^{2}}{768} + \frac{7 \zeta\left(3\right)}{384} + \frac{1}{24}$ |

### 颜色 CA^2*CF

| ε幂次 | 完整颜色系数 |
|---:|---|
| -3 | $- \frac{25}{1152}$ |
| -2 | $- \frac{25 M_{01}}{384} - \frac{283}{6912}$ |
| -1 | $- \frac{583 M_{01}}{2304} - \frac{25 M_{11}}{256} - \frac{10915}{41472} - \frac{25 \pi^{2}}{4608} + \frac{\zeta\left(3\right)}{96}$ |
| 0 | $\frac{75 L^{2} M_{01}}{1024} - \frac{14311 M_{01}}{13824} - \frac{25 \pi^{2} M_{01}}{1536} + \frac{M_{01} \zeta\left(3\right)}{32} - \frac{25 M_{03}}{1024} - \frac{583 M_{11}}{1536} - \frac{75 M_{21}}{1024} - \frac{177067}{248832} - \frac{283 \pi^{2}}{27648} + \frac{\pi^{4}}{5760} + \frac{55 \zeta\left(3\right)}{1152}$ |

### 颜色 CA*CF*TF*nf

| ε幂次 | 完整颜色系数 |
|---:|---|
| -3 | $\frac{5}{288}$ |
| -2 | $\frac{5 M_{01}}{96} + \frac{41}{1728}$ |
| -1 | $\frac{101 M_{01}}{576} + \frac{5 M_{11}}{64} + \frac{5 \pi^{2}}{1152} + \frac{\zeta\left(3\right)}{12} + \frac{1469}{10368}$ |
| 0 | $- \frac{15 L^{2} M_{01}}{256} + \frac{5 \pi^{2} M_{01}}{384} + \frac{M_{01} \zeta\left(3\right)}{4} + \frac{1961 M_{01}}{3456} + \frac{5 M_{03}}{256} + \frac{101 M_{11}}{384} + \frac{15 M_{21}}{256} + \frac{\zeta\left(3\right)}{32} + \frac{41 \pi^{2}}{6912} + \frac{\pi^{4}}{720} + \frac{14597}{62208}$ |

### 颜色 CF^2*TF*nf

| ε幂次 | 完整颜色系数 |
|---:|---|
| -2 | $\frac{1}{48}$ |
| -1 | $\frac{M_{01}}{16} - \frac{\zeta\left(3\right)}{6} + \frac{43}{288}$ |
| 0 | $- \frac{M_{01} \zeta\left(3\right)}{2} + \frac{55 M_{01}}{96} + \frac{3 M_{11}}{32} - \frac{\pi^{4}}{360} - \frac{7 \zeta\left(3\right)}{36} + \frac{\pi^{2}}{192} + \frac{1195}{1728}$ |

由旧定义到现行定义的线性变换、颜色和系数连接已核对；未重新积分，不继承新形式证明。formal=not_run。
