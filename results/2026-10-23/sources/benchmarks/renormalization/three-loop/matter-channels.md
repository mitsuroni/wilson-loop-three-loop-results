# 三圈提前闭合的两个颜色通道

2026-10-06更新：C_F C_A T_F n_f也已完成求和、重整化和RG，现有三个物质颜色通道闭合于现行输入范围；纯C_F C_A²仍缺双TGV10图。[新通道结果及连续推导](workspaces/CF-CA-TF-nf/derivation.md) · [当前选择](current.json)。新腿Π₁／Γ₁恢复物理P时需明确的−1色序桥，原I与旧两通道结果保持。下文较早的26缺项或“另两颜色仍缺项”是2026-10-01范围。

这里的“闭合”限定于现行图集和输入约定：未完成顶点腿自能及中央Γ₁的颜色排除已在[颜色覆盖说明](color-coverage.md)中逐步给出。RG通过本身不证明图集完备，也不能排除遗漏的尺度无关有限项。

本页只发布 $C_F(T_Fn_f)^2$ 与 $C_F^2T_Fn_f$，另外两颜色仍缺项。取 $T=T_Fn_f$、$\lambda=\alpha_s/\pi$、$L=(L_s+L_t)/2$、$\rho=\log z$。

$$
\log W_R\big|_{\lambda^3,\mathrm{closed}}
=\lambda^3[C_FT^2F_{3,ff}+C_F^2TF_{3,Ff}].
$$

ff来自双Pi1两图与反项 $2\widehat w_{2,f}/(3\epsilon)+\widehat w_1/(9\epsilon^2)$；低阶乘积至多一个nf，故这一通道W3与w3相同。Ff来自Pi2两图及 $\widehat w_1/(8\epsilon)$；采用已核对的web路线去掉4e-Pi1可约的外部CF平方，保留Pi2内部夸克CF。不能将整个普通CF²nf的W3与w3视为天然相同。

原单图表若以L指log z，已映射为rho；此处L是物理尺度。现行自能2.0及唯一C/I保持。反项后单位颜色极点为

| epsilon幂 | ff | Ff |
|---:|---|---|
| -4 | $- \frac{11}{81}$ | $0$ |
| -3 | $\frac{25}{243} - \frac{2 L}{27}$ | $- \frac{2}{9}$ |
| -2 | $\frac{10 L}{81} - \frac{\pi^{2}}{162} + \frac{10}{81}$ | $- \frac{L}{6} - \frac{2 \zeta\left(3\right)}{9} + \frac{55}{216}$ |
| -1 | $\frac{2 L}{81} - \frac{14 \zeta\left(3\right)}{81} + \frac{5 \pi^{2}}{486} + \frac{260}{2187}$ | $- \frac{2 L \zeta\left(3\right)}{3} + \frac{55 L}{72} - \frac{19 \zeta\left(3\right)}{27} - \frac{\pi^{4}}{270} - \frac{\pi^{2}}{72} + \frac{1711}{1296}$ |

负六次与负五次均为零，Ff负四次也为零。有限结果为

$$
\begin{aligned}
F_{3,ff}={}&\frac{32}{6561} - \frac{232 L}{729}\\
&- \frac{110 \zeta\left(3\right)}{243} - \frac{25 L^{2}}{81}\\
&- \frac{10 L^{3}}{81} - \frac{7 \pi^{4}}{19440}\\
&- \frac{5 M_{1,1}}{27} - \frac{M_{2,1}}{36}\\
&- \frac{L^{4}}{54} + \frac{\rho^{4}}{864}\\
&+ \frac{25 \rho^{2}}{324} + \frac{38 \pi^{2}}{243}\\
&- \frac{4 L \zeta\left(3\right)}{9} - \frac{L M_{1,1}}{9}\\
&+ \frac{\pi^{2} L^{2}}{18} + \frac{L^{2} \rho^{2}}{36}\\
&+ \frac{\pi^{2} \rho^{2}}{216} + \frac{5 \pi^{2} L}{27}\\
&+ \frac{5 L \rho^{2}}{54}
\end{aligned}
$$

$$
\begin{aligned}
F_{3,Ff}={}&\frac{42727}{7776} - \frac{355 \zeta\left(3\right)}{162}\\
&- \frac{55 \pi^{2}}{108} - \frac{55 \rho^{2}}{192}\\
&- \frac{19 \pi^{4}}{1620} - \frac{14 \zeta\left(5\right)}{9}\\
&+ \frac{M_{1,1}}{8} + \frac{L^{3}}{12}\\
&+ \frac{55 L^{2}}{48} + \frac{1711 L}{432}\\
&- L^{2} \zeta\left(3\right) - \frac{19 L \zeta\left(3\right)}{9}\\
&- \frac{\pi^{2} L}{8} - \frac{L \rho^{2}}{16}\\
&- \frac{\pi^{4} L}{90} + \frac{\rho^{2} \zeta\left(3\right)}{4}\\
&+ \frac{4 \pi^{2} \zeta\left(3\right)}{9}
\end{aligned}
$$

$f_{3,c}(\rho)=F_{3,c}(0,\rho)$，机器数据同时存储。$h=\gamma_4/2$，两路提取为

| 通道 | k3 | h3 |
|---|---|---|
| ff | $- \frac{1}{27}$ | $- \frac{260}{729} - \frac{5 \pi^{2}}{162} + \frac{14 \zeta\left(3\right)}{27}$ |
| Ff | $- \frac{55}{48} + \zeta\left(3\right)$ | $- \frac{1711}{432} + \frac{\pi^{2}}{24} + \frac{\pi^{4}}{90} + \frac{19 \zeta\left(3\right)}{9}$ |

有限RG残差严格为零：

$$
\partial_LF_{3,ff}+\frac23F_{2,f}=-2k_{3,ff}L-h_{3,ff},
$$

$$
\partial_LF_{3,Ff}+\frac14F_1=-2k_{3,Ff}L-h_{3,Ff}.
$$

从极点与有限式独立提取的k3/h3一致，crossing通过。cusp两项匹配[GHKM §6.3式(134)](https://arxiv.org/html/1510.07803)；h3和有限形状采用项目坐标方案，未作完整外部对照。裸系数、反项后系数和依赖见[matter-channels-r001.json](matter-channels-r001.json)。范围是两个完整子集，不是完整三圈总和。未新求单图积分或Lean证明，formal=not_run。
