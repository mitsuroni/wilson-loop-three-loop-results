# 单个二圈 1PI 胶子自能的完整参数化

日期：2026-09-24。本文为 `2e-Pi2-01/02` 建立 montopar 理论基础。原动量图贡献记为 $\mathcal A=C^{\mathrm{src}}J$，归一化后的 $I$ 只按第7节的同一个系数提取规则定义。本文按三个颜色通道共享一个 I 的约定统一；当前C与r001相反，J保持不变。生产实现和实际输出见 [C/J交付](cj-r002.md)，现行prepared与积分结果均已发布为r002。

范围为普通协变 Feynman 规范、无质量夸克、裸二圈 **1PI** 自能；总 Wilson 圈数为3。两个 Π₁ 串联是另一个分支，不放进 Π₂。自能内部的非阿贝尔相互作用完整保留，不增加外接三胶子顶点骨架，也不加入重整化反项。

现行数学约定为[自能2.0](../conventions.md)，输入和结果见[当前清单](../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi2/theory.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。

## 1. 两条规范化源与唯一动量贡献

以[两条源记录及原公式摘录](sources/2e-Pi2-r002.json)为起点。令 $D=4-2e$、$e=\varepsilon$，采用 $(+---)$ 度规、归一化颜色迹 $\mathrm{Tr}/N_c$：

$$
\begin{aligned}
\mathcal A_{ij}={}&\int_0^1ds_i\int_0^1ds_j\,
(igv_i^{\mu_1}t^{a_1})(igv_j^{\mu_2}t^{a_2})\\
&\times\int\frac{d^D\ell}{(2\pi)^D}
\frac{-ig_{\mu_1\nu_1}\delta^{a_1b_1}}{\ell^2+i0}
\Pi^{(2)\nu_1\nu_2}_{b_1b_2}(\ell)
\frac{-ig_{\nu_2\mu_2}\delta^{b_2a_2}}{\ell^2+i0}
e^{-i\ell\cdot X_{ij}},\qquad X_{ij}=x_i(s_i)-x_j(s_j).
\end{aligned}
\tag{P2.1}
$$

当前Π₂是含i的1PI插入，图式不再额外乘外负号；两侧颜色 δ 在源 Markdown 中省写，但由结构化记录明确给出，式(P2.1)补全它们。源公式原先的图标签 $I_{(\cdot)}$ 在本文统一改记 $\mathcal A$。末尾沿用一圈模板的 SE 图名注释不作为二圈子图定义；实际身份由 `vacuum-pol-2loop` 和下述文献算符确定。

| 图 | 原外到内顺序 | $X$ | 总阶数 | C4轨道大小 $m$ |
|---|---|---|---|---|
| 01，相邻3/4 | $s_{31}\to s_{41}$ | $-[(1-s_{31})v_3+s_{41}v_4]$ | $g^6$ | 4 |
| 02，对边2/4 | $s_{21}\to s_{41}$ | $(1-s_{41})v_4+v_1+s_{21}v_2$ | $g^6$ | 2 |

每个参数上下限为0、1，Jacobian为1。JSON 中 $\Phi=-X$、指数为 $e^{+i\ell\cdot\Phi}$。一个有效插入、两个传播子、同一个交换动量；算符内部已有两圈及 $g^4$，故显式只剩一个动量积分并不表示总一圈。无额外 $1/2!$。源元数据 `n_diagrams=1` 不是 C4 重数。

## 2. 从文献标量到插入算符

采用[已核对的来源与约定字典](../../../workflow/self-energy-literature.md)：Ruijl–Ueda–Vermaseren–Vogt 的[官方 gluonB.h](https://arxiv.org/src/1703.08532v2/anc/gluonB.h)中 $a^2$ 块，Feynman 规范取作者的 `xi=1`。现有30项来源核验包括与1998年裸展开的跨文献匹配；本轮复用该记录，不声称独立重新计算了二圈自能。

设 $T_f=T_Fn_f$，$T_F=1/2$；$T^{\mu\nu}=\ell^2g^{\mu\nu}-\ell^\mu\ell^\nu$，$q=-\ell^2-i0$。文献的逆传播子标量与项目的插入约定连接为

$$
\begin{aligned}
\eta(e)&=\frac{\Gamma(1+e)\Gamma(1-e)^2}{\Gamma(1-2e)},\\
\mathcal P_2(\ell^2)&=\frac{g^4}{(4\pi)^D}
\left[\frac{\eta(e)}{1-2e}\right]^2P_2(e)q^{-2e},\\
\Pi^{(2)\mu\nu}_{ab}&=-i\delta_{ab}T^{\mu\nu}\mathcal P_2
=i\delta_{ab}\mathcal A_2(e)T^{\mu\nu}q^{-2e},\\
\mathcal A_2(e)&=-\frac{g^4}{(4\pi)^D}
\left[\frac{\eta(e)}{1-2e}\right]^2P_2(e).
\end{aligned}
\tag{P2.2}
$$

同一字典在一圈给出 $P_1=-\mathcal B_1/(2e)$、$\Pi^{(1)}=-iT\mathcal P_1$，精确恢复[Π₁主推导](../pi1/theory.md)的正 $i\mathcal A_1T$。这固定的是统一插入约定，不能在二圈另加一个符号来拟合答案。

定义

$$
\begin{aligned}
R(e)&=\frac{\Gamma(1+2e)\Gamma(1-2e)^2}
{\Gamma(1+e)^2\Gamma(1-e)\Gamma(1-3e)},\qquad U(e)=R(e)^{-1},\\
P_2(e)&=C_A^2p_{AA}(e)+C_AT_fp_{Af}(e)+C_FT_fp_{Ff}(e),\\
p_c(e)&=A_c(e)+B_c(e)R(e),\qquad c\in\{AA,Af,Ff\}.
\end{aligned}
\tag{P2.3}
$$

六个有理函数完整列为

$$
\begin{aligned}
A_{AA}&=\frac{2e^2-4e-1}{2e^3(2e-3)},\\
B_{AA}&=\frac{(2e-1)(6e^5-e^4-44e^3+79e^2-30e-2)}
{2e^3(e+1)(2e-3)(3e-2)(3e-1)},\\
A_{Af}&=\frac{2(e-1)(2e^2-3e+2)}{e^3(2e-3)},\\
B_{Af}&=-\frac{2(e-1)(2e-1)(2e^4-11e^3-10e^2+17e-4)}
{e^3(e+1)(2e-3)(3e-2)(3e-1)},\\
A_{Ff}&=-\frac{4(e-1)(2e^2-e+2)}{e^3(2e-3)},\\
B_{Ff}&=\frac{8(e-1)(2e-1)^2(e^2-2e+2)}
{e^3(2e-3)(3e-2)(3e-1)}.
\end{aligned}
\tag{P2.4}
$$

$R$ 是官方 `G10` 的精确 Gamma 表示，恰为串联分支的 $U$ 的倒数。它不是未求主积分。本式无 $(T_Fn_f)^2$，不能从串联的 $\mathcal B_1^2$ 猜出 $P_2$。

## 3. 完整收缩与 Fourier 变换

颜色收缩给出 $C_F$；Lorentz 收缩为 $\ell^2(v_i\cdot v_j)-(v_i\cdot\ell)(v_j\cdot\ell)$。式(P2.1)的显式相位乘积为 $(ig)^2(-i)^2i=+ig^2$，再将 $\ell^2$ 化为 $-q$，得到

$$
\mathcal A_{ij}=-ig^2C_F\mathcal A_2
\int ds_i\,ds_j\int\frac{d^D\ell}{(2\pi)^D}e^{-i\ell\cdot X}
\left[(v_i\cdot v_j)q^{-1-2e}+(v_i\cdot\ell)(v_j\cdot\ell)q^{-2-2e}\right].
\tag{P2.5}
$$

传播子相消按共同 Feynman 边界值及维数解析延拓解释，不把有限 $i\eta$ 的分母直接当成无处方有理式。沿用 $Q=-X^2+i0$ 的母式：

$$
\begin{aligned}
\mathcal F_a(X)&=\int\frac{d^D\ell}{(2\pi)^D}e^{-i\ell\cdot X}q^{-a}
=i\frac{\Gamma(D/2-a)}{4^a\pi^{D/2}\Gamma(a)}Q^{a-D/2},\\
\mathcal F_{1+2e}&=ih_1Q^{-1+3e},\quad
h_1=\frac{\Gamma(1-3e)}{4^{1+2e}\pi^{2-e}\Gamma(1+2e)},\\
\mathcal F_{2+2e}&=ih_2Q^{3e},\quad
h_2=\frac{\Gamma(-3e)}{4^{2+2e}\pi^{2-e}\Gamma(2+2e)},\quad
\frac{h_2}{h_1}=-\frac1{12e(1+2e)}.
\end{aligned}
\tag{P2.6}
$$

动量双乘积对应负的二阶坐标导数。直接使用

$$
\partial_\mu\partial_\nu Q^{3e}
=-6e\,g_{\mu\nu}Q^{-1+3e}
+12e(3e-1)X_\mu X_\nu Q^{-2+3e}
\tag{P2.7}
$$

得到完整坐标贡献

$$
\boxed{\mathcal A_{ij}=+g^2C_F\mathcal A_2h_1\int ds_i\,ds_j
\left[\frac{1+4e}{2(1+2e)}(v_i\cdot v_j)Q^{-1+3e}
-\frac{1-3e}{1+2e}(v_i\cdot X)(v_j\cdot X)Q^{-2+3e}\right].}
\tag{P2.8}
$$

这里的正号来自式(P2.5)的 $-i$ 与 Fourier 的 $i$。两个张量系数完整保留，后续不再补任何“横向修正”。D维迹检查为 $(1+4e)D/[2(1+2e)]+(1-3e)/(1+2e)=D-1$；不提前取 $D=4$。

## 4. 无量纲核与源 C

取 $r=v_{14}<0$、$z=v_{12}/r>0$、$a=-2r>0$，定义 $f=X^2/(2r)$、$d=(v_i\cdot v_j)/r$、$u=(v_i\cdot X)/r$、$w=(v_j\cdot X)/r$。原几何仍由 $v_i^2=0$、$\sum_i v_i=0$ 给出：

$$
\mathsf G=\left(\frac{v_i\cdot v_j}{r}\right)=
\begin{pmatrix}0&z&-(1+z)&1\\z&0&1&-(1+z)\\-(1+z)&1&0&z\\1&-(1+z)&z&0\end{pmatrix},\qquad
f=\tfrac12\boldsymbol c^T\mathsf G\boldsymbol c\quad(X=\sum c_iv_i).
\tag{P2.9}
$$

由 $Q=af+i0$、$r/a=-1/2$，式(P2.8)成为

$$
\begin{aligned}
\mathcal A_{ij}&=C^{\rm src}_{\Pi_2}J_{ij},\qquad
C^{\rm src}_{\Pi_2}=-\frac{g^2C_F\mathcal A_2h_1}{2}(a+i0)^{3e},\\
J_{ij}&=\int_0^1ds_i\int_0^1ds_j\,\mathcal K_2,\\
\mathcal K_2&=\boxed{\frac{(1+4e)df+(1-3e)uw}{2(1+2e)}(f+i0)^{-2+3e}.}
\end{aligned}
\tag{P2.10}
$$

此核与串联的核完全相同，原因是两者都含 $q^{-2e}$ 乘一个完整 $T$；前因子从各自算符推导，不能以核相同替代符号检查。

用 $g^2=4\pi\alpha_s$（裸耦合的记号转换）代入所有常数：

$$
\begin{aligned}
C^{\rm src}_{\Pi_2}
&=+\frac{C_FP_2}{32}\left(\frac{\alpha_s}{\pi}\right)^3\pi^{3e}
\frac{\Gamma(1+e)^2\Gamma(1-e)^4\Gamma(1-3e)}
{(1-2e)^2\Gamma(1-2e)^2\Gamma(1+2e)}(-2r+i0)^{3e}\\
&=\boxed{+\frac{C_F}{32(1-2e)^2}\left(\frac{\alpha_s}{\pi}\right)^3
[\pi^e\Gamma(1-e)]^3(-2r+i0)^{3e}\,U(e)P_2(e).}
\end{aligned}
\tag{P2.11}
$$

特别地，$UP_2=\sum_c\kappa_c(A_cU+B_c)$，其中 $\kappa=(C_A^2,C_AT_f,C_FT_f)$。这消去了每个通道中的 $UR$，适合作为实际 C 账本；不能把 $U$ 整体删去。符号来自正号Dyson连接、顶点/传播子/算符的 $i$、Fourier 的 $i$ 和 $r/a$，没有按邻边论文答案反调。

raw 不含 Euler 补偿、C4重数、颜色投影或反项。若另定义无量纲耦合 $g_0=\mu^e g$，须显式补 $\mu^{6e}$；本文以有维裸 $g$ 工作。

## 5. 两个原参数 J 与端点控制

相邻01有 $f=z(1-s)t$、$d=z$、$u=-zt$、$w=-z(1-s)$，故 $uw=df$：

$$
\boxed{J_{01}=\frac{2+e}{2(1+2e)}\int_0^1ds_{31}\int_0^1ds_{41}\,
z[z(1-s_{31})s_{41}+i0]^{-1+3e}.}
\tag{P2.12}
$$

在 $\Re e>0$ 按原次序逐个积分，每个端点贡献 $1/(3e)$，给出理论控制式

$$
J_{01}=\frac{(2+e)z^{3e}}{18e^2(1+2e)}.
\tag{P2.13}
$$

对边02设 $s=s_{21}$、$t=s_{41}$，$F=(1-s)(1-t)+zst$、$d=-(1+z)$、$u=-1+(1+z)t$、$w=1-(1+z)s$，有 $uw=dF+z$：

$$
\boxed{J_{02}=\int_0^1ds_{21}\int_0^1ds_{41}\,
\frac{-(2+e)(1+z)F+(1-3e)z}{2(1+2e)}(F+i0)^{-2+3e}.}
\tag{P2.14}
$$

正域中的核恒等式与原端点边界为

$$
\begin{aligned}
\mathcal K_{02}&=-(1+z)F^{-1+3e}
+\frac{\partial_s\partial_t F^{3e}}{6e(1+2e)},\\
\int_0^1ds\int_0^1dt\,\partial_s\partial_tF^{3e}&=1+z^{3e}\qquad(\Re e>0).
\end{aligned}
\tag{P2.15}
$$

两个零角点的边界值按收敛域取0，不丢掉另外两个角点。相邻两个端点及对边零角点缩放均要求 $\Re e>0$；可共同取 $0<\Re e<1/3$ 后向 $e=0$ 亚纯延拓。Fourier 母式单独按解析延拓解释。一般复 $z$ 时须联合延拓尺度与距离处方，本文未证明该延拓。

## 6. 三个完整颜色通道与首坐标的边界

乘 Wilson 迹后的三种颜色为

$$
\boldsymbol{\mathcal T}=(C_A^2C_F,\ C_AC_FT_f,\ C_F^2T_f),\qquad
C^{\rm src}_{\Pi_2}=\left(\frac{\alpha_s}{\pi}\right)^3
(-2\pi r+i0)^{3e}\sum_c\mathcal T_c\,r_c^{\rm src},\quad
r_c^{\rm src}=+\frac{\Gamma(1-e)^3(A_cU+B_c)}{32(1-2e)^2}.
\tag{P2.16}
$$

按项目固定有序基底完整重构：

| 味道次数 | 有序基底 | $C_FP_2$ 的完整坐标 |
|---|---|---|
| 0 | $(C_A^2C_F,C_AC_F^2,C_F^3)$ | $(p_{AA},0,0)$ |
| 1 | $T_f(C_AC_F,C_F^2)$ | $(p_{Af},p_{Ff})$ |
| 2 | $T_f^2(C_F)$ | $(0)$ |

所以既有“重构后取 $a_1$”对无味道通道选 $p_{AA}$，对一次味道通道选 $p_{Af}$；$p_{Ff}$ 是真实的第二坐标，完整图必须另行保留。不能先把整个 $P_2$ 当作不参与颜色代数的黑箱，再宣称全部三通道都是 $a_1=1$。这也不能照搬单Π₁四端点的 $a_1=-1/2$。

因此 raw 与后续 prepared 都保存完整颜色坐标；需要 $a_1$ 目标时另行选择。颜色投影和是否共享一个 I 是独立问题。以下将完整算符因子提出，三个通道共享每张几何图的唯一 I，不再分别命名三个 $I_c$。此前逐通道记法留在[历史稿](../../../workflow/checks/self-energy-pi2-montopar/history/benchmarks/self-energy/pi2/theory.md)，不作为第二份现行定义。

## 7. 完整算符前因子与唯一的共享 I

为保持原算符提取因子，显式定义 $\mathcal B_2=4e^2P_2$。这只是将同一 C 的因素重新组织，不改变文献算符；当前C/J已吸收物理符号，旧r001另存。一次三圈 Euler 准备为

$$
\begin{aligned}
\mathcal B_2(e)&=4e^2P_2(e),\qquad G(e)=e^{-\gamma_Ee}\Gamma(1-e),\\
\mathcal N_{\Pi_2}&=\left(\frac{\alpha_s}{\pi}\right)^3 C_F\mathcal B_2(e),\\
S_{14}&=(-2\pi r+i0)^{3e},\qquad S_{12}=z^{3e}S_{14},\\
R_*(e)&=+\frac{G(e)^3U(e)}{128e^2(1-2e)^2},\qquad
C=e^{-3\gamma_Ee}C^{\rm src}_{\Pi_2}=\mathcal N_{\Pi_2}S_{14}R_*.
\end{aligned}
\tag{P2.17}
$$

$\log U=6\zeta_3e^3+9\zeta_4e^4+O(e^5)$，没有线性 Euler 项；总补偿取3而不是算符内部圈数2。raw 中的 Gamma 因子保持完整，montopar 不提前补偿。

最终 I 的**唯一定义**为

$$
\boxed{\begin{aligned}
\frac m2 C[J(z,e)+z^{3e}J(z^{-1},e)]
&=\mathcal N_{\Pi_2}(S_{14}+S_{12})I(z,e),\\
I(z,e)&=\frac m2 R_*(e)
\frac{J(z,e)+z^{3e}J(z^{-1},e)}{1+z^{3e}}.
\end{aligned}}
\tag{P2.18}
$$

三个颜色通道只改变外部系数 $4e^2p_c$；每个几何图只有一个 I。I 不含颜色、味数、耦合或参考尺度，也不含 $p_c$。无需除以可能在特殊颜色/ε处为零的 $P_2$：第二行直接定义 I，第一行是其系数装配关系。

原几何、次序、C4重数及完整核与[串联分支](../pi1-chain/theory.md)逐项相同，且 $R_*=-R_{\rm chain}$，故

$$
\boxed{I^{\Pi_2}(z,e)=-I_{\rm chain}(z,e),\qquad
C_F\mathcal B_2=4e^2\left[C_A^2C_Fp_{AA}+C_AC_FT_fp_{Af}+C_F^2T_fp_{Ff}\right].}
\tag{P2.19}
$$

这不是 $\Pi_2=\Pi_1^2$。相同的是几何J；保留原外部算符因子后，两者的I相反；Π₂ 的 $\mathcal B_2$ 与串联的 $\mathcal B_1^2$ 不同。若要各物理颜色系数，对同一个 I 乘 $4e^2p_c$ 并按所需阶数卷积即可。若只取既有首坐标目标，则在 $\mathcal N_{\Pi_2}$ 中选 AA/Af；Ff 仍在 full-color 记录中，I 不变。

相邻理论控制式为

$$
I_{01}=+\frac{G^3U(2+e)}{1152e^4(1-2e)^2(1+2e)}.
\tag{P2.20}
$$

源系数合并后 $p_{AA}\sim-25/(12e^2)$、$p_{Af}\sim5/(3e^2)$、$p_{Ff}\sim2/e$。共享 I 与串联一样，相邻四阶、对边三阶；乘回 $4e^2p_c$ 后 AA/Af 维持这些极点阶数，Ff 分别降为三阶、二阶。三个颜色有限项的现行r002已同步发布；几何核复用与颜色卷积分别核验。

求完整颜色有限项时，先合并 $A_cU+B_c$。采用式(P2.16)的通道前因子，相邻其 $e^2r_c$ 展开至 $e^4$、对边至 $e^3$，J均至 $e^2$ 足够。若分别展开有表观 $e^{-3}$ 的 A/B 项，相邻 U 至少需 $e^5$、对边至少需 $e^4$，再检查消去；不要过早截断。

## 8. montopar 实现边界

| 层 | 本推导要求 |
|---|---|
| 源校验 | `2e-Pi2`、总三圈、$g^6$、恰一个 `vacuum-pol-2loop`、两侧传播子同动量且内部指标完整 |
| 算符 | 式(P2.2)的 $-i$、$g^4$、$q^{-2e}$；`1PI`、三颜色系数、文献约定与来源哈希 |
| 原 C/J | 式(P2.11–14)，完整分子、原边界和外到内次序；raw 不补 Euler、不投影 |
| full color | 保留 AA/Af/Ff，尤其 $C_F^2T_f$ 第二坐标；实际字段直接用 $C_F^2T_f$ |
| prepared | 一次 $e^{-3\gamma_Ee}$；全向量与首坐标分开；未选 Ff 不标记为零 |
| 核验 | 一圈字典回归、相位账本、D维迹、Gamma前因子、两图核/域、C4、颜色投影及错误输入拒绝 |
| 拒绝 | Π₂误标串联、额外传播子/插入、混合动量、缺 δ 桥、未知规范/重整化约定、漏 Ff 或伪造 $n_f^2$ |

六边对求和相消不是逐图参数化的前提；本轮未核验 Π₂ 的闭合围道总和。当前默认47条 Kernel 根保持不变。形式状态为 `formal=not_run`，辅助代数和数值检查见[理论核验](../../../workflow/checks/self-energy-pi2-theory/README.md)。
