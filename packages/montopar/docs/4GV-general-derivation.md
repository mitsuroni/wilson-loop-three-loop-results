# 通用4GV因子分解：有序Wilson圈上的四胶子星型积分

版本：2026-09-18修订版。适用于标准QCD/Yang–Mills的一个裸四胶子顶点，四条传播子分别连接Wilson圈的四个有序插入。理论几何可取一般N边形；当前montopar的`fourv`实现限于类光四边形的`4e-4g`。本文是修正后的推导，不宣称物理分布、复解析延拓或Lean整链已经获证。

修正要点：点间距方向号、无量纲Lorentz因子、同边排序和Jacobian、四个α的单纯形、退化表、颜色零通道证明，以及C各因子的完整账本。原修订前文档和输出已保存于Wilson审计快照；不以修改后的正文替代旧版本证据。

## 1. 几何、参数与所计算对象

令边向量满足 $v_i^2=0$、$\sum_{i=1}^N v_i=0$，取

$$
V_1=0,\qquad V_i=\sum_{k<i}v_k,\qquad X_i(s)=V_i+s v_i.
$$

四条腿记为 $r=1,2,3,4$，其附着边依次为 $(a,b,c,d)$。每条腿有自己的参数 $s_r$，即使边相同也不能合并参数。腿顺序必须对应Wilson闭合颜色迹的顺序；重标号时同时重标号颜色、Lorentz指标、位置、α和排序域。

同边插入的参数域为

$$0\le s_{e,1}\le s_{e,2}\le\cdots\le s_{e,m_e}\le1.$$

各边这些有序单纯形的乘积记为 $\mathcal D$。这里使用有序积分，不再附加重复的 $1/m_e!$。输出是给定代表图的原始贡献 $\mathcal W_g=C_{\rm raw}J_g$；orbit求和重数、crossing、固定颜色投影及Euler补偿均另行明确处理。它不是已经完成求和/重整化的 $\log W$，也不直接等于Wilson项目规范后的 $I$。

采用 $v_{ij}=v_i\cdot v_j$、$D=4-2\epsilon$ 和归一化基本表示迹。对四边形选 $v_{\rm ref}=v_{14}\ne0$，

$$z=\frac{v_{12}}{v_{14}}=\frac{v_{34}}{v_{14}},\quad
v_{23}=v_{14},\quad v_{13}=v_{24}=-(1+z)v_{14}.$$

$z$是运动学尺度比，不是类光四边形的独立对偶共形交比。一般N的Mandelstam计数需与固定维数Gram约束区分，不能据此断言所有变量都是独立交比。

## 2. 动量空间结构与顶点

在 $e^{iS}$、Wilson插入 $ig v^mu T^A$、传播子 $-ig_{\mu\nu}\delta^{AB}/(p^2+i0)$ 的约定下，四腿贡献包含

$$
\frac1{N_c}\operatorname{tr}\!\int_{\mathcal D}d^4s\,
\prod_{r=1}^4(ig v_{e_r}^{\mu_r}T^{A_r})
\prod_{r=1}^4\left[\int\frac{d^D\ell_r}{(2\pi)^D}
\frac{-ig_{\mu_r\nu_r}\delta^{A_rB_r}}{\ell_r^2+i0}
 e^{-i\ell_r\cdot X_r}\right]
\mathcal V_4^{B_1B_2B_3B_4;\nu_1\nu_2\nu_3\nu_4}
(2\pi)^D\delta^{(D)}\!\left(\sum_r\ell_r\right).
$$

上式的生成元乘积按既定有序字理解；不交换非对易因子。所有腿取入射，四胶子顶点为

$$
\begin{aligned}
\mathcal V_4^{ABCD;\mu\nu\rho\sigma}=-ig^2\big[&
f^{ABX}f^{CDX}(g_{\mu\rho}g_{\nu\sigma}-g_{\mu\sigma}g_{\nu\rho})\\
+&f^{ACX}f^{BDX}(g_{\mu\nu}g_{\rho\sigma}-g_{\mu\sigma}g_{\nu\rho})\\
+&f^{ADX}f^{BCX}(g_{\mu\nu}g_{\rho\sigma}-g_{\mu\rho}g_{\nu\sigma})\big].
\end{aligned}
$$

若第二项改写成 $f^{ACX}f^{DBX}$，其Lorentz括号也应反号。顶点没有动量，不产生 $\ell\to i\partial_X$ 分子。此顶点已包含拉氏量四场项的24种标记指派，不再额外乘 $4!$ 或 $1/4!$。[Weinzierl讲义，第114页](https://particlephysics.uni-mainz.de/weinzierl/download/script_qft1.pdf#page=114)。

## 3. 固定有序迹的颜色收缩

定义

$$
c_1=\frac1{N_c}\operatorname{tr}(T^AT^BT^CT^D)f^{ABX}f^{CDX},\quad
c_2=\frac1{N_c}\operatorname{tr}(T^AT^BT^CT^D)f^{ACX}f^{BDX},\quad
c_3=\frac1{N_c}\operatorname{tr}(T^AT^BT^CT^D)f^{ADX}f^{BCX}.
$$

利用 $f^{ABX}T^AT^B=iC_AT^X/2$ 得

$$c_1=-C_A^2C_F/4.$$

对 $c_3$，先收缩相邻的 $B,C$，再用迹的循环性把 $T^D$ 移到前面；$f^{ADX}T^DT^A=-iC_AT^X/2$，故

$$c_3=+C_A^2C_F/4.$$

Jacobi恒等式给 $c_1-c_2+c_3=0$，于是

$$\boxed{(c_1,c_2,c_3)=C_A^2C_F(-1/4,0,1/4).}$$

这才是第二颜色通道为零的论证。“两个生成元被隔开所以不能反对称配对”不是证明。结论依赖上述颜色字；任意重新排列腿而保持旧通道系数不变是不合法的。独立程序还通过 $f^{ABC}=-2i\operatorname{tr}([T^A,T^B]T^C)$ 和基本表示Fierz恒等式重建了一般 $N_c$ 的结果。

## 4. Lorentz收缩与代数零

依上述三种配对，

$$
L_1=v_{ac}v_{bd}-v_{ad}v_{bc},\qquad
L_2=v_{ab}v_{cd}-v_{ad}v_{bc},\qquad
L_3=v_{ab}v_{cd}-v_{ac}v_{bd}.
$$

必须逐通道将颜色和Lorentz相乘：

$$\boxed{\sum_jc_jL_j=-\frac{C_A^2C_F}{4}
(2v_{ac}v_{bd}-v_{ad}v_{bc}-v_{ab}v_{cd}).}$$

对两条类光方向 $u,v$，令 $h=u\cdot v$：

| 腿方向 | $L_1$ | $L_2$ | $L_3$ |
|---|---|---|---|
| $(u,u,v,v)$ | $0$ | $-h^2$ | $-h^2$ |
| $(u,v,u,v)$ | $-h^2$ | $0$ | $h^2$ |
| $(u,v,v,u)$ | $h^2$ | $h^2$ | $0$ |
| 三腿同一方向 | $0$ | $0$ | $0$ |

后三种腿置换若用于颜色合并，也须同步改变颜色字。当前4e-4g的01/03/04因三腿同边而为零；在一般非退化正 $z$ 四边形、按轮廓排序的边多重集合内，其余最多两腿同边的模式非恒等零。对任意N、特殊运动学或单个数值点，不把这个充分条件写成无条件的“当且仅当”。

## 5. 四传播子的Fourier变换与星型积分

按声明的Feynman处方，单传播子为

$$\int\frac{d^Dp}{(2\pi)^D}e^{-ip\cdot x}\frac{-i}{p^2+i0}
=-\frac{\Gamma(1-\epsilon)}{4\pi^{2-\epsilon}}
(-x^2+i0)^{-1+\epsilon}.$$

动量守恒δ引入一个共同位置 $y$。定义几何多项式

$$Q_0=-\sum_{i<j}\alpha_i\alpha_j(X_i-X_j)^2.$$

四点标量星核参数化为

$$
\begin{aligned}
&\int d^Dy\prod_{i=1}^4[-(X_i-y)^2+i0]^{-1+\epsilon}\\
&\quad=\frac{i\pi^{2-\epsilon}\Gamma(2-3\epsilon)}{\Gamma(1-\epsilon)^4}
\int_{\alpha_i\ge0}d^4\alpha\,\delta(1-\sum_i\alpha_i)
\frac{(\alpha_1\alpha_2\alpha_3\alpha_4)^{-\epsilon}}
{(Q_0+i0)^{2-3\epsilon}}.
\end{aligned}
$$

这是在收敛/解析起始区推导后作亚纯及Feynman边界值延拓的公式；在类光端点奇异处不能仅凭形式代换宣称Fubini或逐项ε展开。四条传播子的Γ幂抵消后，公共Fourier系数为

$$\frac{i\Gamma(2-3\epsilon)}{256\pi^{6-3\epsilon}}.$$

星核幂次来自 $4(1-\epsilon)-D/2=2-3\epsilon$。相位账本为

$$i^4\times(-i)\times(-1)^4\times i=1.$$

依次对应四个Wilson顶点、4GV、四个坐标传播子及位置积分。没有TGV的导数因子。相位账本为实数不意味着尺度复幂在所有物理区也为实数。

## 6. 正确距离：起点差与首尾边

对两条腿附着在 $a,b$ 边，令 $V_{ba}=V_b-V_a$。正确点差为

$$X_a-X_b=-V_{ba}+s_av_a-s_bv_b,$$

所以

$$
\boxed{\Delta_{ab}=-[V_{ba}-s_av_a+s_bv_b]^2+i0
=-[V_{ba}^2-2s_a v_a\cdot V_{ba}+2s_b v_b\cdot V_{ba}-2s_as_bv_{ab}]+i0.}
$$

不可将其写成 $-[V_{ba}+s_av_a-s_bv_b]^2+i0$。对通常相邻边 $b=a+1$，有

$$\Delta_{ab}=-2v_{ab}(1-s_a)s_b+i0.$$

对首尾边 $N\to1$，闭合条件给 $V_N=-v_N$，因此

$$\Delta_{1N}=-2v_{1N}s_1(1-s_N)+i0.$$

同边两点的几何间距平方恒等为零；这里指几何多项式，不删除公共Feynman处方。不同边也可能在端点或特殊运动学处距离为零，不能写成逐点的“距离零当且仅当同边”。

## 7. 无量纲核与分支

选 $A=-2v_{\rm ref}+i0$，在解析起始区定义

$$d_{ij}=\frac{(X_i-X_j)^2}{2v_{\rm ref}},\qquad
S=\sum_{i<j}\alpha_i\alpha_jd_{ij},\qquad
\ell_{abcd}=\frac{-L_1+L_3}{4v_{\rm ref}^2}.$$

$\ell$无量纲。为避免不合法的复幂拆分，首先在 $z>0,v_{\rm ref}<0$ 起始区使用实正尺度，然后沿一致的Feynman延拓定义尺度分解。完整物理母式始终是 $(Q_0+i0)^{-2+3\epsilon}$；不把有限的独立 $i0$ 符号当成可任意相乘/相除的实数。

$\sum c_jL_j=C_A^2C_F v_{\rm ref}^2\ell$。尺度的代数账本为

$$v_{\rm ref}^2(-2v_{\rm ref})^{-2+3\epsilon}
=\frac14(-2v_{\rm ref})^{3\epsilon}.$$

因此净尺度之外应留下 $\ell$，而不是仍含 $v_{\rm ref}^2$ 的原始 $L_j$。

## 8. 完整C与原参数J

公共因子采用可逐项相乘的定义：

$$
\begin{aligned}
C_{\rm raw}&=C_A^2C_F\,(4\pi\alpha_s)^3
\frac{\Gamma(2-3\epsilon)}{256\pi^{6-3\epsilon}}
(-2v_{\rm ref}+i0)^{3\epsilon}\,\frac14\\
&=\boxed{\frac{C_A^2C_F}{16}\left(\frac{\alpha_s}{\pi}\right)^3
\Gamma(2-3\epsilon)\pi^{3\epsilon}(-2v_{\rm ref}+i0)^{3\epsilon}.}
\end{aligned}
$$

$$\boxed{J_{abcd}=\int_{\mathcal D}d^4s
\int_{\alpha_i\ge0}d^4\alpha\,\delta(1-\sum_i\alpha_i)
\ell_{abcd}(z)(\alpha_1\alpha_2\alpha_3\alpha_4)^{-\epsilon}
S^{-2+3\epsilon}.}$$

两个 $1/4$ 来源不同：$\ell$中的一个来自颜色通道，C中的一个来自尺度消去，不能删掉或重复。零图可保留同一个公共C并令J=0，零原因记录在Lorentz层。

raw输出未加入Euler补偿。若接入Wilson三圈规范，须在独立预处理步骤中恰好乘一次 $e^{-3\gamma_E\epsilon}$，保留原始C与变换记录；此步骤不能与Γ因子相混。

## 9. 显式α单纯形和有序域换元

消去 $\alpha_4$ 后采用普通δ测度：

$$\alpha_4=1-\alpha_1-\alpha_2-\alpha_3,\quad
0\le\alpha_1\le1,\quad0\le\alpha_2\le1-\alpha_1,\quad
0\le\alpha_3\le1-\alpha_1-\alpha_2.$$

消去δ的Jacobian为1；不能把四个α分别积在独立的 $[0,1]$ 中。权重 $(\prod\alpha)^{-\epsilon}$ 已在integrand内，不再在measure字段重复乘一次。

对于 $s_1\le\cdots\le s_m$，正确乘法映射为

$$s_m=x_1,\quad s_{m-1}=x_1x_2,\quad\ldots,\quad s_1=x_1\cdots x_m,$$

$$\left|\frac{\partial s}{\partial x}\right|=x_1^{m-1}x_2^{m-2}\cdots x_{m-1}.$$

在非退化内部映射可逆；边界按测度/可积性解释。各边映射合成后，cube J的分子是 $\operatorname{Jac}\times\ell$，不能因4GV无导数分子而将它写成1。

## 10. 四边形九个代表

当前四条腿按各自α序排列；下表不乘orbit重数。

| 图 | 边 $(a,b,c,d)$ | $\ell(z)$ | 换元Jacobian | 有效颜色×Lorentz通道 |
|---|---|---|---|---|
| 01 | $(3,4,4,4)$ | $0$ | $x_2^2x_3$ | 无 |
| 02 | $(3,3,4,4)$ | $-z^2/4$ | $x_1x_2$ | 3 |
| 03 | $(3,3,3,4)$ | $0$ | $x_1^2x_3$ | 无 |
| 04 | $(2,4,4,4)$ | $0$ | $x_2^2x_3$ | 无 |
| 05 | $(2,3,4,4)$ | $z(1+z)/4$ | $x_3$ | 3 |
| 06 | $(2,3,3,4)$ | $-z/4$ | $x_2$ | 1 |
| 07 | $(2,2,4,4)$ | $-(1+z)^2/4$ | $x_1x_2$ | 3 |
| 08 | $(2,2,3,4)$ | $(1+z)/4$ | $x_1$ | 3 |
| 09 | $(1,2,3,4)$ | $-(1+4z+z^2)/4$ | $1$ | 1、3 |

06同边配对是腿2、3，即 $d_{23}=0$，不是 $d_{12}=0$。每条非零轨道的S始终由本腿序的六对几何距离计算，不从另一图的多项式猜测。

对09，$s_{11}=x_1,s_{21}=x_2,s_{31}=x_3,s_{41}=x_4$，

$$\begin{aligned}
S_{09}={}&\alpha_1\alpha_2 z(1-x_1)x_2
+\alpha_1\alpha_3[z(1-x_1)(1-x_3)+x_1x_3]\\
&+\alpha_1\alpha_4 x_1(1-x_4)
+\alpha_2\alpha_3(1-x_2)x_3\\
&+\alpha_2\alpha_4[(1-x_2)(1-x_4)+zx_2x_4]
+\alpha_3\alpha_4 z(1-x_3)x_4.
\end{aligned}$$

旧09导出把第三项改成 $\alpha_1\alpha_4x_4(1-x_1)$，其差为 $\alpha_1\alpha_4(x_4-x_1)$，不为零。应重新生成旧文件，不能只修改说明文字后沿用它们。

## 11. montopar输出契约

当前修订使用`schema_version=1.2`，每个orbit只有一个规范代表：

- `C.components`和`factor_order`是可乘账本，乘积等于`C.combined`；相位来源另列`phase_provenance`。
- `J_raw_s`保留四个独立腿参数、原始同边依赖界、外到内次序、S和Lorentz因子。
- `J_raw`保留cube式、原s到x映射和Jacobian。
- 两者均显式记录四个α的单纯形、非负条件、消去变量、积分界和δ约束。
- `representative`保存原始拓扑及腿序，`normalization`声明未补Euler、未投影颜色、未乘重数、未组装crossing。
- `branch_prescription`保存物理母式和解析起始区；不以裸字符串`-i0`替代完整定义。
- `final.json`、`state.json`与`export.dimreduce.json`只从已验证的新主输出派生，不复制历史错误核。

实现为纯Python代数，不声称逐图调用FeynCalc。程序正确性检查、积分求值、形式整链分别标记。4GV不含动量导数，不要求为了模仿TGV而先输出某个特定IBP形式。

## 12. 适用边界

标准重整化Yang–Mills作用量只有三、四胶子裸局部顶点。一般n点标量星积分的参数公式可以推广，但不意味着存在任意 $n\ge5$ 且无导数、分子为1的QCD裸胶子顶点；有效高点顶点需单独定义。

TGV求导前标量星核幂次为 $D-3=1-2\epsilon$；求导后的原TGV分母幂次为 $D-2=2-2\epsilon$。4GV的 $2-3\epsilon$ 来自四传播子星核，与上述两个幂次不要混用。

本修订得到自包含的参数积分C/J及代数检查。图枚举完备性、一般分布/Fubini、物理复分支、实际Laurent系数和Lean证明均需后续工作，不因格式修复自动完成。