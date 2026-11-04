# 双一圈胶子自能串联的完整参数化

本文是`2e-Pi1+Pi1`分支的理论与实现契约。以[项目现行规范化源](sources/2e-Pi1+Pi1-r002.json)为起点，原动量图贡献统一记作 $\mathcal A=C^{\mathrm{src}}J$；$I$只指第7节的归一化系数。两个自能均采用[单插入主稿](../pi1/theory.md)中的同一个裸一圈算符。本文不改变旧13图的定义，也不把串联项并入二圈1PI算符 $\Pi_2$。

当前两图C/J与I均已完成，采用2.0约定；数值相对r001保持。

现行数学约定为[自能2.0](../conventions.md)，输入和结果见[当前清单](../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1-chain/theory.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。

## 1. 源积分、指标与圈数

设两个Wilson端点位于边i、j，$X=x_i(s)-x_j(t)$。采用 $D=4-2\varepsilon$、度规 $(+---)$，记 $e=\varepsilon$。源积分为

$$
\begin{aligned}
\mathcal A_{ij}={}&\int_0^1ds\int_0^1dt\,
(igv_i^{\mu_1}t^{a_1})(igv_j^{\mu_2}t^{a_2})\\
&\times\int\frac{d^D\ell}{(2\pi)^D}
\frac{-ig_{\mu_1\nu_1}\delta^{a_1b_1}}{\ell^2+i0}
\Pi_{b_1b_2}^{\nu_1\nu_2}(\ell)
\frac{-ig_{\nu_2\nu_3}\delta^{b_2b_3}}{\ell^2+i0}
\Pi_{b_3b_4}^{\nu_3\nu_4}(\ell)
\frac{-ig_{\nu_4\mu_2}\delta^{b_4a_2}}{\ell^2+i0}
e^{-i\ell\cdot X}.
\end{aligned}
\tag{SC1}
$$

两个Π均为含i的1PI插入，按统一Dyson规则以+1连接。旧r001的两个额外负号相乘也为+1，故本分支数学表达保持；现行r002与其他分支使用同一算符契约。

每个Π中已经包含一个g²、一个i及其内部圈积分。式(SC1)的总阶数为g⁶、总圈数为3，但只剩一个显式交换动量积分；三段传播子共享同一ℓ。链系数为1，没有额外1/2!。两个相同算符的标签置换不另算一个图。

仍用归一化颜色迹Tr/Nc提取单位矩阵系数。源JSON的相位若写为`exp(+i l1·Phi)`，则 $\Phi=-X$。

## 2. 两算符的完整收缩

记 $q=-\ell^2-i0$，一圈算符为

$$
\begin{aligned}
\Pi_{ab}^{\mu\nu}(\ell)&=i\delta_{ab}\,\mathcal A_1(e)
T^{\mu\nu}(\ell)q^{-e},\qquad
T^{\mu\nu}=\ell^2g^{\mu\nu}-\ell^\mu\ell^\nu,\\
\mathcal A_1(e)&=\frac{g^2}{2(4\pi)^{2-e}}
\frac{\Gamma(e)\Gamma(1-e)^2}{\Gamma(2-2e)}\,\mathcal B_1(e),\\
\mathcal B_1(e)&=\frac{(10-6e)C_A-8(1-e)T_Fn_f}{3-2e}.
\end{aligned}
\tag{SC2}
$$

此处 $T_F=1/2$时与原NF记号一致；不设 $n_f=0$。因为

$$
T^{\mu\alpha}g_{\alpha\beta}T^{\beta\nu}=\ell^2T^{\mu\nu},
\qquad \ell_\mu T^{\mu\nu}=0,
\tag{SC3}
$$

两个自能收缩后仍保留一个完整T，而不是两个独立坐标核的乘积。颜色指标的五个δ串联后给出 $t^at^a=C_F\mathbf1$；中间传播子不会再产生一个 $C_F$或额外 $C_A$。

按式(SC1)逐项相乘，$(ig)^2(-i)^3i^2=+ig^2$。结合(SC3)，在共同Feynman边界值与维数解析延拓意义下整理传播子幂，得到

$$
\mathcal A_{ij}=-ig^2C_F\mathcal A_1(e)^2
\int_0^1ds\int_0^1dt\int\frac{d^D\ell}{(2\pi)^D}e^{-i\ell\cdot X}
\left[(v_i\cdot v_j)q^{-1-2e}
+(v_i\cdot\ell)(v_j\cdot\ell)q^{-2-2e}\right].
\tag{SC4}
$$

这里使用的是边界值分布的乘法/解析延拓，不将有限iη的分母直接当作无处方的普通有理式。唯一被平方的是标量系数 $\mathcal A_1$，并不是单插入的整个J或I。

## 3. Fourier变换与D维坐标核

沿用同一个Fourier母式，$Q=-X^2+i0$：

$$
\mathcal F_a(X)=\int\frac{d^D\ell}{(2\pi)^D}e^{-i\ell\cdot X}q^{-a}
=i\frac{\Gamma(D/2-a)}{4^a\pi^{D/2}\Gamma(a)}Q^{a-D/2}.
\tag{SC5}
$$

它按维数解析延拓解释；不要求拆出的每个Fourier积分都在同一实ε区间绝对收敛。动量二次项对应负的二阶坐标导数，故

$$
\mathcal A_{ij}=-ig^2C_F\mathcal A_1^2\int ds\,dt\,
v_i^\mu v_j^\nu\left[g_{\mu\nu}\mathcal F_{1+2e}-\partial_\mu\partial_\nu\mathcal F_{2+2e}\right].
\tag{SC6}
$$

具体系数为

$$
\begin{aligned}
\mathcal F_{1+2e}&=ih_1Q^{-1+3e},&
h_1&=\frac{\Gamma(1-3e)}{4^{1+2e}\pi^{2-e}\Gamma(1+2e)},\\
\mathcal F_{2+2e}&=ih_2Q^{3e},&
h_2&=\frac{\Gamma(-3e)}{4^{2+2e}\pi^{2-e}\Gamma(2+2e)},\\
\frac{h_2}{h_1}&=-\frac1{12e(1+2e)}.
\end{aligned}
\tag{SC7}
$$

直接作用坐标导数：

$$
\partial_\mu\partial_\nu Q^{3e}
=-6e\,g_{\mu\nu}Q^{-1+3e}
+12e(3e-1)X_\mu X_\nu Q^{-2+3e}.
\tag{SC8}
$$

于是单一完整坐标表达式为

$$
\boxed{\mathcal A_{ij}=g^2C_F\mathcal A_1^2h_1\int ds\,dt\,
\left[\frac{1+4e}{2(1+2e)}(v_i\cdot v_j)Q^{-1+3e}
-\frac{1-3e}{1+2e}(v_i\cdot X)(v_j\cdot X)Q^{-2+3e}\right].}
\tag{SC9}
$$

这一步的实系数为正，因为式(SC6)的−i与Fourier的i相乘。D维迹为 $(D-1)h_1Q^{-1+3e}$：用 $X^2=-Q$，可验证 $(1+4e)D/[2(1+2e)]+(1-3e)/(1+2e)=D-1$。不能把张量迹中的D提前设为4。

## 4. C/J拆分与完整前因子

与原montopar几何一致，取

$$
r=v_{14},\quad z=v_{12}/v_{14},\quad A=-2r,\qquad
f=\frac{X^2}{2r},\quad d=\frac{v_i\cdot v_j}{r},\quad
u=\frac{v_i\cdot X}{r},\quad w=\frac{v_j\cdot X}{r}.
\tag{SC10}
$$

在 $A>0,z>0$ 起始域中 $Q=Af+i0$，$r/A=-1/2$。统一提出的源C为负：

$$
\boxed{\begin{aligned}
\mathcal A_{ij}&=C^{\mathrm{src}}_{\rm chain}J_{ij},\\
C^{\mathrm{src}}_{\rm chain}&=-\frac{g^2C_F\mathcal A_1^2h_1}{2}(A+i0)^{3e},\\
J_{ij}&=\int_0^1ds\int_0^1dt\,\mathcal K_{\rm chain}(f,d,u,w;e),\\
\mathcal K_{\rm chain}&=
\frac{(1+4e)df+(1-3e)uw}{2(1+2e)}(f+i0)^{-2+3e}.
\end{aligned}}
\tag{SC11}
$$

本式是首版montopar应直接输出的完整核。符号链为：源外部+、两个Wilson顶点、三段−i传播子、两个Π中的i，给出+i；转成q幂后为−i；Fourier给出实正；最终 $r/A=-1/2$使C为负。

代入g²=4παs（裸耦合的记号转换），未约化的Gamma形式为

$$
C^{\mathrm{src}}_{\rm chain}=-\frac{C_F\mathcal B_1^2}{128}
\left(\frac{\alpha_s}{\pi}\right)^3\pi^{3e}
\frac{\Gamma(e)^2\Gamma(1-e)^4\Gamma(1-3e)}
{\Gamma(2-2e)^2\Gamma(1+2e)}(-2v_{14}+i0)^{3e}.
\tag{SC12}
$$

定义只含Gamma函数的明确比值

$$
\mathcal U(e)=\frac{\Gamma(1+e)^2\Gamma(1-e)\Gamma(1-3e)}
{\Gamma(1-2e)^2\Gamma(1+2e)}.
\tag{SC13}
$$

则适合代码账本的等价形式是

$$
\boxed{C^{\mathrm{src}}_{\rm chain}=
-\frac{C_F\mathcal B_1(e)^2\mathcal U(e)}{128e^2(1-2e)^2}
\left(\frac{\alpha_s}{\pi}\right)^3
[\pi^e\Gamma(1-e)]^3(-2v_{14}+i0)^{3e}.}
\tag{SC14}
$$

$\mathcal U$不是新的未知积分；它不能因单插入公式没有该比值而省略。特别地，$\log\mathcal U=6\zeta_3e^3+9\zeta_4e^4+O(e^5)$，无线性Euler项，也无e²项。把单插入公式中的e整体替换为2e会错误地改变泡算符本身的系数，不能得到(SC14)。

raw阶段C不含Euler补偿、轨道重数、额外阶乘或反项。若后续改用 $g_0=\mu^e g$，则必须显式补总 $\mu^{6e}$，不能在本账本中隐含加入。

## 5. 两张原参数C/J

两图共用(SC14)，只改变几何。几何由类光闭合四边形的原Gram表计算，$v_i^2=0$、$\sum_i v_i=0$，不引入新坐标。

**01：相邻边3、4。** $s=s_{31},t=s_{41}$，$X=-[(1-s)v_3+tv_4]$；

$$
f=z(1-s)t,\quad d=z,\quad u=-zt,\quad w=-z(1-s),\qquad uw=df.
\tag{SC15}
$$

所以

$$
\boxed{J_{01}=\frac{2+e}{2(1+2e)}\int_0^1ds_{31}\int_0^1ds_{41}\,
z[z(1-s_{31})s_{41}+i0]^{-1+3e}.}
\tag{SC16}
$$

相邻核的系数是 $(2+e)/[2(1+2e)]$。它来自完整分子 $[(1+4e)+(1-3e)]df$，不能照搬单插入的 $1/(1+e)$，也不等于其平方。原外到内顺序为s31→s41，Jacobian为1。

**02：对边2、4。** $s=s_{21},t=s_{41}$，$X=(1-t)v_4+v_1+sv_2$，与冻结源的 $\Phi=-X$ 一致；

$$
F=(1-s)(1-t)+zst,\quad d=-(1+z),\quad
u=-1+(1+z)t,\quad w=1-(1+z)s,\quad uw=dF+z.
\tag{SC17}
$$

因此

$$
\boxed{J_{02}=\int_0^1ds_{21}\int_0^1ds_{41}\,
\frac{-(2+e)(1+z)F+(1-3e)z}{2(1+2e)}(F+i0)^{-2+3e}.}
\tag{SC18}
$$

原外到内顺序为s21→s41，Jacobian为1。在正域有辅助恒等式

$$
\mathcal K_{02}=-(1+z)F^{-1+3e}
+\frac{\partial_s\partial_tF^{3e}}{6e(1+2e)}.
\tag{SC19}
$$

它可用于独立核验及后续求解；生产主输出仍使用(SC18)的完整显式核，不把边界项删除。

## 6. 颜色通道与一次Euler准备

设 $a=(10-6e)/(3-2e)$、$b=8(1-e)/(3-2e)$。完整颜色为

$$
C_F\mathcal B_1^2=C_F\left[a^2C_A^2-2abC_AT_Fn_f+b^2(T_Fn_f)^2\right].
\tag{SC20}
$$

必须保留三个通道；混合项的2来自平方展开，不能漏掉。对相应完整有序颜色基底，坐标分别为

| 通道基底 | 完整坐标 |
|---|---|
| $(C_A^2C_F,C_AC_F^2,C_F^3)$ | $(a^2,0,0)$ |
| $T_Fn_f(C_AC_F,C_F^2)$ | $(-2ab,0)$ |
| $(T_Fn_f)^2(C_F)$ | $(b^2)$ |

两端点Wilson迹就是CF；将共同算符因子 $\mathcal B_1^2$提出后，剩余首坐标为1。这不同于单Π₁加裸交换的 $a_1=-1/2$。raw必须保留全颜色，不用纯胶子三圈的 $C_A^2C_F$将味道部分丢掉。

准备阶段只加一次 $e^{-3\gamma_Ee}$。定义 $G=e^{-\gamma_Ee}\Gamma(1-e)$，则(SC14)中的 $\Gamma(1-e)^3$变为G³，$\mathcal U$保持。由于U内部仍含Gamma函数，旧的“只数Gamma(1−e)字面总幂”等检查不能直接套用；应验证(SC12)与(SC14)等价，并验证U的线性Euler系数为零。

## 7. 唯一I的归一化与控制例

定义

$$
\begin{aligned}
\mathcal N_{\rm chain}&=\left(\frac{\alpha_s}{\pi}\right)^3C_F\mathcal B_1(e)^2,\\
S_{14}&=(-2\pi v_{14}+i0)^{3e},\qquad S_{12}=z^{3e}S_{14},\\
R_{\rm chain}&=-\frac{G^3\mathcal U(e)}{128e^2(1-2e)^2},\qquad
C=e^{-3\gamma_Ee}C^{\rm src}_{\rm chain}=\mathcal N_{\rm chain}S_{14}R_{\rm chain}.
\end{aligned}
\tag{SC21}
$$

带插入结构的C4轨道大小为相邻m=4、对边m=2。链反向时两个相同Π₁只是交换位置，不再产生额外标签重数。与其他族相同的系数提取关系定义唯一I：

$$
\frac m2 C[J(z,e)+z^{3e}J(1/z,e)]
=\mathcal N_{\rm chain}(S_{14}+S_{12})I(z,e),
\tag{SC22}
$$
$$
\boxed{I(z,e)=\frac m2 R_{\rm chain}
\frac{J(z,e)+z^{3e}J(1/z,e)}{1+z^{3e}}.}
\tag{SC23}
$$

三个颜色通道共享这一标量I。色因子、耦合和参考尺度从I中提出，不能将 $C^{\rm src}J$另称为第二个I。

作为实现时的端点控制，(SC16)在Re(e)>0给出

$$
J_{01}=\frac{(2+e)z^{3e}}{18e^2(1+2e)},\qquad
I_{01}=-\frac{G^3\mathcal U(e)(2+e)}{1152e^4(1-2e)^2(1+2e)}.
\tag{SC24}
$$

其首项是 $-1/(576e^4)$，只用作理论与将来代码的控制例；当前尚无串联生产输出、prepared清单或求解结果修订。对边图在本阶段保持原参数式(SC18)，不提前标记已求解。

## 8. 后续montopar实现契约

| 对象 | 必须输出/检查的内容 |
|---|---|
| 族与来源 | `2e-Pi1+Pi1`，nloops=3，g_power=6，与`2e-Pi2`区分 |
| 图连接 | 两个`vacuum-pol`、三段同动量传播子；按内部指标连接核验一条无分叉路径 |
| 算符元数据 | insertion_count=2，子圈数[1,1]，完整算符列表，Dyson系数1，源外部符号+1 |
| `C.components` | Wilson颜色CF、完整B₁²、耦合三次、$[\pi^e\Gamma(1-e)]^3$、U、$-1/[128e^2(1-2e)^2]$、三圈尺度 |
| `C.combined` | 上述分量严格相乘，与(SC12)/(SC14)匹配；raw不含Euler或m |
| `J_raw` | 完整多项式分子、原f及其−2+3e幂、原s域、Jacobian1、外到内次序 |
| 颜色准备 | 三味道通道及完整坐标，不沿用4e-Pi1的−1/2 |
| C4与符号 | 轨道4/2，算符标签交换不加重数；完整相位账本可复验 |
| 拒绝条件 | 不同动量、缺中间传播子、分叉/断链、指标不成对、误标Π₂、错误g阶或额外插入 |

可复用旧几何、Gram多项式和原域输出器；需要新的链源验证、算符收缩/前因子分支及Gamma准备检查，不能仅把当前单插入的`insertion_count`改成2。

本推导的共同起始区域取 $z>0,v_{14}<0,0<\Re e<1/3$；Fourier母式本身仍以解析延拓解释。一般复运动学下必须联合保留尺度和各距离的Feynman处方，不独立选分支。本文没有验证一般复z的物理延拓。

34项精确核验及3组高精度Gamma前因子对照见[理论检查](../../../workflow/checks/self-energy-pi1-chain-theory/README.md)。CAS辅助检查不等于形式整链证明，`formal=not_run`。本轮没有修改外部montopar或qgraf-pipeline生产代码。
