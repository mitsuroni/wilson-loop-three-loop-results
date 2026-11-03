# 单个一圈胶子自能的完整参数化

本文是Π₁的唯一主推导。动量图贡献记为$\mathcal A=C^{\rm src}J$；最终求解的$I(z,\varepsilon)$只由第6节的颜色、尺度及轨道系数提取关系定义。这样与已有2ne/3g流程一致，同时保留完整$C_A/N_F$自能前因子。

## 1. 积分和自能算符

相邻边3、4上的积分定义为

$$
\boxed{
\begin{aligned}
\mathcal A={}&-\int_0^1ds_{3,1}\int_0^1ds_{4,1}\,
(igv_3^{\mu_1}t^{a_1})(igv_4^{\mu_2}t^{a_2})\\
&\times\int\frac{d^D\ell}{(2\pi)^D}
\frac{-ig_{\mu_1\nu_1}\delta^{a_1b_1}}{\ell^2+i0}
\Pi_{b_1b_2}^{\nu_1\nu_2}(\ell)
\frac{-ig_{\nu_2\mu_2}\delta^{b_2a_2}}{\ell^2+i0}
e^{-i\ell\cdot X}.
\end{aligned}}
\tag{1}
$$

采用D=4−2ε、度规 `(+−−−)` 和类光闭合围道。记v_ij=v_i·v_j，s=s31、t=s41，则

$$
\begin{aligned}
v_i^2&=0,\qquad \sum_{i=1}^4v_i=0,\\
x_3(s)&=v_1+v_2+s v_3,\\
x_4(t)&=v_1+v_2+v_3+t v_4,\\
X&=x_3(s)-x_4(t)=-[(1-s)v_3+t v_4].
\end{aligned}
\tag{2}
$$

$\mathcal A$取式(1)的单位颜色矩阵系数，等价于取归一化迹Tr/Nc。最外层负号是式(1)的一部分，始终保留。

Feynman规范下，裸一圈自能写成

$$
\begin{aligned}
\Pi_{b_1b_2}^{\nu_1\nu_2}(\ell)
&=i\delta_{b_1b_2}\mathcal A_1(\varepsilon)
(\ell^2g^{\nu_1\nu_2}-\ell^{\nu_1}\ell^{\nu_2})
(-\ell^2-i0)^{-\varepsilon},\\
\mathcal A_1(\varepsilon)
&=\frac{g^2}{2(4\pi)^{2-\varepsilon}}
\frac{\Gamma(\varepsilon)\Gamma^2(1-\varepsilon)}{\Gamma(2-2\varepsilon)}
\mathcal B_1(\varepsilon),\\
\mathcal B_1(\varepsilon)
&=\frac{(10-6\varepsilon)C_A-8(1-\varepsilon)T_Fn_f}{3-2\varepsilon}.
\end{aligned}
\tag{3}
$$

TF=1/2，nf为无质量夸克味数。算符已经包含胶子、鬼和夸克贡献及其内部符号，且已经包含一个i和一个g²。这里g²=4παs为裸耦合的记号转换；式(1)连同Π₁总阶数为g⁴。

算符系数依据[Davydychev–Osland–Tarasov，附录C](https://arxiv.org/pdf/hep-ph/9605348)；本文以式(3)明确插入约定。

## 2. 收缩并执行 Fourier 变换

颜色与Lorentz指标分别收缩为

$$
\begin{aligned}
t^{a_1}t^{a_2}\delta^{a_1b_1}\delta_{b_1b_2}\delta^{b_2a_2}
&=C_F\mathbf1,\\
v_3^{\mu_1}g_{\mu_1\nu_1}
(\ell^2g^{\nu_1\nu_2}-\ell^{\nu_1}\ell^{\nu_2})
g_{\nu_2\mu_2}v_4^{\mu_2}
&=\ell^2v_{34}-(v_3\cdot\ell)(v_4\cdot\ell).
\end{aligned}
\tag{4}
$$

式(1)在Π之外的相位是−(ig)²(−i)²=−g²。令q=−ℓ²−i0，代入式(3)并整理分母，得到

$$
\mathcal A=ig^2C_F\mathcal A_1
\int_0^1ds\int_0^1dt\int\frac{d^D\ell}{(2\pi)^D}
e^{-i\ell\cdot X}
\left[v_{34}q^{-1-\varepsilon}
+(v_3\cdot\ell)(v_4\cdot\ell)q^{-2-\varepsilon}\right].
\tag{5}
$$

使用同一个Fourier母式

$$
\mathcal F_a(X):=\int\frac{d^D\ell}{(2\pi)^D}
e^{-i\ell\cdot X}(-\ell^2-i0)^{-a}
=i\frac{\Gamma(D/2-a)}{4^a\pi^{D/2}\Gamma(a)}
(-X^2+i0)^{a-D/2}.
\tag{6}
$$

该式由欧氏Schwinger参数高斯积分及Wick转动得到，再按维数作解析延拓。全文的传播子幂均指共同Feynman边界值；坐标底数和i0随这一母式保持一致。

因为两个动量作用于指数等于两个坐标导数乘−1，式(5)化为

$$
\mathcal A=ig^2C_F\mathcal A_1\int_0^1ds\int_0^1dt\,
v_3^\mu v_4^\nu
\left[g_{\mu\nu}\mathcal F_{1+\varepsilon}
-\partial_\mu\partial_\nu\mathcal F_{2+\varepsilon}\right].
\tag{7}
$$

记Q=−X²+i0。两个Fourier系数为

$$
\begin{aligned}
\mathcal F_{1+\varepsilon}&=ih_1Q^{-1+2\varepsilon},
&h_1&=\frac{\Gamma(1-2\varepsilon)}{4^{1+\varepsilon}\pi^{2-\varepsilon}\Gamma(1+\varepsilon)},\\
\mathcal F_{2+\varepsilon}&=ih_2Q^{2\varepsilon},
&h_2&=\frac{\Gamma(-2\varepsilon)}{4^{2+\varepsilon}\pi^{2-\varepsilon}\Gamma(2+\varepsilon)},\\
\frac{h_2}{h_1}&=-\frac1{8\varepsilon(1+\varepsilon)}.
\end{aligned}
\tag{8}
$$

系数比由Gamma递推确定。完整张量求导使用

$$
\partial_\mu\partial_\nu Q^{2\varepsilon}
=-4\varepsilon g_{\mu\nu}Q^{-1+2\varepsilon}
+8\varepsilon(2\varepsilon-1)X_\mu X_\nu Q^{-2+2\varepsilon}.
\tag{9}
$$

将式(8)、(9)直接代入式(7)，得到单一的坐标空间表达式

$$
\boxed{
\mathcal A=-g^2C_F\mathcal A_1h_1
\int_0^1ds\int_0^1dt\,
\left[
\frac{1+2\varepsilon}{2(1+\varepsilon)}v_{34}Q^{-1+2\varepsilon}
-\frac{1-2\varepsilon}{1+\varepsilon}
(v_3\cdot X)(v_4\cdot X)Q^{-2+2\varepsilon}
\right].}
\tag{10}
$$

式(10)包含式(1)的全部张量贡献。其整体负号来自式(7)中的i与Fourier系数中的i相乘。D维迹保持D=4−2ε；由□F_(2+ε)=F_(1+ε)也能检验这里两个坐标系数的连接。

## 3. 无量纲化并得到原参数表示

统一取

$$
r=v_{14},\qquad z=\frac{v_{12}}{v_{14}}=\frac{v_{34}}{v_{14}},
\qquad A=-2r.
\tag{11}
$$

先在A>0、z>0的起始区域提取尺度，再共同解析延拓。对任意SE交换线X=x_i(s_i)−x_j(s_j)，定义

$$
f=\frac{X^2}{2r},\qquad d=\frac{v_i\cdot v_j}{r},\qquad
u=\frac{v_i\cdot X}{r},\qquad w=\frac{v_j\cdot X}{r}.
\tag{12}
$$

它们由类光四边形的无量纲Gram表计算：

$$
G_{ij}=\frac{v_i\cdot v_j}{r}=
\begin{pmatrix}
0&z&-(1+z)&1\\ z&0&1&-(1+z)\\
-(1+z)&1&0&z\\ 1&-(1+z)&z&0
\end{pmatrix}.
\tag{13}
$$

若X=Σk ck vk，则f=cᵀGc/2、u=(Gc)i、w=(Gc)j、d=Gij。四个边向量系数用于表达围道几何，时空张量仍采用D维运算。

在式(10)中代入Q=Af+i0以及r/A=−1/2，整体提出的系数为正。完整无量纲核是

$$
\boxed{\mathcal K(f,d,u,w;\varepsilon)=
\frac{(1+2\varepsilon)df+(1-2\varepsilon)uw}{2(1+\varepsilon)}
(f+i0)^{-2+2\varepsilon}.}
\tag{14}
$$

两个端点的积分因此满足

$$
\mathcal A=C_2^{\rm src}J,\qquad
C_2^{\rm src}=\frac{g^2C_F\mathcal A_1h_1}{2}(A+i0)^{2\varepsilon},
\qquad J=\int_0^1ds_i\int_0^1ds_j\,\mathcal K.
\tag{15}
$$

将式(3)、(8)完整代入C₂：

$$
C_2^{\rm src}=\frac{g^4C_F\mathcal B_1}{4(4\pi)^{2-\varepsilon}}
\frac{\Gamma(\varepsilon)\Gamma^2(1-\varepsilon)}{\Gamma(2-2\varepsilon)}
\frac{\Gamma(1-2\varepsilon)}{4^{1+\varepsilon}\pi^{2-\varepsilon}\Gamma(1+\varepsilon)}
(A+i0)^{2\varepsilon}.
\tag{16}
$$

使用Γ(ε)/Γ(1+ε)=1/ε、Γ(1−2ε)/Γ(2−2ε)=1/(1−2ε)以及g²=4παs，得到

$$
\boxed{C_2^{\rm src}=
\frac{C_F\mathcal B_1(\varepsilon)}{16\varepsilon(1-2\varepsilon)}
\left(\frac{\alpha_s}{\pi}\right)^2
[\pi^\varepsilon\Gamma(1-\varepsilon)]^2
(-2v_{14}+i0)^{2\varepsilon}.}
\tag{17}
$$

这一定义中，C包含所有颜色、耦合、Fourier和尺度因子，J只包含原参数积分核。C保留完整颜色，不包含Euler补偿、轨道重数或反项。

## 4. 相邻图的原参数积分

对式(2)的端点，直接有

$$
X^2=2v_{34}(1-s)t,\quad
v_3\cdot X=-tv_{34},\quad v_4\cdot X=-(1-s)v_{34}.
\tag{18}
$$

于是f=z(1−s)t、d=z、u=−zt、w=−z(1−s)，故uw=zf。式(14)的分子整体简化为

$$
(1+2\varepsilon)df+(1-2\varepsilon)uw=2zf.
\tag{19}
$$

同一个$\mathcal A$的C/J表达式为

$$
\boxed{\mathcal A=C_2^{\rm src}J,\qquad
J=\frac{z}{1+\varepsilon}\int_0^1ds\int_0^1dt\,
[z(1-s)t+i0]^{-1+2\varepsilon}.}
\tag{20}
$$

Re ε>0时，沿原次序执行两个端点积分，每个给出1/(2ε)，因此

$$
J=\frac{z^{2\varepsilon}}{4\varepsilon^2(1+\varepsilon)}.
\tag{21}
$$

代入式(17)，并将尺度合并为同一个边界值，最终得到

$$
\boxed{\mathcal A=
\frac1{64}\left(\frac{\alpha_s}{\pi}\right)^2 C_F\mathcal B_1(\varepsilon)
\frac{\Gamma^2(1-\varepsilon)}{\varepsilon^3(1-2\varepsilon)(1+\varepsilon)}
(-2\pi v_{34}+i0)^{2\varepsilon}.}
\tag{22}
$$

式(1)、(10)、(20)、(22)是同一个积分的不同表达。式(21)先在收敛区域得到，再按亚纯函数延拓至ε=0附近。这里的1/(1+ε)已经包含在完整核的化简中，计算时不得再额外乘一次。

## 5. 对边与三圈图沿用同一规则

对边2/4的端点s=s21、t=s41给出

$$
f=(1-s)(1-t)+zst,\quad d=-(1+z),\quad
u=-1+(1+z)t,\quad w=1-(1+z)s.
\tag{23}
$$

把式(23)直接代入完整核(14)，仍使用C₂。相邻图的分子化简(19)是特定几何关系；一般端点均按完整分子计算。

总三圈的4e-Pi1包含一条SE交换和一条裸交换。裸线连同两个Wilson顶点，剥离颜色后为

$$
-\frac12\left(\frac{\alpha_s}{\pi}\right)\pi^\varepsilon\Gamma(1-\varepsilon)
(-2v_{14}+i0)^\varepsilon\,d_b(f_b+i0)^{-1+\varepsilon}.
\tag{24}
$$

令T_W为所有Wilson生成元按原顺序求得的归一化颜色迹。对两端点T_W=CF；当前十一张四端点来源都是交错配对，T_W=CF(CF−CA/2)。总圈数L=2、3的统一前因子为

$$
\boxed{C_L^{\rm src}(T_W)=\frac{(-1)^{L-2}T_W\mathcal B_1(\varepsilon)}
{2^{L+2}\varepsilon(1-2\varepsilon)}
\left(\frac{\alpha_s}{\pi}\right)^L
[\pi^\varepsilon\Gamma(1-\varepsilon)]^L
(-2v_{14}+i0)^{L\varepsilon}.}
\tag{25}
$$

每张图的$\mathcal A$始终由其原动量规则定义，参数化统一为

$$
\mathcal A=C_L^{\rm src}(T_W)J,\qquad
J=\int_{\Omega}^{\rm ordered}d\boldsymbol s\,
\mathcal K(f_{\rm SE},d_{\rm SE},u_{\rm SE},w_{\rm SE};\varepsilon)
\prod_{b\in\text{裸交换}}d_b(f_b+i0)^{-1+\varepsilon}.
\tag{26}
$$

本批裸交换数为0或1。Ω保留原端点顺序；同边的第k个参数满足0≤s_ik≤s_i,k+1，最外层上限为1。J中的各变量必须按下表顺序解释。

| 图 | SE端点 | 裸线端点 | 原积分外到内顺序 |
|---|---|---|---|
| 2e-Pi1-01 | s31, s41 | — | s31 → s41 |
| 2e-Pi1-02 | s21, s41 | — | s21 → s41 |
| 4e-Pi1-01 | s31, s41 | s32, s42 | s32 → s31 → s42 → s41 |
| 4e-Pi1-02 | s32, s42 | s31, s41 | s32 → s31 → s42 → s41 |
| 4e-Pi1-03 | s21, s41 | s31, s42 | s21 → s31 → s42 → s41 |
| 4e-Pi1-04 | s31, s42 | s21, s41 | s21 → s31 → s42 → s41 |
| 4e-Pi1-05 | s21, s32 | s31, s41 | s21 → s32 → s31 → s41 |
| 4e-Pi1-06 | s31, s41 | s21, s32 | s21 → s32 → s31 → s41 |
| 4e-Pi1-07 | s21, s41 | s22, s42 | s22 → s21 → s42 → s41 |
| 4e-Pi1-08 | s22, s42 | s21, s41 | s22 → s21 → s42 → s41 |
| 4e-Pi1-09 | s21, s31 | s22, s41 | s22 → s21 → s31 → s41 |
| 4e-Pi1-10 | s22, s41 | s21, s31 | s22 → s21 → s31 → s41 |
| 4e-Pi1-11 | s11, s31 | s21, s41 | s11 → s21 → s31 → s41 |

两条交换线在同一有序域上，四端点积分始终按Ω计算。端点、配对和上下限的完整结构化输入见[二圈包](sources/2e-Pi1-r001.json)、[三圈包](sources/4e-Pi1-r001.json)。

## 6. 颜色、Euler准备与唯一的I

上文$C^{\rm src}J$是原动量图贡献。求解沿用2ne/3g的固定基底、一次Euler补偿和C4系数提取。保留原完整颜色向量；对$C_A$与$T_Fn_f$两通道分别重构后选第一坐标，等价于先对Wilson迹作

$$
T_W=\sum_{k=1}^{n-1}a_k C_A^{n-1-k}C_F^k,
\qquad T_W\longmapsto a_1 C_A^{n-2}C_F.
\tag{27}
$$

本批二圈$a_1=1$；三圈完整Wilson坐标为$(−1/2,1)$，所以$a_1=−1/2$。投影不是大Nc近似；未选通道仍保存在输入的完整颜色记录中。定义

$$
\begin{aligned}
G(\varepsilon)&=e^{-\gamma_E\varepsilon}\Gamma(1-\varepsilon),\\
\mathcal N_n^{\Pi_1}&=\left(\frac{\alpha_s}{\pi}\right)^n
C_A^{n-2}C_F\,\mathcal B_1(\varepsilon),\\
S_{14}&=(-2\pi v_{14}+i0)^{n\varepsilon},\qquad S_{12}=z^{n\varepsilon}S_{14}.
\end{aligned}
\tag{28}
$$

其中$\mathcal B_1$是式(3)的整个自能颜色系数，含$C_A$和$T_Fn_f$。在$T_F=1/2$时，它等于$[(3D-2)C_A-2(D-2)N_F]/(D-1)$。因此两种颜色通道共用同一个$I$，没有删去夸克贡献。

经过一次$e^{-n\gamma_E\varepsilon}$准备的C为

$$
C=a_1\mathcal N_n^{\Pi_1}S_{14}R_n,\qquad
R_n=\frac{(-1)^{n-2}G(\varepsilon)^n}{2^{n+2}\varepsilon(1-2\varepsilon)}.
\tag{29}
$$

令m为带SE标记的C4轨道大小，最终$I$的唯一定义是

$$
\boxed{\frac m2 C[J(z,\varepsilon)+z^{n\varepsilon}J(z^{-1},\varepsilon)]
=\mathcal N_n^{\Pi_1}(S_{14}+S_{12})I(z,\varepsilon).}
\tag{30}
$$

等价地，

$$
\boxed{I(z,\varepsilon)=\frac{m a_1}{2}R_n
\frac{J(z,\varepsilon)+z^{n\varepsilon}J(z^{-1},\varepsilon)}{1+z^{n\varepsilon}}.}
\tag{31}
$$

这是形式颜色及微扰系数提取关系。I不再带颜色、耦合、参考尺度或EulerGamma；它们分别按上述公式提出或补偿。不能把源文件的$C^{\rm src}J$直接当成I，也不能对准备好的C再补一次Euler因子。

相邻基准的m=4。把式(21)代入，得到

$$
\boxed{I_{B_A}(\varepsilon)=
\frac{G(\varepsilon)^2}{32\varepsilon^3(1-2\varepsilon)(1+\varepsilon)}.}
\tag{32}
$$

两份[专用准备输入](prepared/README.md)保留原J及完整颜色记录；[B_A](families/b-adjacent/solutions.md)、[B_O](families/b-opposite/solutions.md)给出按式(31)得到的完整I和系数。原始图的符号没有在此修改；历史文献比较记录保留原有对象，不把它们冒充本归一化后的结果。

## 7. montopar与求解输入对应

| 推导结果 | 代码对象 | 计算内容 |
|---|---|---|
| 式(1)–(3) | `validate_source`、`operator`、`phase_ledger` | 一次Π₁、两个内部δ连接、动量相位、总圈数、原负号及Π所含i |
| 式(4)、T_W | `color_data` | 全部Wilson生成元的归一化颜色迹 |
| 式(12)、(13) | `geometry_for_line` | 原端点位置及f、d、u、w |
| 式(14)、(26) | `J_raw.integrand` | 完整分子与距离幂，乘上额外裸交换核 |
| 式(25) | `C.components`、`C.combined` | 完整颜色、耦合、Gamma、有理系数和尺度的乘积 |
| Ω与上表 | `bounds`、`iterated_order_outer_to_inner` | 实际有序域及迭代次序 |

montopar原输出为源C和J；专用准备后以下游C和原J进入式(31)。符号、圈阶和尺度由上述推导确定；轨道重数另行保存，只有进行明确的图求和时才装配。源JSON若采用Φ记相位，则Φ=−X，其指数`exp(+iℓ·Φ)`与式(1)完全一致。

本文的幂运算在规定的起始区域进行，物理复参数下继承统一Feynman处方。Γ函数和Fourier母式使用维数解析延拓；逐图参数化不等于所有剩余积分都已求值。

本归一化的准备、求值和检查见[当前记录](../../../workflow/checks/self-energy-pi1-baselines/README.md)。旧稿、修订和参考比较仅作[历史证据](history.md)，不作为另一份现行推导。当前形式验证状态为`formal=not_run`。
