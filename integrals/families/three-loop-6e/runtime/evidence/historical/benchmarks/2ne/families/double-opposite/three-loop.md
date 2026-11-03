# [1,1]双通道对边族：完整生成公式与任意阶M/N约化

日期：2026-09-28。本文是本族的当前统一阅读入口，覆盖三圈6e-32、6e-33、6e-34，并收录两圈4e-06作对照，共四图。完整I、原排序的因子化依据、生成源、任意阶系数规则与已有有限项均放在同页。

**四图均可由同一个一圈M矩源的乘积生成，不需要新增函数族。** 这里的闭包允许M的乘积；不是声称所有系数都属于M的线性张成。显式对数可通过已知恒等式换成M/N多项式。全阶生成清单此前已经把这四图计为达标，本轮补齐自然生成式与逐图推导，不增加38/47或36/47的计数。formal=not_run。

## 1. 先看完整I需要哪些源

取正实z，令

$$
\ell=\log z,\qquad x=\frac{\epsilon\ell}{2},\qquad
Z=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),\qquad
c_n=\cosh(nx).
\tag{D1}
$$

等式先在实 $0<\epsilon<1$ 建立，再对固定z在ε=0附近展开。Z保留当前归一化中唯一的Euler补偿；复z/i0不在本次范围。

共同M源与[统一语言](../../theory/generation/unified-language.md)完全相同：

$$
\begin{gathered}
r=\frac{z}{1+z},\quad W=\log[\beta(1-\beta)],\quad
L=\ell+\log\frac{1-\beta}{\beta},\\
\mathscr U(a,b;z)=\int_0^1\frac{e^{aW}(e^{bL}-1)}{\beta-r}\,d\beta
=\sum_{k\ge0,\,j\ge1}\frac{a^kb^j}{k!j!}M_{k,j}(z),\\
M_{k,j}(z)=\int_0^1\frac{W^kL^j}{\beta-r}\,d\beta.
\end{gathered}
\tag{D2}
$$

内部点β=r的分子为零，按可去点处理。先取a、b足够小，端点也可积。

本族只需要两种由U明确给出的组合：无插入的一圈源A，以及一个真实端点差的源D。定义

$$
\boxed{\displaystyle
\mathcal A(\epsilon,z)=
\frac{\mathscr U(\epsilon/2,\epsilon/2;z)
-\mathscr U(\epsilon/2,-\epsilon/2;z)}{\epsilon}.}
\tag{D3}
$$

它就是一圈的 $\mathcal V_M(\epsilon,z)$；如果采用双参数记号，则是 $\mathcal V_M(\epsilon/2,\epsilon/2;z)$。随后定义σ=±1的镜像对：

$$
\boxed{\displaystyle
\mathcal D_\sigma(\epsilon,z)=\frac1\epsilon
\left[\mathcal A(\epsilon,z)
-\frac{\sigma e^{-\sigma x}}{\epsilon}
\mathscr U(\epsilon,\sigma\epsilon;z)\right].}
\tag{D4}
$$

两个D不是新基函数，也不是未求解积分的改名：(D4)已全部用原U给出，它们满足
$\mathcal D_-(\epsilon,z)=\mathcal D_+(\epsilon,1/z)$。方括号是O(ε)，因此D在ε=0正则。这个相消直接来自原相邻端点，见§3。

## 2. 四个完整I的生成解

下列A、D均在(ε,z)取值，σ求和只有+1与−1两项。

两圈对照图：

$$
\boxed{\displaystyle
I_{4e-06}=-\frac{Z^2}{16c_2}\mathcal A^2.}
\tag{D5}
$$

三圈两张含一个未减除端点的图：

$$
\boxed{\displaystyle
I_{6e-32}=I_{6e-34}
=-\frac{Z^3}{32\epsilon c_3}
\sum_{\sigma=\pm1}e^{\sigma x}
(\mathcal A-\epsilon\mathcal D_\sigma)\mathcal D_\sigma.}
\tag{D6}
$$

三圈两个端点都带差的图：

$$
\boxed{\displaystyle
I_{6e-33}=-\frac{Z^3}{16c_3}
\sum_{\sigma=\pm1}e^{\sigma x}\mathcal D_\sigma^2.}
\tag{D7}
$$

现在可以在未展开时直接看出：4e-06与6e-33在ε=0正则；6e-32/34至多为简单极点。后面的首项核对确认32/34的简单极点确实存在。公式中没有HD、GD、τ或新的待求母函数，全部输入来自(D2)。

## 3. 为什么原积分会生成这些乘积

### 3.1 因子化来自实际变量与排序

两个对边核是

$$
D_z(u,v)=(1-u)(1-v)+zuv,\qquad
\widetilde D_z(u,v)=z(1-u)(1-v)+uv.
$$

本族每条输入包含一个D和一个反向核，且它们的四个端点互不相同。三圈还有一个相邻传播子 $z(1-b)a$。这里必须先依据各图原bounds积分a、b，不能仅凭[1,1]标签假定分离。

例如实际出现的端点积分有

$$
\int_0^s a^{\epsilon-1}da=\frac{s^\epsilon}{\epsilon},\qquad
\int_s^1 a^{\epsilon-1}da=\frac{1-s^\epsilon}{\epsilon}.
\tag{D8}
$$

另一端把s替换为1−s。完成端点积分后，每个权重只依赖其中一个对边核的端点，因此剩余四维积分分离。对实ε>0，放宽到立方体后的绝对值积分可积；这给出本族所需的分组/换序依据，原prepared的外到内排序仍用于输入核对。

### 3.2 单个端点差为什么就是D

定义单核被积函数

$$
k_\epsilon(t,s;z)=-(1+z)z^{-\epsilon/2}D_z(t,s)^{\epsilon-1}.
$$

一圈已有推导给出

$$
\int_0^1\!\!\int_0^1k_\epsilon\,dt\,ds=\mathcal A,
\qquad
\int_0^1\!\!\int_0^1k_\epsilon s^\epsilon\,dt\,ds
=\frac{e^{-x}}{\epsilon}\mathscr U(\epsilon,\epsilon;z).
\tag{D9}
$$

于是从(D8)留下的真实权重直接得到

$$
\boxed{\displaystyle
\int_0^1\!\!\int_0^1k_\epsilon(t,s;z)
\frac{1-s^\epsilon}{\epsilon}\,dt\,ds
=\mathcal D_+(\epsilon,z).}
\tag{D10}
$$

这解释了(D4)的来源：先做原相邻端点积分，再用已经求出的单核源相减。这里没有从Hlog拟合D。由于 $(1-s^\epsilon)/\epsilon\to-\log s$，D在零点是带一个明确对数插入的单核源。

反向核满足 $\widetilde D_z(1-t,1-s)=D_z(t,s)$。反射时同时变换端点权重，或者使用 $\widetilde D_z=zD_{1/z}$ 并保留尺度，都得到(D4)的负号镜像。

### 3.3 从原J回到每个完整I

下表直接用A和D写出完成相邻端点积分后的两核因子；原来的ε分母已经合并到D或表中的显式分母。

| 图 | 原排序的非平凡端点链 | 反射后两核的因子 | 原J |
|---|---|---|---|
| 4e-06 | 无相邻端点 | A，A | $z^\epsilon\mathcal A^2$ |
| 6e-32 | OB，OA | A−εD+，εD+ | $z^{2\epsilon}(\mathcal A-\epsilon\mathcal D_+)\mathcal D_+/\epsilon$ |
| 6e-33 | BO，OA | εD+，εD+ | $z^{2\epsilon}\mathcal D_+^2$ |
| 6e-34 | BO，AO | εD+，A−εD+ | $z^{2\epsilon}(\mathcal A-\epsilon\mathcal D_+)\mathcal D_+/\epsilon$ |

O表示对边端点，A表示 $s^{\epsilon-1}$ 型相邻端点，B表示 $(1-s)^{\epsilon-1}$ 型相邻端点；链按参数从小到大排列。

32与34的两核因子交换，因而J相等。其颜色与轨道重数还必须独立核对：

| 图 | n | a1 | sym | κ |
|---|---:|---:|---:|---:|
| 4e-06 | 2 | −1/2 | 1 | −1/16 |
| 6e-32 | 3 | 1/4 | 4 | −1/16 |
| 6e-33 | 3 | 1/2 | 4 | −1/8 |
| 6e-34 | 3 | 1/4 | 4 | −1/16 |

将表中各自的J代入同一当前归一化

$$
I=\kappa Z^n\frac{J(z)+z^{n\epsilon}J(1/z)}{1+z^{n\epsilon}},
\qquad\kappa=\frac{\mathrm{sym}\,a_1}{2(-2)^n},
\tag{D11}
$$

便逐条得到(D5)–(D7)。因此32与34是完整I的全ε相等，不是仅仅有限项相同。

## 4. 任意阶怎样生成，为什么不需要新函数族

写成普通Taylor级数

$$
\mathcal A=\sum_{q\ge0}a_q\epsilon^q,\qquad
\mathcal D_\sigma=\sum_{q\ge0}d_{\sigma,q}\epsilon^q.
$$

由(D2)直接展开，得到固定M矩中的有限公式：

$$
\boxed{\displaystyle
a_q=\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{M_{q-2j,2j+1}}{2^q(q-2j)!(2j+1)!}.}
\tag{D12}
$$

为书写D的系数，令 $p_{\sigma,q}$ 表示(D4)方括号中第二项的q阶系数，它也完全显式：

$$
\begin{aligned}
p_{\sigma,q}
&=\sum_{\substack{m,k\ge0,\ j\ge1\\m+k+j=q+1}}
\frac{(-1)^m\sigma^{m+j+1}\ell^m}{2^m m!k!j!}M_{k,j},\\
d_{\sigma,q}&=a_{q+1}-p_{\sigma,q+1}.
\end{aligned}
\tag{D13}
$$

$p_{\sigma,0}=a_0=M_{0,1}$，所以没有ε⁻¹残留在D中。这是任意q的推导，不是从有限几阶外推闭包。

所有源系数都是M与显式ℓ的有限线性组合。完整I只对它们做乘法、有限卷积，再乘Gamma/cosh前因子。因此任意Laurent系数属于

$$
\mathbb Q[\ell,\zeta_2,\zeta_3,\ldots]
[M_{k,j}:k\ge0,j\ge1].
\tag{D14}
$$

若要求基函数外不保留ℓ，使用

$$
\ell^2=-2M_{0,1}-6\zeta_2,\qquad
\ell M_{k,j}=M_{k,j+1}-N_{k,j},\qquad
N_{k,j}=\int_0^1\frac{W^k\log[(1-\beta)/\beta]L^j}{\beta-r}\,d\beta.
\tag{D15}
$$

先消去ℓ的偶次幂，再把剩余一次ℓ吸收到一个M因子。这里所有非零物理项均含源矩因子，故得到M/N与ζ常数的有限多项式。**不需要E、τD或新的高核源；也不声称这些乘积代数独立或最小。**

### 完整I的有限卷积公式

令 $P_n=Z^n/c_n=\sum_{r\ge0}p_r^{(n)}\epsilon^r$；其系数由
$\log Z=\sum_{j\ge2}\zeta_j\epsilon^j/j$ 与普通cosh级数直接获得。定义

$$
\begin{aligned}
T_q&=\sum_{\sigma=\pm1}\sum_{i+j+k=q}
\frac{(\sigma\ell/2)^i}{i!}a_jd_{\sigma,k},\\
S_q&=\sum_{\sigma=\pm1}\sum_{i+j+k=q}
\frac{(\sigma\ell/2)^i}{i!}d_{\sigma,j}d_{\sigma,k}.
\end{aligned}
$$

全部求和指标非负；负下标系数记0。完整公式为

$$
\begin{aligned}
[\epsilon^q]I_{4e-06}
&=-\frac1{16}\sum_{r+j+k=q}p_r^{(2)}a_ja_k,\\
[\epsilon^q]I_{6e-32}=[\epsilon^q]I_{6e-34}
&=-\frac1{32}\sum_{r+j=q+1}p_r^{(3)}T_j
+\frac1{32}\sum_{r+j=q}p_r^{(3)}S_j,\\
[\epsilon^q]I_{6e-33}
&=-\frac1{16}\sum_{r+j=q}p_r^{(3)}S_j.
\end{aligned}
\tag{D16}
$$

于是指定任意q，就有有限个明确的M/N乘积；无需先算Hlog再拟合。

## 5. 与低圈图的全ε关系

本族还满足两条保留物理前因子的关系：

$$
\boxed{\displaystyle
I_{4e-06}=-\frac{c_1^2}{4c_2}\,I_{2e-02}^2,}
\tag{D17}
$$

$$
\boxed{\displaystyle
I_{6e-32}+I_{6e-34}+I_{6e-33}
=-\frac{c_1c_2}{c_3}\,I_{2e-02}I_{4e-02},
\qquad I_{4e-02}=I_{4e-05}.}
\tag{D18}
$$

这里所用低圈解是

$$
I_{2e-02}=-\frac{Z\mathcal A(\epsilon)}{2c_1},\qquad
I_{4e-02}=\frac{Z^2}{4\epsilon^2c_2}
[\mathcal A(2\epsilon)-c_1\mathcal A(\epsilon)].
$$

证明只需把(D6)中 $\mathcal A-\epsilon\mathcal D_\sigma$ 展开，平方项与(D7)相消，再用

$$
\sum_{\sigma=\pm1}e^{\sigma x}\mathcal D_\sigma
=\frac{2[c_1\mathcal A(\epsilon)-\mathcal A(2\epsilon)]}{\epsilon}.
$$

这些关系体现双通道乘积结构。系数包含当前C/I对称化的cosh因子，不能把它们省略成裸J的关系。

## 6. 极点与有限项的实际核对

记 $M_0=M_{0,1}$、$M_1=M_{1,1}$、$M_2=M_{2,1}$。使用已知精确恒等式 $M_{0,2}=2\ell M_0/3$，(D13)的首项是

$$
d_{\sigma,0}=-\frac{M_1}{2}+\frac{\sigma\ell M_0}{6}.
$$

于是32/34的简单极点与33的有限首项直接由(D6)/(D7)得到。全表为

| 完整积分 | ε⁻¹系数 | ε⁰系数 |
|---|---|---|
| 4e-06 | 0 | $-M_0^2/16$ |
| 6e-32 | $M_0M_1/32$ | $5M_0^3/1152+7\zeta_2M_0^2/192+3M_0M_2/128+M_1^2/32$ |
| 6e-33 | 0 | $M_0^3/144+\zeta_2M_0^2/48-M_1^2/32$ |
| 6e-34 | $M_0M_1/32$ | 与6e-32相同 |

四图的ε⁻⁶…ε⁻²均为零，表后余项为O(ε)。本页给出任意阶规则，但不把未逐条列出的正ε系数称为已发布的新系数表。

## 7. 验证、旧结果与当前状态

本轮[生成检查](../../../../workflow/checks/double-opposite-generators/README.md)完成：

- 四个完整生成式与当前prepared输入降维式精确连接；原顺序J在ε=1、2、3重新做12个多项式控制。
- A及两个带权源的q=0…4与已有Hlog逐项精确比较；四图全部28个ε⁻⁶…ε⁰位置由新生成式重建。
- U的一维生成积分与独立端点Euler积分比较，并检查四个完整I在三个非整数参数点的归一化，共21个数值控制。
- 旧原域Sobol检查、旧结果和证据保持原字节；本轮没有重新运行旧QMC或HyperInt。数值比较不当作严格区间或Lean证明。

旧结果分别为[4e-06/r002](../../results/4e-06/r002.json)、[6e-32/r001](../../results/6e-32/r001.json)、[6e-33/r001](../../results/6e-33/r001.json)、[6e-34/r001](../../results/6e-34/r001.json)。旧[方法](method.md)、[解表](solutions.md)和[33端点有限性证明](endpoint-finiteness.md)保留作来源；它们当时“只核验至有限项”的措辞不覆盖本页从全ε源得到的任意阶结论。

当前本族4/4已达到自然生成公式与任意阶M/N乘积表示的目标。全项目原有计数保持；更复杂[1]的36/39矩族压缩以及[2]/[3]剩余问题并未因此解决。未提交Git。
