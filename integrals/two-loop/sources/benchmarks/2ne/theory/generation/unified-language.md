# 统一生成语言与一、两圈2ne完整解

日期：2026-09-28。本文同时作为生成语言的说明和低圈解表：收录2ne中的一圈2图、两圈6图，共8个完整物理积分I。全ε结果统一写成生成函数的组合，展开系数统一使用M/N及ζ常数；τ作为来源明确的系数组合简称保留。

本文整理已有推导并补齐表达之间的代数连接，没有重新定义原积分、改变归一化或重新求解8张图。主线研究的完成状态仍见[全阶生成盘点](all-epsilon-status.md)：32/47达标，双端源及高核的剩余问题不因统一记号而自动解决。

阅读顺序：§1–4说明语言；§5–6列出全部8图的全阶解；§7给极点至有限项；§8给任意阶取系数规则。所有主解等式左边都写完整图编号。

## 1. 固定约定

本文取正实 $z$、$\ell=\log z$，先在小正 $\epsilon$ 的共同收敛域建立等式，再在 $\epsilon=0$ 的小复邻域作局部展开。定义

$$
r=\frac{z}{1+z},\qquad
Z(\epsilon)=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),
\qquad
c_n(\epsilon,z)=\cosh\!\left(\frac{n\epsilon\ell}{2}\right).
\tag{U1}
$$

其中 $Z$ 精确保留项目约定的一次Euler补偿。后文 $c_1,c_2$ 的下标记录圈阶对应的物理双曲前因子，不是Taylor系数。

空集族还需要

$$
R_2(\epsilon)=\frac{\Gamma(1+\epsilon)^2}{\Gamma(1+2\epsilon)},
\qquad
\chi(\epsilon)=\gamma_E+\psi(1+2\epsilon),
\qquad \psi=\frac{\Gamma'}{\Gamma}.
\tag{U2}
$$

$R_2$ 是Gamma比值；它与下面的对数权重 $R$ 不同。$\chi$ 不依赖z，来自4e-04求值时的端点常数。

一维积分中的权重固定为

$$
W(\beta)=\log[\beta(1-\beta)],\qquad
R(\beta)=\log\frac{1-\beta}{\beta},\qquad
L(\beta,z)=\ell+R(\beta).
\tag{U3}
$$

$W$ 在本文始终表示这个对数权重，不表示三核稿中的花体 $\mathcal W$。在内部点 $\beta=r$ 有 $L=0$，后面所有带分母 $\beta-r$ 的式子都按可去极限理解，不引入主值处方。

## 2. 一个共同的矩源：M、N从哪里来

### 2.1 生成函数与M

取两个独立的小源参数a、b，定义

$$
\boxed{\displaystyle
\mathscr U(a,b;z)
=\int_0^1\frac{e^{aW}(e^{bL}-1)}{\beta-r}\,d\beta.}
\tag{U4}
$$

这里 $e^{bL}-1$ 消去内部点的表观极点。端点收敛可先取 $\Re(a-b)>-1$、$\Re(a+b)>-1$、$\Re a>-1$ 的小参数邻域；逐项展开由稍强的端点幂支配。

展开两个源，而不是先拟合某张图的答案：

$$
\boxed{\displaystyle
\mathscr U(a,b;z)
=\sum_{\substack{k\ge0\\j\ge1}}
\frac{a^kb^j}{k!j!}M_{k,j}(z),\qquad
M_{k,j}(z)=\int_0^1\frac{W^kL^j}{\beta-r}\,d\beta.}
\tag{U5}
$$

因此 $M_{k,j}=\left.\partial_a^k\partial_b^j\mathscr U\right|_{a=b=0}$。M的定义不随图号或展开阶数改变；只有指标范围随所求阶数增长。

### 2.2 N是同一矩源的带权导数

保留原定义

$$
N_{k,j}(z)=\int_0^1\frac{W^kRL^j}{\beta-r}\,d\beta,
\qquad k\ge0,\quad j\ge1.
\tag{U6}
$$

由 $R=L-\ell$，直接得到

$$
\boxed{\displaystyle
N_{k,j}=M_{k,j+1}-\ell M_{k,j}
=\left.
\partial_a^k\partial_b^j(\partial_b-\ell)\mathscr U
\right|_{a=b=0}.}
\tag{U7}
$$

所以允许显式ℓ时，N可由M表达；希望最终系数只含有理数、ζ常数而不把ℓ留在基函数外时，保留N很方便。这里不把M与N宣称为代数独立的两套最小基底。

常用的对数消去关系是

$$
M_{0,1}=-\frac{\ell^2+\pi^2}{2},\qquad
\ell^2=-2M_{0,1}-6\zeta_2,\qquad
\ell M_{k,j}=M_{k,j+1}-N_{k,j},
\tag{U8}
$$

其中 $\zeta_j=\zeta(j)$。高次ℓ先消去偶次幂，再处理剩余一次；若出现 $\ell N$，先用(U7)恢复成M与ℓ后重复这一步。

## 3. 低圈实际使用的两个生成函数

### 3.1 一圈生成函数

定义

$$
\boxed{\begin{aligned}
\mathcal V_M(\epsilon,z)
&=\int_0^1\frac{e^{\epsilon W/2}}{\beta-r}
\frac{\sinh(\epsilon L/2)}{\epsilon/2}\,d\beta\\
&=\frac{\mathscr U(\epsilon/2,\epsilon/2;z)
-\mathscr U(\epsilon/2,-\epsilon/2;z)}{\epsilon}.
\end{aligned}}
\tag{U9}
$$

它就是旧单核 $H_0$，与裸2e-02的生成核一致。令 $\mathcal V_M=\sum_{q\ge0}h_q\epsilon^q$，则

$$
\boxed{\displaystyle
h_q=
\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{M_{q-2j,\,2j+1}}
{2^q(q-2j)!(2j+1)!}.}
\tag{U10}
$$

例如

$$
h_0=M_{0,1},\qquad h_1=\frac{M_{1,1}}2,\qquad
h_2=\frac{M_{2,1}}8+\frac{M_{0,3}}{24}.
$$

### 3.2 4e-04中自然出现的带R权重生成函数

定义

$$
\boxed{\displaystyle
\mathcal V_N(\epsilon,z)
=\int_0^1\frac{e^{\epsilon W}R}{\beta-r}
\frac{\cosh(\epsilon L)-1}{\epsilon^2}\,d\beta.}
\tag{U11}
$$

它来自4e-04未展开积分的实际变形，详见[连续Q推导](4e04-direct-mn.md)。内部点可去，$\epsilon=0$ 也可去。

这仍是(U4)的组合。为明确它没有引入另一套任意定义，记 $\mathscr U_b=\partial_b\mathscr U$，则

$$
\mathcal V_N=
\frac{
\mathscr U_b(\epsilon,\epsilon)+\mathscr U_b(\epsilon,-\epsilon)
-2\mathscr U_b(\epsilon,0)
-\ell[\mathscr U(\epsilon,\epsilon)+\mathscr U(\epsilon,-\epsilon)]
}{2\epsilon^2}.
\tag{U12}
$$

各项均在相同z处取值。由 $\mathcal V_N=\sum_{q\ge0}n_q\epsilon^q$，

$$
\boxed{\displaystyle
n_q=\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{N_{q-2j,\,2j+2}}{(q-2j)!(2j+2)!}.}
\tag{U13}
$$

例如 $n_0=N_{0,2}/2$、$n_1=N_{1,2}/2$。$\mathcal V_M,\mathcal V_N$ 是读图方便的生成组合，基础矩源仍是同一个 $\mathscr U$。

## 4. τ在统一语言中的位置

τ的下标记录“哪种端点权重、取第几阶”，M/N的下标记录“对数权重的幂次”。因此τ通常是多个M/N的组合。

令 $H_m$ 是单端插入 $s^{m\epsilon}$ 的母函数，并定义 $H_m^{(+1)}=H_m(\epsilon,z)$、$H_m^{(-1)}=H_m(\epsilon,1/z)$。括号中的符号在此选择自变量方向，不表示偶/奇投影。对本项目使用的 $m=0,1,2$，

$$
\boxed{\displaystyle
H_m^{(\sigma)}
=\frac{e^{-\sigma m\epsilon\ell/2}}{\epsilon}
\left[
\mathscr U\!\left(\frac{(m+1)\epsilon}{2},
                 \frac{(\sigma m+1)\epsilon}{2};z\right)
-\mathscr U\!\left(\frac{(m+1)\epsilon}{2},
                  \frac{(\sigma m-1)\epsilon}{2};z\right)
\right],\quad \sigma=\pm1.}
\tag{U14}
$$

这是原单核完成一个端点积分后得到的权重公式。[逐步连接与任意阶通项](all-epsilon-status.md)已核对两个自变量方向。特别地，

$$
H_0=\mathcal V_M,\qquad
H_1(\epsilon,z)=\frac{e^{-\epsilon\ell/2}}{\epsilon}\mathscr U(\epsilon,\epsilon;z),
\qquad
H_1(\epsilon,1/z)=-\frac{e^{\epsilon\ell/2}}{\epsilon}\mathscr U(\epsilon,-\epsilon;z).
\tag{U15}
$$

旧系数记号就是

$$
\tau_{m,q}=[\epsilon^q]\frac{H_m(\epsilon,z)+H_m(\epsilon,1/z)}2,\qquad
\eta_{m,q}=\ell[\epsilon^q]\frac{H_m(\epsilon,z)-H_m(\epsilon,1/z)}2.
\tag{U16}
$$

这里不额外乘q!。原有序字母核 $T_m$ 满足 $H_m=m!T_m$，因此有序阶乘已经计入。

例如

$$
\tau_{0,0}=M_{0,1},\qquad
\tau_{0,1}=\frac{M_{1,1}}2,\qquad
\tau_{0,2}=\frac{M_{2,1}}8-\frac{M_{0,1}^2}{24}.
$$

这些普通τ可通过(U14)在任意阶转换，不需要逐阶Hlog拟合。扩展的双端标签D/X仍按各自的母函数定义解释，不能套用(U14)。本文8图的最终公式无需D/X，也无需保留τ。

全局采用以下规则：

| 层次 | 主表达 | 用途 |
|---|---|---|
| 原图与归一化 | 完整 $I_g$ | 保持图号、尺度、颜色投影与边界连接 |
| 未展开生成式 | $\mathscr U$及其明确组合 $\mathcal V_M,\mathcal V_N$ | 解释函数族来源，给出全ε表达 |
| Laurent系数 | M/N、ζ及其乘积 | 明确每阶结果及约化范围 |
| 旧τ/η | 带来源的系数组合简称 | 追溯旧解表，并附转换规则 |

## 5. 一圈：全部2个完整解

### 2e-01：空集族

原端点积分与Gamma字积分给出

$$
\boxed{\displaystyle
I_{2e-01}(\epsilon,z)=-\frac{Z(\epsilon)}{\epsilon^2}.}
\tag{I-2e01}
$$

本图完整I不依赖z。

### 2e-02：单对边族

原裸单核等于 $\mathcal V_M$，恢复物理前因子：

$$
\boxed{\displaystyle
I_{2e-02}(\epsilon,z)
=-\frac{Z(\epsilon)}{2c_1(\epsilon,z)}
\mathcal V_M(\epsilon,z).}
\tag{I-2e02}
$$

右边每个因子均已在前文固定，任意阶由(U10)及§8的前因子递推直接得到。

## 6. 两圈：全部6个完整解

### 4e-01：空集族

$$
\boxed{\displaystyle
I_{4e-01}=-\frac{Z^2}{16\epsilon^4}.}
\tag{I-4e01}
$$

### 4e-02：单对边族

旧解表使用 $H_0,H_1$ 及镜像。将(U15)直接代入完整I：

$$
\begin{aligned}
I_{4e-02}
&=-\frac{Z^2}{4\epsilon^2c_2}
\left[
c_1H_0-\frac{e^{\epsilon\ell/2}H_1(\epsilon,z)
+e^{-\epsilon\ell/2}H_1(\epsilon,1/z)}2
\right]\\
&=\boxed{\displaystyle
\frac{Z^2}{4\epsilon^2c_2}
\left[\mathcal V_M(2\epsilon,z)-c_1\mathcal V_M(\epsilon,z)\right].}
\end{aligned}
\tag{I-4e02}
$$

中间的端点权重组合恰好成为 $\mathcal V_M(2\epsilon,z)$。这里仅把裸生成函数的参数变为 $2\epsilon$，没有把外部 $Z(\epsilon)^2$ 换成 $Z(2\epsilon)$。

### 4e-03：空集族的有序字差

$$
\boxed{\displaystyle
I_{4e-03}=-\frac{Z^2(1-R_2)}{4\epsilon^4c_2}.}
\tag{I-4e03}
$$

$1-R_2$ 从二阶开始，所以本图实际只有二阶极点。

### 4e-04：同通道双核

已有[全阶直接推导](4e04-direct-mn.md)将裸同序双核Q求为

$$
Q(\epsilon,z)
=\frac12\mathcal V_M(\epsilon,z)^2
-2\mathcal V_N(\epsilon,z)
+\frac{2\chi(\epsilon)}{\epsilon}\mathcal V_M(2\epsilon,z).
\tag{U17}
$$

代回本图归一化，完整结果是

$$
\boxed{\displaystyle
I_{4e-04}
=-\frac{Z^2}{8c_2}
\left[
\frac12\mathcal V_M(\epsilon,z)^2
-2\mathcal V_N(\epsilon,z)
+\frac{2\chi(\epsilon)}{\epsilon}\mathcal V_M(2\epsilon,z)
\right].}
\tag{I-4e04}
$$

$\chi/\epsilon$ 在零点正则：

$$
\frac{\chi(\epsilon)}{\epsilon}
=\sum_{q\ge0}(-1)^q2^{q+1}\zeta_{q+2}\epsilon^q.
\tag{U18}
$$

因此本图无负Laurent幂。N权重由双核的未展开求值产生，系数直接由(U13)读出。

### 4e-05：与4e-02全阶相等

该相等关系来自完整I映射，不是由有限项相同猜测：

$$
\boxed{\displaystyle
I_{4e-05}=I_{4e-02}
=\frac{Z^2}{4\epsilon^2c_2}
\left[\mathcal V_M(2\epsilon,z)-c_1\mathcal V_M(\epsilon,z)\right].}
\tag{I-4e05}
$$

### 4e-06：两个通道的单核乘积

原域分组与物理归一化给出

$$
\boxed{\displaystyle
I_{4e-06}
=-\frac{Z^2}{16c_2}\mathcal V_M(\epsilon,z)^2.}
\tag{I-4e06}
$$

八张图至此都写成了同一生成语言。没有以新τ符号代替待处理的全ε对象。

## 7. 八图的极点与有限项

以下全部来自上述完整公式。所有未列出的更低负幂为零，各行整体余项为 $O(\epsilon)$。采用 $M_{0,3}=-M_{0,1}^2$ 和(U8)消去冗余矩及外部ℓ。

| 积分 | $\epsilon^{-4}$ | $\epsilon^{-3}$ | $\epsilon^{-2}$ | $\epsilon^{-1}$ | $\epsilon^0$ |
|---|---:|---:|---:|---:|---|
| 2e-01 | 0 | 0 | $-1$ | 0 | $-\zeta_2/2$ |
| 2e-02 | 0 | 0 | 0 | 0 | $-M_{0,1}/2$ |
| 4e-01 | $-1/16$ | 0 | $-\zeta_2/16$ | $-\zeta_3/24$ | $-7\zeta_2^2/160$ |
| 4e-02 | 0 | 0 | 0 | $M_{1,1}/8$ | $(M_{0,1}^2+3M_{2,1}+6\zeta_2M_{0,1})/32$ |
| 4e-03 | 0 | 0 | $-\zeta_2/4$ | $\zeta_3/2$ | $-\zeta_2M_{0,1}/4-49\zeta_2^2/40$ |
| 4e-04 | 0 | 0 | 0 | 0 | $-M_{0,1}^2/48$ |
| 4e-05 | 0 | 0 | 0 | $M_{1,1}/8$ | $(M_{0,1}^2+3M_{2,1}+6\zeta_2M_{0,1})/32$ |
| 4e-06 | 0 | 0 | 0 | 0 | $-M_{0,1}^2/16$ |

例如4e-03旧式的有限项为 $\pi^2\ell^2/48-19\pi^4/1440$，代入(U8)即得到表中形式，未改变归一化。

为了看清N怎样进入，4e-04的一阶是

$$
[\epsilon^1]I_{4e-04}
=-\frac{M_{0,1}M_{1,1}}{16}
+\zeta_3M_{0,1}
-\frac{\zeta_2M_{1,1}}2
+\frac{N_{1,2}}8.
$$

其已有至 $\epsilon^4$ 的完整结果见[高阶解表](4e04-high-order.md)。本页§8给的是任意阶生成规则，不把无限多系数逐项列成表。

## 8. 不再积分，怎样生成任意Laurent阶

### 8.1 所有标量前因子的递推

定义

$$
Z^n=\sum_{r\ge0}z_{n,r}\epsilon^r,\qquad
\frac{Z^n}{c_n}=\sum_{r\ge0}p_{n,r}\epsilon^r,\qquad
R_2=\sum_{r\ge0}r_r\epsilon^r.
$$

$z_{n,r}$ 是系数符号，不是运动学变量z的幂。由Gamma对数级数，

$$
\log Z=\sum_{k\ge2}\frac{\zeta_k}{k}\epsilon^k,\qquad
\log R_2=\sum_{k\ge2}
\frac{(-1)^k(2-2^k)\zeta_k}{k}\epsilon^k.
$$

因此 $z_{n,0}=p_{n,0}=r_0=1$，对 $r\ge1$，

$$
\begin{aligned}
z_{n,r}&=\frac nr\sum_{k=2}^{r}\zeta_k z_{n,r-k},\\
p_{n,r}&=z_{n,r}
-\sum_{k=1}^{\lfloor r/2\rfloor}
\frac{(n\ell/2)^{2k}}{(2k)!}\,p_{n,r-2k},\\
r_r&=\frac1r\sum_{k=2}^{r}(-1)^k(2-2^k)\zeta_k r_{r-k}.
\end{aligned}
\tag{U19}
$$

空和为零。另记 $b_r=\delta_{r,0}-r_r$，即 $1-R_2=\sum b_r\epsilon^r$。

### 8.2 两种组合的系数

使用(U10)、(U13)的 $h_q,n_q$，令

$$
\mathcal V_M(2\epsilon,z)-c_1\mathcal V_M(\epsilon,z)
=\sum_{q\ge0}d_q\epsilon^q .
$$

则

$$
d_q=(2^q-1)h_q
-\sum_{k=1}^{\lfloor q/2\rfloor}
\frac{(\ell/2)^{2k}}{(2k)!}h_{q-2k}.
\tag{U20}
$$

特别地 $d_0=0$、$d_1=h_1$，解释4e-02/05为何只有简单极点。

由(U17)，$Q=\sum_{q\ge0}Q_q\epsilon^q$，其中

$$
Q_q=\frac12\sum_{i=0}^{q}h_i h_{q-i}
-2n_q
+2^{q+2}\sum_{k=0}^{q}(-1)^k\zeta_{k+2}h_{q-k}.
\tag{U21}
$$

### 8.3 八个完整I的通项

对任意整数p，下列所有求和指标均非负；不满足求和约束时取零，负下标系数也取零：

$$
\boxed{\begin{aligned}
[\epsilon^p]I_{2e-01}
&=-z_{1,p+2},\\
[\epsilon^p]I_{2e-02}
&=-\frac12\sum_{r+q=p}p_{1,r}h_q,\\
[\epsilon^p]I_{4e-01}
&=-\frac1{16}z_{2,p+4},\\
[\epsilon^p]I_{4e-02}
=[\epsilon^p]I_{4e-05}
&=\frac14\sum_{r+q=p+2}p_{2,r}d_q,\\
[\epsilon^p]I_{4e-03}
&=-\frac14\sum_{r+s=p+4}p_{2,r}b_s,\\
[\epsilon^p]I_{4e-04}
&=-\frac18\sum_{r+q=p}p_{2,r}Q_q,\\
[\epsilon^p]I_{4e-06}
&=-\frac1{16}\sum_{r+i+j=p}p_{2,r}h_i h_j.
\end{aligned}}
\tag{U22}
$$

给定p，所有求和都有限。所需的h、n由明确的M/N指标公式产生，前因子由(U19)产生，最后用(U8)消去外部ℓ。这就是本页“全阶生成公式”的具体含义。

## 9. 来源、检查及后续扩展

| 图 | 全ε来源与当前原结果 |
|---|---|
| 2e-01 | [空集方法](../../families/empty/method.md)、[r002](../../results/2e-01/r002.json) |
| 2e-02 | [单核生成桥](all-epsilon-status.md)、[r002](../../results/2e-02/r002.json) |
| 4e-01 | [空集解表](../../families/empty/solutions.md)、[r002](../../results/4e-01/r002.json) |
| 4e-02 | [单对边解表](../../families/single-opposite/solutions.md)、[r002](../../results/4e-02/r002.json) |
| 4e-03 | [空集方法](../../families/empty/method.md)、[r002](../../results/4e-03/r002.json) |
| 4e-04 | [直接M/N推导](4e04-direct-mn.md)、[原有限项r002](../../results/4e-04/r002.json)及[高阶补充](4e04-high-order.md) |
| 4e-05 | [单对边解表](../../families/single-opposite/solutions.md)、[r001](../../results/4e-05/r001.json) |
| 4e-06 | [双通道乘积推导](../../families/double-opposite/method.md)、[r002](../../results/4e-06/r002.json) |

本次[连接与展开检查](../../../../workflow/checks/2ne-unified-language/README.md)核对八条完整I、任意阶递推的有限阶控制、八图极点至有限项及数学渲染。旧积分求值与数值证据复用，不重跑全部CAS；formal=not_run，不改变32/47与56/329两种完成数。未提交Git。

后续统一语言的扩展需从原积分导出新权重、边界和任意阶规则，再决定是否增加矩族。HX/HD及[2]/[3]剩余对象仍按[当前路线](all-epsilon-status.md)处理；给Meijer G导数或原多维积分改名，不能代替这一工作。
