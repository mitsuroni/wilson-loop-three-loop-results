# 3e-3g+Pi1-07：对边核心的极点至有限项

## 当前物理积分的入口（2026-10-15，r003）

从共同作用量出发，U′=igAU确定真实Wilson词；FT/Dyson规则在提取C前生效。当前[原生C/J](../../sources/montopar-r002.json)给出下式，J就是本文实际积分的原核，域、分支和IBP全保留：

$$R_P=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad I_P=2R_P\frac{J(z)+z^{3\epsilon}J(1/z)}{1+z^{3\epsilon}}.$$

用下文已经求出的J代入此定义，得到[当前物理I](../../results/3e-3g+Pi1-07/r003.json)；组装系数为+1，Euler只补一次。以下带日期的求值过程与r001公式属于旧根历史，计算J的步骤可按同域积分线性复用，旧形式状态不继承。

2026-10-15表示修订：主I、原J和主Laurent系数保持r002的值；当前源哈希、配对引用及附属公式统一到新修订。当前[原M表示](../../results/3e-3g+Pi1-07/mn-r002.json)已按物理符号重发；旧mn-r001仅供历史复现。[修订说明](../../revisions/r003.md)。以下2026-10-01/02的r001系数表保留历史符号，不能直接当作当前结果。

2026-10-02。07的附着边为244，自能标记在边4的后槽（腿号3）。本次得到指定唯一 $I$ 的全部极点与有限项，余项 $O(\epsilon)$；有限项中的参数积分全部积成字母0、1的Hlog，最高词长5。全 $\epsilon$ 仍保留一个收敛的一维积分，尚未得到全阶积完闭式。结果状态是 `laurent_complete`，formal=`not_run`。

[结果r001](../../results/3e-3g+Pi1-07/r001.json) · [实际检查及复算入口](../../../../../workflow/checks/tgv-pi1-07/README.md)。源是不可变的 [Montopar raw/IBP包](../../sources/montopar-r001.json)，SHA256=`cc9a45b9a3aecd1437fea5fbe222680f67001f8137c09cb38f46ba6fff337e90`。01–06的结果与源原字节保留。

## 1. 同一个I与原排序域

写 $s=s_{21}$、$t=s_{41}$、$u=s_{42}$。原路径顺序是先 $s$、再 $u$、最后 $t$，$0\le t\le u\le1$；其后是 $\alpha_1$ 与 $\alpha_2$，$\alpha_3=1-\alpha_1-\alpha_2$。原权重为 $\alpha_1^{-\epsilon}\alpha_2^{-\epsilon}\alpha_3^{-2\epsilon}$，保留完整两项张量核。

唯一归一化沿用

$$
R(\epsilon)=\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad
I(\epsilon,z)=\frac{2R(\epsilon)}{1+z^{3\epsilon}}\big[J(\epsilon,z)+z^{3\epsilon}J(\epsilon,1/z)\big].
\tag{1}
$$

这里 $z=z_1>0$；C4重数为4，一次三圈Euler补偿已在 $R$。颜色和 $\lambda_b^3(S_{14}+S_{12})$ 留在 $I$ 外面。没有对原raw再做颜色投影。

原核精确等于三个轮廓通量的散度，recoil为零。若 $\nu=1-3\epsilon$，通量为 $b_{0,i}(q+i\eta)^{-\nu}/\nu+b_{1,i}(q+i\eta)^{-\nu-1}$，则两项原分子分别是

$$K_1=\sum_i(\partial_i b_{1,i}-b_{0,i}\partial_iq),\qquad
K_2=-(\nu+1)\sum_i b_{1,i}\partial_iq.$$

本次从实际源重新核验这两条恒等式及 $\sum_i\partial_i b_{0,i}=0$。因此IBP并未丢失：原五维积分变成下边界 $t=0$、上边界 $u=1$ 和同边接触面 $t=u$ 的和，法向符号由原有序域决定。

## 2. 先积α尺度，保留三个完整边界

取

$$\alpha_1=a,\quad\alpha_2=(1-a)v,\quad\alpha_3=(1-a)(1-v),\quad
w(v)=v^{-\epsilon}(1-v)^{-2\epsilon},\quad\delta=3\epsilon-1.$$

Jacobian为 $1-a$。令 $D(s,b)=(1-s)(1-b)+zsb$，原 $q=a(1-a)D(s,vt+(1-v)u)$。每个通量的两幂提出 $a(1-a)$ 后，测度同为 $a^{2\epsilon-1}w(v)\,da\,dv$；剩余分子在 $a$ 上至多一次。因此 $\int_0^1a^{2\epsilon-1}(p_0+p_1a)da=p_0/(2\epsilon)+p_1/(1+2\epsilon)$ 精确积掉尺度。

记 $m=2-3\epsilon$。三个边界在积掉 $a$ 后均为 $\int_0^1dv\,w(v)\int_0^1ds\int_0^1dx\,[N_0D^{-\nu}+N_1D^{-m}]$，其中

| 面 | $b$ | $N_0$ | $N_1$ |
|---|---|---|---|
| $t=0$，$x=u$ | $(1-v)u$ | $-(\epsilon+1)(1+z)/(4\epsilon\delta)$ | $0$ |
| $u=1$，$x=t$ | $1-v(1-t)$ | $-(1+z)[(2\epsilon-3)v+12\epsilon+2]/[8\epsilon(1+2\epsilon)\delta]$ | $A D+B_s h_s+C z$ |
| $t=u$，$x=u$ | $u$ | $(1+z)[4\epsilon^2+(2\epsilon-3)v+18\epsilon+4]/[8\epsilon(1+2\epsilon)\delta]$ | $d[(1+z)D+z]$ |

这里

$$A=\frac{(1+z)(v-2\epsilon)}{8\epsilon(1+2\epsilon)},\quad
B_s=-\frac{zv}{8\epsilon},\quad h_s=(1+z)s-1,$$
$$C=-\frac{(1+4\epsilon)v+2\epsilon}{8\epsilon(1+2\epsilon)},\qquad
d=\frac{v+1}{4(1+2\epsilon)}.$$

共同起始区域取 $0<\Re\epsilon<1/3$、$z>0$。先保留共同 $i\eta$ 和α内部截断作普通IBP，然后移除调节。各通量按 $|(q+i\eta)^{-p}|\le q^{-\Re p}$ 控制（固定复ε时再乘有界相位）。尺度端点为 $a^{2\Re\epsilon-1}$。上边界在 $s=0,b=1$ 附近，积 $s$ 后的最强项由 $v^{2\Re\epsilon-1}(1-t)^{3\Re\epsilon-1}$ 控制；在 $s=1,b=0$ 是 $(1-v)+t$ 的可积角点。接触面的二维角点要求 $m<2$，下边界只含 $\nu<1$。这些条件及 $w$ 的端点条件在该条带共同成立。

因此三个边界均可分别取调节极限；本页后面的换序使用这一实际绝对可积条带。原ordered输入没有被默认为任意乘积测度。结果最后按ε作亚纯延拓。

## 3. 通过混合偏导再降到一维

置

$$G(b)=(zb)^{3\epsilon}-(1-b)^{3\epsilon},\quad h_b=(1+z)b-1,\quad
H(b)=\int_0^1ds\,D(s,b)^{-\nu}=\frac{G(b)}{3\epsilon h_b}.$$

$h_b=0$ 处是可去点，极限为 $(1-b)^{-\nu}$，无需主值处方。第二个分母幂仍可改写为偏导：

$$
zD^{-m}=\frac{3\epsilon(1+z)D^{-\nu}-(3\epsilon)^{-1}\partial_s\partial_bD^{3\epsilon}}{\delta},\qquad
h_sD^{-m}=\frac1\delta\partial_bD^{-\nu}.
\tag{2}
$$

(2)针对已经移除调节的正 $D$；有限 $\eta$ 时不能直接将其照搬为同一个式子。这一步把上边界和接触面进一步分成边界值加余项。相应的 $D^{-\nu}$ 系数化为

$$U_c(v)=-\frac{(1+z)[\delta v+3\epsilon+1]}{4\epsilon\delta},\qquad
C_c(v)=\frac{(1+z)[3(2\epsilon-1)v+4(2\epsilon+1)]}{8\epsilon\delta}.$$

令 $y=1-b$，并用 $\mathrm B$ 表示Euler Beta的亚纯函数：

$$B_*=\mathrm B(1-\epsilon,1-2\epsilon),\quad M_*=\mathrm B(2-\epsilon,1-2\epsilon),$$
$$F(y)=\int_0^yw(v)dv,\quad L(y)=\int_0^y\frac{w(v)}{1-v}dv,\quad
U(y)=\int_y^1\frac{w(v)}v dv.$$

下边界换 $b=(1-v)u$ 后给出 $L(y)$；上边界换 $b=1-v(1-t)$ 后给出 $B_*-F(y)$ 与 $U(y)$；接触面直接给出 $B_*,M_*$。这些都是原域决定的上下限。合并得

$$
J=J_{\rm elem}+J_{\rm reg},\qquad
J_{\rm reg}=\frac1{8\epsilon\delta}\int_0^1db\,H(b)\big[(1+z)\rho(b,\epsilon)+z b^{-2\epsilon}(1-b)^{-\epsilon}\big],
\tag{3}
$$
$$\rho=3(2\epsilon-1)M_*+4(2\epsilon+1)B_*-2(\epsilon+1)L(y)
-2\delta[B_*-F(y)]-2(3\epsilon+1)U(y).$$

上边界的 $B_s h_s$ 项给出 $-z^{3\epsilon}B_*/(24\epsilon^2\delta)$ 和(3)中的额外权重项。混合导数的上边界值是 $[G(1)-G(1-v)]/v$；接触面是 $G(1)-G(0)=1+z^{3\epsilon}$。因此

$$\begin{aligned}
J_{\rm elem}={}&-\frac{z^{3\epsilon}B_*}{24\epsilon^2\delta}
-\frac{(1+z^{3\epsilon})(M_*+B_*)}{12\epsilon\delta(1+2\epsilon)}\\
&+\frac{(1+4\epsilon)\{z^{3\epsilon}[B_*-\mathrm B(1-\epsilon,1+\epsilon)]+\mathrm B(1+2\epsilon,1-2\epsilon)\}
+2\epsilon\{z^{3\epsilon}[\mathrm B(-\epsilon,1-2\epsilon)-\mathrm B(-\epsilon,1+\epsilon)]+\mathrm B(2\epsilon,1-2\epsilon)\}}
{24\epsilon^2(1+2\epsilon)\delta}.
\end{aligned}
\tag{4}
$$

含 $\mathrm B(-\epsilon,\cdot)$ 的两个项来自同一个可积差分 $\int_0^1dv\,w(v)[G(1)-G(1-v)]/v$。可先在共同收敛的辅助端点指数下相减再延拓；不能将这两个Beta值各自解释为当前ε条带中的收敛Euler积分。(3)–(4)是全ε精确一维表示。

## 4. 为什么只需三个可积矩

由于 $R$ 带 $1/\epsilon$，求 $I$ 的有限项必须求 $J$ 至 $\epsilon^1$。(4)是Gamma/Beta，可直接展开；(3)只有显式 $1/\epsilon$，所以其内部积分需要到 $\epsilon^2$。

本节用 $r=1/(1+z)$ 表示运动学位置，**它不是自能腿号3**。令 $\ell=\log z$、$A_b=\ell+\log b$、$B_b=\log(1-b)$，展开

$$\rho(b,\epsilon)+(1-r)b^{-2\epsilon}(1-b)^{-\epsilon}
=\sum_{k=0}^2q_k(b,r)\epsilon^k+O(\epsilon^3).$$

内部一维积分的三个系数为

$$\mathcal F_n(r)=\int_0^1\frac{db}{b-r}\sum_{k=0}^n
\frac{3^{n-k}q_k(b,r)}{(n-k+1)!}
\big[A_b^{n-k+1}-B_b^{n-k+1}\big],\qquad n=0,1,2.
\tag{5}
$$

$A_b-B_b=0$ 在 $b=r$，所以(5)内部没有真实极点。对固定 $0<r<1$，在 $|\Re\epsilon|<\kappa<1/6$ 的邻域，(3)内部 integrand 及有限次ε导数可由端点幂至多 $b^{-5\kappa},(1-b)^{-5\kappa}$ 乘有限次对数控制，均可积。因此可以在显式因子之外逐项展开，余项在该固定运动学处为解析的Taylor余项。

例如不完全Beta $F(y)=y+\epsilon f_1(y)+\epsilon^2 f_2(y)+\cdots$ 的实际系数是

$$f_1=3y-y\log y+2(1-y)\log(1-y),$$
$$f_2=\frac y2\log^2y-3y\log y+9y-2(1-y)\log y\log(1-y)
-2(1-y)\log^2(1-y)+6(1-y)\log(1-y)-2\operatorname{Li}_2(y).$$

其导数分别是 $-[\log y+2\log(1-y)]$ 和 $[\log y+2\log(1-y)]^2/2$，下端均为0。$L,U$ 由字母1的前缀积分生成；$B_*=1+3\epsilon+(9-2\zeta_2)\epsilon^2+\cdots$、$M_*=1/2+7\epsilon/4+(45/8-\zeta_2)\epsilon^2+\cdots$。

HyperInt已经实际积完(5)的三个矩，并转换到仅依赖 $r$ 的Hlog；[F0](../../../../../workflow/checks/tgv-pi1-07/F0.txt)、[F1](../../../../../workflow/checks/tgv-pi1-07/F1.txt)、[F2](../../../../../workflow/checks/tgv-pi1-07/F2.txt)及结果JSON保留完整表达式。开启发散检查，三矩用时约0.69 CPU秒；这不是将未知矩改名。各矩在三个运动学点与直接一维对数积分核对。

它们在原 $J$ 中的贡献分别为 $-\mathcal F_0/(8\epsilon)$、$-(\mathcal F_1+3\mathcal F_0)/8$ 和 $-\epsilon(\mathcal F_2+3\mathcal F_1+9\mathcal F_0)/8$。若写 $J_{\rm elem}=\sum j_n^{\rm elem}\epsilon^n$，实际所需项为

$$\begin{aligned}
j_{-2}^{\rm elem}&=-1/24,\qquad j_{-1}^{\rm elem}=(\ell+1)/8,\\
j_0^{\rm elem}&=3\ell^2/16+3\ell/4+5/12-13\pi^2/144,\\
j_1^{\rm elem}&=3\ell^3/16+9\ell^2/8+(31/16-5\pi^2/48)\ell+137/48-11\pi^2/48-\zeta_3/2.
\end{aligned}$$

## 5. I的四个解析系数

用(1)组装后

$$I_{07}(\epsilon,z)=-\frac1{768\epsilon^3}
+\frac{c_{-2}(z)}{\epsilon^2}+\frac{c_{-1}(z)}\epsilon+c_0(z)+O(\epsilon),\qquad
c_{-2}=\frac5{768}+\frac{M_{0,1}(z)}{64}+\frac{M_{1,1}(z)}{128}.
\tag{6}$$

$\epsilon^{-6},\epsilon^{-5},\epsilon^{-4}$ 的系数均为零。这里 $M$ 沿用项目原定义，而非另建同名对象：

$$M_{k,j}(z)=\int_0^1d\beta\,
\frac{[\log(\beta(1-\beta))]^k[\log z+\log((1-\beta)/\beta)]^j}{\beta-z/(1+z)},\quad j\ge1.$$

内点分子消零。换 $\beta=1-b$ 后，已用新的HyperInt积分精确核验(6)和实际Hlog系数相等。

### 5.1 M形式（2026-10-02补充）

$c_{-1},c_0$ 都可以写回原有M/N函数，不需要额外特殊函数；进一步使用已有矩恒等式后，N也可在这两个式子中消去。这里允许明确的有理 $z$ 系数。记

$$\chi(z)=\frac{1-z}{1+z}=2r-1,\qquad \zeta_k=\zeta(k),$$

所有M均在同一个 $z$ 处取值，则完整表达为

$$\boxed{\displaystyle
c_{-1}=\frac1{192}-\frac{7\zeta_2}{512}
+\frac{(87-12\zeta_2)M_{0,1}+43M_{1,1}}{1024}
+\frac{3M_{2,1}}{512}+\frac{\chi M_{0,2}}{2048}.}
\tag{7}
$$

$$\boxed{\begin{aligned}
c_0={}&\frac{35}{512}+\frac{11\zeta_2}{512}-\frac{19\zeta_3}{768}
+\frac{(667+207\zeta_2)M_{0,1}+331M_{1,1}}{2048}
+\frac{\chi M_{0,2}}{2048}\\
&+\frac{129M_{2,1}+115M_{0,1}^2}{4096}
+\frac{3M_{3,1}}{1024}+\frac{5M_{1,3}}{1536}
+\frac{33\zeta_2M_{1,1}}{512}
-\frac{\zeta_3M_{0,1}}{32}+\frac{9M_{0,1}M_{1,1}}{512}.
\end{aligned}}
\tag{8}
$$

M/N的识别是在实际积分完成后进行的等价表示转换。重新按原定义积分权重至5的十个M矩，再用精确shuffle消元；(7)–(8)已在新HyperInt会话完整代回旧Hlog，差均为零。另外直接数值积分原M定义，在 $z=1/2,1,2$ 对照两系数，最大差约 $3.71\times10^{-29}$。[等价表示JSON](../../results/3e-3g+Pi1-07/mn-r001.json) · [核验记录](../../../../../workflow/checks/tgv-pi1-07-mn/README.md)。原r001结果和Hlog表没有改写。

中间M/N式使用 $N_{0,2}$，由 $N_{k,j}=M_{k,j+1}-\ell M_{k,j}$ 及

$$M_{0,3}=-M_{0,1}^2,\qquad M_{0,2}=\frac23\ell M_{0,1},\qquad
N_{0,2}=\frac13M_{0,1}^2+4\zeta_2M_{0,1}$$

得到上述纯M式。三个矩恒等式也在本次精确代回核验。$M_{0,2}$ 与 $\chi$ 都在crossing下变号，所以它们的乘积仍是偶函数。

严格区分两种范围：(7)–(8)属于允许有理 $z$ 系数的M/N多项式表示；若将系数限定为有理常数与ζ，仅枚举权重≤5的原M/N乘积，本次固定Hlog/周期坐标检查不能覆盖共同项 $\chi M_{0,2}/2048$。不把这个范围区别误称为需要新的特殊函数，也不由两阶约化推断全εM/N闭包。

### 5.2 Hlog系数表

原积分的Hlog结果继续保留作另一套精确坐标。定义 $H_w(r)=\operatorname{Hlog}(r,w)$，递归核为 $dt/(t-a)$，$H_0=\log r$、$H_1=\log(1-r)$，$H_{0^n}=\log^nr/n!$。取正实区间上的切向基点正则化；特别注意字母1的符号约定。令 $\zeta_k=\zeta(k)$，则

$$c_p(z)=K_p+\sum_w A_{p,w}(r)H_w(r),\qquad p=-1,0,\quad r=\frac1{1+z}.$$

下表是计算结果的全部非零词系数（列内0表示该词不出现），无剩余积分。

| 字 $w$ | $A_{-1,w}(r)$ | $A_{0,w}(r)$ |
|---|---|---|
| 常数 $K_p$ | $- \frac{9 \zeta_{2}^{2}}{2560} - \frac{275 \zeta_{2}}{1024} + \frac{43 \zeta_{3}}{512} + \frac{1}{192}$ | $- \frac{1323 \zeta_{2}^{2}}{5120} + \frac{31 \zeta_{2} \zeta_{3}}{256} - \frac{1957 \zeta_{2}}{2048} + \frac{917 \zeta_{3}}{3072} + \frac{37 \zeta_{5}}{256} + \frac{35}{512}$ |
| $0$ | $\frac{\zeta_{2} r - 65 \zeta_{2} + 12 \zeta_{3}}{512}$ | $\frac{- 1254 \zeta_{2}^{2} + 20 \zeta_{2} r - 4975 \zeta_{2} + 1290 \zeta_{3}}{10240}$ |
| $1$ | $- \frac{\zeta_{2} r + 64 \zeta_{2} - 12 \zeta_{3}}{512}$ | $- \frac{1254 \zeta_{2}^{2} + 20 \zeta_{2} r + 4955 \zeta_{2} - 1290 \zeta_{3}}{10240}$ |
| $0,0$ | $- \frac{3 \left(4 \zeta_{2} + 29\right)}{1024}$ | $- \frac{120 \zeta_{2} - 72 \zeta_{3} + 667}{2048}$ |
| $0,1$ | $- \frac{3 \left(20 \zeta_{2} - 29\right)}{1024}$ | $- \frac{654 \zeta_{2} - 72 \zeta_{3} - 667}{2048}$ |
| $1,0$ | $- \frac{3 \left(20 \zeta_{2} - 29\right)}{1024}$ | $- \frac{654 \zeta_{2} - 72 \zeta_{3} - 667}{2048}$ |
| $1,1$ | $- \frac{3 \left(4 \zeta_{2} + 29\right)}{1024}$ | $- \frac{120 \zeta_{2} - 72 \zeta_{3} + 667}{2048}$ |
| $0,0,0$ | $\frac{r - 22}{512}$ | $\frac{188 \zeta_{2} + 4 r - 333}{2048}$ |
| $0,0,1$ | $- \frac{r - 22}{512}$ | $- \frac{188 \zeta_{2} + 4 r - 333}{2048}$ |
| $0,1,0$ | $- \frac{r - 22}{512}$ | $- \frac{188 \zeta_{2} + 4 r - 333}{2048}$ |
| $0,1,1$ | $\frac{r - 22}{512}$ | $\frac{- 244 \zeta_{2} + 4 r - 333}{2048}$ |
| $1,0,0$ | $- \frac{r + 21}{512}$ | $- \frac{244 \zeta_{2} + 4 r + 329}{2048}$ |
| $1,0,1$ | $\frac{r + 21}{512}$ | $\frac{- 188 \zeta_{2} + 4 r + 329}{2048}$ |
| $1,1,0$ | $\frac{r + 21}{512}$ | $\frac{- 188 \zeta_{2} + 4 r + 329}{2048}$ |
| $1,1,1$ | $- \frac{r + 21}{512}$ | $- \frac{- 188 \zeta_{2} + 4 r + 329}{2048}$ |
| $0,0,0,0$ | $- \frac{3}{256}$ | $\frac{27}{256}$ |
| $0,0,0,1$ | $\frac{3}{256}$ | $- \frac{27}{256}$ |
| $0,0,1,0$ | $\frac{3}{256}$ | $- \frac{27}{256}$ |
| $0,0,1,1$ | $- \frac{3}{256}$ | $\frac{27}{256}$ |
| $0,1,0,0$ | $- \frac{3}{256}$ | $- \frac{237}{1024}$ |
| $0,1,0,1$ | $\frac{3}{256}$ | $\frac{237}{1024}$ |
| $0,1,1,0$ | $\frac{3}{256}$ | $\frac{237}{1024}$ |
| $0,1,1,1$ | $- \frac{3}{256}$ | $- \frac{237}{1024}$ |
| $1,0,0,0$ | $- \frac{3}{256}$ | $- \frac{237}{1024}$ |
| $1,0,0,1$ | $\frac{3}{256}$ | $\frac{237}{1024}$ |
| $1,0,1,0$ | $\frac{3}{256}$ | $\frac{237}{1024}$ |
| $1,0,1,1$ | $- \frac{3}{256}$ | $- \frac{237}{1024}$ |
| $1,1,0,0$ | $- \frac{3}{256}$ | $\frac{27}{256}$ |
| $1,1,0,1$ | $\frac{3}{256}$ | $- \frac{27}{256}$ |
| $1,1,1,0$ | $\frac{3}{256}$ | $- \frac{27}{256}$ |
| $1,1,1,1$ | $- \frac{3}{256}$ | $\frac{27}{256}$ |
| $0,0,0,0,0$ | $0$ | $\frac{71}{512}$ |
| $0,0,0,0,1$ | $0$ | $- \frac{71}{512}$ |
| $0,0,0,1,0$ | $0$ | $- \frac{71}{512}$ |
| $0,0,0,1,1$ | $0$ | $\frac{71}{512}$ |
| $0,0,1,0,0$ | $0$ | $- \frac{71}{512}$ |
| $0,0,1,0,1$ | $0$ | $\frac{71}{512}$ |
| $0,0,1,1,0$ | $0$ | $\frac{71}{512}$ |
| $0,0,1,1,1$ | $0$ | $- \frac{71}{512}$ |
| $0,1,0,0,0$ | $0$ | $- \frac{35}{512}$ |
| $0,1,0,0,1$ | $0$ | $\frac{35}{512}$ |
| $0,1,0,1,0$ | $0$ | $\frac{35}{512}$ |
| $0,1,0,1,1$ | $0$ | $- \frac{35}{512}$ |
| $0,1,1,0,0$ | $0$ | $- \frac{1}{512}$ |
| $0,1,1,0,1$ | $0$ | $\frac{1}{512}$ |
| $0,1,1,1,0$ | $0$ | $\frac{1}{512}$ |
| $0,1,1,1,1$ | $0$ | $- \frac{1}{512}$ |
| $1,0,0,0,0$ | $0$ | $- \frac{1}{512}$ |
| $1,0,0,0,1$ | $0$ | $\frac{1}{512}$ |
| $1,0,0,1,0$ | $0$ | $\frac{1}{512}$ |
| $1,0,0,1,1$ | $0$ | $- \frac{1}{512}$ |
| $1,0,1,0,0$ | $0$ | $- \frac{35}{512}$ |
| $1,0,1,0,1$ | $0$ | $\frac{35}{512}$ |
| $1,0,1,1,0$ | $0$ | $\frac{35}{512}$ |
| $1,0,1,1,1$ | $0$ | $- \frac{35}{512}$ |
| $1,1,0,0,0$ | $0$ | $- \frac{71}{512}$ |
| $1,1,0,0,1$ | $0$ | $\frac{71}{512}$ |
| $1,1,0,1,0$ | $0$ | $\frac{71}{512}$ |
| $1,1,0,1,1$ | $0$ | $- \frac{71}{512}$ |
| $1,1,1,0,0$ | $0$ | $\frac{71}{512}$ |
| $1,1,1,0,1$ | $0$ | $- \frac{71}{512}$ |
| $1,1,1,1,0$ | $0$ | $- \frac{71}{512}$ |
| $1,1,1,1,1$ | $0$ | $\frac{71}{512}$ |


四个 $c_p$ 均已在Hlog基底精确检验 $r\leftrightarrow1-r$，即 $z\leftrightarrow1/z$。数值例：

| $z$ | $c_{-2}$ | $c_{-1}$ | $c_0$ |
|---:|---:|---:|---:|
| 1 | 0.01571887935280993 | −0.05752769538672054 | 0.1945672459628129 |
| 2 | 0.01661605460751705 | −0.05894646374058299 | 0.1969399530454971 |

这些数值用于控制解析结果，不作为数值拟合的来源。

## 6. 两个颜色通道与完成范围

完整颜色是 $C_FC_A B_1$，其中 $B_1=[(10-6\epsilon)C_A-8(1-\epsilon)T_Fn_f]/(3-2\epsilon)$。因此只有这两个独立颜色通道：

$$g(\epsilon,z)=\frac{10-6\epsilon}{3-2\epsilon}I_{07},\qquad
f(\epsilon,z)=-\frac{8(1-\epsilon)}{3-2\epsilon}I_{07}.$$

如果 $I=\sum_{n=-3}^0c_n\epsilon^n+O(\epsilon)$，两通道所需的全部卷积可简洁写成

$$\begin{array}{c|c|c}
n&g_n&f_n\\\hline
-3&\frac{10}3c_{-3}&-\frac83c_{-3}\\
-2&\frac{10}3c_{-2}+\frac29c_{-3}&-\frac83c_{-2}+\frac89c_{-3}\\
-1&\frac{10}3c_{-1}+\frac29c_{-2}+\frac4{27}c_{-3}&-\frac83c_{-1}+\frac89c_{-2}+\frac{16}{27}c_{-3}\\
0&\frac{10}3c_0+\frac29c_{-1}+\frac4{27}c_{-2}+\frac8{81}c_{-3}&-\frac83c_0+\frac89c_{-1}+\frac{16}{27}c_{-2}+\frac{32}{81}c_{-3}
\end{array}$$

最高留数为 $g_{-3}=-5/1152$、$f_{-3}=1/288$。八个实际展开系数已保存到结果JSON，并由SymPy独立核验；乘回未展开的颜色因子后，两个运动学点的余项也随ε按一阶缩小。尺度因子仍留在外面，不能把 $c_0$ 直接当作三圈重整化后的物理有限项。

本次另外核验了原raw/flux连接13项、降维12项、四个crossing及一个M基底等式；9个系数矩数值点、两组原IBP边界的逐级Gauss积分、6个全ε/Laurent对照及12个颜色余项点支持结果。原边界对照是在精确积掉a和s后的二维数值积分，不是独立五维raw积分；末次差约 $7\times10^{-5}$，没有严格误差包围。

原源、01–06旧结果、生产模块和Kernel/Lean均保持。当前无可调用Wilson数学MCP，所以没有形式收据，formal=`not_run`。本分支完成极点至有限项7/12，其中01–06另有全ε闭式，07只有(3)–(4)的全ε一维表示。下一图08仍须核对实际标记权重、分子和排序域，不能自动继承07。
