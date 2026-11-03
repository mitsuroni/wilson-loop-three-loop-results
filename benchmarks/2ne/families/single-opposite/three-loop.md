# [1]单对边族完整解：统一生成函数与18图解表

日期：2026-09-29。本页给出本族18图的已知解及当前表示边界。16图已有简洁M/N生成式；36/39有全阶计算公式，但只用VM、VN、VP的紧凑表示尚未完成。这里先把原Gamma因子和实际Laurent矩组合写清，再讨论还需要什么约化。

阅读顺序：[统一语言](#definitions) → [完整I](#complete-i) → [任意阶规则](#all-orders) → [有限阶解表](#laurent) → [推导](#derivations)。早期矩阵、Volterra和Green路线放在[历史索引](#history)，不再作为正文中的平行解法。上一版原字节保存为[r009](three-loop-r009.md)，更早的完整研究记录仍在r008。

## 1. 范围与“平等对待”的含义

M、N、P都可以作为定义明确的矩族，使用各自的生成函数参与解的组合。这里的平等是**同样明确地定义、生成、展开和使用**，不表示它们已经被证明是彼此独立的最小基底。

| 矩族 | 完整生成源 | 已知关系与用途 |
|---|---|---|
| M | $\mathscr U(a,b;z)$ | 基本的一维对数矩源 |
| N | $\mathscr U_N(a,b;z)$ | 由U的参数导数与减除项得到；保留N使表达简洁 |
| P | $\mathcal V_P(a,b;z)$ | 减除Euler源，与Gauss函数显式相连 |
| 扩展P | $\widehat{\mathcal V}_P(a,b,c;z)$ | 对VP做一次Mellin积分；36/39实际需要 |

所以本页使用M/N与P/扩展P两组生成入口，明确它们的依赖关系，不把扩展P藏在未定义的系数符号中。P与M/N的一般转换关系尚未建立。

| 圈数 | 本族图号 | 全阶生成 | 已证全阶M/N |
|---|---|---:|---:|
| 1 | 2e-02 | 1/1 | 1/1 |
| 2 | 4e-02、05 | 2/2 | 2/2 |
| 3 | 6e-04/05/06、13/14/15/16、26/27/28/29/30/31、36/39 | 15/15 | 13/15 |

固定 $z>0$、$\ell=\log z$，采用实对数，先在小正ε域连接原积分，再在ε=0邻域展开。一般复z/i0未在本页完成。解析计算与辅助核验不等于Lean整链，当前formal=not_run。

<a id="definitions"></a>
## 2. 统一生成语言

### 2.1 归一化与前因子

$$
\begin{gathered}
c_n=\cosh(n\epsilon\ell/2),\qquad
s_1=\frac{\sinh(\epsilon\ell/2)}{\epsilon\ell/2},\\
R_2=\frac{\Gamma(1+\epsilon)^2}{\Gamma(1+2\epsilon)},\qquad
R_{12}=\frac{\Gamma(1+\epsilon)\Gamma(1+2\epsilon)}{\Gamma(1+3\epsilon)}.
\end{gathered}
\tag{C1}
$$

$s_1$在εℓ=0取极限1。完整I中的Gamma因子均直接写出：$e^{-n\gamma_E\epsilon}\Gamma(1-\epsilon)^n$；它就是过去Zⁿ的展开写法，不是新物理因子，不可省略。

### 2.2 M、N及其完整生成源

$$
\begin{gathered}
r=\frac{z}{1+z},\qquad
W(\beta)=\log[\beta(1-\beta)],\qquad
R(\beta)=\log\frac{1-\beta}{\beta},\qquad L=\ell+R,\\
M_{k,j}=\int_0^1\frac{W^kL^j}{\beta-r}\,d\beta,\qquad
N_{k,j}=\int_0^1\frac{W^kRL^j}{\beta-r}\,d\beta,
\quad k\ge0,\ j\ge1.
\end{gathered}
\tag{C2}
$$

内部点β=r处L=0，按连续极限处理，不另加主值处方。R(β)是对数权重，与Gamma比值R₂及整体源$\mathcal R$不同。

$$
N_{k,j}=M_{k,j+1}-\ell M_{k,j},\qquad
\ell^2=-2M_{0,1}-6\zeta_2,\qquad
M_{0,3}=-M_{0,1}^2,\qquad \zeta_j=\zeta(j).
\tag{C3}
$$

$$
\mathscr U(a,b;z)=\int_0^1\frac{e^{aW}(e^{bL}-1)}{\beta-r}\,d\beta
=\sum_{k\ge0,\ j\ge1}\frac{a^kb^j}{k!j!}M_{k,j}.
\tag{C5}
$$

N同样可以配备完整生成函数：

$$
\boxed{\begin{aligned}
\mathscr U_N(a,b;z)
&=\int_0^1\frac{e^{aW}R(e^{bL}-1)}{\beta-r}\,d\beta
=\sum_{k\ge0,j\ge1}\frac{a^kb^j}{k!\,j!}N_{k,j},\\
\mathscr U_N(a,b;z)
&=(\partial_b-\ell)\mathscr U(a,b;z)
-\left.\partial_b\mathscr U(a,b;z)\right|_{b=0}.
\end{aligned}}
\tag{C8}
$$

第二行由R=L−ℓ直接得到。最后的零点项不可省略，否则生成源在b=0不为零。若允许显式log z，N已由M导出；“给N一个生成函数”不等于增加独立函数族。

完整I常用下面两个宇称组合：

$$
\boxed{\begin{aligned}
\mathcal V_M(a,b;z)
&=\int_0^1\frac{e^{aW}}{\beta-r}\frac{\sinh(bL)}b\,d\beta
=\sum_{k,j\ge0}\frac{a^kb^{2j}}{k!(2j+1)!}M_{k,2j+1},\\
\mathcal V_N(a,b;z)
&=\int_0^1\frac{e^{aW}R}{\beta-r}\frac{\cosh(bL)-1}{b^2}\,d\beta
=\sum_{k,j\ge0}\frac{a^kb^{2j}}{k!(2j+2)!}N_{k,2j+2}.
\end{aligned}}
\tag{C4}
$$

它们与完整源的关系为

$$
\mathcal V_M=\frac{\mathscr U(a,b)-\mathscr U(a,-b)}{2b},
\qquad
\mathcal V_N=\frac{\mathscr U_N(a,b)+\mathscr U_N(a,-b)}{2b^2}.
\tag{C9}
$$

b=0均取解析极限。本页V采用双参数约定；旧4e-04单参数VM、VN分别对应本页的VM(ε/2,ε/2)、VN(ε,ε)。为缩短解表，只另设一个完全已知的组合：

$$
\boxed{\displaystyle
\mathcal C_{MN}(\epsilon,z)=
\left.\partial_b[b\mathcal V_M(\epsilon,b;z)]\right|_{b=\epsilon}
-\mathcal V_M(\epsilon,0;z)
-\epsilon^2\mathcal V_N(\epsilon,\epsilon;z).}
\tag{C6}
$$

$$
\mathcal C_{MN}
=\frac{\ell}{2}\big[\mathscr U(\epsilon,\epsilon;z)+\mathscr U(\epsilon,-\epsilon;z)\big],
\quad
[\epsilon^q]\mathcal C_{MN}
=\sum_{j=1}^{\lfloor q/2\rfloor}
\frac{M_{q-2j,2j+1}-N_{q-2j,2j}}{(q-2j)!(2j)!}.
\tag{C7}
$$

先对独立b求导再取b=ε；$\mathcal C_{MN}=O(\epsilon^2)$，它不是新矩族。

### 2.3 P：由Gauss函数自然生成的一维矩

令 $A(t)=\log[t/(1-t)]$、$B(t)=\log(1-t)$、$C_z(t)=\log(1+zt)$。则

$$
\boxed{\displaystyle
\mathsf P_{k,j}(z)=\int_0^1\frac{B(t)^j}{t}
\left[(A(t)-C_z(t))^k-A(t)^k\right]dt,}
\qquad k\ge1,\ j\ge0.
\tag{RG8}
$$

$$
\boxed{\begin{aligned}
\mathcal V_P(a,b;z)
&=\int_0^1\frac{e^{aA(t)+bB(t)}
                  (e^{-aC_z(t)}-1)}{t}\,dt\\
&=\sum_{k\ge1,j\ge0}\frac{a^k b^j}{k!\,j!}\mathsf P_{k,j}(z),\\
{}_2F_1(a,a;1+b;-z)&=1+a\,g(a,b)\mathcal V_P(a,b;z),\\
g(a,b)&=\frac{\Gamma(1+b)}{\Gamma(1+a)\Gamma(1+b-a)}.
\end{aligned}}
\tag{RG9}
$$

t=0处的减除保证可积，t=1只剩可积对数幂。Gauss等式从Euler积分中减去常数1得到；初始收敛条件可取Re a>−1、Re(b−a)>−1的零点邻域。a=0按可去极限解释。最简单的矩为

$$
\begin{aligned}
\mathsf P_{1,0}(z)&=\operatorname{Li}_2(-z),\\
\mathsf P_{2,0}(z)&=2\int_0^z\frac{\log^2(1+x)}x\,dx
                  =4G(0,-1,-1;z),\\
\mathsf P_{1,1}(z)&=G(0,0,-1;z)-G(0,-1,-1;z).
\end{aligned}
\tag{RG11}
$$

### 2.4 扩展P：同一源的一次Mellin积分

$$
\boxed{\begin{aligned}
\widehat{\mathcal V}_P(a,b,c;z)
&=\int_0^1v^{c-1}\mathcal V_P(a,b;zv)\,dv
=\frac{{}_3F_2(a,a,c;1+b,1+c;-z)-1}{a c\,g(a,b)},\\
\widehat{\mathsf P}_{k,j,m}(z)
&=\int_0^1\frac{\log^m v}{v}\mathsf P_{k,j}(zv)\,dv,
\qquad k\ge1,\quad j,m\ge0,\\
\widehat{\mathcal V}_P
&=\sum_{k\ge1,j,m\ge0}
\frac{a^k b^j c^m}{k!\,j!\,m!}\widehat{\mathsf P}_{k,j,m}(z).
\end{aligned}}
\tag{PM6}
$$

这同时给出明确矩定义与已求值的生成函数。a或c为零时取可去极限；c=0可直接使用
$$
\widehat{\mathcal V}_P(a,b,0;z)
=-\frac{az}{(1+b)g(a,b)}
{}_4F_3(1,1,1+a,1+a;2,2,2+b;-z).
\tag{DG4}
$$

先在|z|<1逐项积分，再沿超几何自变量负实轴延拓到z>0；减除积分在Re c>−1的参数零点邻域收敛。

$$
\boxed{\widehat{\mathsf P}_{1,0,m}(z)
=(-1)^m m!\operatorname{Li}_{m+3}(-z).}
\tag{PM7}
$$

$$
\theta\widehat{\mathsf P}_{k,j,0}=\mathsf P_{k,j},
\qquad
\theta\widehat{\mathsf P}_{k,j,m}
=-m\,\widehat{\mathsf P}_{k,j,m-1}\quad(m\ge1),
\tag{PM8}
$$

同时有 $(\theta+c)\widehat{\mathcal V}_P=\mathcal V_P$，其中θ=z∂z。扩展P是固定二重对数矩：把(RG8)代入(PM6)即可，阶数增加不会增加嵌套深度。该明确表示并不证明它能约回一维P。

## 3. 全部18图的完整I
<a id="complete-i"></a>

下列相等与比例均为全ε关系。所有生成源已在§2定义；同一图只保留一套最终公式。

### 3.1 一圈与两圈：本族3图

$$
\boxed{\displaystyle
I_{2e-02}=-\frac{e^{-\gamma_E\epsilon}\Gamma(1-\epsilon)}{2c_1}\mathcal V_M(\epsilon/2,\epsilon/2;z).}
\tag{F1}
$$

$$
\boxed{\displaystyle
I_{4e-02}=I_{4e-05}
=\frac{e^{-2\gamma_E\epsilon}\Gamma(1-\epsilon)^2}{4\epsilon^2c_2}
\left[\mathcal V_M(\epsilon,\epsilon;z)
-c_1\mathcal V_M(\epsilon/2,\epsilon/2;z)\right].}
\tag{F2}
$$

### 3.2 三圈：单端权重与幂差组合

$$
\boxed{\begin{aligned}
I_{6e-04}=I_{6e-26}=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3}{64\epsilon^4c_3}\Big[
&2c_2\mathcal V_M(\epsilon/2,\epsilon/2;z)
-4c_1\mathcal V_M(\epsilon,\epsilon;z)\\
&+3\mathcal V_M(3\epsilon/2,3\epsilon/2;z)
-\mathcal V_M(3\epsilon/2,\epsilon/2;z)
-2s_1\mathcal C_{MN}\Big],\\
I_{6e-05}=I_{6e-28}&=\frac12I_{6e-04}.
\end{aligned}}
\tag{F3}
$$

$$
\boxed{\begin{aligned}
I_{6e-06}=I_{6e-27}=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3}{64\epsilon^4c_3}\Big[
&2c_1\mathcal V_M(\epsilon,\epsilon;z)
-3\mathcal V_M(3\epsilon/2,3\epsilon/2;z)\\
&+\mathcal V_M(3\epsilon/2,\epsilon/2;z)
+s_1\mathcal C_{MN}\Big].
\end{aligned}}
\tag{F4}
$$

### 3.3 三圈：两端异侧权重的四图

$$
\boxed{\begin{aligned}
I_{6e-13}=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3R_2}{16\epsilon^4c_3}\Big[
&\mathcal V_M(\epsilon/2,\epsilon/2;z)
-2c_1\mathcal V_M(\epsilon,\epsilon;z)\\
&+\frac{R_2}{R_{12}}\mathcal V_M(3\epsilon/2,\epsilon/2;z)
+s_1\mathcal C_{MN}\Big],\\
I_{6e-14}&=\frac{2(1-R_2)}{R_2}I_{6e-13}.
\end{aligned}}
\tag{F5}
$$

$$
\boxed{\displaystyle
I_{6e-15}=I_{6e-16}=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3(1-R_2)}{32\epsilon^4c_3}
\left[2c_1\mathcal V_M(\epsilon,\epsilon;z)
-\frac{2R_2}{R_{12}}\mathcal V_M(3\epsilon/2,\epsilon/2;z)
-s_1\mathcal C_{MN}\right].}
\tag{F6}
$$

### 3.4 三圈：同端两个幂权重

$$
\boxed{\displaystyle
I_{6e-29}=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3}{8\epsilon^4c_3}
\left[\mathcal V_M(\epsilon/2,\epsilon/2;z)
-2c_1\mathcal V_M(\epsilon,\epsilon;z)
+\mathcal V_M(3\epsilon/2,\epsilon/2;z)
+s_1\mathcal C_{MN}\right].}
\tag{F7}
$$

$$
\boxed{\displaystyle
I_{6e-30}=I_{6e-31}=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3}{32\epsilon^4c_3}
\left[2c_1\mathcal V_M(\epsilon,\epsilon;z)
-\mathcal V_M(3\epsilon/2,\epsilon/2;z)
-R_2\mathcal V_M(\epsilon/2,\epsilon/2;z)
-s_1\mathcal C_{MN}\right].}
\tag{F8}
$$

<a id="mellin-moments"></a>
### 3.5 三圈：36/39当前还没有简洁的三源形式

这两图的积分已求值，但目前不能写成已证明的、只含VM、VN、VP的简洁组合。上一版把多层组装式放在此处，容易让人误以为这项目标已经完成。**P只来自Gauss核；完整结果当前仍需其Mellin扩展。**

先看可以直接核对的结果：
$$
I_{6e-36}=I_{6e-39}
=\frac{1}{\epsilon^2}
\left(\frac{M_{0,1}\zeta_2}{64}-\frac{M_{2,1}}{128}\right)
+\frac{1}{\epsilon}
\left(-\frac{M_{0,1}M_{1,1}}{64}
+\frac{M_{0,1}\zeta_3}{8}
-\frac{7M_{1,1}\zeta_2}{64}
+\frac{M_{1,3}}{768}
-\frac{M_{3,1}}{128}
+\frac{3N_{1,2}}{128}\right)
+[\epsilon^0]I_{6e-36}+O(\epsilon).
\tag{DG2}
$$

两项极点完全落在M/N。有限项的明确M/N/P/扩展P组合列在[§5.10](#explicit-36-coefficients)，不再用τD4代替待说明的内容。[附录A](#direct-full-source)给出直接代入全部中间量后的全阶公式，保留原Gamma前因子；它仍较长，故不把它称为已经达到目标的紧凑表示。

为什么还不能只写M/N/P？在权重六，将P本身、镜像、全部允许乘积及两个基本对数都加入候选后，τD4仍使精确坐标秩从48增到49。这个结论只针对明确列出的齐次多项式空间；没有穷尽一般有理z系数或变换自变量。它足以说明目前缺的是进一步约化，而不只是排版缩写。见[实际检查](../../../../workflow/checks/6e36-39-direct-display/README.md)。

<a id="all-orders"></a>
## 4. 任意Laurent阶怎样生成

### 4.1 固定矩族的有限和

$$
\begin{aligned}
[\epsilon^q]\mathcal V_M(\lambda\epsilon,\mu\epsilon;z)
&=\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{\lambda^{q-2j}\mu^{2j}}{(q-2j)!(2j+1)!}M_{q-2j,2j+1},\\
[\epsilon^q]\mathcal V_N(\lambda\epsilon,\mu\epsilon;z)
&=\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{\lambda^{q-2j}\mu^{2j}}{(q-2j)!(2j+2)!}N_{q-2j,2j+2}.
\end{aligned}
\tag{R1}
$$

P及扩展P按同样方式展开：

$$
\begin{aligned}
[\epsilon^q]\mathcal V_P(\alpha\epsilon,\beta\epsilon;z)
&=\sum_{\substack{k\ge1,j\ge0\\k+j=q}}
\frac{\alpha^k\beta^j}{k!\,j!}\mathsf P_{k,j}(z),\\
[\epsilon^q]\widehat{\mathcal V}_P(\alpha\epsilon,\beta\epsilon,\gamma\epsilon;z)
&=\sum_{\substack{k\ge1,j,m\ge0\\k+j+m=q}}
\frac{\alpha^k\beta^j\gamma^m}{k!\,j!\,m!}
\widehat{\mathsf P}_{k,j,m}(z).
\end{aligned}
\tag{R5}
$$

零次幂按常数项取1。U、UN和C_MN的通项已在(C5)、(C8)、(C7)给出。Gamma因子g来自
$\log g=\sum_{n\ge2}(-1)^n\zeta_n[b^n-a^n-(b-a)^n]/n$。
M、N、P、扩展P的权重分别为k+j+1、k+j+2、k+j+1、k+j+m+2；权重计数不代替独立性证明。

### 4.2 前因子、卷积与截断预算

$$
\begin{aligned}
\log\!\left[e^{-\gamma_E\epsilon}\Gamma(1-\epsilon)\right]&=\sum_{k\ge2}\frac{\zeta_k}{k}\epsilon^k,\\
\log R_2&=\sum_{k\ge2}\frac{(-1)^k(2-2^k)\zeta_k}{k}\epsilon^k,\\
\log R_{12}&=\sum_{k\ge2}\frac{(-1)^k(1+2^k-3^k)\zeta_k}{k}\epsilon^k,\\
\epsilon d(\epsilon)&=\sum_{k\ge2}[1+(-1)^k2^k]\zeta_k\epsilon^k.
\end{aligned}
\tag{R2}
$$

若 $F=\exp(\sum_{k\ge1}a_k\epsilon^k)=\sum f_r\epsilon^r$，则
$f_0=1$、$f_r=r^{-1}\sum_{k=1}^r k a_kf_{r-k}$。显式幂函数和双曲函数按普通Taylor展开，乘积取有限Cauchy卷积。

求三圈完整I的εᵖ系数时，(F3)–(F8)保守保留正则源至p+4阶；直接使用(DG1)的分块公式可保守保留至p+6阶，再整体相消。不能根据最终只有二阶极点，提前截断中间含ε分母的各项。

36/39已有M、P、扩展P中的有限阶生成规则；但这仍不是只用VM、VN、VP的简洁通解。扩展矩的进一步约化是当前研究缺口。

<a id="laurent"></a>
## 5. 极点至有限项：全部非零系数

记 $I_g=\sum_{p=-6}^0c_{g,p}(z)\epsilon^p+O(\epsilon)$。以下列出全部非零系数；未列出的ε⁻⁶…ε⁰位置为零。相等或常数比例的图按§3复用，不重复长多项式。

以下每项都写出所用矩组合。36/39的有限项目前还需要扩展P；该局限在§5.10明确显示，不再用τ记号遮住。

### 5.1 2e-02

$$
[\epsilon^{0}]I_{\mathrm{2e-02}}=- \frac{M_{0,1}}{2}.
$$

### 5.2 4e-02、4e-05

$$
[\epsilon^{-1}]I_{\mathrm{4e-02}}=\frac{M_{1,1}}{8}.
$$

$$
[\epsilon^{0}]I_{\mathrm{4e-02}}=\frac{M_{0,1}^{2}}{32} + \frac{3 M_{0,1} \zeta(2)}{16} + \frac{3 M_{2,1}}{32}.
$$

### 5.3 6e-04、6e-26、6e-05、6e-28

其中6e-05、28的各阶系数均为6e-04的一半；6e-26与6e-04相同。

$$
[\epsilon^{-2}]I_{\mathrm{6e-04}}=\frac{M_{0,1}^{2}}{384} - \frac{M_{0,1} \zeta(2)}{64} - \frac{M_{2,1}}{128}.
$$

$$
[\epsilon^{-1}]I_{\mathrm{6e-04}}=- \frac{M_{3,1}}{128} - \frac{N_{1,2}}{64}.
$$

$$
\begin{aligned}
[\epsilon^{0}]I_{\mathrm{6e-04}}={}&\frac{103 M_{0,1}^{3}}{23040} - \frac{131 M_{0,1}^{2} \zeta(2)}{7680} - \frac{11 M_{0,1} M_{2,1}}{512}\\
&- \frac{459 M_{0,1} \zeta(2)^{2}}{2560} - \frac{39 M_{2,1} \zeta(2)}{512} - \frac{19 M_{2,3}}{3072}\\
&- \frac{25 M_{4,1}}{6144} - \frac{N_{2,2}}{128}
\end{aligned}
$$

### 5.4 6e-06、6e-27

$$
[\epsilon^{-3}]I_{\mathrm{6e-06}}=\frac{M_{1,1}}{64}.
$$

$$
[\epsilon^{-2}]I_{\mathrm{6e-06}}=\frac{5 M_{0,1}^{2}}{768} + \frac{7 M_{0,1} \zeta(2)}{128} + \frac{5 M_{2,1}}{256}.
$$

$$
\begin{aligned}
[\epsilon^{-1}]I_{\mathrm{6e-06}}={}&\frac{11 M_{0,1} M_{1,1}}{256} + \frac{39 M_{1,1} \zeta(2)}{256} + \frac{19 M_{1,3}}{1536}\\
&+ \frac{19 M_{3,1}}{1536} + \frac{N_{1,2}}{128}
\end{aligned}
$$

$$
\begin{aligned}
[\epsilon^{0}]I_{\mathrm{6e-06}}={}&\frac{617 M_{0,1}^{3}}{46080} + \frac{2351 M_{0,1}^{2} \zeta(2)}{15360} + \frac{49 M_{0,1} M_{2,1}}{1024}\\
&+ \frac{2319 M_{0,1} \zeta(2)^{2}}{5120} + \frac{M_{1,1} \zeta(3)}{64} + \frac{177 M_{2,1} \zeta(2)}{1024}\\
&+ \frac{77 M_{2,3}}{6144} + \frac{65 M_{4,1}}{12288} + \frac{N_{2,2}}{256}
\end{aligned}
$$

### 5.5 6e-13

$$
[\epsilon^{-2}]I_{\mathrm{6e-13}}=- \frac{M_{0,1}^{2}}{192} - \frac{M_{0,1} \zeta(2)}{32} - \frac{M_{2,1}}{64}.
$$

$$
\begin{aligned}
[\epsilon^{-1}]I_{\mathrm{6e-13}}={}&- \frac{M_{0,1} M_{1,1}}{32} + \frac{M_{0,1} \zeta(3)}{4} - \frac{3 M_{1,1} \zeta(2)}{16}\\
&- \frac{M_{1,3}}{64} - \frac{M_{3,1}}{64} + \frac{N_{1,2}}{32}
\end{aligned}
$$

$$
\begin{aligned}
[\epsilon^{0}]I_{\mathrm{6e-13}}={}&- \frac{133 M_{0,1}^{3}}{11520} - \frac{409 M_{0,1}^{2} \zeta(2)}{3840} - \frac{13 M_{0,1} M_{2,1}}{256}\\
&- \frac{691 M_{0,1} \zeta(2)^{2}}{1280} + \frac{3 M_{1,1} \zeta(3)}{8} - \frac{59 M_{2,1} \zeta(2)}{256}\\
&- \frac{13 M_{2,3}}{1536} - \frac{25 M_{4,1}}{3072} + \frac{N_{2,2}}{64}
\end{aligned}
$$

### 5.6 6e-14

$$
[\epsilon^{0}]I_{\mathrm{6e-14}}=- \frac{M_{0,1}^{2} \zeta(2)}{96} - \frac{M_{0,1} \zeta(2)^{2}}{16} - \frac{M_{2,1} \zeta(2)}{32}.
$$

### 5.7 6e-15、6e-16

$$
[\epsilon^{-1}]I_{\mathrm{6e-15}}=\frac{M_{1,1} \zeta(2)}{32}.
$$

$$
\begin{aligned}
[\epsilon^{0}]I_{\mathrm{6e-15}}={}&\frac{M_{0,1}^{2} \zeta(2)}{384} + \frac{3 M_{0,1} \zeta(2)^{2}}{64} - \frac{M_{1,1} \zeta(3)}{16}\\
&+ \frac{5 M_{2,1} \zeta(2)}{128}
\end{aligned}
$$

### 5.8 6e-29

$$
[\epsilon^{-2}]I_{\mathrm{6e-29}}=- \frac{M_{0,1}^{2}}{96} + \frac{M_{0,1} \zeta(2)}{16} - \frac{M_{2,1}}{32}.
$$

$$
\begin{aligned}
[\epsilon^{-1}]I_{\mathrm{6e-29}}={}&- \frac{M_{0,1} M_{1,1}}{16} - \frac{3 M_{1,1} \zeta(2)}{16} - \frac{M_{1,3}}{32}\\
&- \frac{M_{3,1}}{32} + \frac{N_{1,2}}{16}
\end{aligned}
$$

$$
\begin{aligned}
[\epsilon^{0}]I_{\mathrm{6e-29}}={}&- \frac{133 M_{0,1}^{3}}{5760} + \frac{101 M_{0,1}^{2} \zeta(2)}{1920} - \frac{13 M_{0,1} M_{2,1}}{128}\\
&+ \frac{369 M_{0,1} \zeta(2)^{2}}{640} - \frac{45 M_{2,1} \zeta(2)}{128} - \frac{13 M_{2,3}}{768}\\
&- \frac{25 M_{4,1}}{1536} + \frac{N_{2,2}}{32}
\end{aligned}
$$

### 5.9 6e-30、6e-31

$$
[\epsilon^{-2}]I_{\mathrm{6e-30}}=\frac{M_{0,1}^{2}}{384} - \frac{3 M_{0,1} \zeta(2)}{64} + \frac{M_{2,1}}{128}.
$$

$$
\begin{aligned}
[\epsilon^{-1}]I_{\mathrm{6e-30}}={}&\frac{M_{0,1} M_{1,1}}{64} + \frac{M_{0,1} \zeta(3)}{16} + \frac{M_{1,1} \zeta(2)}{32}\\
&+ \frac{M_{1,3}}{128} + \frac{M_{3,1}}{128} - \frac{N_{1,2}}{64}
\end{aligned}
$$

$$
\begin{aligned}
[\epsilon^{0}]I_{\mathrm{6e-30}}={}&\frac{133 M_{0,1}^{3}}{23040} - \frac{631 M_{0,1}^{2} \zeta(2)}{7680} + \frac{13 M_{0,1} M_{2,1}}{512}\\
&- \frac{1101 M_{0,1} \zeta(2)^{2}}{2560} + \frac{M_{1,1} \zeta(3)}{32} + \frac{43 M_{2,1} \zeta(2)}{512}\\
&+ \frac{13 M_{2,3}}{3072} + \frac{25 M_{4,1}}{6144} - \frac{N_{2,2}}{128}
\end{aligned}
$$

<a id="explicit-36-coefficients"></a>
### 5.10 6e-36、6e-39：逐项看到新矩在哪里出现

两图相同，极点项仍是
$$
[\epsilon^{-2}]I_{\mathrm{6e-36}}=\frac{M_{0,1} \zeta(2)}{64} - \frac{M_{2,1}}{128}.
$$

$$
\begin{aligned}
[\epsilon^{-1}]I_{\mathrm{6e-36}}={}&- \frac{M_{0,1} M_{1,1}}{64} + \frac{M_{0,1} \zeta(3)}{8} - \frac{7 M_{1,1} \zeta(2)}{64}\\
&+ \frac{M_{1,3}}{768} - \frac{M_{3,1}}{128} + \frac{3 N_{1,2}}{128}
\end{aligned}
$$

有限项已经展开。这里只用一个标准的偶投影记号：
$$
\mathsf P^+_{k,j}(z)=\frac{\mathsf P_{k,j}(z)+\mathsf P_{k,j}(1/z)}2,\qquad
\widehat{\mathsf P}^{+}_{k,j,m}(z)
=\frac{\widehat{\mathsf P}_{k,j,m}(z)+\widehat{\mathsf P}_{k,j,m}(1/z)}2.
$$

于是同一完整I的有限项为
$$
\begin{aligned}
[\epsilon^0]I_{6e-36}={}& - \frac{43 M_{0,1}^{3}}{46080} - \frac{13 M_{0,5}}{4096} - \frac{25 M_{4,1}}{6144}\\
&-  \frac{M_{2,3}}{192} + \frac{11 N_{2,2}}{512} + \frac{7 \zeta_{3}^{2}}{24}\\
&-  \frac{17 M_{0,1} M_{2,1}}{512} + \frac{3 M_{1,1} \zeta_{3}}{16} - \frac{81 M_{2,1} \zeta_{2}}{512}\\
&+ \frac{91 M_{0,1}^{2} \zeta_{2}}{2560} + \frac{401 M_{0,1} \zeta_{2}^{2}}{2560}\\
&-  \frac{5 \mathsf P^+_{2,3}}{16} - \frac{11 \mathsf P^+_{3,2}}{48} - \frac{29 \mathsf P^+_{4,1}}{384}\\
&-  \frac{\mathsf P^+_{1,4}}{6} - \frac{3 \mathsf P^+_{5,0}}{320}\\
&+ \mathsf P^+_{1,2} \left(- \frac{7 M_{0,1}}{48} + \frac{277 \zeta_{2}}{120}\right)\\
&+ \mathsf P^+_{2,1} \left(- \frac{7 M_{0,1}}{96} + \frac{121 \zeta_{2}}{80}\right) + \frac{449 \mathsf P^+_{3,0} \zeta_{2}}{1440}\\
&+ \frac{17717 \widehat{\mathsf P}^+_{1,0,3}}{3840} - \frac{69 \widehat{\mathsf P}^+_{1,1,0} \zeta_{2}}{40} - \frac{339 \widehat{\mathsf P}^+_{1,1,2}}{80}\\
&-  \frac{917 \widehat{\mathsf P}^+_{1,2,1}}{120} - \frac{23 \widehat{\mathsf P}^+_{1,3,0}}{24} - \frac{2711 \widehat{\mathsf P}^+_{2,1,1}}{480}\\
&-  \frac{37 \widehat{\mathsf P}^+_{2,2,0}}{32} - \frac{431 \widehat{\mathsf P}^+_{3,0,1}}{360} - \frac{47 \widehat{\mathsf P}^+_{3,1,0}}{72}\\
&-  \frac{25 \widehat{\mathsf P}^+_{4,0,0}}{192}
\end{aligned}
\tag{DG3}
$$

这条式子把旧τD4完整消去了，所有成员均有§2的积分定义。它由已导出的源展开后作精确坐标转换得到，并逐项回代原Hlog为零；不是拿转换结果反过来定义P。系数较长本身说明现有P/扩展P语言还没有达到理想的简洁程度。未声称该表达最短，也未声称这些扩展矩分别独立。

<a id="derivations"></a>
## 6. 原积分怎样导出这些生成函数

### 6.1 单对边核与端点权重

原核是 $D_z(t,s)=(1-t)(1-s)+zts$。按原排序积掉相邻传播子的端点后，得到

$$
\mathcal H_{a,b;c,d}=-(1+z)z^{-\epsilon/2}
\int_0^1\!\!\int_0^1D_z^{\epsilon-1}
t^{a\epsilon}(1-t)^{b\epsilon}s^{c\epsilon}(1-s)^{d\epsilon}\,dt\,ds.
$$

| 旧源 | 端点指数(a,b;c,d) | 当前处理 |
|---|---|---|
| H0、H1、H2 | (0,0;m,0)，m=0,1,2 | U的参数取值与物理镜像 |
| HW | (0,0;1,1) | VM(3ε/2,ε/2) |
| HX | (1,0;0,1) | (R₂/R₁₂)HW |
| HD | (1,0;1,0) | (DG1)的P/Mellin生成式 |

这些源的原Meijer G表示与每个prepared根的连接保留在[原公共推导](basis.md)和历史核验中。当前正文给出可直接展开的生成式。

### 6.2 单端H如何消去

对不带插入的端点做线性幂积分，换元 $\beta=1-s$，可得

$$
\begin{aligned}
H_m(\epsilon,z)
&=\frac{e^{-m\epsilon\ell/2}}{\epsilon}
\left[\mathscr U\!\left(\frac{(m+1)\epsilon}{2},\frac{(m+1)\epsilon}{2};z\right)
-\mathscr U\!\left(\frac{(m+1)\epsilon}{2},\frac{(m-1)\epsilon}{2};z\right)\right],\\
H_m(\epsilon,1/z)
&=\frac{e^{m\epsilon\ell/2}}{\epsilon}
\left[\mathscr U\!\left(\frac{(m+1)\epsilon}{2},\frac{(1-m)\epsilon}{2};z\right)
-\mathscr U\!\left(\frac{(m+1)\epsilon}{2},-\frac{(m+1)\epsilon}{2};z\right)\right],
\quad m=0,1,2.
\end{aligned}
$$

其中关键合并为

$$
\begin{aligned}
e^{\epsilon\ell}H_1+e^{-\epsilon\ell}\bar H_1
&=2c_1\mathcal V_M(\epsilon,\epsilon;z)+s_1\mathcal C_{MN},\\
H_1+\bar H_1
&=2c_1\mathcal V_M(\epsilon,\epsilon;z)-s_1\mathcal C_{MN},\\
e^{\epsilon\ell}H_2+e^{-\epsilon\ell}\bar H_2
&=3\mathcal V_M(3\epsilon/2,3\epsilon/2;z)
-\mathcal V_M(3\epsilon/2,\epsilon/2;z).
\end{aligned}
$$

bar表示z→1/z。这些恒等式连同 $H_0=\mathcal V_M(\epsilon/2,\epsilon/2;z)$、
$H_W=\mathcal V_M(3\epsilon/2,\epsilon/2;z)$，实际代入每张完整I，得到§3的单端部分；没有丢掉奇部乘奇指数的贡献。

### 6.3 HX的Gamma比例

先让端点指数独立。对

$$
B_{\eta;\alpha,\beta;\gamma,\delta}
=\int_0^1\!\!\int_0^1D_z^{\eta-1}
t^\alpha(1-t)^\beta s^\gamma(1-s)^\delta\,dt\,ds
$$

做一次Mellin–Barnes分解与两个Beta积分，得到

$$
\begin{aligned}
B_{\eta;\alpha,\beta;\gamma,\delta}
={}&\frac1{\Gamma(1-\eta)\Gamma(1+\eta+\alpha+\beta)\Gamma(1+\eta+\gamma+\delta)}
\int_{\mathcal C}\frac{dw}{2\pi i}\,z^w\\
&\times\Gamma(-w)\Gamma(1-\eta+w)
\Gamma(1+\alpha+w)\Gamma(\eta+\beta-w)
\Gamma(1+\gamma+w)\Gamma(\eta+\delta-w).
\end{aligned}
$$

围道满足

$$
\max(\Re\eta-1,-1-\Re\alpha,-1-\Re\gamma)
<\Re w<\min(0,\Re(\eta+\beta),\Re(\eta+\delta)).
$$

在本文小正参数域可固定实部−1/2。交换α与γ保持Mellin分子和围道，只改变两个外部Beta分母。因此

$$
\boxed{\displaystyle
H_X=\frac{\Gamma(1+\epsilon)\Gamma(1+3\epsilon)}
{\Gamma(1+2\epsilon)^2}H_W=\frac{R_2}{R_{12}}H_W.}
$$

这条未展开恒等式解释了(F5)–(F6)中的比例，不是由有限阶系数拟合得出。

<a id="regular-source"></a>
### 6.4 36/39：端点差 → 简单源方程 → 两个Gauss原函数

由两个相邻端点积分得到 $u_\epsilon(t)=(1-t^\epsilon)/\epsilon$。完整I首先写为

$$
\begin{aligned}
I_{6e-36}=I_{6e-39}
&=\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3(1+z)}{64\epsilon^2c_3}
\left[z^{\epsilon/2}\mathcal R(\epsilon,z)
+z^{-1-\epsilon/2}\mathcal R(\epsilon,1/z)\right],\\
\mathcal R(\epsilon,z)
&=\int_0^1\!\!\int_0^1D_z^{\epsilon-1}u_\epsilon(t)u_\epsilon(s)\,dt\,ds.
\end{aligned}
\tag{S1}
$$

展开两个端点差，$\epsilon^2\mathcal R=F_0-2F_1+F_D$，其中F0、F1是已知M源：

$$
\begin{aligned}
F_0(\epsilon,z)
&=-\frac{z^{\epsilon/2}}{1+z}
\mathcal V_M(\epsilon/2,\epsilon/2;z),\\
F_1(\epsilon,z)
&=-\frac{\mathscr U(\epsilon,\epsilon;z)}{\epsilon(1+z)}.
\end{aligned}
\tag{Q1}
$$

双端源为

$$
F_D=\int_0^1\!\!\int_0^1D_z^{\epsilon-1}t^\epsilon s^\epsilon\,dt\,ds.
$$

§6.3的Mellin核取对应参数，其移位比为
$-(w+1+\epsilon)^2/[(w+1)(w+1-\epsilon)]$。移动围道时保留w=ε残数，得到带源方程

$$
[\theta(\theta-\epsilon)+z(\theta+1+\epsilon)^2]F_D=z^\epsilon.
\tag{S2}
$$

令 $\mathcal Y_D=(1+z)^{1+3\epsilon}F_D$，方程变成

$$
[\theta(\theta-\epsilon)+z(\theta-2\epsilon)^2]\mathcal Y_D
=z^\epsilon(1+z)^{1+3\epsilon}.
\tag{S3}
$$

为了求解这个方程，在推导中引入两个Gauss齐次解；它们最终都会代回VP，不作为完整I的辅助层保留：

$$
\Phi_\epsilon(z)={}_2F_1(-2\epsilon,-2\epsilon;1-\epsilon;-z),
\qquad
\Xi_\epsilon(z)=z^\epsilon{}_2F_1(-\epsilon,-\epsilon;1+\epsilon;-z),
\tag{RG1}
$$

$$
\begin{aligned}
\Phi_\epsilon(z)&=1-2\epsilon\,g(-2\epsilon,-\epsilon)
                         \mathcal V_P(-2\epsilon,-\epsilon;z),\\
\Xi_\epsilon(z)&=z^\epsilon
\left[1-\epsilon\,g(-\epsilon,\epsilon)
                         \mathcal V_P(-\epsilon,\epsilon;z)\right].
\end{aligned}
\tag{RG10}
$$

其Wronskian为
$\Phi\theta\Xi-\Xi\theta\Phi=\epsilon z^\epsilon(1+z)^{3\epsilon}$。
常数变易中，简单源与Wronskian的幂相消，故只剩两个原函数：

$$
V_\epsilon(z)=\int_0^z\frac{\Xi_\epsilon(u)}u\,du
=z^\epsilon\left[\frac1\epsilon+
\mathfrak J(-\epsilon,\epsilon,\epsilon;z)\right].
\tag{PM11}
$$

U积分中的减除1把零端点的对数发散显式分离为log z，不能将这项丢弃。记下面这个由Gamma函数给定的边界常数：

$$
R_2=\frac{\Gamma(1+\epsilon)^2}{\Gamma(1+2\epsilon)},\qquad
d_\epsilon=\gamma_E-\psi(1-\epsilon)+2\psi(1+2\epsilon).
\tag{PM1}
$$

原Mellin表示的w=0简单极点与w=ε双极点固定边界：

$$
\mathcal Y_D=
\frac{R_2^2}{\epsilon^2}
+\frac{z^\epsilon}{\epsilon}
\left(\log z-\frac1\epsilon+d_\epsilon\right)
+O(z)+O(z^{1+\epsilon}\log z).
\tag{PM12}
$$

这固定了两个齐次解的系数；将FD代回端点差，再由(Q1)替换F0/F1并代回VP，就得到附录(DG1)的完整I。全过程没有从旧Hlog拟合生成函数。

(RG9)来自减除Euler积分；(PM6)把Gauss级数第n项乘以1/(n+c)，(DG4)是c=0的可去极限。标准约定参考[NIST Gauss方程](https://dlmf.nist.gov/15.10.E2)、[Euler积分](https://dlmf.nist.gov/15.6.E1)及[广义超几何积分](https://dlmf.nist.gov/16.5.E2)。原先复杂的P/M卷积因此已被解析消去。

<a id="hd-structure"></a>
### 6.5 性质与旧τ/η记号

整体R在Re ε>−1全纯，ε=0的Taylor半径为1；单独HD在ε=−1/2的表观奇点在完整端点差中相消。R还具有已求值的一维正谱表示；这些性质及端点证明见[保留的结构推导](three-loop-r008.md#hd-structure)，不把HD单项的奇点误作完整I的性质。

普通旧系数只是源的偶/奇组合：

$$
\tau_{m,q}=[\epsilon^q]\frac{H_m(z)+H_m(1/z)}2,\qquad
\eta_{m,q}=\ell[\epsilon^q]\frac{H_m(z)-H_m(1/z)}2.
$$

这里不额外乘q!；有序重复字核的阶乘已计入Hm。普通τ可由§6.2直接转换M/N，例如
$\tau_{0,0}=M_{0,1}$、$\tau_{0,1}=M_{1,1}/2$、
$\tau_{0,2}=M_{2,1}/8-M_{0,1}^2/24$。旧双端τD4是HD的四阶偶系数；本页§5.10已将它完整展开，不再将它作为结果中的独立缩写。

<a id="history"></a>
## 7. 来源、验证与历史索引

以下原结果路径保持不变。首项及相等/比例规则帮助核对完整解表。

| 图 | 圈数 | 首项 | 当前原结果 |
|---|---:|---:|---|
| 2e-02 | 1 | $\epsilon^{0}$ | [r002](../../results/2e-02/r002.json) |
| 4e-02 | 2 | $\epsilon^{-1}$ | [r002](../../results/4e-02/r002.json) |
| 4e-05 | 2 | $\epsilon^{-1}$ | [r001](../../results/4e-05/r001.json) |
| 6e-04 | 3 | $\epsilon^{-2}$ | [r002](../../results/6e-04/r002.json) |
| 6e-05 | 3 | $\epsilon^{-2}$ | [r001](../../results/6e-05/r001.json) |
| 6e-06 | 3 | $\epsilon^{-3}$ | [r001](../../results/6e-06/r001.json) |
| 6e-13 | 3 | $\epsilon^{-2}$ | [r001](../../results/6e-13/r001.json) |
| 6e-14 | 3 | $\epsilon^{0}$ | [r001](../../results/6e-14/r001.json) |
| 6e-15 | 3 | $\epsilon^{-1}$ | [r001](../../results/6e-15/r001.json) |
| 6e-16 | 3 | $\epsilon^{-1}$ | [r001](../../results/6e-16/r001.json) |
| 6e-26 | 3 | $\epsilon^{-2}$ | [r001](../../results/6e-26/r001.json) |
| 6e-27 | 3 | $\epsilon^{-3}$ | [r001](../../results/6e-27/r001.json) |
| 6e-28 | 3 | $\epsilon^{-2}$ | [r001](../../results/6e-28/r001.json) |
| 6e-29 | 3 | $\epsilon^{-2}$ | [r001](../../results/6e-29/r001.json) |
| 6e-30 | 3 | $\epsilon^{-2}$ | [r001](../../results/6e-30/r001.json) |
| 6e-31 | 3 | $\epsilon^{-2}$ | [r001](../../results/6e-31/r001.json) |
| 6e-36 | 3 | $\epsilon^{-2}$ | [r002](../../results/6e-36/r002.json) |
| 6e-39 | 3 | $\epsilon^{-2}$ | [r001](../../results/6e-39/r001.json) |

上一轮重编见[历史编辑报告](../../../../workflow/checks/single-family-compact-edit/README.md)。本轮[直接表达与有限项检查](../../../../workflow/checks/6e36-39-direct-display/README.md)核验显式Gamma前因子、直接生成式、展开后的有限项及P空间缺口；没有新增物理积分或改变旧结果。

既有数学证据：[18图整合与126位置](../../../../workflow/checks/single-family-consolidation/README.md)、[单端V转换](../../../../workflow/checks/6e-single-generators/README.md)、[HX/HD](../../../../workflow/checks/6e-double-endpoint/README.md)、[Gauss核](../../../../workflow/checks/6e36-39-gauss-resolvent/README.md)、[完整Mellin矩生成](../../../../workflow/checks/6e36-39-mellin-moments/README.md)。各报告保留自身的适用范围与原始绑定。

<a id="volterra-generator"></a>
<a id="volterra-derivation"></a>
<a id="gauss-resolvent"></a>

旧链接保留在此。以下路线是推导与独立核验的历史材料，不是当前需要逐一代入的平行定义：

| 历史路线 | 原文位置 |
|---|---|
| HD常矩阵与词生成，旧D0–D2/R4 | [r008 §2.3](three-loop-r008.md#definitions) |
| 正则源、正谱与半整数控制，旧P1–P15 | [r008结构](three-loop-r008.md#hd-structure) |
| Volterra递推及边界，旧Q2–Q17 | [r008生成](three-loop-r008.md#volterra-generator)、[推导](three-loop-r008.md#volterra-derivation) |
| Green积分与混合卷积，旧RG2–RG7/RG13–RG14 | [r008 §3.8](three-loop-r008.md#gauss-resolvent) |
| P/M卷积消去的完整研究记录 | [r008 §3.9](three-loop-r008.md#mellin-moments) |

旧[τ解表](solutions.md)、[M/N解表](mn-solutions.md)、[公共推导](basis.md)、[Hlog系数](basis-functions.md)及r001–r007保留原字节。本页的公式编号沿用已保留原式，编号跳跃对应移入历史的过程，不表示缺失定义。

当前本族18/18全阶生成、16/18为M/N；36/39不计为简洁三源形式已完成。全项目44/47算法生成、36/47 M/N及56/329新框架有限阶验收不变。进一步一维化、最小性和形式整链保持开放。下一步先处理36/39的紧凑生成缺口；[2]的10/23顺延。

<a id="direct-full-source"></a>
## 附录A：36/39已知全阶公式的直接代入版

以下把旧式中的Z、YR、YD、Phi、Xi、J及C_MN全部代回，只留下已定义的生成函数与显式Gamma/对数因子。它仍使用扩展源 $\widehat{\mathcal V}_P$，所以不是已经得到“只用VM、VN、VP”的简洁通解。写出它是为了显示当前确切算到了哪里，不以辅助函数改名代替约化。

$$
\begin{aligned}
I_{6e-36}=I_{6e-39}
={}&\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3}
{32\epsilon^4\cosh(3\epsilon\ell/2)}
\Bigg\{
2\cosh(\epsilon\ell/2)\mathcal V_M(\epsilon,\epsilon;z)\\
&\quad-\cosh(\epsilon\ell)\mathcal V_M(\epsilon/2,\epsilon/2;z)\\
&\quad+\frac{2\sinh(\epsilon\ell/2)}{\epsilon\ell}
\left[
\left.\partial_b\!\left(b\mathcal V_M(\epsilon,b;z)\right)\right|_{b=\epsilon}
-\mathcal V_M(\epsilon,0;z)-\epsilon^2\mathcal V_N(\epsilon,\epsilon;z)
\right]\Bigg\}\\
&+\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3}
{64\epsilon^6\cosh(3\epsilon\ell/2)}
\sum_{\xi\in\{z,1/z\}}\frac{\xi^{\epsilon/2}}{(1+\xi)^{3\epsilon}}
\Bigg\{
\frac{\Gamma(1+\epsilon)^4}{\Gamma(1+2\epsilon)^2}
\left[1-2\epsilon g(-2\epsilon,-\epsilon)\mathcal V_P(-2\epsilon,-\epsilon;\xi)\right]\\
&\qquad+\xi^\epsilon
\left[1-\epsilon g(-\epsilon,\epsilon)\mathcal V_P(-\epsilon,\epsilon;\xi)\right]\\
&\qquad\quad\times\left[
\epsilon\log\xi+\epsilon\gamma_E-\epsilon\psi(1-\epsilon)
+2\epsilon\psi(1+2\epsilon)
-2\epsilon^2g(-2\epsilon,-\epsilon)\widehat{\mathcal V}_P(-2\epsilon,-\epsilon,0;\xi)
\right]\\
&\qquad-\xi^\epsilon
\left[1-2\epsilon g(-2\epsilon,-\epsilon)\mathcal V_P(-2\epsilon,-\epsilon;\xi)\right]
\left[1-\epsilon^2g(-\epsilon,\epsilon)\widehat{\mathcal V}_P(-\epsilon,\epsilon,\epsilon;\xi)\right]
\Bigg\}.
\end{aligned}
\tag{DG1}
$$

这里 $\ell=\log z$，g仅为(RG9)给定的Gamma比值。$z=1$处的sinh差商取连续极限。各项的表观ε⁻⁶至ε⁻³在整体中相消；完整I从ε⁻²开始。每个生成函数按其明确矩定义展开即可，但目前这种长式仍不符合所要求的简洁三源形式。

旧组装式与推导中间量保留在[r009](three-loop-r009.md#mellin-moments)，对应等价核验见[本轮报告](../../../../workflow/checks/6e36-39-direct-display/README.md)。
