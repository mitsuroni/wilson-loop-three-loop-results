# T1[1]：两个通道的原核归约

对象为T1[1]₂₄的5e-3g-11、12，以及T1[1]₁₃的5e-3g-32。这个目录联合展示两组计算，不更改既有分类。以下图号均指5e-3g。唯一I采用prepared r001；固定z>0，先在共同绝对收敛条带1/2<Re(e)<1换元和换序，再作e=ε的亚纯延拓。formal=not_run。

## 1. 从原始S、P和排序域消元

TGV孤立端点t在边3，两同边端点0<x<y<1在边4。令u=1-t。三图的源多项式均为
$$
S=zu\alpha_1(\alpha_2x+\alpha_3y),\qquad
P=\frac{z^2u}{4}\alpha_1(\alpha_2-\alpha_3).
$$
取alpha1=a、alpha2=(1-a)b、alpha3=(1-a)(1-b)，Jacobian=1-a。a的总幂为e-1，积分给出1/e。剩余核是
$$
\frac{z^{2e}}{4e}\,u^{2e-1}
b^{-e}(1-b)^{-e}(2b-1)[bx+(1-b)y]^{-2+2e}.
$$
这一步与[T1[]](../t1-empty/method.md)相同，但新pure线连接对边，其Lorentz因子为−(1+z)，不能继续套邻边幂积分。

定义正实域上的对边距离
$$D_z(s,w)=(1-s)(1-w)+zsw.$$
11/12的s在边2、w在边4；32令w=1-s_p3、s=1-s_p1后也得到同一D。原排序为：

| 图 | 与TGV相连边的额外排序 | 可先积的TGV孤立端点 |
|---|---|---|
| 11 | x<w<y | u∈(0,1)，给出1/(2e) |
| 12 | w<x<y | 同11 |
| 32 | 0<w<u<1，x<y | u∈(w,1)，给出(1-w^(2e))/(2e) |

记
$$
h_0(e)=\int_0^1dr\,db\,
b^{-e}(1-b)^{-e}(2b-1)(1-b+br)^{-2+2e}
=\frac{U(e)}{(1-2e)e},
$$
$$
U(e)=2\Gamma(1+e)\Gamma(1-e)-B(1-e,1-e).
$$
完整xy三角域的积分为h0/(2e)。把它按固定w分为x<y<w、x<w<y、w<x<y三块；第一块经(x,y)→w(x,y)带出w^(2e)。因此后两块之和恰好是
$$\frac{h_0}{2e}(1-w^{2e}).$$
与32的u积分比较得到一个原始有序域恒等式：
$$\boxed{J_{11}+J_{12}=J_{32}.}$$
实际三图C、颜色坐标(−1/4,1/2,0)、sym=4及一次Euler补偿相同，故同样有完整ε关系
$$\boxed{I_{11}+I_{12}=I_{32}.}$$
13/24通道仍保留分类意义；上述等式来自变量换元和排序分区。

## 2. 12的同边积分与两个结构函数

令λ=1-2e、B_e=B(1-e,1-e)，并令F12(e,w)表示12中b,x,y的积分，即w<x<y<1。

原核的联合原函数为
$$
\frac{[bx+(1-b)y]^{2e}}{2e(2e-1)b(1-b)}.
$$
代回三角边界，再用b上的Beta分部积分，得到
$$
F_{12}(e,w)=\frac{\mathcal N(e,w)}{2\lambda e^2},
$$
$$
\mathcal N(e,w)=B_e\left[
 w^{2e}-1+
 2e(1-w)\,{}_2F_1(1-2e,1-e;2-2e;1-w)\right].
$$
Beta分项可先在固定0<w<1、Re(e)<0处理，再延拓回完整差分；不能把分项Beta的发散当成原F12的发散。e=0附近括号的常数项和一次项精确相消，故N从e²起，而F12本身有限。

定义两个有明确积分来源的结构函数：
$$
\mathcal T(e,z)=(1+z)\int_0^1ds\,dw\,D_z(s,w)^{-1+e}(1-w^{2e}),
$$
$$
\mathcal K(e,z)=(1+z)\int_0^1ds\,dw\,D_z(s,w)^{-1+e}\mathcal N(e,w).
$$
T是带端点权重的既有单对边母积分差；K来自TGV同边排序的精确积分。这些定义用于说明结构，下一节给出无剩余积分的完整解析表示。
$$
J_{32}=-\frac{z^{2e}U(e)}{16\lambda e^4}\mathcal T(e,z),\qquad
J_{12}=-\frac{z^{2e}}{16\lambda e^4}\mathcal K(e,z).
$$

设ℓ=log z，并定义
$$
G=e^{-3\gamma_Ee}\Gamma(1-e)\Gamma(1-2e),\qquad
\mathcal E_X(e,z)=\frac{z^{e/2}X(e,z)+z^{-e/2}X(e,1/z)}2.
$$
prepared C中的Euler补偿仅出现一次，通用规范为
$$I=\frac{\lambda G}{8}\frac{J(z,e)+z^{3e}J(1/z,e)}{1+z^{3e}}.$$
代入上面的J得到
$$
\boxed{
I_{32}=-\frac{GU\,\mathcal E_{\mathcal T}}{128e^4\cosh(3e\ell/2)},\qquad
I_{12}=-\frac{G\,\mathcal E_{\mathcal K}}{128e^4\cosh(3e\ell/2)},\qquad
I_{11}=I_{32}-I_{12}.
}
$$
实际最高极点为32三阶、12二阶、11三阶，不能从显式1/e⁴直接读出。

## 3. 完整ε：已逐项积完的收敛级数

令Λ=1+z、c_n=(1-e)_n/n!，并将以下有限多项式写成
$$
S_n(w)=\frac1{n+1}\sum_{k=0}^n
 \left(\frac{z+w}{\Lambda}\right)^{n-k}
 \left(1-\frac{zw}{\Lambda}\right)^k
=\sum_{j=0}^n a_{n,j}(z)w^j.
$$
它就是s积分后的[1-D/Λ]^n；因此a由一个显式有限多项式确定，无待求积分。

定义
$$
R_j(e)=\frac1{j+1}-\frac1{j+1+2e},
$$
$$
N_j(e)=B_e\left[
-R_j(e)+\frac{2e}{(j+1)(j+2)}
{}_3F_2(1-2e,1-e,2;2-2e,j+3;1)\right].
$$
则
$$
\boxed{
\mathcal T(e,z)=\Lambda^e\sum_{n=0}^{\infty}c_n
 \sum_{j=0}^{n}a_{n,j}(z)R_j(e),\qquad
\mathcal K(e,z)=\Lambda^e\sum_{n=0}^{\infty}c_n
 \sum_{j=0}^{n}a_{n,j}(z)N_j(e).
}
$$
R_j和N_j分别是w^j(1-w^(2e))及w^j N(e,w)的已计算矩。所有n项已积完。内层j和有限，必须先求内层和，不任意交换成独立双重无穷和。

在实1/2<e<1、z>0，D/Λ∈(0,1)，binomial系数非负。b↔1-b配对说明F12≥0，且F12(e,w)≤F12(e,0)=h0/(2e)，所以|N(e,w)|≤|U(e)|且N符号固定；由Tonelli得到以上级数，再作ε亚纯延拓。对截断到n=N，有一个显式、虽不紧的尾界：
$$
|\mathcal T-\mathcal T_N|
\le \Lambda^e\,\frac{(1-e)_{N+1}}{(N+1)!}
\frac{\ell^2+\pi^2}{2},
$$
$$
|\mathcal K-\mathcal K_N|
\le |U(e)|\Lambda^e\,\frac{(1-e)_{N+1}}{(N+1)!}
\frac{\ell^2+\pi^2}{2}.
$$
利用c_n递减、几何尾和及
$(1+z)\int D^{-1}ds\,dw=(\ell^2+\pi^2)/2$。界随N趋于零；没有声称级数高效或是有限Gamma闭式。e=1/2等中间可去点取完整组合的极限。

## 4. 有限展开的实际积分与既有基底

s先积掉，给出
$$
P(e,w,z)=\frac{(zw)^e-(1-w)^e}{e[w-1/(1+z)]}.
$$
内点w=1/(1+z)的奇点可去，其值是(1+z)(1-w)^(e-1)。

令
$$
L_b=-\log b-\log(1-b)+2\log(1-b+bw),\qquad
H_j(w)=\int_0^1db\,\frac{1-w}{1-b+bw}\frac{L_b^j}{j!}.
$$
这四个临时矩H0…H3只用于计算展开，H0=−log w，不等于2ne的同名单对边母函数。新鲜HyperInt先积b，再积w，所有发散检查开启。b消元后字母只有0、1；w阶段再加1/(1+z)，最终fibration到z的字母为0、−1。这给出该计算的线性可约次序。

实际求出T的e¹…e⁴、K的e¹…e⁴（K¹=0），再做上述E投影和Gamma展开。四个T矩另用此前2ne的H0/H2源矩独立精确重构。21个Laurent位置全部回写为现有权重≤6基底，并保留完整Hlog坐标。

11图的ε⁻³、ε⁻²、ε⁻¹、ε⁰系数分别纯权重3、4、5、6；12和32含较低权重项，它们在I32−I12中相消。有限项用到已有N、τDq4与E，但本轮没有新增基底生成元，也不声称任意高阶闭合或这些表示全局最短。当前固定Q-Hlog坐标中，M/N的权重6空间秩为12，加入本批有限项后为13；11/12需要的额外方向仅为既有组合τDq4+10E（分别系数−1/192和+1/192），32无需该方向。

展开可在固定z>0、足够小的|e|<δ中逐项进行。例如取δ<1/8，对边角点使s积分至多产生w^(-δ)，b角点至多再产生w^(-3δ)，并带可积对数；总上界w^(-4δ)可积。另一端点同样受控。N的前两阶精确相消，T从一次起；保留至e⁴，代入显式e⁻⁴后余项为O(e)。Gamma、U及1/cosh在零点邻域解析。

[最终解表](solutions.md)给出全部21项；[核验记录](../../../../workflow/checks/tgv-t1-opposite/README.md)区分新算矩、旧结果复验、高精度降维求积及原七维QMC。没有新增Kernel/Lean收据，物理复z/i0延拓不在本轮结论内。
