# 4e-04直接导出M/N：未展开恒等式与任意阶约化

日期：2026-09-27。本轮回答“能否从4e-04直接推出M/N，而不是先求Hlog再拟合”。**可以：下面给出未展开ε的积分恒等式。它直接产生原N的对数权重，并给出4e-04任意阶的M/N生成规则。**

本轮不以此前有限阶偶源拟合作为前提。推导使用有序Beta全ε式、同一矩源的Cauchy变换方程及解析边界常数；Hlog只在最后对已有五阶做回归，不参与寻找系数。分析推导及辅助检查完成，formal=not_run。

## 1. 直接结果：N的权重在积分层出现

沿用原定义，固定z>0：

$$
\begin{gathered}
r=\frac{z}{1+z},\quad \rho=r(1-r),\quad \ell=\log z,\\
W(x)=\log[x(1-x)],\qquad
R(x)=\log\frac{1-x}{x},\qquad L(x;z)=\ell+R(x),\\
M_{k,j}=\int_0^1\frac{W^kL^j}{x-r}\,dx,\qquad
N_{k,j}=\int_0^1\frac{W^kRL^j}{x-r}\,dx \quad(j\ge1).
\end{gathered}
\tag{R1}
$$

L(r;z)=0，故这些积分的内部点连续可去。下面的N生成函数没有新基底定义：

$$
\boxed{\displaystyle
\mathcal N(\epsilon,z)=
\int_0^1\frac{e^{\epsilon W(x)}R(x)
[\cosh(\epsilon L(x;z))-1]}{x-r}\,dx .}
\tag{R2}
$$

它从ε²开始；x=r处被积函数趋于0，不取主值。令H为原一圈正规化母核，χ(ε)=γ_E+ψ(1+2ε)，则本轮得到

$$
\boxed{\begin{aligned}
Q(\epsilon,z)
&=\frac12H(\epsilon,z)^2
-\frac{2}{\epsilon^2}\mathcal N(\epsilon,z)
+\frac{2\chi(\epsilon)}{\epsilon}H(2\epsilon,z),\\
\Delta(\epsilon,z)
&=\frac12H(\epsilon,z)^2
-\frac{4}{\epsilon^2}\mathcal N(\epsilon,z)
+\frac{4\chi(\epsilon)}{\epsilon}H(2\epsilon,z).
\end{aligned}}
\tag{R3}
$$

Q就是4e-04的裸双核，Δ=Q−S。H(2ε,z)是同一一圈母核的调节参数替换，不是另一个二圈图。完整物理I仍为

$$
I_{4e-04}=-\frac{e^{-2\gamma_E\epsilon}\Gamma(1-\epsilon)^2}
{8\cosh(\epsilon\ell)}\,Q.
\tag{R4}
$$

χ是已求值的Gamma边界项，不是待定常数。χ/ε在0处趋于2ζ₂，𝓝/ε²也正则。因此(R3)没有真实的ε负幂。

**N出现的原因现在明确了：转换过程中产生的两项是L插入与外部ℓ乘同一偶幂源，相减后恰好留下R=L−ℓ。** 这个组合在未展开积分上成立，先于任何Hlog或有限阶坐标比较。

## 2. 任意阶的M/N系数公式

令h_n=[εⁿ]H。原一圈生成式给出

$$
h_n=\sum_{j=0}^{\lfloor n/2\rfloor}
\frac{M_{n-2j,\,2j+1}}{2^n(n-2j)!(2j+1)!}.
\tag{R5}
$$

直接展开(R2)而非拟合，有

$$
\frac{\mathcal N(\epsilon,z)}{\epsilon^2}
=\sum_{q\ge0}\epsilon^q
\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{N_{q-2j,\,2j+2}}{(q-2j)!(2j+2)!}.
\tag{R6}
$$

于是

$$
\boxed{\begin{aligned}
Q_q={}&\frac12\sum_{i=0}^{q}h_i h_{q-i}
-2\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{N_{q-2j,\,2j+2}}{(q-2j)!(2j+2)!}\\
&+2^{q+2}\sum_{k=0}^{q}(-1)^k\zeta_{k+2}h_{q-k},
\qquad q\ge0 .
\end{aligned}}
\tag{R7}
$$

这不是只验到q=4的猜想：它由(R3)的解析生成函数逐项展开得到。它也把此前行列式稿(D11)的偶源商关系提升为任意阶解析关系；历史稿仍保留其当时的有限阶状态。局部一致可积性可取|ε|<1/4；端点的最差幂为x^{2 Re ε}或(1−x)^{2 Re ε}，有限次对数插入仍可积。

物理I再按(R4)卷积。外部sech只产生偶次ℓ，可用ℓ²=−2M₀,₁−π²消去。因此**4e-04任意阶系数属于M/N多项式代数，标量系数只需有理数和ζ常数**。这项结论只针对本图及其裸交换配对，不推广为全部三圈图或全部双源f_mn的闭包。

例如(R7)在一次幂立即给出

$$
\boxed{\displaystyle
[\,\epsilon\,]I_{4e-04}
=-\frac{M_{0,1}M_{1,1}}{16}
+\zeta_3M_{0,1}-\frac{\pi^2M_{1,1}}{12}
+\frac{N_{1,2}}8 .}
\tag{R8}
$$

这与旧简洁式相同，但本轮的系数由(R3)决定，没有进行候选拟合。

## 3. 更统一的对象是同一个矩源

仍使用一圈的双源

$$
\mathscr U(a,b;z)=
\int_0^1\frac{e^{aW}(e^{bL}-1)}{x-r}\,dx .
\tag{R9}
$$

M由参数导数产生；N是同一源上的另一个算子：

$$
\begin{aligned}
M_{k,j}&=\left.\partial_a^k\partial_b^j\mathscr U\right|_{a=b=0},\\
N_{k,j}&=\left.\partial_a^k\partial_b^j
(\partial_b-\ell)\mathscr U\right|_{a=b=0}.
\end{aligned}
\tag{R10}
$$

这是R=L−ℓ的积分恒等式。对4e-04，所需的是b=±ε的偶组合，并减去b=0的值。具体记

$$
\begin{aligned}
Y&=\partial_b\mathscr U(\epsilon,0)
=\int_0^1\frac{e^{\epsilon W}L}{x-r}\,dx,\\
Y_c&=\frac12[\partial_b\mathscr U(\epsilon,\epsilon)
+\partial_b\mathscr U(\epsilon,-\epsilon)]
=\int_0^1\frac{e^{\epsilon W}L\cosh(\epsilon L)}{x-r}\,dx,\\
E&=\mathscr U(\epsilon,\epsilon)+\mathscr U(\epsilon,-\epsilon).
\end{aligned}
$$

则

$$
\boxed{\mathcal N=Y_c-Y-\frac{\ell E}{2}.}
\tag{R11}
$$

因此本图不必另造一个比M/N更大的最终表示族。更统一的是𝓤及其源导数，而不是宣称M和N分别属于不同圈数。最小独立性仍是另一项问题。

## 4. 连接有序Beta函数与矩源的两条解析恒等式

复用[有序Beta推导](ordered-beta-generator.md)的J、T、Ψ，其中
J_a(x)=−ax ₃F₂(1,1,1−a;2,2;x)，
T_a=𝔉(a,0;x)，
Ψ(a,b;z)=𝔉(a,b;−z)+𝔉(a,b;−1/z)−2𝔉(a,b;1)。

只在中间引入一个已减除的零指标矩

$$
V(\epsilon,z)=\int_0^1
\frac{e^{\epsilon W(x)}-\rho^\epsilon}{x-r}\,dx .
$$

V的内部点同样可去；它最终会消去，不是最终增加的基底。令d=J_{2ε}(1)=−χ(ε)，则两条关键转换是

$$
\boxed{\rho^\epsilon\Psi(\epsilon,\epsilon;z)
=Y-\ell V-\rho^\epsilon M_{0,1}.}
\tag{R12}
$$

$$
\boxed{\begin{aligned}
\rho^\epsilon\big[&
\Psi(2\epsilon,0;z)
-\ell(J_{2\epsilon}(-z)-J_{2\epsilon}(-1/z))
+2M_{0,1}\big]\\
&=2Y_c+2\epsilon d\,H(2\epsilon,z).
\end{aligned}}
\tag{R13}
$$

另对两个单边幂作线性换元，得到

$$
E=\rho^\epsilon[J_{2\epsilon}(-1/z)-J_{2\epsilon}(-z)]-2V.
\tag{R14}
$$

将(R12)–(R14)代入此前
Q=H²/2+ρ^ε[2Ψ(ε,ε)−Ψ(2ε,0)]/ε²，所有V、显式边界J差和M₀,₁项按原系数消去，留下
Q=H²/2+[2(Y−Y_c)+ℓE−2εdH(2ε)]/ε²。
再用(R11)，即得(R3)。

下面给出(R12)/(R13)的证明和常数固定，避免把它们当作数值猜测。

## 5. 矩源方程及反射

定义

$$
C(a,b)=\frac{\Gamma(1+a-b)\Gamma(1+a+b)}{\Gamma(1+2a)},
\qquad \theta=z\partial_z.
$$

对矩源作一次分部积分，得到

$$
\boxed{\theta(\rho^{-a}\mathscr U(a,b))
=\rho^{-a}[C(a,0)-C(a,b)z^b].}
\tag{R15}
$$

推导可通过辅助Cauchy变换
𝓒_{s,t}(r)=PV∫₀¹x^s(1−x)^t/(x−r)dx组织。先令r在积分线外，积分
∂_x[x(1−x)x^s(1−x)^t/(x−r)]，得到

$$
r(1-r)\mathcal C'_{s,t}
=[s-(s+t)r]\mathcal C_{s,t}
-\frac{\Gamma(1+s)\Gamma(1+t)}{\Gamma(1+s+t)}.
$$

再取上下边界平均。𝓤的两个Cauchy留数相同，差中完全抵消，所以最终(R15)是普通可去积分的恒等式；没有给原图另加主值处方。此路线先在端点绝对收敛参数域成立，再延拓到小参数邻域。

记A=C(ε,0)ρ^{-ε}、D=ρ^{-ε}V、B=ρ^{-ε}Y。由(R15)及相应零指标减除式，

$$
\theta D=1-A,\qquad \theta B=-\ell A.
$$

另一方面，由有序Beta源方程及Gauss反射，

$$
G_\epsilon(z)+G_\epsilon(1/z)-1=C(\epsilon,0)\rho^{-\epsilon},
\quad
G_\epsilon(z)={}_2F_1(1,-\epsilon;1+\epsilon;-z),
$$

因此θ²Ψ(ε,ε)=A−1。D在z→1/z下为奇，Ψ为偶，故θΨ=−D。继续积分可知
Ψ−ρ^{-ε}(Y−ℓV)−ℓ²/2仅差z无关常数。附录证明该常数为π²/2，给出(R12)。

对(R15)分别在b=ε、−ε求b导数，使用C(ε,±ε)=1、
∂_bC(ε,ε)=−d、∂_bC(ε,−ε)=d，即得

$$
\theta\{\rho^{-\epsilon}[2Y_c+2\epsilon dH(2\epsilon)]\}
=-\ell\rho^{-\epsilon}(z^\epsilon+z^{-\epsilon}).
$$

这也等于θ{Ψ(2ε,0)−ℓ[J_{2ε}(−z)−J_{2ε}(−1/z)]+2M₀,₁}。附录固定其常数差为0，得到(R13)。

## 6. 证据与适用范围

- 上述转换先取0<ε<1/4推导边界渐近，再由解析延拓覆盖z>0、实0<ε<1及ε=0的足够小复邻域。
- [辅助检查](../../../../workflow/checks/4e04-direct-mn/README.md)核对Cauchy分部积分、两条边界常数、N源展开和全部重新组合。遗漏digamma边界的负控制被检出。
- 13个ε/z点直接积分普通一维Y、Y_c、V、𝓝，与₄F₃比较，包括负ε与复ε样点；不使用数值PV积分，也不使用Hlog数值来寻找公式。
- 最后才按(R5)–(R7)生成五阶M/N表达，回代旧Hlog：Q、Δ、完整I共15项精确相同。没有rref、PSLQ或候选拟合。

本轮给出4e-04的全阶M/N解析约化，不证明任意双源f_mn都可约，也不推广到带端点插入、6e或3g。M/N的最小独立性仍未证明。旧结果、旧阶段“尚未证明”的历史报告及catalog保持原字节；formal=not_run，原329位置覆盖不重复增加。

## 附录：边界常数如何解析固定

以下只用Gamma/Beta渐近，不读取任何Laurent答案。令u=log z→−∞，先取0<ε<1/4，保证两次使用的幂a=ε、2ε都在(0,1)。

### A. 两个统一的渐近公式

由Euler积分及[Gauss连接公式](https://dlmf.nist.gov/15.8#E2)，辅助Cauchy变换有精确分解

$$
\mathcal C_{s,t}(r)
=\mathrm B(s,t+1)\,
{}_2F_1(1,-s-t;1-s;r)
-\pi\cot(\pi s)\,r^s(1-r)^t .
\tag{R16}
$$

先取非整数0<Re s<1，再按消去可去极点的方式延拓。故
𝓒_{s,t}/[r^s(1−r)^t]
=B(s,t+1)r^{-s}(1−r)^{-t}−πcot(πs)+o(1)；
在s→0时必须把Beta项和cot项合起来展开，不能丢掉其中一项。

对𝔉(a,b;−y)，Gauss函数的Mellin核为

$$
\int_0^\infty y^{\sigma-1}
{}_2F_1(1,-a;1+b;-y)\,dy
=\frac{\Gamma(\sigma)\Gamma(1-\sigma)
\Gamma(-a-\sigma)\Gamma(1+b)}
{\Gamma(-a)\Gamma(1+b-\sigma)}.
\tag{R17}
$$

先在Re a<0、Re b>0及0<Re σ<min(1,−Re a)由两次Euler/Beta积分推出，再作参数延拓。由于𝔉是Gauss函数减1的两次Euler原函数，σ=0处的三重极点固定log²、log和常数；延拓到0<a<1时，另保留移动幂极点的贡献。得到

$$
\begin{gathered}
\mathfrak F(a,b;-y)
=\frac{\Gamma(1+a)\Gamma(1+b)}{\Gamma(1+a+b)}
\frac{y^a}{a^2}
-\frac12\log^2y+c_{ab}\log y+d_{ab}+o(1),\\
c_{ab}=\psi(1+b)-\psi(-a),\\
d_{ab}=-\frac12\left[c_{ab}^2+\psi_1(-a)-\psi_1(1+b)+2\zeta_2\right].
\end{gathered}
\tag{R18}
$$

这些系数也可直接由(R17)在σ=0的Gamma展开读取。主分支始终取−y的负实超几何自变量；没有穿越其[1,∞)割线。

### B. 第一条转换的常数是π²/2

记P=πcot(πε)、E₂=π²csc²(πε)、A=C(ε,0)ρ^{-ε}。由(R16)及b=0处的参数导数，

$$
D=\frac A\epsilon+u-P+o(1),\qquad
B=A\left(\frac u\epsilon+\frac1{\epsilon^2}\right)-E_2+o(1).
$$

由(R18)、𝔉(ε,ε;1)=[ψ₁(1+ε)−ζ₂]/2及[polygamma反射](https://dlmf.nist.gov/5.15)，

$$
\Psi(\epsilon,\epsilon)
=\frac A{\epsilon^2}-\frac{u^2}{2}+Pu-E_2+\frac{\pi^2}{2}+o(1).
$$

所以Ψ−(B−uD)−u²/2→π²/2，固定(R12)的常数。不是通过某个有限阶拟合补常数。

### C. 第二条转换的常数是0

记a=2ε、d=J_a(1)、P=πcot(πa)、E₂=π²csc²(πa)、T₁=T_a(1)，并设

$$
d_T=-\frac12[(d-P)^2+E_2-\psi_1(1+a)+\zeta_2].
$$

由(R18)，(R13)左侧除以ρ^ε后的渐近为

$$
(1+1/z)^a\left(\frac u a+\frac1{a^2}\right)
-\frac{u^2}{2}+d_T-2T_1-\pi^2+o(1).
$$

由(R16)分别对b=ε、−ε求导，右侧除以ρ^ε后的渐近具有相同幂项和−u²/2，常数为

$$
\frac{d^2-\psi_1(1+a)-5\zeta_2}{2}-E_2+dP.
$$

在b=ε一端，先对(R16)的完整Gauss项消去s=0的可去极点并作源导数，再取r→0；高阶r项给出的余项仍趋零。保留到O(1)时，所需合并极点的局部式为
B(s,1+a-s)=1/s+d+s[d²+ζ₂+ψ₁(1+a)]/2+O(s²)，
再乘e^{-su}并加入−πcot(πs)。这说明上述常数中的5ζ₂不能省略。
用T₁=[ψ₁(1+a)−ζ₂−d²]/2及E₂−P²=π²，两常数精确相同。因此(R13)的常数差为0。所有渐近余项在此参数范围及所用参数导数下趋零；随后以原可去积分的解析性延拓。
