# 4e-04的通解推导：从原积分到M/N生成函数

本文按“原积分 → 降维 → 整理剩余排序 → 一维生成函数 → Laurent展开”的顺序推导4e-04。组织方式参考[GF cal中的2e-02推导](<../../../../GF cal 推导草稿/2e-02.md>)，归一化采用项目当前约定。

与2e-02相比，多出来的工作是处理两条对边传播子的相对排序。最后它会变成三个明确部分：**一圈生成函数的平方、N的生成函数、端点产生的一圈边界项**。正文逐步说明它们怎样出现；涉及Cauchy变换与边界常数的证明放在文末附录。

## 1. 起点：原参数积分与物理归一化

### 1.1 原四重有序积分

令z>0，先取实0<ε<1。原参数对应关系为
t₁=s₂₁、t₂=s₂₂、s₁=s₄₁、s₂=s₄₂，原外到内顺序为t₂、s₂、t₁、s₁。记

$$
D_z(t,s)=(1-t)(1-s)+zts,\qquad
\Omega=\{0<t_1<t_2<1,\ 0<s_1<s_2<1\}.
$$

当前输入的参数部分是

$$
\boxed{
J_{4e-04}(\epsilon,z)
=(1+z)^2\int_\Omega
D_z(t_1,s_1)^{\epsilon-1}D_z(t_2,s_2)^{\epsilon-1}
\,dt_1dt_2ds_1ds_2 .}
\tag{J}
$$

两个低端点相配，两个高端点相配。原输入、颜色和顺序见[逐图根推导](../../workspaces/4e-04/derivation.md)。以下换序是下游操作：在0<ε<1时，把Ω放大成两个独立核的正方形乘积，即得到可积上界，允许所用的Fubini换序。

### 1.2 去掉简单尺度幂

统一使用

$$
\alpha=\frac1{1+z},\qquad
\bar\alpha=\frac{z}{1+z},\qquad
\rho=\alpha\bar\alpha=\frac{z}{(1+z)^2},\qquad
\ell=\log z.
$$

定义正规化母核Q=z^{-ε}J₄ₑ₋₀₄。两条边同时反射并反转排序，利用
D_z(1−t,1−s)=zD_{1/z}(t,s)，可得
J(1/z)=z^{-2ε}J(z)，所以Q(1/z)=Q(z)。

本图sym=2、a₁=−1/2，因此κ=sym·a₁/[2(−2)²]=−1/8。记

$$
Z(\epsilon)=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon).
$$

当前唯一I约定给出

$$
\begin{aligned}
I_{4e-04}
&=-\frac{Z^2}{8}\,
\frac{J(z)+z^{2\epsilon}J(1/z)}{1+z^{2\epsilon}}\\
&=\boxed{-\frac{Z^2}{8\cosh(\epsilon\ell)}\,Q(\epsilon,z).}
\end{aligned}
\tag{I}
$$

**从现在起只需求Q，最后再乘这一前因子。** Z已精确包含Euler补偿，不另行“丢弃γ_E项”或再乘一次补偿。

## 2. 先回顾一圈：M的生成函数怎样出现

这一步既与2e-02对照，也定义两圈通解中将出现的一圈子结构。

### 2.1 积去一个线性变量

定义不含Gamma前因子的一圈母核

$$
\mathcal V_M(\epsilon,z)
=-(1+z)z^{-\epsilon/2}
\int_0^1dt\int_0^1ds\,D_z(t,s)^{\epsilon-1}.
$$

对s用基本原函数

$$
\int_a^b(A+Bx)^{\epsilon-1}dx
=\frac{(A+Bb)^\epsilon-(A+Ba)^\epsilon}{\epsilon B}
$$

积分，再令t=1−β，得到

$$
\mathcal V_M
=\frac1\epsilon\int_0^1
\frac{z^{-\epsilon/2}
[(z(1-\beta))^\epsilon-\beta^\epsilon]}
{\beta-\bar\alpha}\,d\beta .
$$

### 2.2 用双曲函数整理幂差

定义

$$
W(\beta)=\log[\beta(1-\beta)],\qquad
R(\beta)=\log\frac{1-\beta}{\beta},\qquad
L(\beta;z)=\ell+R(\beta).
$$

分子恰好是

$$
z^{-\epsilon/2}
[(z(1-\beta))^\epsilon-\beta^\epsilon]
=2e^{\epsilon W/2}\sinh(\epsilon L/2).
$$

因此

$$
\boxed{\displaystyle
\mathcal V_M(\epsilon,z)=
\int_0^1\frac{d\beta}{\beta-\bar\alpha}\,
e^{\epsilon W/2}
\frac{\sinh(\epsilon L/2)}{\epsilon/2}.}
\tag{VM}
$$

展开指数与sinh，就自然产生

$$
M_{k,j}(z)=\int_0^1\frac{W^kL^j}{\beta-\bar\alpha}\,d\beta,
\qquad k\ge0,\ j\ge1.
$$

这就是GF cal的2e-02文档中作为积分部分的生成函数。本文记作𝓥_M，以免与物理前因子混淆；此前研究稿的裸单核H也等于它。

## 3. 两圈第一次降维：4D → 2D

令u=t₂、v=s₁。先积t₁，再积s₂；两步都使用§2.1的线性幂原函数：

$$
\begin{aligned}
\int_0^uD_z(t_1,v)^{\epsilon-1}dt_1
&=\frac{D_z(u,v)^\epsilon-(1-v)^\epsilon}
{\epsilon[(1+z)v-1]},\\
\int_v^1D_z(u,s_2)^{\epsilon-1}ds_2
&=\frac{(zu)^\epsilon-D_z(u,v)^\epsilon}
{\epsilon[(1+z)u-1]}.
\end{aligned}
$$

记D=D_z(u,v)，则

$$
\boxed{\displaystyle
Q=\frac{(1+z)^2z^{-\epsilon}}{\epsilon^2}
\int_0^1du\int_0^1dv\,
\frac{[D^\epsilon-(1-v)^\epsilon][(zu)^\epsilon-D^\epsilon]}
{[(1+z)v-1][(1+z)u-1]}.}
\tag{Q2D}
$$

**与一圈的差别在这里：现在得到两个幂差的乘积，它们共同依赖D(u,v)。** 不能把它们当成两个独立一圈积分。

分母为零的点由原函数连续延拓。例如第一个差商在(1+z)v−1=0时趋于u(1−v)^{ε−1}。必须先保留差的整体，不能把四个展开项各自当成收敛积分。

## 4. 把剩余排序求值：为什么出现有序Beta函数

### 4.1 中心换元与可积分项的整理

令p=(1+z)u−1、q=(1+z)v−1，则p,q都在(−1,z)，并且
D=(z+pq)/(1+z)。暂记

$$
A(p)=(1+p)^\epsilon,\quad B(p)=(1-p/z)^\epsilon,\quad
C(p,q)=(1+pq/z)^\epsilon .
$$

式(Q2D)变成

$$
Q=\frac{\rho^\epsilon}{\epsilon^2}
\int_{-1}^z\frac{dp}{p}\int_{-1}^z\frac{dq}{q}\,
[C-B(q)][A(p)-C].
$$

利用p、q交换下的反对称项积分为零，可将分子替换为

$$
[A(p)+B(p)](C-1)-(C^2-1)-[A(p)-1][B(q)-1].
$$

两分子之差为$[B(p)-B(q)](C-1)$，确实反对称。替换后的各项在p=0或q=0处分别可去，所以现在可以拆项积分。

定义端点原函数

$$
\mathcal J_a(x)=\int_0^x\frac{(1-t)^a-1}{t}\,dt,\qquad
\mathcal T_a(x)=\int_0^x\frac{\mathcal J_a(t)}t\,dt .
$$

积掉q所用的关系只是

$$
\int_{-1}^z\frac{(1+pq/z)^a-1}{q}\,dq
=\mathcal J_a(-p)-\mathcal J_a(p/z).
$$

因此四重积分已变成一维原函数的乘积。把p的正、负半区分别缩放到(0,1)，剩余排序统一为下一小节的K。

### 4.2 对排序积分求z导数，再作一次分部积分

两种边界幂分别记为A_a(t;z)=(1+zt)^a−1、B_b(t)=(1−t)^b−1。定义

$$
K(a,b;z)=\int_{0<u<v<1}
\frac{A_a(u;z)B_b(v)-B_b(u)A_a(v;z)}{uv}\,du\,dv .
$$

令f(v)=𝓙_a(−zv)、g(v)=𝓙_b(v)，则K=∫₀¹(fg'−gf')dv。
利用θ_zf=A_a和θ_zA_a=v∂_vA_a，分部积分得到

$$
\begin{aligned}
\theta_zK
&=2P(a,b;z)-\mathcal J_b(1)[(1+z)^a-1],\\
P(a,b;z)&=\int_0^1
\frac{[(1+zv)^a-1][(1-v)^b-1]}v\,dv .
\end{aligned}
$$

嵌套原函数消去了，v=1的边界项仍保留。随后展开第一项幂，并使用普通Beta积分：

$$
P=\sum_{n\ge1}\frac{(-a)_n(-z)^n}{n!}
\left[\mathrm B(n,1+b)-\frac1n\right].
$$

由此自然得到

$$
\boxed{\displaystyle
\mathfrak F(a,b;x)
=\sum_{n\ge1}\frac{(-a)_n}{(1+b)_n}\frac{x^n}{n^2}
=-\frac{ax}{1+b}\,
{}_4F_3\!\left(\begin{matrix}1,1,1,1-a\\2,2,2+b\end{matrix};x\right).}
\tag{F}
$$

其中T_a=𝔉(a,0)。再沿z积分，用K(a,b;0)=0固定常数：

$$
K(a,b;z)=2\mathfrak F(a,b;-z)
-2\mathfrak F(a,0;-z)-\mathcal J_b(1)\mathcal J_a(-z).
$$

所以排序积分已经实际求值，不再保留未知的二维函数。上面的级数先在|z|<1使用；标准₄F₃沿负实自变量延拓，覆盖z>0。

### 4.3 所有边界项重新合并

记同一函数的crossing组合

$$
\Psi(a,b;z)=\mathfrak F(a,b;-z)+\mathfrak F(a,b;-1/z)
-2\mathfrak F(a,b;1).
$$

为把边界项也明确列出，暂记j₋=𝓙_ε(−z)、j₀=𝓙_ε(1)、j₊=𝓙_ε(−1/z)，以及

$$
\mathcal C_a=\mathcal T_a(-z)+\mathcal T_a(-1/z)-2\mathcal T_a(1).
$$

§4.1的一维原函数积分经一次分部积分，给出

$$
Q=\frac{\rho^\epsilon}{\epsilon^2}\left[
K(\epsilon,\epsilon;z)+K(\epsilon,\epsilon;1/z)
+\frac{(j_-+j_+)^2}{2}-j_0(j_-+j_+)
+2\mathcal C_\epsilon-\mathcal C_{2\epsilon}\right].
$$

代入§4.2的K求值式，并用
𝓥_M=−ρ^{ε/2}(j₋+j₊−2j₀)/ε。
两个权重相同时交换子为零，另给出端点关系
𝔉(ε,ε;1)=𝓣_ε(1)+j₀²/2。这样T_ε的边界项和j₀的平方均按原系数合并，得到

$$
\boxed{\displaystyle
Q=\frac12\mathcal V_M(\epsilon,z)^2
+\frac{\rho^\epsilon}{\epsilon^2}
[2\Psi(\epsilon,\epsilon;z)-\Psi(2\epsilon,0;z)].}
\tag{QF}
$$

这里的一圈平方有明确来源：两种配对的有序区域拼起来给
Q+S=𝓥_M²/2；S为交换配对的裸母核。该关系属于母积分层，不直接套用于带不同颜色和前因子的物理I。

式(QF)已经是有限标准函数通解。下一步是把它整理成与一圈M生成函数相同的矩语言。

## 5. 把Beta表达变成同一个一维矩核

### 5.1 使用一圈已有的双源

仍使用§2的W、R、L和分母$\beta-\bar\alpha$，定义

$$
\mathscr U(a,b;z)=\int_0^1
\frac{e^{aW}(e^{bL}-1)}{\beta-\bar\alpha}\,d\beta .
$$

对它求源参数导数，产生L插入。此处需要三个组合：

$$
\begin{aligned}
Y&=\int_0^1\frac{e^{\epsilon W}L}{\beta-\bar\alpha}\,d\beta,\\
Y_c&=\int_0^1\frac{e^{\epsilon W}L\cosh(\epsilon L)}
{\beta-\bar\alpha}\,d\beta,\\
E&=2\int_0^1\frac{e^{\epsilon W}[\cosh(\epsilon L)-1]}
{\beta-\bar\alpha}\,d\beta .
\end{aligned}
$$

Y是b=0处的源导数，Y_c是b=±ε处源导数的平均，E是b=±ε处源本身的和。

### 5.2 关键约化：保留边界后，整个组合怎样变化

记d=𝓙_{2ε}(1)。有如下未展开恒等式：

$$
\boxed{\displaystyle
\rho^\epsilon[2\Psi(\epsilon,\epsilon)-\Psi(2\epsilon,0)]
=2Y-2Y_c+\ell E-2\epsilon d\,\mathcal V_M(2\epsilon,z).}
\tag{BR}
$$

这是从有序Beta表达到M/N表达的关键步骤。它由矩源的一阶Cauchy方程及解析边界常数得到，不是系数拟合；完整证明在附录A、B。

具体消去过程也可以逐项看清。暂记
V=∫₀¹[e^{εW}−ρ^ε]dβ/($\beta-\bar\alpha$)，并记J_-=𝓙_{2ε}(−z)、J_+=𝓙_{2ε}(−1/z)。三条转换为

$$
\begin{aligned}
\rho^\epsilon\Psi(\epsilon,\epsilon)&=Y-\ell V-\rho^\epsilon M_{0,1},\\
\rho^\epsilon[\Psi(2\epsilon,0)-\ell(J_--J_+)+2M_{0,1}]
&=2Y_c+2\epsilon d\,\mathcal V_M(2\epsilon,z),\\
E&=\rho^\epsilon(J_+-J_-)-2V.
\end{aligned}
$$

把前两式代入(QF)，再用第三式，V、J差及M₀,₁边界项全部消去，正好留下(BR)。V只是中间减除量，不加入最终基底。

## 6. 双曲分解后，N自然出现

现在把Y、Y_c、E的积分放到同一分母上：

$$
\begin{aligned}
2Y-2Y_c+\ell E
&=-2\int_0^1
\frac{e^{\epsilon W}(L-\ell)[\cosh(\epsilon L)-1]}
{\beta-\bar\alpha}\,d\beta\\
&=-2\int_0^1
\frac{e^{\epsilon W}R[\cosh(\epsilon L)-1]}
{\beta-\bar\alpha}\,d\beta .
\end{aligned}
$$

**额外的R不是猜出来的：它来自L插入减去外部的ℓ插入。** 于是定义正则的N生成函数

$$
\boxed{\displaystyle
\mathcal V_N(\epsilon,z)=
\int_0^1\frac{d\beta}{\beta-\bar\alpha}\,
e^{\epsilon W}R\,
\frac{\cosh(\epsilon L)-1}{\epsilon^2}.}
\tag{N}
$$

展开它，就自然产生原来的

$$
N_{k,j}(z)=\int_0^1
\frac{W^kRL^j}{\beta-\bar\alpha}\,d\beta
=M_{k,j+1}-\ell M_{k,j},\qquad k\ge0,\ j\ge1.
$$

两类生成函数现在可以直接对照：

| 来源 | 指数生成的权重 | 双曲函数生成的权重 |
|---|---|---|
| 一圈𝓥_M | W^k | sinh产生L^{2j+1} |
| 两圈新增𝓥_N | W^k | cosh−1产生RL^{2j+2} |

$\beta=\bar\alpha$处L=0。𝓥_M的被积函数取连续值；𝓥_N的被积函数趋于0。它们都不需要主值处方。

## 7. 端点项χ是什么，为什么要保留

剩下的d来自一个已经积完的端点：

$$
\mathcal J_{2\epsilon}(1)
=\int_0^1\frac{(1-t)^{2\epsilon}-1}{t}\,dt
=-\gamma_E-\psi(1+2\epsilon).
$$

因此记

$$
\boxed{\chi(\epsilon)=\gamma_E+\psi(1+2\epsilon)=-d,\qquad
\psi(x)=\frac{d}{dx}\log\Gamma(x).}
\tag{chi}
$$

它只依赖ε，与z无关。零点附近

$$
\frac{\chi(\epsilon)}{\epsilon}
=2\zeta_2-4\zeta_3\epsilon+8\zeta_4\epsilon^2
-16\zeta_5\epsilon^3+\cdots .
$$

所以χ/ε没有极点。这个项来自§4.2保留下来的边界，不能为了简化公式而删除。

## 8. 合并：4e-04的完整母函数

由(QF)、(BR)、(N)及d=−χ，

$$
Q=\frac12\mathcal V_M(\epsilon,z)^2
-2\mathcal V_N(\epsilon,z)
+\frac{2\chi(\epsilon)}{\epsilon}\mathcal V_M(2\epsilon,z).
$$

最后乘回式(I)的物理前因子：

$$
\boxed{\begin{aligned}
I_{4e-04}(\epsilon,z)
={}&-\frac18\,
\exp\!\left(2\sum_{p=2}^\infty\frac{\zeta_p}{p}\epsilon^p\right)
\operatorname{sech}(\epsilon\ell)\\
&\times\left[
\frac12\mathcal V_M(\epsilon,z)^2
-2\mathcal V_N(\epsilon,z)
+\frac{2\chi(\epsilon)}{\epsilon}\mathcal V_M(2\epsilon,z)
\right].
\end{aligned}}
\tag{GF}
$$

这就是与GF cal一圈文档对应的母函数形式。指数中的级数在ε=0附近使用；未展开时它就是Z²。标准函数形式(QF)及上述积分恒等式先在z>0、实0<ε<1成立，并给出零附近的解析延拓。

| 因子或项 | 作用 |
|---|---|
| exp(2Σζ_p ε^p/p) | 两圈Gamma及唯一Euler补偿 |
| sech(εℓ) | 两圈crossing尺度组合，只产生偶次ℓ |
| 𝓥_M²/2 | 两个已知一圈子结构的乘积 |
| −2𝓥_N | 排序产生的额外R权重 |
| 2χ𝓥_M(2ε)/ε | 原函数端点留下的已知边界项 |

𝓥_M(2ε,z)只是同一一圈母核的参数替换，不是另一个二圈图；它不额外携带Gamma前因子。

## 9. 展开生成函数：M、N的系数逐项出现

由(VM)，

$$
\begin{aligned}
\mathcal V_M(\epsilon,z)
={}&M_{0,1}+\frac{\epsilon}{2}M_{1,1}\\
&+\epsilon^2\left(\frac{M_{2,1}}8+\frac{M_{0,3}}{24}\right)
+O(\epsilon^3).
\end{aligned}
$$

由(N)，

$$
\mathcal V_N(\epsilon,z)
=\frac{N_{0,2}}2+\frac{\epsilon N_{1,2}}2
+\epsilon^2\left(\frac{N_{2,2}}4+\frac{N_{0,4}}{24}\right)
+O(\epsilon^3).
$$

一般地，令𝓥_M=Σh_q ε^q、𝓥_N=Σn_q ε^q，则

$$
\boxed{\begin{aligned}
h_q&=\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{M_{q-2j,\,2j+1}}{2^q(q-2j)!(2j+1)!},\\
n_q&=\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{N_{q-2j,\,2j+2}}{(q-2j)!(2j+2)!}.
\end{aligned}}
\tag{series}
$$

各项在ε=0处均正则，因此4e-04没有负幂Laurent系数。对固定z，在足够小的复ε邻域内，端点幂和有限次对数插入有共同可积上界，允许上述逐项展开。

## 10. 通项公式与前两阶检验

将三部分按ε阶数卷积，Q=ΣQ_q ε^q的通项为

$$
\boxed{\displaystyle
Q_q=\frac12\sum_{i=0}^q h_i h_{q-i}
-2n_q
+2^{q+2}\sum_{k=0}^q(-1)^k\zeta_{k+2}h_{q-k}.}
\tag{coeff}
$$

最后一项的2^{q+2}来自χ的参数2ε与𝓥_M(2ε)共同缩放。

物理前因子的前几项是

$$
Z^2\operatorname{sech}(\epsilon\ell)
=1+\left(\zeta_2-\frac{\ell^2}{2}\right)\epsilon^2
+\frac{2\zeta_3}{3}\epsilon^3+O(\epsilon^4).
$$

### 10.1 零阶

由(coeff)，

$$
Q_0=\frac12M_{0,1}^2-N_{0,2}+4\zeta_2M_{0,1}.
$$

使用一圈低权重恒等式

$$
N_{0,2}=\frac13M_{0,1}^2+4\zeta_2M_{0,1},
\qquad M_{0,1}=-\frac{\ell^2+\pi^2}{2},
$$

便得

$$
\boxed{I_{4e-04}^{(0)}
=-\frac{M_{0,1}^2}{48}
=-\frac{(\ell^2+\pi^2)^2}{192}.}
$$

### 10.2 一阶

此时Gamma/cosh前因子没有一次项，所以I^(1)=−Q₁/8。代入
h₀=M₀,₁、h₁=M₁,₁/2、n₁=N₁,₂/2，

$$
Q_1=\frac12M_{0,1}M_{1,1}-N_{1,2}
+4\zeta_2M_{1,1}-8\zeta_3M_{0,1}.
$$

于是

$$
\boxed{
I_{4e-04}^{(1)}
=-\frac{M_{0,1}M_{1,1}}{16}
+\zeta_3M_{0,1}-\frac{\pi^2M_{1,1}}{12}
+\frac{N_{1,2}}8 .}
$$

这明确展示了第一项N怎样从通解展开得到，而不是先选N再匹配答案。ε²起还要卷积Gamma/cosh的高阶项，不能直接把−Q_q/8当作完整I的系数。

若写Z²sech(εℓ)=Σp_r ε^r，则任意阶最终系数就是

$$
I_{4e-04}^{(q)}=-\frac18\sum_{r=0}^q p_r Q_{q-r}.
$$

p_r由(GF)两个显式因子直接产生；其中偶次ℓ可以用ℓ²=−2M₀,₁−6ζ₂消去。故4e-04任意阶均可写成M/N及ζ常数的多项式。这里不再额外声称所有高j矩都能约化为j=1的子族。

## 11. 阅读范围与已有验证

正文给出了从原J到通解的完整路线；[有序Beta主稿](ordered-beta-generator.md)保存相应标准函数和词递推，[M/N高阶表](4e04-high-order.md)及[有序矩系数表](4e04-ordered-coefficients.md)列出完整I至ε⁴。

直接约化的13个源积分点和15个Q/Δ/I系数回归已经通过，见[原计算记录](../../../../workflow/checks/4e04-direct-mn/README.md)。本文本次重写的记号等价、源到通解的衔接及渲染检查见[编排检查](../../../../workflow/checks/4e04-derivation-edit/README.md)。[上一版原文](4e04-direct-mn-r001.md)保留原字节和R编号，旧检查按该版本解释。

当前结论限于正实z的4e-04及其裸交换配对；物理复z、M/N最小独立性和其他三圈图的闭包不在此结论内。解析计算与推导不等于Kernel/Lean形式整链，formal=not_run。本轮改写不新增数学结果或改变覆盖数。

---

## 附录A：关键转换式(BR)的Cauchy证明

本附录为正文§5的三条转换固定依据。为与既有计算逐项对应，以下临时记H=𝓥_M、$r=\bar\alpha$、d=𝓙_{2ε}(1)，其余Y、Y_c、V、E均沿用正文定义；附录中的J_a、T_a分别指正文的𝓙_a、𝓣_a。

定义

$$
C(a,b)=\frac{\Gamma(1+a-b)\Gamma(1+a+b)}{\Gamma(1+2a)},
\qquad \theta=z\partial_z.
$$

对矩源作一次分部积分，得到

$$
\boxed{\theta(\rho^{-a}\mathscr U(a,b))
=\rho^{-a}[C(a,0)-C(a,b)z^b].}
\tag{A1}
$$

推导可通过辅助Cauchy变换
𝓒_{s,t}(r)=PV∫₀¹x^s(1−x)^t/(x−r)dx组织。先令r在积分线外，积分
∂_x[x(1−x)x^s(1−x)^t/(x−r)]，得到

$$
r(1-r)\mathcal C'_{s,t}
=[s-(s+t)r]\mathcal C_{s,t}
-\frac{\Gamma(1+s)\Gamma(1+t)}{\Gamma(1+s+t)}.
$$

再取上下边界平均。𝓤的两个Cauchy留数相同，差中完全抵消，所以最终(A1)是普通可去积分的恒等式；没有给原图另加主值处方。此路线先在端点绝对收敛参数域成立，再延拓到小参数邻域。

记A=C(ε,0)ρ^{-ε}、D=ρ^{-ε}V、B=ρ^{-ε}Y。由(A1)及相应零指标减除式，

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
Ψ−ρ^{-ε}(Y−ℓV)−ℓ²/2仅差z无关常数。附录B证明该常数为π²/2，给出正文§5的第一条转换。

对(A1)分别在b=ε、−ε求b导数，使用C(ε,±ε)=1、
∂_bC(ε,ε)=−d、∂_bC(ε,−ε)=d，即得

$$
\theta\{\rho^{-\epsilon}[2Y_c+2\epsilon dH(2\epsilon)]\}
=-\ell\rho^{-\epsilon}(z^\epsilon+z^{-\epsilon}).
$$

这也等于θ{Ψ(2ε,0)−ℓ[J_{2ε}(−z)−J_{2ε}(−1/z)]+2M₀,₁}。附录B固定其常数差为0，得到正文§5的第二条转换。

## 附录B：边界常数如何解析固定

以下只用Gamma/Beta渐近，不读取任何Laurent答案。令u=log z→−∞，先取0<ε<1/4，保证两次使用的幂a=ε、2ε都在(0,1)。

### B.1 两个统一的渐近公式

由Euler积分及[Gauss连接公式](https://dlmf.nist.gov/15.8#E2)，辅助Cauchy变换有精确分解

$$
\mathcal C_{s,t}(r)
=\mathrm B(s,t+1)\,
{}_2F_1(1,-s-t;1-s;r)
-\pi\cot(\pi s)\,r^s(1-r)^t .
\tag{B1}
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
\tag{B2}
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
\tag{B3}
$$

这些系数也可直接由(B2)在σ=0的Gamma展开读取。主分支始终取−y的负实超几何自变量；没有穿越其[1,∞)割线。

### B.2 第一条转换的常数是π²/2

记P=πcot(πε)、E₂=π²csc²(πε)、A=C(ε,0)ρ^{-ε}。由(B1)及b=0处的参数导数，

$$
D=\frac A\epsilon+u-P+o(1),\qquad
B=A\left(\frac u\epsilon+\frac1{\epsilon^2}\right)-E_2+o(1).
$$

由(B3)、𝔉(ε,ε;1)=[ψ₁(1+ε)−ζ₂]/2及[polygamma反射](https://dlmf.nist.gov/5.15)，

$$
\Psi(\epsilon,\epsilon)
=\frac A{\epsilon^2}-\frac{u^2}{2}+Pu-E_2+\frac{\pi^2}{2}+o(1).
$$

所以Ψ−(B−uD)−u²/2→π²/2，固定正文§5的第一条转换的常数。不是通过某个有限阶拟合补常数。

### B.3 第二条转换的常数是0

记a=2ε、d=J_a(1)、P=πcot(πa)、E₂=π²csc²(πa)、T₁=T_a(1)，并设

$$
d_T=-\frac12[(d-P)^2+E_2-\psi_1(1+a)+\zeta_2].
$$

由(B3)，正文§5的第二条转换左侧除以ρ^ε后的渐近为

$$
(1+1/z)^a\left(\frac u a+\frac1{a^2}\right)
-\frac{u^2}{2}+d_T-2T_1-\pi^2+o(1).
$$

由(B1)分别对b=ε、−ε求导，右侧除以ρ^ε后的渐近具有相同幂项和−u²/2，常数为

$$
\frac{d^2-\psi_1(1+a)-5\zeta_2}{2}-E_2+dP.
$$

在b=ε一端，先对(B1)的完整Gauss项消去s=0的可去极点并作源导数，再取r→0；高阶r项给出的余项仍趋零。保留到O(1)时，所需合并极点的局部式为
B(s,1+a-s)=1/s+d+s[d²+ζ₂+ψ₁(1+a)]/2+O(s²)，
再乘e^{-su}并加入−πcot(πs)。这说明上述常数中的5ζ₂不能省略。
用T₁=[ψ₁(1+a)−ζ₂−d²]/2及E₂−P²=π²，两常数精确相同。因此正文§5的第二条转换的常数差为0。所有渐近余项在此参数范围及所用参数导数下趋零；随后以原可去积分的解析性延拓。
