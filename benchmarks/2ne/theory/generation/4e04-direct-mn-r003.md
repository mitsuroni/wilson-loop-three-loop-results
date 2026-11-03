# 4e-04的通解推导：从原积分到M/N生成函数

本文按“原积分 → 降维 → 整理剩余排序 → 一维生成函数 → Laurent展开”的顺序推导4e-04。组织方式参考[GF cal中的2e-02推导](<../../../../GF cal 推导草稿/2e-02.md>)，归一化采用项目当前约定。

与2e-02相比，多出来的工作是处理两条对边传播子的相对排序。最后它会变成三个明确部分：**一圈生成函数的平方、N的生成函数、端点产生的一圈边界项**。正文逐步说明它们怎样出现；涉及Cauchy变换与边界常数的证明放在文末附录。

## 1. 起点：原参数积分与物理归一化

### 1.1 原四重有序积分

令 $z>0$，先取实 $0<\epsilon<1$。原参数对应关系为
$t_1=s_{21}$、$t_2=s_{22}$、$s_1=s_{41}$、$s_2=s_{42}$，原外到内顺序为 $t_2,s_2,t_1,s_1$。记

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

两个低端点相配，两个高端点相配。原输入、颜色和顺序见[逐图根推导](../../workspaces/4e-04/derivation.md)。以下换序是下游操作：在 $0<\epsilon<1$ 时，把 $\Omega$ 放大成两个独立核的正方形乘积，即得到可积上界，允许所用的Fubini换序。

### 1.2 去掉简单尺度幂

统一使用

$$
\alpha=\frac1{1+z},\qquad
\bar\alpha=\frac{z}{1+z},\qquad
\rho=\alpha\bar\alpha=\frac{z}{(1+z)^2},\qquad
\ell=\log z.
$$

定义正规化母核 $Q=z^{-\epsilon}J_{4e-04}$。两条边同时反射并反转排序，利用
$D_z(1-t,1-s)=zD_{1/z}(t,s)$，可得
$J(1/z)=z^{-2\epsilon}J(z)$，所以 $Q(1/z)=Q(z)$。

本图 $\mathrm{sym}=2$、$a_1=-1/2$，因此 $\kappa=\mathrm{sym}\,a_1/[2(-2)^2]=-1/8$。记

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

**从现在起只需求 $Q$，最后再乘这一前因子。** $Z$ 已精确包含Euler补偿，不另行“丢弃 $\gamma_E$ 项”或再乘一次补偿。

## 2. 先回顾一圈：M的生成函数怎样出现

这一步既与2e-02对照，也定义两圈通解中将出现的一圈子结构。

### 2.1 积去一个线性变量

定义不含Gamma前因子的一圈母核

$$
\mathcal V_M(\epsilon,z)
=-(1+z)z^{-\epsilon/2}
\int_0^1dt\int_0^1ds\,D_z(t,s)^{\epsilon-1}.
$$

对 $s$ 用基本原函数

$$
\int_a^b(A+Bx)^{\epsilon-1}dx
=\frac{(A+Bb)^\epsilon-(A+Ba)^\epsilon}{\epsilon B}
$$

积分，再令 $t=1-\beta$，得到

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

这就是GF cal的2e-02文档中作为积分部分的生成函数。本文记作 $\mathcal V_M$，以免与物理前因子混淆；此前研究稿的裸单核 $H$ 也等于它。

## 3. 两圈第一次降维：4D → 2D

令 $u=t_2$、$v=s_1$。先积 $t_1$，再积 $s_2$；两步都使用§2.1的线性幂原函数：

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

记 $D=D_z(u,v)$，则

$$
\boxed{\displaystyle
Q=\frac{(1+z)^2z^{-\epsilon}}{\epsilon^2}
\int_0^1du\int_0^1dv\,
\frac{[D^\epsilon-(1-v)^\epsilon][(zu)^\epsilon-D^\epsilon]}
{[(1+z)v-1][(1+z)u-1]}.}
\tag{Q2D}
$$

**与一圈的差别在这里：现在得到两个幂差的乘积，它们共同依赖 $D_z(u,v)$。** 不能把它们当成两个独立一圈积分。

分母为零的点由原函数连续延拓。例如第一个差商在 $(1+z)v-1=0$ 时趋于 $u(1-v)^{\epsilon-1}$。必须先保留差的整体，不能把四个展开项各自当成收敛积分。

## 4. 算出剩余排序积分：先求导降层，再积分还原

第三节已得到二维差商式(Q2D)。这一节完成它的求值，分成三个动作：

1. 把两个分母变成 $p$、$q$，合法地拆开可去差商。
2. 对剩余排序积分求 $z$ 导数，消去一层嵌套原函数。
3. 求完普通 Beta 积分后，再对 $z$ 积分还原原函数，并保留边界项。

### 4.1 中心换元：让两个分母具有相同形式

令

$$
p=(1+z)u-1,\qquad q=(1+z)v-1.
$$

于是 $p,q\in(-1,z)$，两个差商的分母正好变成 $pq$，而共同核变成

$$
D_z(u,v)=\frac{z+pq}{1+z}.
$$

为缩短这一小节的代数，暂记

$$
A(p)=(1+p)^\epsilon,\quad
B(p)=(1-p/z)^\epsilon,\quad
C(p,q)=(1+pq/z)^\epsilon .
$$

代入(Q2D)，得到

$$
Q=\frac{\rho^\epsilon}{\epsilon^2}
\int_{-1}^z\frac{dp}{p}\int_{-1}^z\frac{dq}{q}\,
[C(p,q)-B(q)][A(p)-C(p,q)].
$$

这一步只是换元，排序耦合仍保留在 $C(p,q)$ 中。

### 4.2 加入一个积分为零的项，使每个拆项都可积

我们希望先积掉 $q$，但不能直接把上式的四个幂项分开。做法是将分子改写为

$$
[A(p)+B(p)](C-1)-(C^2-1)-[A(p)-1][B(q)-1].
$$

**新旧分子并不逐点相等。** 它们之差是

$$
[B(p)-B(q)][C(p,q)-1].
$$

交换 $p,q$ 后，这个差变号，而积分区域和因子 $1/(pq)$ 不变，所以它的积分为零。因而替换前后的 $Q$ 相同。

替换的好处是：每个拆项在 $p=0$、$q=0$ 处都可去。现在可以逐项积分，而不引入主值或假极点。

定义两个普通原函数：

$$
\mathcal J_a(x)=\int_0^x\frac{(1-t)^a-1}{t}\,dt,\qquad
\mathcal T_a(x)=\int_0^x\frac{\mathcal J_a(t)}t\,dt.
$$

积掉 $q$ 时使用

$$
\int_{-1}^z\frac{(1+pq/z)^a-1}{q}\,dq
=\mathcal J_a(-p)-\mathcal J_a(p/z).
$$

因此剩下的是单变量 $p$ 上的原函数乘积。将负半区用 $p=-v$、正半区用 $p=zv$ 缩放到 $(0,1)$，这些乘积可以统一为下面的排序积分。

### 4.3 写出需要求值的排序积分

两个边界幂的指数分别记为 $a,b$，定义带减除的权重

$$
\mathcal A_a(v;z)=(1+zv)^a-1,\qquad
\mathcal B_b(v)=(1-v)^b-1.
$$

剩余的反对称排序部分是

$$
K(a,b;z)=
\int_{0<u<v<1}
\frac{\mathcal A_a(u;z)\mathcal B_b(v)
-\mathcal B_b(u)\mathcal A_a(v;z)}{uv}\,du\,dv .
$$

把内层积分分别写成原函数：

$$
f(v)=\mathcal J_a(-zv),\qquad
g(v)=\mathcal J_b(v).
$$

因为

$$
f'(v)=\frac{\mathcal A_a(v;z)}v,\qquad
g'(v)=\frac{\mathcal B_b(v)}v,
$$

所以上面的二重积分就是

$$
K=\int_0^1(fg'-gf')\,dv.
$$

此时仍然是“积分套积分”。下一步的目的是消去这一层原函数。

### 4.4 为什么对 $z$ 求导：它恰好消去 $f$ 的原函数

定义 $\theta_z=z\partial_z$。由 $\mathcal J_a$ 的定义，

$$
\begin{aligned}
\theta_z f(v)
&=(-zv)\mathcal J_a'(-zv)\\
&=(1+zv)^a-1
=\mathcal A_a(v;z).
\end{aligned}
$$

**右边已经不含任何积分。** 而 $g(v)$ 不依赖 $z$，所以

$$
\begin{aligned}
\theta_z K
&=\int_0^1
\left[\mathcal A_a(v;z)g'(v)
-g(v)\partial_v\mathcal A_a(v;z)\right]dv\\
&=P(a,b;z)-\int_0^1g(v)\partial_v\mathcal A_a(v;z)\,dv,
\end{aligned}
$$

其中

$$
P(a,b;z)=
\int_0^1\frac{\mathcal A_a(v;z)\mathcal B_b(v)}v\,dv .
$$

接着只对第二项分部积分：

$$
-\int_0^1g\,\partial_v\mathcal A_a\,dv
=-[g\mathcal A_a]_0^1+\int_0^1g'\mathcal A_a\,dv .
$$

最后一个积分又等于 $P$。在 $v=0$，边界乘积为零；在 $v=1$，它等于 $\mathcal J_b(1)[(1+z)^a-1]$。因此

$$
\boxed{\theta_z K
=2P(a,b;z)-\mathcal J_b(1)[(1+z)^a-1].}
$$

现在嵌套原函数已经消失，右边只剩普通一维幂积分和明确的端点项。端点项必须保留。

### 4.5 用 Beta 积分求出 $P$，再积回 $K$

先在 $|z|<1$ 展开第一个幂：

$$
(1+zv)^a-1
=\sum_{n\ge1}\frac{(-a)_n(-z)^n}{n!}v^n .
$$

每一项都能直接积分：

$$
\int_0^1v^{n-1}[(1-v)^b-1]\,dv
=\mathrm B(n,1+b)-\frac1n.
$$

于是

$$
\begin{aligned}
P
&=\sum_{n\ge1}\frac{(-a)_n(-z)^n}{n!}
\left[\mathrm B(n,1+b)-\frac1n\right]\\
&=\sum_{n\ge1}\frac{(-a)_n}{(1+b)_n}\frac{(-z)^n}{n}
-\mathcal J_a(-z).
\end{aligned}
$$

要从 $\theta_zK$ 恢复 $K$，还需作一次 $dz/z$ 积分。它将每个 $z^n/n$ 变成 $z^n/n^2$，所以自然定义

$$
\boxed{\displaystyle
\mathfrak F(a,b;x)
=\sum_{n\ge1}\frac{(-a)_n}{(1+b)_n}\frac{x^n}{n^2}
=-\frac{ax}{1+b}\,
{}_4F_3\!\left(\begin{matrix}1,1,1,1-a\\2,2,2+b\end{matrix};x\right).}
\tag{F}
$$

这里的超几何函数符号只是已求得级数的标准名称。利用 $\mathcal T_a=\mathfrak F(a,0)$，逐项积回去：

$$
K(a,b;z)
=2\mathfrak F(a,b;-z)-2\mathfrak F(a,0;-z)
-\mathcal J_b(1)\mathcal J_a(-z).
$$

积分常数由 $K(a,b;0)=0$ 固定。先得到的小正 $z$ 等式沿负实超几何自变量延拓，即覆盖所有 $z>0$。

### 4.6 把结果代回原来的 $Q$

本步只做边界账本的代回，不再引入新的求解方法。记

$$
\begin{gathered}
j_-=\mathcal J_\epsilon(-z),\quad
j_0=\mathcal J_\epsilon(1),\quad
j_+=\mathcal J_\epsilon(-1/z),\\
\mathcal C_a=\mathcal T_a(-z)+\mathcal T_a(-1/z)-2\mathcal T_a(1).
\end{gathered}
$$

§4.2的剩余一维积分整理为

$$
\begin{aligned}
Q=\frac{\rho^\epsilon}{\epsilon^2}\bigg[
&K(\epsilon,\epsilon;z)+K(\epsilon,\epsilon;1/z)\\
&+\frac{(j_-+j_+)^2}{2}-j_0(j_-+j_+)
+2\mathcal C_\epsilon-\mathcal C_{2\epsilon}\bigg].
\end{aligned}
$$

将刚求出的 $K$ 代入，并用

$$
\mathcal V_M=-\frac{\rho^{\epsilon/2}}{\epsilon}(j_-+j_+-2j_0),
\qquad
\mathfrak F(\epsilon,\epsilon;1)=\mathcal T_\epsilon(1)+\frac{j_0^2}{2},
$$

边界平方与 $\mathcal T_\epsilon$ 项便能合并。第二个端点关系可在辅助函数中取 $K(\epsilon,\epsilon;-1)=0$ 得到：此时两种权重都变成 $(1-v)^\epsilon-1$，交换子逐点为零。这只用于求辅助函数的端点常数，原图仍取 $z>0$。

为缩短结果，定义同一函数的反射组合

$$
\Psi(a,b;z)=\mathfrak F(a,b;-z)+\mathfrak F(a,b;-1/z)
-2\mathfrak F(a,b;1).
$$

最终得到

$$
\boxed{\displaystyle
Q=\frac12\mathcal V_M(\epsilon,z)^2
+\frac{\rho^\epsilon}{\epsilon^2}
[2\Psi(\epsilon,\epsilon;z)-\Psi(2\epsilon,0;z)].}
\tag{QF}
$$

这是第四节的终点：原积分已写成有限个标准函数。下一节不再计算同一个排序积分，而是证明它与一维矩源之间的转换关系。

## 5. 转成M/N的积分语言：比较导数，再固定常数

从(QF)可见，已知的一圈平方之外，只需处理

$$
\rho^\epsilon[2\Psi(\epsilon,\epsilon;z)-\Psi(2\epsilon,0;z)].
$$

本节的目标是把它写到共同分母 $\beta-\bar\alpha$ 下。方法是：**分别求两种表示对 $z$ 的导数，证明差只剩常数，再用解析边界固定常数。** 这一步不是普通换元，也不是把标准函数直接认作某个基底。

### 5.1 准备同一个矩源中的三个插入

使用一圈已经有的双源

$$
\mathscr U(a,b;z)=\int_0^1
\frac{e^{aW}(e^{bL}-1)}{\beta-\bar\alpha}\,d\beta.
$$

对 $b$ 求导，积分里就多一个 $L$。取 $b=0$ 或 $b=\pm\epsilon$，得到

$$
\begin{aligned}
Y
&=\left.\partial_b\mathscr U(\epsilon,b)\right|_{b=0}
=\int_0^1\frac{e^{\epsilon W}L}{\beta-\bar\alpha}\,d\beta,\\
Y_c
&=\frac12\left[
\left.\partial_b\mathscr U(\epsilon,b)\right|_{b=\epsilon}
+\left.\partial_b\mathscr U(\epsilon,b)\right|_{b=-\epsilon}
\right]\\
&=\int_0^1\frac{e^{\epsilon W}L\cosh(\epsilon L)}
{\beta-\bar\alpha}\,d\beta,\\
E
&=\mathscr U(\epsilon,\epsilon)+\mathscr U(\epsilon,-\epsilon)
=2\int_0^1\frac{e^{\epsilon W}[\cosh(\epsilon L)-1]}
{\beta-\bar\alpha}\,d\beta.
\end{aligned}
$$

它们的作用可以直接从被积函数读出：$Y$ 有一个 $L$，$Y_c$ 多了 $\cosh(\epsilon L)$，$E$ 则只有偶幂差。它们都是同一个源的取值或导数。

证明中还暂时需要一个减除量

$$
V=\int_0^1
\frac{e^{\epsilon W}-\rho^\epsilon}{\beta-\bar\alpha}\,d\beta.
$$

它的分子在 $\beta=\bar\alpha$ 处为零，内部点可去。$V$ 最后会消去，不加入最终基底。

### 5.2 第一条转换：把 $\Psi(\epsilon,\epsilon)$ 接到 $Y-\ell V$

令

$$
A_\epsilon(z)=
\frac{\Gamma(1+\epsilon)^2}{\Gamma(1+2\epsilon)}\rho^{-\epsilon}.
$$

矩源的分部积分方程给出前两条导数关系，有序Beta源的方程给出第三条：

$$
\begin{aligned}
\theta_z(\rho^{-\epsilon}V)&=1-A_\epsilon,\\
\theta_z(\rho^{-\epsilon}Y)&=-\ell A_\epsilon,\\
\theta_z^2\Psi(\epsilon,\epsilon)&=A_\epsilon-1.
\end{aligned}
$$

这些方程的来源见附录A。先把第一、第三式相加：

$$
\theta_z\left[
\theta_z\Psi(\epsilon,\epsilon)+\rho^{-\epsilon}V
\right]=0.
$$

方括号只可能是常数。在 $z=1$，$\Psi$ 的反射偶性使其一阶导数为零，而 $V=0$，所以这个常数为零：

$$
\theta_z\Psi(\epsilon,\epsilon)=-\rho^{-\epsilon}V.
$$

现在寻找具有同样导数的矩组合。逐项求导：

$$
\begin{aligned}
&\theta_z\left[
\rho^{-\epsilon}(Y-\ell V)+\frac{\ell^2}{2}
\right]\\
&=-\ell A_\epsilon-\rho^{-\epsilon}V
-\ell(1-A_\epsilon)+\ell\\
&=-\rho^{-\epsilon}V.
\end{aligned}
$$

因此这两个表达式只差一个与 $z$ 无关的常数。附录B.2解析算得该常数为 $\pi^2/2$。使用
$M_{0,1}=-(\ell^2+\pi^2)/2$，得到

$$
\boxed{
\rho^\epsilon\Psi(\epsilon,\epsilon)
=Y-\ell V-\rho^\epsilon M_{0,1}.}
\tag{BR1}
$$

这才是第一条转换的依据：导数相同，加上已固定的常数。

### 5.3 第二条转换：处理 $\Psi(2\epsilon,0)$

为书写边界值，记

$$
J_-=\mathcal J_{2\epsilon}(-z),\qquad
J_+=\mathcal J_{2\epsilon}(-1/z),\qquad
d=\mathcal J_{2\epsilon}(1).
$$

$d$ 的显式 digamma 值在第7节给出。由原函数定义，

$$
\begin{aligned}
\theta_z\Psi(2\epsilon,0)&=J_--J_+,\\
\theta_z(J_--J_+)&=(1+z)^{2\epsilon}
+(1+1/z)^{2\epsilon}-2.
\end{aligned}
$$

于是乘积求导中的 $J_--J_+$ 抵消，且 $\theta_zM_{0,1}=-\ell$，所以

$$
\begin{aligned}
&\theta_z\left[
\Psi(2\epsilon,0)-\ell(J_--J_+)+2M_{0,1}
\right]\\
&=-\ell\rho^{-\epsilon}(z^\epsilon+z^{-\epsilon}).
\end{aligned}
$$

另一方面，把矩源方程分别在 $b=\epsilon$、$b=-\epsilon$ 处求源导数并相加，可得

$$
\theta_z\left\{
\rho^{-\epsilon}
[2Y_c+2\epsilon d\,\mathcal V_M(2\epsilon,z)]
\right\}
=-\ell\rho^{-\epsilon}(z^\epsilon+z^{-\epsilon}).
$$

两边又有相同的导数；附录B.3证明这次的常数差为零。因此

$$
\boxed{\begin{aligned}
\rho^\epsilon\big[
\Psi(2\epsilon,0)-\ell(J_--J_+)+2M_{0,1}
\big]
=2Y_c+2\epsilon d\,\mathcal V_M(2\epsilon,z).
\end{aligned}}
\tag{BR2}
$$

### 5.4 第三条转换：用一次换元消去 $V$

前两条式子含有中间量 $V$。要消去它，先展开 $E$ 的两个幂：

$$
e^{\epsilon W}e^{\epsilon L}
=[z(1-\beta)^2]^\epsilon,\qquad
e^{\epsilon W}e^{-\epsilon L}
=[\beta^2/z]^\epsilon.
$$

各自加减 $\rho^\epsilon$ 后，可以直接用 $\mathcal J_{2\epsilon}$ 积分。第一项用
$t=(\beta-\bar\alpha)/(1-\bar\alpha)$，第二项用 $t=1-\beta/\bar\alpha$，得到

$$
\begin{aligned}
\mathscr U(\epsilon,\epsilon)&=\rho^\epsilon(d-J_-)-V,\\
\mathscr U(\epsilon,-\epsilon)&=\rho^\epsilon(J_+-d)-V.
\end{aligned}
$$

两式相加，端点值 $d$ 消去：

$$
\boxed{E=\rho^\epsilon(J_+-J_-)-2V.}
\tag{BR3}
$$

这条关系不需要新的特殊函数或常数。

### 5.5 把三条转换逐项代回去

将(BR1)、(BR2)代入本节要处理的组合：

$$
\begin{aligned}
&\rho^\epsilon[2\Psi(\epsilon,\epsilon)-\Psi(2\epsilon,0)]\\
&=2(Y-\ell V-\rho^\epsilon M_{0,1})\\
&\quad-\left[
2Y_c+2\epsilon d\,\mathcal V_M(2\epsilon,z)
+\rho^\epsilon\ell(J_--J_+)-2\rho^\epsilon M_{0,1}
\right]\\
&=2Y-2Y_c-2\epsilon d\,\mathcal V_M(2\epsilon,z)
+\ell\left[\rho^\epsilon(J_+-J_-)-2V\right].
\end{aligned}
$$

这里 $M_{0,1}$ 项已经抵消。最后的方括号由(BR3)正好等于 $E$，因此

$$
\boxed{\displaystyle
\rho^\epsilon[2\Psi(\epsilon,\epsilon)-\Psi(2\epsilon,0)]
=2Y-2Y_c+\ell E
-2\epsilon d\,\mathcal V_M(2\epsilon,z).}
\tag{BR}
$$

现在中间量 $V$ 和两个端点原函数的差都已消去。留下的 $Y,Y_c,E$ 具有同一积分分母，可以直接合并被积函数；第六节将由此得到 $R=L-\ell$，也就是 N 的权重。

## 6. 双曲分解后，N自然出现

现在把 $Y,Y_c,E$ 的积分放到同一分母上：

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

**额外的 $R$ 不是猜出来的：它来自 $L$ 插入减去外部的 $\ell$ 插入。** 于是定义正则的N生成函数

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
| 一圈 $\mathcal V_M$ | $W^k$ | $\sinh$ 产生 $L^{2j+1}$ |
| 两圈新增 $\mathcal V_N$ | $W^k$ | $\cosh-1$ 产生 $RL^{2j+2}$ |

在 $\beta=\bar\alpha$ 处有 $L=0$。$\mathcal V_M$ 的被积函数取连续值；$\mathcal V_N$ 的被积函数趋于零。它们都不需要主值处方。

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

它只依赖 $\epsilon$，与 $z$ 无关。零点附近

$$
\frac{\chi(\epsilon)}{\epsilon}
=2\zeta_2-4\zeta_3\epsilon+8\zeta_4\epsilon^2
-16\zeta_5\epsilon^3+\cdots .
$$

所以 $\chi(\epsilon)/\epsilon$ 没有极点。这个项来自§4.2保留下来的边界，不能为了简化公式而删除。

## 8. 合并：4e-04的完整母函数

由(QF)、(BR)、(N)及 $d=-\chi(\epsilon)$，

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

这就是与GF cal一圈文档对应的母函数形式。指数中的级数在 $\epsilon=0$ 附近使用；未展开时它就是 $Z^2$。标准函数形式(QF)及上述积分恒等式先在 $z>0$、实 $0<\epsilon<1$ 成立，并给出零附近的解析延拓。

| 因子或项 | 作用 |
|---|---|
| $\exp(2\sum_{p\ge2}\zeta_p\epsilon^p/p)$ | 两圈Gamma及唯一Euler补偿 |
| $\operatorname{sech}(\epsilon\ell)$ | 两圈反射尺度组合，只产生偶次 $\ell$ |
| $\mathcal V_M^2/2$ | 两个已知一圈子结构的乘积 |
| $-2\mathcal V_N$ | 排序产生的额外 $R$ 权重 |
| $2\chi(\epsilon)\mathcal V_M(2\epsilon)/\epsilon$ | 原函数端点留下的已知边界项 |

$\mathcal V_M(2\epsilon,z)$ 只是同一一圈母核的参数替换，不是另一个二圈图；它不额外携带Gamma前因子。

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

一般地，令 $\mathcal V_M=\sum_{q\ge0}h_q\epsilon^q$、$\mathcal V_N=\sum_{q\ge0}n_q\epsilon^q$，则

$$
\boxed{\begin{aligned}
h_q&=\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{M_{q-2j,\,2j+1}}{2^q(q-2j)!(2j+1)!},\\
n_q&=\sum_{j=0}^{\lfloor q/2\rfloor}
\frac{N_{q-2j,\,2j+2}}{(q-2j)!(2j+2)!}.
\end{aligned}}
\tag{series}
$$

各项在 $\epsilon=0$ 处均正则，因此4e-04没有负幂Laurent系数。对固定 $z$，在足够小的复 $\epsilon$ 邻域内，端点幂和有限次对数插入有共同可积上界，允许上述逐项展开。

## 10. 通项公式与前两阶检验

将三部分按 $\epsilon$ 阶数卷积，$Q=\sum_{q\ge0}Q_q\epsilon^q$ 的通项为

$$
\boxed{\displaystyle
Q_q=\frac12\sum_{i=0}^q h_i h_{q-i}
-2n_q
+2^{q+2}\sum_{k=0}^q(-1)^k\zeta_{k+2}h_{q-k}.}
\tag{coeff}
$$

最后一项的 $2^{q+2}$ 来自 $\chi$ 的参数 $2\epsilon$ 与 $\mathcal V_M(2\epsilon)$ 共同缩放。

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

此时Gamma/cosh前因子没有一次项，所以 $I^{(1)}=-Q_1/8$。代入
$h_0=M_{0,1}$、$h_1=M_{1,1}/2$、$n_1=N_{1,2}/2$，

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

这明确展示了第一项N怎样从通解展开得到，而不是先选N再匹配答案。从 $\epsilon^2$ 起还要卷积Gamma/cosh的高阶项，不能直接把 $-Q_q/8$ 当作完整 $I$ 的系数。

若写 $Z^2\operatorname{sech}(\epsilon\ell)=\sum_{r\ge0}p_r\epsilon^r$，则任意阶最终系数就是

$$
I_{4e-04}^{(q)}=-\frac18\sum_{r=0}^q p_r Q_{q-r}.
$$

$p_r$ 由(GF)两个显式因子直接产生；其中偶次 $\ell$ 可以用 $\ell^2=-2M_{0,1}-6\zeta_2$ 消去。故4e-04任意阶均可写成M/N及ζ常数的多项式。这里不再额外声称所有高 $j$ 矩都能约化为 $j=1$ 的子族。

## 11. 阅读范围与已有验证

正文给出了从原J到通解的完整路线；[有序Beta主稿](ordered-beta-generator.md)保存相应标准函数和词递推，[M/N高阶表](4e04-high-order.md)及[有序矩系数表](4e04-ordered-coefficients.md)列出完整 $I$ 至 $\epsilon^4$。

直接约化的13个源积分点和15个 $Q,\Delta,I$ 系数回归已经通过，见[原计算记录](../../../../workflow/checks/4e04-direct-mn/README.md)。[r001原稿](4e04-direct-mn-r001.md)和[r002讲解稿](4e04-direct-mn-r002.md)均保留原字节，历史检查按相应版本解释。当前第四、第五节的逐步说明及全稿行内数学修复，见[本次检查](../../../../workflow/checks/4e04-readable-sections/README.md)；前一次记号编排另见[原编排检查](../../../../workflow/checks/4e04-derivation-edit/README.md)。

当前结论限于正实 $z$ 的4e-04及其裸交换配对；物理复 $z$、M/N最小独立性和其他三圈图的闭包不在此结论内。解析计算与推导不等于Kernel/Lean形式整链，形式验证未运行。本轮改写不新增数学结果或改变覆盖数。

---

## 附录A：关键转换式(BR)的Cauchy证明

本附录为正文§5的三条转换固定依据。为与既有计算逐项对应，以下临时记 $H=\mathcal V_M$、$r=\bar\alpha$、$d=\mathcal J_{2\epsilon}(1)$；其余 $Y,Y_c,V,E$ 均沿用正文定义。附录中的 $J_a,T_a$ 分别指正文的 $\mathcal J_a,\mathcal T_a$。

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

为推导矩源方程，先引入辅助Cauchy变换：

$$
\mathcal C_{s,t}(r)=\operatorname{PV}\int_0^1
\frac{x^s(1-x)^t}{x-r}\,dx.
$$

先令 $r$ 在积分线外，对下面的全导数积分：

$$
\partial_x\left[\frac{x(1-x)x^s(1-x)^t}{x-r}\right].
$$

得到

$$
r(1-r)\mathcal C'_{s,t}
=[s-(s+t)r]\mathcal C_{s,t}
-\frac{\Gamma(1+s)\Gamma(1+t)}{\Gamma(1+s+t)}.
$$

再取上下边界平均。$\mathscr U$ 的两个Cauchy留数相同，差中完全抵消，所以最终(A1)是普通可去积分的恒等式；没有给原图另加主值处方。此路线先在端点绝对收敛参数域成立，再延拓到小参数邻域。

记 $A=C(\epsilon,0)\rho^{-\epsilon}$、$D=\rho^{-\epsilon}V$、$B=\rho^{-\epsilon}Y$。由(A1)及相应零指标减除式，

$$
\theta D=1-A,\qquad \theta B=-\ell A.
$$

另一方面，由有序Beta源方程及Gauss反射，

$$
G_\epsilon(z)+G_\epsilon(1/z)-1=C(\epsilon,0)\rho^{-\epsilon},
\quad
G_\epsilon(z)={}_2F_1(1,-\epsilon;1+\epsilon;-z),
$$

因此 $\theta^2\Psi(\epsilon,\epsilon)=A-1$。$D$ 在 $z\mapsto1/z$ 下为奇，$\Psi$ 为偶，故 $\theta\Psi=-D$。继续积分可知

$$
\Psi-\rho^{-\epsilon}(Y-\ell V)-\frac{\ell^2}{2}
$$

只差与 $z$ 无关的常数。附录B证明该常数为 $\pi^2/2$，给出正文§5的第一条转换。

对(A1)分别在 $b=\epsilon$、$b=-\epsilon$ 处求 $b$ 导数，使用

$$
C(\epsilon,\pm\epsilon)=1,\qquad
\partial_bC(\epsilon,\epsilon)=-d,\qquad
\partial_bC(\epsilon,-\epsilon)=d,
$$

即得

$$
\theta\{\rho^{-\epsilon}[2Y_c+2\epsilon dH(2\epsilon)]\}
=-\ell\rho^{-\epsilon}(z^\epsilon+z^{-\epsilon}).
$$

这也等于

$$
\theta\left\{\Psi(2\epsilon,0)
-\ell[J_{2\epsilon}(-z)-J_{2\epsilon}(-1/z)]+2M_{0,1}\right\}.
$$

附录B固定其常数差为零，得到正文§5的第二条转换。

## 附录B：边界常数如何解析固定

以下只用Gamma/Beta渐近，不读取任何Laurent答案。令 $u=\log z\to-\infty$，先取 $0<\epsilon<1/4$，保证两次使用的幂 $a=\epsilon,2\epsilon$ 都在 $(0,1)$。

### B.1 两个统一的渐近公式

由Euler积分及[Gauss连接公式](https://dlmf.nist.gov/15.8#E2)，辅助Cauchy变换有精确分解

$$
\mathcal C_{s,t}(r)
=\mathrm B(s,t+1)\,
{}_2F_1(1,-s-t;1-s;r)
-\pi\cot(\pi s)\,r^s(1-r)^t .
\tag{B1}
$$

先取非整数 $0<\operatorname{Re}s<1$，再按消去可去极点的方式延拓。因此

$$
\frac{\mathcal C_{s,t}}{r^s(1-r)^t}
=\mathrm B(s,t+1)r^{-s}(1-r)^{-t}-\pi\cot(\pi s)+o(1).
$$

在 $s\to0$ 时，必须把Beta项和cot项合起来展开，不能丢掉其中一项。

对 $\mathfrak F(a,b;-y)$，Gauss函数的Mellin核为

$$
\int_0^\infty y^{\sigma-1}
{}_2F_1(1,-a;1+b;-y)\,dy
=\frac{\Gamma(\sigma)\Gamma(1-\sigma)
\Gamma(-a-\sigma)\Gamma(1+b)}
{\Gamma(-a)\Gamma(1+b-\sigma)}.
\tag{B2}
$$

先在 $\operatorname{Re}a<0$、$\operatorname{Re}b>0$ 及 $0<\operatorname{Re}\sigma<\min(1,-\operatorname{Re}a)$ 由两次Euler/Beta积分推出，再作参数延拓。由于 $\mathfrak F$ 是Gauss函数减一的两次Euler原函数，$\sigma=0$ 处的三重极点固定 $\log^2 y$、$\log y$ 和常数；延拓到 $0<a<1$ 时，另保留移动幂极点的贡献。得到

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

这些系数也可直接由(B2)在 $\sigma=0$ 的Gamma展开读取。主分支始终取 $-y$ 的负实超几何自变量；没有穿越其 $[1,\infty)$ 割线。

### B.2 第一条转换的常数是 $\pi^2/2$

记 $P=\pi\cot(\pi\epsilon)$、$E_2=\pi^2\csc^2(\pi\epsilon)$、$A=C(\epsilon,0)\rho^{-\epsilon}$。由(B1)及 $b=0$ 处的参数导数，

$$
D=\frac A\epsilon+u-P+o(1),\qquad
B=A\left(\frac u\epsilon+\frac1{\epsilon^2}\right)-E_2+o(1).
$$

由(B3)、$\mathfrak F(\epsilon,\epsilon;1)=[\psi_1(1+\epsilon)-\zeta_2]/2$ 及[polygamma反射](https://dlmf.nist.gov/5.15)，

$$
\Psi(\epsilon,\epsilon)
=\frac A{\epsilon^2}-\frac{u^2}{2}+Pu-E_2+\frac{\pi^2}{2}+o(1).
$$

所以 $\Psi-(B-uD)-u^2/2\to\pi^2/2$，固定正文§5的第一条转换的常数。不是通过某个有限阶拟合补常数。

### B.3 第二条转换的常数是0

记 $a=2\epsilon$、$d=J_a(1)$、$P=\pi\cot(\pi a)$、$E_2=\pi^2\csc^2(\pi a)$、$T_1=T_a(1)$，并设

$$
d_T=-\frac12[(d-P)^2+E_2-\psi_1(1+a)+\zeta_2].
$$

由(B3)，正文§5的第二条转换左侧除以 $\rho^\epsilon$ 后的渐近为

$$
(1+1/z)^a\left(\frac u a+\frac1{a^2}\right)
-\frac{u^2}{2}+d_T-2T_1-\pi^2+o(1).
$$

由(B1)分别对 $b=\epsilon,-\epsilon$ 求导，右侧除以 $\rho^\epsilon$ 后的渐近具有相同幂项和 $-u^2/2$，常数为

$$
\frac{d^2-\psi_1(1+a)-5\zeta_2}{2}-E_2+dP.
$$

在 $b=\epsilon$ 一端，先对(B1)的完整Gauss项消去 $s=0$ 的可去极点并作源导数，再取 $r\to0$；高阶 $r$ 项给出的余项仍趋零。保留到 $O(1)$ 时，所需合并极点的局部式为

$$
\mathrm B(s,1+a-s)=\frac1s+d
+\frac{s}{2}[d^2+\zeta_2+\psi_1(1+a)]+O(s^2).
$$

再乘 $e^{-su}$ 并加入 $-\pi\cot(\pi s)$。这说明上述常数中的 $5\zeta_2$ 不能省略。
用 $T_1=[\psi_1(1+a)-\zeta_2-d^2]/2$ 及 $E_2-P^2=\pi^2$，两常数精确相同。因此正文§5的第二条转换的常数差为0。所有渐近余项在此参数范围及所用参数导数下趋零；随后以原可去积分的解析性延拓。
