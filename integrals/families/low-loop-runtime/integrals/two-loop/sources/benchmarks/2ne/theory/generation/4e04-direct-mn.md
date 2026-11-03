# 4e-04的通解推导：沿同一个Q逐步变形

本文以 $Q$ 为唯一推导主线。正文式(Q-1)至(Q-13)的左边始终是同一个 $Q$：每一步只说明右边发生了什么变化、哪些变量已积掉、还有哪一项未处理。

**辅助函数的定义不等于推导目标改变。** $K$ 的求值和转换恒等式的证明放在附录，并标明服务于正文哪一步；正文使用后立即代回 $Q$。这些Q编号是等式编号，不是Laurent系数 $Q_q$。

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

## 2. 主线：始终追踪同一个Q

下面每个带“Q-”编号的方框，都写出**当前同一个 $Q$ 的完整表达式**。框前说明本步操作，框后说明还剩什么。辅助函数只用于计算右边的某一项；其证明放在附录，不替换推导目标。

### 步骤1：先写出要求解的四重积分

由 $Q=z^{-\epsilon}J_{4e-04}$，

$$
\boxed{\displaystyle
Q=(1+z)^2z^{-\epsilon}
\int_\Omega
D_z(t_1,s_1)^{\epsilon-1}D_z(t_2,s_2)^{\epsilon-1}
\,dt_1dt_2ds_1ds_2.}
\tag{Q-1}
$$

**当前状态：** 四个积分变量都还在，两个排序条件已经包含在 $\Omega$ 中。

### 步骤2：积掉 $t_1$ 和 $s_2$

令 $u=t_2$、$v=s_1$。分别对 $t_1\in(0,u)$、$s_2\in(v,1)$ 使用线性幂原函数，记 $D=D_z(u,v)$，得到

$$
\boxed{\displaystyle
Q=\frac{(1+z)^2z^{-\epsilon}}{\epsilon^2}
\int_0^1du\int_0^1dv\,
\frac{[D^\epsilon-(1-v)^\epsilon][(zu)^\epsilon-D^\epsilon]}
{[(1+z)v-1][(1+z)u-1]}.}
\tag{Q-2}
$$

**当前状态：** 已积掉两个变量，还剩相互独立的 $u,v\in(0,1)$。没有额外的 $u<v$ 条件。两条差商都按连续极限取值，不能先拆成带内部极点的单项。

### 步骤3：只换元，把两个分母变成 $pq$

令

$$
p=(1+z)u-1,\qquad q=(1+z)v-1,
$$

并暂记

$$
A(p)=(1+p)^\epsilon,\quad
B(p)=(1-p/z)^\epsilon,\quad
C(p,q)=(1+pq/z)^\epsilon .
$$

因为 $D=(z+pq)/(1+z)$，式(Q-2)变成

$$
\boxed{\displaystyle
Q=\frac{\rho^\epsilon}{\epsilon^2}
\int_{-1}^z\frac{dp}{p}\int_{-1}^z\frac{dq}{q}\,
[C(p,q)-B(q)][A(p)-C(p,q)].}
\tag{Q-3}
$$

**当前状态：** 仍是二维积分，只是两个变量现在具有相同区间。

### 步骤4：加一个积分为零的反对称项

新旧分子的差为 $[B(p)-B(q)][C(p,q)-1]$。它交换 $p,q$ 后变号，而区域和 $1/(pq)$ 不变；整个差商又可积，所以其积分为零。

因此可以改写为

$$
\boxed{\begin{aligned}
Q=\frac{\rho^\epsilon}{\epsilon^2}
\int_{-1}^z\frac{dp}{p}\int_{-1}^z\frac{dq}{q}\,
\big\{&
[A(p)+B(p)](C-1)\\
&-(C^2-1)-[A(p)-1][B(q)-1]\big\}.
\end{aligned}}
\tag{Q-4}
$$

**当前状态：** 积分值没有改变；右边每个拆项都已可积，可以真正进行下一次积分。

### 步骤5：固定 $p$，实际积掉 $q$

本步使用普通原函数

$$
\mathcal J_a(x)=\int_0^x\frac{(1-t)^a-1}{t}\,dt.
$$

换元 $t=-pq/z$ 给出

$$
\int_{-1}^z\frac{(1+pq/z)^a-1}{q}\,dq
=\mathcal J_a(-p)-\mathcal J_a(p/z).
$$

为缩短表达，记

$$
\begin{gathered}
f(p)=\mathcal J_\epsilon(-p),\qquad
g(p)=\mathcal J_\epsilon(p/z),\\
j_-=\mathcal J_\epsilon(-z),\quad
j_0=\mathcal J_\epsilon(1),\quad
j_+=\mathcal J_\epsilon(-1/z).
\end{gathered}
$$

对(Q-4)的三项分别积 $q$，完整结果是

$$
\boxed{\begin{aligned}
Q=\frac{\rho^\epsilon}{\epsilon^2}\bigg\{&
\int_{-1}^z\frac{dp}{p}
\big[(A(p)+B(p))(f-g)\\
&\hspace{35mm}-\mathcal J_{2\epsilon}(-p)
+\mathcal J_{2\epsilon}(p/z)\big]\\
&-(j_--j_0)(j_0-j_+)\bigg\}.
\end{aligned}}
\tag{Q-5}
$$

**当前状态：** 只剩一个 $p$ 积分。最后一项已经是两个已知端点值的乘积。

### 步骤6：分部积分，分出边界项和唯一的反对称余项

由定义，$f'=(A-1)/p$、$g'=(B-1)/p$，因此

$$
\frac{(A+B)(f-g)}p
=(f'+g')(f-g)+\frac{2(f-g)}p,
$$

而

$$
(f'+g')(f-g)=\frac{d}{dp}\frac{f^2-g^2}{2}+fg'-gf'.
$$

再记

$$
\mathcal T_a(x)=\int_0^x\frac{\mathcal J_a(t)}t\,dt,\qquad
\mathcal C_a=\mathcal T_a(-z)+\mathcal T_a(-1/z)-2\mathcal T_a(1).
$$

将这些原函数代回(Q-5)：

$$
\boxed{\begin{aligned}
Q=\frac{\rho^\epsilon}{\epsilon^2}\bigg[
&\frac12[f^2-g^2]_{p=-1}^{p=z}
+\int_{-1}^z(fg'-gf')\,dp\\
&+2\mathcal C_\epsilon-\mathcal C_{2\epsilon}
-(j_--j_0)(j_0-j_+)\bigg].
\end{aligned}}
\tag{Q-6}
$$

**当前状态：** 唯一尚未求值的是方框中的反对称积分。其余各项都已有端点表达式。

### 步骤7：把这个余项统一到标准有序域

将剩余积分在 $p=0$ 分开，正半区令 $p=zy$，负半区令 $p=-y$。为表示两段，使用辅助积分

$$
K(a,b;z)=\int_{0<x<y<1}
\frac{[(1+zx)^a-1][(1-y)^b-1]
-[(1-x)^b-1][(1+zy)^a-1]}{xy}\,dx\,dy.
$$

于是反对称余项为 $K(\epsilon,\epsilon;z)+K(\epsilon,\epsilon;1/z)$。同时代入(Q-6)的端点值，得到

$$
\boxed{\begin{aligned}
Q=\frac{\rho^\epsilon}{\epsilon^2}\bigg[
&K(\epsilon,\epsilon;z)+K(\epsilon,\epsilon;1/z)\\
&+\frac{(j_-+j_+)^2}{2}-j_0(j_-+j_+)
+2\mathcal C_\epsilon-\mathcal C_{2\epsilon}\bigg].
\end{aligned}}
\tag{Q-7}
$$

**当前状态：** 这里只是统一了剩余积分的形式，还没有把给 $K$ 命名当作完成求值。下一步将代入它的实际求值结果。

### 步骤8：求出 $K$，立即代回 $Q$

辅助求值只做一件事：对 $K$ 求 $z$ 导数并分部积分，把嵌套原函数消去，再用Beta积分求和。详细过程放在附录C；这里使用其结果

$$
\begin{gathered}
\mathfrak F(a,b;x)=
-\frac{ax}{1+b}\,
{}_4F_3\!\left(\begin{matrix}1,1,1,1-a\\2,2,2+b\end{matrix};x\right),\\
K(a,b;z)=2\mathfrak F(a,b;-z)
-2\mathfrak F(a,0;-z)-\mathcal J_b(1)\mathcal J_a(-z).
\end{gathered}
$$

定义反射组合

$$
\Psi(a,b;z)=\mathfrak F(a,b;-z)+\mathfrak F(a,b;-1/z)
-2\mathfrak F(a,b;1).
$$

一圈裸母核记为 $\mathcal V_M$，其端点表示为

$$
\mathcal V_M(\epsilon,z)
=-\frac{\rho^{\epsilon/2}}{\epsilon}(j_-+j_+-2j_0).
$$

代入(Q-7)，并用 $\mathfrak F(\epsilon,\epsilon;1)=\mathcal T_\epsilon(1)+j_0^2/2$ 合并边界项：

$$
\boxed{\displaystyle
Q=\frac12\mathcal V_M(\epsilon,z)^2
+\frac{\rho^\epsilon}{\epsilon^2}
[2\Psi(\epsilon,\epsilon;z)-\Psi(2\epsilon,0;z)].}
\tag{Q-8}
$$

**当前状态：** $Q$ 已完全求值为标准函数。下面三步是在改写这个结果，使它显式产生M/N，而不是重新计算原四重积分。

### 步骤9：只转换第一个 $\Psi$

引入与一圈相同的矩变量：

$$
\begin{gathered}
W=\log[\beta(1-\beta)],\qquad
R=\log\frac{1-\beta}{\beta},\qquad L=\ell+R,\\
Y=\int_0^1\frac{e^{\epsilon W}L}{\beta-\bar\alpha}\,d\beta,\qquad
V=\int_0^1\frac{e^{\epsilon W}-\rho^\epsilon}{\beta-\bar\alpha}\,d\beta.
\end{gathered}
$$

**本步依据：** 比较两侧的 $z$ 导数并解析固定常数，得到

$$
\rho^\epsilon\Psi(\epsilon,\epsilon)
=Y-\ell V-\rho^\epsilon M_{0,1},
\qquad M_{0,1}=-\frac{\ell^2+\pi^2}{2}.
$$

这条恒等式的证明在附录A、B。这里仅使用该结论，不离开Q主线；将它代入(Q-8)：

$$
\boxed{\displaystyle
Q=\frac12\mathcal V_M(\epsilon,z)^2
+\frac1{\epsilon^2}
\left[2(Y-\ell V-\rho^\epsilon M_{0,1})
-\rho^\epsilon\Psi(2\epsilon,0;z)\right].}
\tag{Q-9}
$$

**当前状态：** 第一个 $\Psi$ 已转换；第二个 $\Psi(2\epsilon,0;z)$ 还没有处理。$V$ 是暂时量，最终会消去。

### 步骤10：转换第二个 $\Psi$，消去 $M_{0,1}$ 边界

本步需要

$$
\begin{gathered}
Y_c=\int_0^1
\frac{e^{\epsilon W}L\cosh(\epsilon L)}{\beta-\bar\alpha}\,d\beta,\\
J_-=\mathcal J_{2\epsilon}(-z),\qquad
J_+=\mathcal J_{2\epsilon}(-1/z),\qquad
d=\mathcal J_{2\epsilon}(1).
\end{gathered}
$$

注意大写 $J_\pm$ 的参数为 $2\epsilon$，小写 $j_\pm$ 的参数为 $\epsilon$。

**本步依据：** 第二条已固定常数的转换是

$$
\rho^\epsilon[\Psi(2\epsilon,0)-\ell(J_--J_+)+2M_{0,1}]
=2Y_c+2\epsilon d\,\mathcal V_M(2\epsilon,z).
$$

代入(Q-9)，两份 $M_{0,1}$ 正好抵消：

$$
\boxed{\begin{aligned}
Q=\frac12\mathcal V_M(\epsilon,z)^2+\frac1{\epsilon^2}\big[
&2Y-2Y_c-2\epsilon d\,\mathcal V_M(2\epsilon,z)\\
&-\ell\{2V+\rho^\epsilon(J_--J_+)\}\big].
\end{aligned}}
\tag{Q-10}
$$

**当前状态：** 两个 $\Psi$ 都已消去。还剩 $Y,Y_c$ 以及花括号中的暂时组合。

### 步骤11：消去 $V$ 和两个端点原函数的差

定义同一个矩核中的偶幂差

$$
E=2\int_0^1
\frac{e^{\epsilon W}[\cosh(\epsilon L)-1]}{\beta-\bar\alpha}\,d\beta.
$$

对两个端点幂作线性换元，可得

$$
E=\rho^\epsilon(J_+-J_-)-2V.
$$

所以(Q-10)的花括号等于 $-E$，立即得到

$$
\boxed{\displaystyle
Q=\frac12\mathcal V_M(\epsilon,z)^2
+\frac1{\epsilon^2}
[2Y-2Y_c+\ell E-2\epsilon d\,\mathcal V_M(2\epsilon,z)].}
\tag{Q-11}
$$

**当前状态：** 暂时量 $V,J_-,J_+$ 已全部消去。三个尚需整理的量 $Y,Y_c,E$ 具有同一个积分分母。

### 步骤12：合并三个积分，直接看见N的权重

将 $Y,Y_c,E$ 放进同一个积分号，有

$$
2Y-2Y_c+\ell E
=-2\int_0^1
\frac{e^{\epsilon W}(L-\ell)[\cosh(\epsilon L)-1]}
{\beta-\bar\alpha}\,d\beta.
$$

因此(Q-11)成为

$$
\boxed{\begin{aligned}
Q={}&\frac12\mathcal V_M(\epsilon,z)^2\\
&-\frac2{\epsilon^2}\int_0^1
\frac{e^{\epsilon W}R[\cosh(\epsilon L)-1]}
{\beta-\bar\alpha}\,d\beta
-\frac{2d}{\epsilon}\mathcal V_M(2\epsilon,z).
\end{aligned}}
\tag{Q-12}
$$

**当前状态：** $R=L-\ell$ 已直接出现在被积函数里。这个式子就是N的来源；还没有借助任何Hlog拟合。

### 步骤13：给已经得到的两种生成函数命名

一圈部分与新出现的部分分别是

$$
\begin{aligned}
\mathcal V_M(\epsilon,z)
&=\int_0^1\frac{e^{\epsilon W/2}}{\beta-\bar\alpha}
\frac{\sinh(\epsilon L/2)}{\epsilon/2}\,d\beta,\\
\mathcal V_N(\epsilon,z)
&=\int_0^1\frac{e^{\epsilon W}R}{\beta-\bar\alpha}
\frac{\cosh(\epsilon L)-1}{\epsilon^2}\,d\beta .
\end{aligned}
$$

两者内部点均可去，且在 $\epsilon=0$ 正则。一圈生成函数的直接降维见附录D。

端点值已知为

$$
d=\int_0^1\frac{(1-t)^{2\epsilon}-1}{t}\,dt
=-\gamma_E-\psi(1+2\epsilon)=-\chi(\epsilon).
$$

于是(Q-12)最终写成

$$
\boxed{\displaystyle
Q=\frac12\mathcal V_M(\epsilon,z)^2
-2\mathcal V_N(\epsilon,z)
+\frac{2\chi(\epsilon)}{\epsilon}\mathcal V_M(2\epsilon,z).}
\tag{Q-13}
$$

**当前状态：** 这就是同一个 $Q$ 的M/N生成函数形式。从(Q-1)到(Q-13)，原积分没有更换；改变的是右边的表示。

## 3. 最后才从Q恢复物理I

将完整的(Q-13)代入开头的归一化式(I)，得到

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

这里 $\chi(\epsilon)=\gamma_E+\psi(1+2\epsilon)$，$\psi=\Gamma'/\Gamma$。它来自端点原函数，并满足

$$
\frac{\chi(\epsilon)}{\epsilon}
=2\zeta_2-4\zeta_3\epsilon+8\zeta_4\epsilon^2-16\zeta_5\epsilon^3+\cdots.
$$

因此 $\chi/\epsilon$ 没有极点。$\mathcal V_M$ 和 $\mathcal V_N$ 也在零点正则。$\mathcal V_M(2\epsilon,z)$ 只表示一圈裸母核的参数替换，不额外携带Gamma前因子。

公式中的指数在零附近展开时等于 $Z^2$，只计一次项目规定的Euler补偿。M/N生成式和标准函数式是同一 $Q$ 的两种表示，不是两次不同的近似求解。

## 4. 展开生成函数：M、N的系数逐项出现

本节使用的矩为

$$
M_{k,j}=\int_0^1\frac{W^kL^j}{\beta-\bar\alpha}\,d\beta,\qquad
N_{k,j}=\int_0^1\frac{W^kRL^j}{\beta-\bar\alpha}\,d\beta,
\quad k\ge0,\ j\ge1.
$$

由(VM)，

$$
\begin{aligned}
\mathcal V_M(\epsilon,z)
={}&M_{0,1}+\frac{\epsilon}{2}M_{1,1}\\
&+\epsilon^2\left(\frac{M_{2,1}}8+\frac{M_{0,3}}{24}\right)
+O(\epsilon^3).
\end{aligned}
$$

由(Q-13)前的 $\mathcal V_N$ 定义，

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

## 5. 通项公式与前两阶检验

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

### 5.1 零阶

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

### 5.2 一阶

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

## 6. 版本与验证

本稿把已有解析推导重新组织为连续的Q等式链，数学结论不变。13个源积分点和15个 $Q,\Delta,I$ 系数回归保留在[原计算记录](../../../../workflow/checks/4e04-direct-mn/README.md)，本次等式衔接与渲染检查见[Q主线检查](../../../../workflow/checks/4e04-q-chain/README.md)。

[r001研究稿](4e04-direct-mn-r001.md)、[r002讲解稿](4e04-direct-mn-r002.md)、[r003分步说明稿](4e04-direct-mn-r003.md)保留各自原字节。旧报告按相应版本解释，不将它们改写为当前文本的验证报告。

结论仍限于正实 $z$ 下4e-04及其裸交换配对；没有新增展开阶数、全族最小性或形式证明。

## 附录A：正文步骤9—11使用的转换恒等式

本附录只证明(Q-9)至(Q-11)代入时使用的恒等式。阅读后返回相应Q等式；本附录不另定义一个Q。为与既有计算逐项对应，以下临时记 $H=\mathcal V_M$、$r=\bar\alpha$、$d=\mathcal J_{2\epsilon}(1)$；其余 $Y,Y_c,V,E$ 均沿用正文定义。附录中的 $J_a,T_a$ 分别指正文的 $\mathcal J_a,\mathcal T_a$。

本附录使用的双源是

$$
\mathscr U(a,b;z)=\int_0^1
\frac{e^{aW}(e^{bL}-1)}{\beta-\bar\alpha}\,d\beta.
$$

三条转换合并后的关系为

$$
\boxed{\displaystyle
\rho^\epsilon[2\Psi(\epsilon,\epsilon)-\Psi(2\epsilon,0)]
=2Y-2Y_c+\ell E
-2\epsilon d\,\mathcal V_M(2\epsilon,z).}
\tag{BR}
$$

下面证明其中的转换。定义

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

只差与 $z$ 无关的常数。附录B证明该常数为 $\pi^2/2$，给出正文步骤9的转换。

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

附录B固定其常数差为零，得到正文步骤10的转换。


### A.2 正文步骤11使用的换元

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


此处的结果代回(Q-10)，便得到(Q-11)。

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

所以 $\Psi-(B-uD)-u^2/2\to\pi^2/2$，固定正文步骤9的转换的常数。不是通过某个有限阶拟合补常数。

### B.3 第二条转换的常数是0

记 $a=2\epsilon$、$d=J_a(1)$、$P=\pi\cot(\pi a)$、$E_2=\pi^2\csc^2(\pi a)$、$T_1=T_a(1)$，并设

$$
d_T=-\frac12[(d-P)^2+E_2-\psi_1(1+a)+\zeta_2].
$$

由(B3)，正文步骤10的转换左侧除以 $\rho^\epsilon$ 后的渐近为

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
用 $T_1=[\psi_1(1+a)-\zeta_2-d^2]/2$ 及 $E_2-P^2=\pi^2$，两常数精确相同。因此正文步骤10的转换的常数差为0。所有渐近余项在此参数范围及所用参数导数下趋零；随后以原可去积分的解析性延拓。
## 附录C：只求K，服务于正文步骤8

正文在(Q-7)处仍保留 $K$。这里独立求出它，然后把结果代回得到(Q-8)。下面的等号可以以辅助函数开头，但推导目标始终是(Q-7)右边的两项 $K$。

### C.1 写出需要求值的排序积分

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

### C.2 为什么对 $z$ 求导：它恰好消去 $f$ 的原函数

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

### C.3 用 Beta 积分求出 $P$，再积回 $K$

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

## 附录D：正文使用的线性幂原函数

D.1说明正文所用的已知一圈子结构；D.2核对(Q-1)到(Q-2)实际积掉的两个变量。阅读后返回原来的Q步骤。

### D.1 一圈M生成函数

这一步既与2e-02对照，也定义两圈通解中将出现的一圈子结构。

#### D.1.1 积去一个线性变量

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

#### D.1.2 用双曲函数整理幂差

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

### D.2 正文步骤2的两个线性原函数

令 $u=t_2$、$v=s_1$。先积 $t_1$，再积 $s_2$；两步都使用附录D.1.1的线性幂原函数：

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

