# [2]族：4e-04与六个三圈衍生积分的自然全阶生成

日期：2026-09-29。本页按“已解母图→相邻端点插入→待求生成函数”分析[2]族，作为当前研究入口。**本族7/7已达到自然全阶生成目标。** 4e-04全阶使用M/N；六个6e衍生图现由明确的排序Beta源、带边界三角系统、一维平方源和已解HD生成，每阶有有限词规则。尚未得到这些新源与M/N同样简洁的少指标矩基底，也没有证明最小性。先前分析版保存在r001。

阅读时可先看第7节的完整I公式(G2-2)、(G2-14)、(G2-18)；第1–6节说明它们从原积分怎样产生。

本族共七图，三圈按完整I分为三对：10=23、11=25、12=24。两圈母图和所有衍生图保留相同的有序双核结构；12/24使用交换配对，因此“都是4e-04衍生”应理解为同一有序双核体系中的插入，不能把配对区别省掉。

## 1. 七图的结构和当前状态

固定

$$
D_z(t,s)=(1-t)(1-s)+zts,\qquad
\Omega_2=\{0<t_1<t_2<1,\ 0<s_1<s_2<1\},
$$

并区分同序配对id与交换配对swap：

$$
K_{\rm id}=D_z(t_1,s_1)^{\epsilon-1}D_z(t_2,s_2)^{\epsilon-1},\qquad
K_{\rm sw}=D_z(t_1,s_2)^{\epsilon-1}D_z(t_2,s_1)^{\epsilon-1}.
$$

先完成新增相邻传播子的两个端点积分。其中一个端点不在双核支撑边上，给出1/ε；另一个落在带有两个对边端点的边上，其位置决定以下权重。

| 图 | 代表图配对 | 支撑边上的相邻端点位置 | 积分后留下的端点权重 | 当前生成状态 |
|---|---|---|---|---|
| 4e-04 | id | 没有新增相邻端点 | 1 | 全阶M/N已完成 |
| 6e-10、23 | id | 两个对边点之后；镜像为之前 | $(1-s_2^\epsilon)/\epsilon$ | 本轮生成完成；少指标矩约化开放 |
| 6e-11、25 | id | 两个对边点之间 | $(s_2^\epsilon-s_1^\epsilon)/\epsilon$ | 本轮生成完成；少指标矩约化开放 |
| 6e-12、24 | swap | 两个对边点之后；镜像为之前 | $(1-s_2^\epsilon)/\epsilon$ | 本轮生成完成；少指标矩约化开放 |

例如三种位置来自

$$
\int_{s_2}^1a^{\epsilon-1}da=\frac{1-s_2^\epsilon}{\epsilon},\qquad
\int_{s_1}^{s_2}a^{\epsilon-1}da=\frac{s_2^\epsilon-s_1^\epsilon}{\epsilon}.
\tag{S2-1}
$$

之前位置及B型端点通过反射处理，必须同时反射配对端点并保留z的尺度幂。实际原输入的边字、颜色与镜像已逐图核对。

## 2. 4e-04给了什么，尚未给什么

沿用共同归一化

$$
g_\epsilon=-(1+z)z^{-\epsilon/2},\quad
Z=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),\quad
c_n=\cosh(n\epsilon\log z/2).
$$

将端点幂暂时作为独立源参数：

$$
\mathcal Q_\pi(\epsilon;\lambda_1,\lambda_2;z)
=g_\epsilon^2\int_{\Omega_2}K_\pi
s_1^{\lambda_1}s_2^{\lambda_2}\,d^2t\,d^2s.
\tag{S2-2}
$$

这里ε改变传播子的幂，λ1、λ2分别改变较小、较大右端点的权重；它们是不同方向。

4e-04已经求完的是零源截面：

$$
\begin{aligned}
I_{4e-04}&=-\frac{Z^2}{8c_2}\mathcal Q_{\rm id}(\epsilon;0,0;z),\\
\mathcal Q_{\rm id}(\epsilon;0,0;z)
&=\frac12\mathcal V_M(\epsilon,z)^2
-2\mathcal V_N(\epsilon,z)
+\frac{2\chi(\epsilon)}{\epsilon}\mathcal V_M(2\epsilon,z).
\end{aligned}
\tag{S2-3}
$$

这里使用[4e-04主稿](../../theory/generation/4e04-direct-mn.md)的单参数VM/VN约定；$\chi=\gamma_E+\psi(1+2\epsilon)$，各可去ε分母取整体极限。区域分割还给出

$$
\mathcal Q_{\rm sw}(\epsilon;0,0;z)
=\frac12\mathcal V_M(\epsilon,z)^2-\mathcal Q_{\rm id}(\epsilon;0,0;z).
\tag{S2-4}
$$

因此**两个未插入的配对都已求解**。衍生图额外需要带端点幂的取值，例如 $\mathcal Q_{\rm id}(\epsilon;0,\epsilon)$。仅知道λ1=λ2=0这条截面，不能确定横向的λ导数。对ε求导插入的是
$\log D_1+\log D_2-\log z$；对λ1求导插入的是 $\log s_1$。若要在两者之间转换，必须额外推导带边界的恒等式，不能直接把ε平移或求导当作端点插入。

## 3. 三个插入源及其完整I

把(S2-1)中的可去ε分母保留在权重内，定义三个正则源：

$$
\begin{aligned}
\mathcal F_A&=\frac{\mathcal Q_{\rm id}(\epsilon;0,0)
-\mathcal Q_{\rm id}(\epsilon;0,\epsilon)}{\epsilon},\\
\mathcal F_B&=\frac{\mathcal Q_{\rm id}(\epsilon;0,\epsilon)
-\mathcal Q_{\rm id}(\epsilon;\epsilon,0)}{\epsilon},\\
\mathcal F_C&=\frac{\mathcal Q_{\rm sw}(\epsilon;0,0)
-\mathcal Q_{\rm sw}(\epsilon;0,\epsilon)}{\epsilon}.
\end{aligned}
\tag{S2-5}
$$

这些正是旧收敛级数对应的三个物理插入源；此处先接回母积分，第7节再给实际求值及任意阶规则，不能把本节定义单独当作完成。令

$$
\langle F\rangle_\epsilon=
\frac{e^{\epsilon\log z/2}F(\epsilon,z)
+e^{-\epsilon\log z/2}F(\epsilon,1/z)}2.
$$

当前唯一C/I归一化给出

$$
\boxed{\begin{aligned}
I_{6e-10}=I_{6e-23}&=-\frac{Z^3}{8\epsilon c_3}\langle\mathcal F_A\rangle_\epsilon,\\
I_{6e-11}=I_{6e-25}&=-\frac{Z^3}{16\epsilon c_3}\langle\mathcal F_B\rangle_\epsilon,\\
I_{6e-12}=I_{6e-24}&=-\frac{Z^3}{16\epsilon c_3}\langle\mathcal F_C\rangle_\epsilon.
\end{aligned}}
\tag{S2-6}
$$

三种插入在ε=0分别趋于 $-\log s_2$、$\log(s_2/s_1)$、$-\log s_2$。它们在原有序域可积且非负，乘双核后积分非零。因此三圈的实际极点是简单极点，来源是未插入支撑边的另一个相邻端点；不是仍有两个未抵消的ε分母。

## 4. 沿4e-04路线降维：三种插入都可保留初等原函数

先固定 $0<\epsilon<1$、z>0。利用原排序域的可积上界作所需换序。定义已求值的线性幂原函数

$$
P_\epsilon(x;y)=\int_0^xD_z(r,y)^{\epsilon-1}dr
=\frac{D_z(x,y)^\epsilon-(1-y)^\epsilon}
{\epsilon[(1+z)y-1]},\qquad T_\epsilon(y)=P_\epsilon(1;y).
\tag{S2-7}
$$

差商内部的零分母按连续极限处理，不能拆成孤立的极点项。令u、v都在(0,1)，写出两个明确的二维核：

$$
\begin{aligned}
\mathscr K_{\rm id}(u,v)
&=g_\epsilon^2P_\epsilon(v;u)
[T_\epsilon(v)-P_\epsilon(u;v)],\\
\mathscr K_{\rm sw}(u,v)
&=g_\epsilon^2P_\epsilon(v;u)P_\epsilon(u;v).
\end{aligned}
\tag{S2-8}
$$

同序的较大端点插入保留u=t1、v=s2，积掉s1和t2；较小端点插入则先保留t2、s1，再交换两个哑变量。由核的t/s对称性，两种选择可以统一到同一个Ki。因此

$$
\boxed{\begin{aligned}
\mathcal F_A&=\int_0^1\!du\int_0^1\!dv\,
\frac{1-v^\epsilon}{\epsilon}\mathscr K_{\rm id}(u,v),\\
\mathcal F_B&=\int_0^1\!du\int_0^1\!dv\,
\frac{v^\epsilon-u^\epsilon}{\epsilon}\mathscr K_{\rm id}(u,v),\\
\mathcal F_C&=\int_0^1\!du\int_0^1\!dv\,
\frac{1-v^\epsilon}{\epsilon}\mathscr K_{\rm sw}(u,v).
\end{aligned}}
\tag{S2-9}
$$

将这三行分别代回(S2-6)，就是三个完整I各自的二维起式。这里没有预先积分带源幂的端点，所以暂时不需要带额外参数的Gauss原函数。

**变量角色的边界：** 第二行是两个分别换元后的积分之差；不能将一般双源Q直接写成Ki乘 $u^{\lambda_1}v^{\lambda_2}$。中间权重在这个正方形上可以变号，其积分与原有序域中非负的中间插入相等。

交换核Ks关于u、v对称；中间插入则只取Ki的反对称部分：

$$
\mathcal F_B=\frac12\int_0^1\!du\int_0^1\!dv\,
\frac{v^\epsilon-u^\epsilon}{\epsilon}
[\mathscr K_{\rm id}(u,v)-\mathscr K_{\rm id}(v,u)].
\tag{S2-10}
$$

因此三种位置不是任意增加三个符号，而是同一个双核的边界差、反对称插入和交换配对插入。

## 5. 为什么4e-04中的消项不能直接继承

4e-04主稿(Q-3)→(Q-4)用到了一个真正积分为零的反对称项。在其变量 $p,q\in[-1,z]$ 下，令

$$
B(p)=(1-p/z)^\epsilon,\qquad
C(p,q)=(1+pq/z)^\epsilon,\qquad
F(p,q)=\frac{[B(p)-B(q)][C(p,q)-1]}{pq}.
$$

F交换p、q后变号，所以无权积分为零。但加入单端权重w(q)后，正确公式变为

$$
\boxed{\displaystyle
\int\!dp\,dq\,w(q)F(p,q)
=\frac12\int\!dp\,dq\,[w(q)-w(p)]F(p,q).}
\tag{S2-11}
$$

右端是“权重差×反对称核”，一般不为零。一个精确控制例子是裸源ε=1、z=1、$w(q)=(1+q)/2$：此时 $F=q-p$，无权积分为0，而带权积分为2/3。这个例子只检查换元后积分的代数结构，不在Gamma极点上评估物理I。

这定位了推广中的具体障碍：要沿4e-04路线继续，必须计算(S2-11)产生的额外项，或者用新的恒等式约化它。它不构成一般M/N不可约性的证明，也不要求增加运动学字母。

## 6. 一个已经得到的一维组合：低核原函数的平方

(S2-8)有精确抵消

$$
\mathscr K_{\rm id}+\mathscr K_{\rm sw}
=g_\epsilon^2P_\epsilon(v;u)T_\epsilon(v).
$$

复用已经求值的一圈矩形边界原函数

$$
f_\epsilon(v;z)=\int_0^1dt\int_0^vds\,k_\epsilon(t,s;z),
\qquad k_\epsilon=g_\epsilon D_z^{\epsilon-1}.
$$

它满足 $f_\epsilon'=g_\epsilon T_\epsilon$、$f_\epsilon(0)=0$。先对u积分，再分部积分，得到

$$
\begin{aligned}
\mathcal F_A+\mathcal F_C
&=\int_0^1\frac{1-v^\epsilon}{\epsilon}f_\epsilon'(v)f_\epsilon(v)dv\\
&=\boxed{\frac12\int_0^1v^{\epsilon-1}f_\epsilon(v;z)^2dv}.
\end{aligned}
\tag{S2-12}
$$

两端分部积分边界均为零：v=1处端点差为0，v=0处f为0，且在初始收敛域足以控制该乘积。相应完整物理组合为

$$
\boxed{\displaystyle
I_{6e-10}+2I_{6e-12}
=-\frac{Z^3}{16\epsilon c_3}
\left\langle\int_0^1v^{\epsilon-1}f_\epsilon(v;z)^2dv\right\rangle_\epsilon.}
\tag{S2-13}
$$

这里f不是新待求母函数。由[已有矩形源闭式](../../theory/generation/three-kernel-generators.md)，令 $q=(1+z)v-1$、$\rho=z/(1+z)^2$，则

$$
\begin{gathered}
J_\epsilon(y)=-\epsilon y\,{}_3F_2(1,1,1-\epsilon;2,2;y),\\
f_\epsilon(v;z)=-\frac{\rho^{\epsilon/2}}{\epsilon}
[J_\epsilon(-q)-J_\epsilon(1)-J_\epsilon(q/z)+J_\epsilon(-1/z)].
\end{gathered}
\tag{S2-14}
$$

$J_\epsilon(1)=-\gamma_E-\psi(1+\epsilon)$。所有自变量不超过1，取与原实积分一致的分支。ε=0按整体可去极限理解。

这条组合已将双核耦合降为一个显式低核原函数平方的一维矩，适合继续寻找自然生成与矩族。**这一降维本身不等于完成；第7.5节进一步把它转成已指定GPL的有限词规则，并接回12/24的完整I。**

## 7. 从二维起式得到六图的自然全阶生成

本节是2026-09-29的实际求解结果，按10/23→12/24→11/25推进。新增源均有明确边界和任意阶有限规则；不是把第3节的Q系数重新命名。采用的固定函数语言是Beta生成源及0、−1迭代积分，暂不宣称它们具有与M/N一样简短的少指标矩基底。

### 7.1 保留带权项，分离出真正新增的源

在第4节的二维式中令 $p=(1+z)u-1$、$q=(1+z)v-1$，再记

$$
A_p=(1+p)^\epsilon,\quad B_p=(1-p/z)^\epsilon,\quad
C_{pq}=(1+pq/z)^\epsilon,\quad \rho=\frac{z}{(1+z)^2}.
$$

同序二维核乘测度变为
$\rho^\epsilon(C_{pq}-B_p)(A_q-C_{pq})\,dp\,dq/(\epsilon^2pq)$。
额外的较大端点幂是 $v^\epsilon=(1+z)^{-\epsilon}A_q$。直接展开并保留全部可去差，得到

$$
\begin{gathered}
a_\lambda(z)=J_\lambda(-z)-J_\lambda(1),\qquad
b_\lambda(z)=J_\lambda(1)-J_\lambda(-1/z),\\
\mathcal R_\epsilon(z)=\mathcal L(2\epsilon,\epsilon;z)
-\mathcal L(\epsilon,2\epsilon;z)+\mathcal S_\epsilon(z)
-b_\epsilon(z)[a_{2\epsilon}(z)-a_\epsilon(z)].
\end{gathered}
\tag{G2-1}
$$

于是第一个完整I已改写为

$$
\boxed{\displaystyle
I_{6e-10}=I_{6e-23}
=-\frac{Z^3}{8\epsilon^2c_3}
\left\langle Q-
\frac{\rho^\epsilon}{\epsilon^2(1+z)^\epsilon}\mathcal R_\epsilon
\right\rangle_\epsilon.}
\tag{G2-2}
$$

Q是第2节已经求解的4e-04母核。L只有一侧带权，可以用已知排序Beta源求值；S才是双侧同时带权产生的新对象。下面分别给出它们的实际解法。

### 7.2 单侧带权部分L已经可以显式求值

沿用已知

$$
\begin{gathered}
J_a(x)=-ax\,{}_3F_2(1,1,1-a;2,2;x),\\
\mathfrak F(a,b;x)=-\frac{ax}{1+b}
{}_4F_3(1,1,1,1-a;2,2,2+b;x),\\
\Psi(a,b;z)=\mathfrak F(a,b;-z)+\mathfrak F(a,b;-1/z)-2\mathfrak F(a,b;1),\\
K(a,b;z)=2\mathfrak F(a,b;-z)-2\mathfrak F(a,0;-z)-J_b(1)J_a(-z).
\end{gathered}
\tag{G2-3}
$$

还需一个普通的一变量排序Beta源：

$$
\mathfrak B(a,b;x)=\int_0^x\frac{(1-t)^a-1}{t}J_b(t)\,dt.
\tag{G2-4}
$$

它有明确的任意阶已积分表达。使用GPL约定
$G(c_1,\ldots,c_n;x)=\int_0^xdt\,G(c_2,\ldots;c_n;t)/(t-c_1)$，空字为1，则

$$
\boxed{\displaystyle
[a^m b^n]\mathfrak B(a,b;x)
=\sum_{w\in (1^m)\mathbin{\amalg}(0,1^n)}G(0,w;x),
\qquad m,n\ge1.}
\tag{G2-5}
$$

某个指标为0时为零，shuffle保留重数。这里没有尚待求值的系数积分；每阶是明确的有限词和，在x=1也取通常收敛端点值。它来自两个幂的直接展开与一次积分，没有读旧Hlog作拟合。

现在单侧带权源的完整求值是

$$
\boxed{\begin{aligned}
\mathcal L(b,c;z)={}&\Psi(c,0;z)+\mathfrak B(b,c;-z)-\mathfrak B(b,c;1)\\
&-\frac12\left[J_b(-z)J_c(1)-J_c(-1/z)J_b(1)
-K(b,c;z)-K(c,b;1/z)\right].
\end{aligned}}
\tag{G2-6}
$$

证明先对p积分，得到
$\int_{-1}^z (1+q)^b[J_c(-q)-J_c(q/z)]\,dq/q$；再将q分为正负两段。相同自变量的部分给B，不同自变量的两项用已解交换子K和乘积shuffle拆开，即得(G2-6)。因此L没有新的未定边界常数。

### 7.3 新双端源S：由边界导出的闭合方程

从真实剩余项定义

$$
\mathcal S_\epsilon(z)=\int_{-1}^z\frac{dp}{p}\int_{-1}^z\frac{dq}{q}
(1-p/z)^\epsilon(1+q)^\epsilon
[(1+pq/z)^\epsilon-1].
\tag{G2-7}
$$

它满足 $\mathcal S_\epsilon(1/z)=\mathcal S_\epsilon(z)$。做 $x=p/z,y=-q$ 后，积分内部只含
$(1-x)^\epsilon(1-y)^\epsilon[(1-xy)^\epsilon-1]$，z只进入互为倒数的两个负端点。故z微分只留下两条边界，不会留下新的二维积分。

令

$$
A(z)=(1+z)^\epsilon,\quad B(z)=(1+1/z)^\epsilon,\quad
R_2=\frac{\Gamma(1+\epsilon)^2}{\Gamma(1+2\epsilon)},\quad\theta=z\partial_z.
$$

对边界Beta/Cauchy源使用与4e-04相同的分部积分，得到三条封闭方程：

$$
\boxed{\begin{aligned}
\theta E_\epsilon&=B(R_2A-1),\\
\theta F_\epsilon&=A(1-R_2B),\\
\theta\mathcal S_\epsilon&=A E_\epsilon-B F_\epsilon.
\end{aligned}}
\tag{G2-8}
$$

E、F是这里的辅助Beta原函数，不是旧有限项表中的E基底。它们可先独立求出，然后右端完全已知。其边界及S的边界为

$$
\boxed{E_\epsilon(-1)=F_\epsilon(-1)=\mathcal S_\epsilon(-1)=0.}
\tag{G2-9}
$$

这里的−1是辅助解析边界：从z=−1沿上半单位圆到z=1，再沿正实轴到目标z，固定相应对数分支。最终z>0的值为实数；这不是宣称已完成物理复z/i0延拓。

边界不是拟合的。对S将p、q各自缩放到单位区间，$z\to-1$时它是
$O((1+z)^{2+2\epsilon})$，所以为零。对E、F可用已知Cauchy源D核对：

$$
E=J_\epsilon(-1/z)-J_\epsilon(1)-D,\qquad
F=J_\epsilon(-z)-J_\epsilon(1)+D,
$$
$$
\theta D=1-R_2\rho^{-\epsilon},\qquad D(1)=0.
$$

在单位圆上 $\rho^{-1}=4\cos^2(\phi/2)$，Gamma倍增公式给出
$R_2\int_0^\pi[4\cos^2(\phi/2)]^\epsilon d\phi=\pi$，因而D(−1)=0。两个J同时趋于J(1)，E、F的边界随之固定。

这便完成S的求值问题：三个显式原函数、固定的零边界，没有待定函数或每阶拟合常数。还有 $E(1/z)=F(z)$。

### 7.4 S的任意阶规则与字母表

记 $\mathcal I f=\int_{-1}^z f(t)dt/t$，采用上述同一辅助路径；A、B、R2的Taylor系数都是明确的对数/Gamma多项式。逐阶使用

$$
\begin{aligned}
E_q&=\mathcal I[\epsilon^q]\{B(R_2A-1)\},\\
F_q&=\mathcal I[\epsilon^q]\{A(1-R_2B)\},\\
S_q&=\mathcal I\sum_{i+j=q}(A_iE_j-B_iF_j).
\end{aligned}
\tag{G2-10}
$$

每阶只有有限个对数乘积及一次已知原函数操作。积分用shuffle后补字首0完成，再减去−1处的已知GPL端点值；收敛端点值属于0/1的MZV，尾零产生的log(−1)=iπ要保留至相消，不能擅自丢弃。

也可取
$\mathbf v=(\epsilon^2S,\epsilon AE,\epsilon BF,R_2A^2B,R_2AB^2,AB)^T$，得到

$$
d\mathbf v=\epsilon\left(A_0\frac{dz}{z}+A_{-1}\frac{dz}{1+z}\right)\mathbf v,
$$
$$
A_0=\begin{pmatrix}
0&1&-1&0&0&0\\0&0&0&1&0&-1\\0&0&-1&0&-1&1\\
0&0&0&-1&0&0\\0&0&0&0&-2&0\\0&0&0&0&0&-1
\end{pmatrix},\qquad A_{-1}=\operatorname{diag}(0,1,1,3,3,2).
\tag{G2-11}
$$

后面三个分量是已写出的幂函数；前面三个的边界由(G2-9)给出，不能在奇点处把六个分量都设成普通零初值。这个上三角系统解释了为何新源仍只需要0、−1字母，但不证明M/N闭包或最小性。

### 7.5 12/24：把一维平方源积分成固定的有限词规则

令 $q=(1+z)v-1$，把第6节已求值的f代回平方积分。定义

$$
\begin{aligned}
h_\epsilon(q;z)&=J_\epsilon(-q)-J_\epsilon(q/z)
-J_\epsilon(1)+J_\epsilon(-1/z),\\
\mathcal C_\epsilon(z)&=\int_{-1}^z\frac{dq}{1+q}(1+q)^\epsilon h_\epsilon(q;z)^2.
\end{aligned}
\tag{G2-12}
$$

f已经是明确的3F2，所以这里没有未知内层原函数。更重要的是，C的任意阶也可以直接积完。令

$$
h_n(q;z)=G(0,(-1)^n;q)-G(0,z^n;q)
-(-1)^n\zeta_{n+1}+G(0,(-1)^n;1/z),\qquad n\ge1.
$$

括号中的 $(-1)^n$、$z^n$ 在G的字槽里表示重复n个字母，外面的 $(-1)^n\zeta$ 才是普通幂。则

$$
[\epsilon^N]\mathcal C_\epsilon
=\sum_{\substack{k\ge0,\ m,n\ge1\\k+m+n=N}}
\left[\operatorname{Prim}_{-1}
\big(G((-1)^k;q)h_m(q;z)h_n(q;z)\big)\right]_{q=-1}^{q=z}.
\tag{G2-13}
$$

Prim的操作完全明确：先作有限shuffle，再给每个词补一个字首−1；它是除以1+q的原函数。边界必须对完整组合取极限，因为 $h_n(-1;z)=0$ 保证发散对数相消。不能分别扔掉不收敛的单词。由此每个系数都是已指定GPL的有限和，没有剩余系数积分；字母碰撞只有z=0、−1，可转换到相同的外部字母表。

第二对完整I为

$$
\boxed{\displaystyle
I_{6e-12}=I_{6e-24}
=-\frac12I_{6e-10}
-\frac{Z^3\rho^{3\epsilon/2}}{64\epsilon^3c_3}
[\mathcal C_\epsilon(z)+\mathcal C_\epsilon(1/z)].}
\tag{G2-14}
$$

C从ε²开始，因此显式ε⁻³只产生实际简单极点。第一对已经由(G2-2)求出，(G2-13)给第二对的任意阶，不需要再求一个加权四维母积分。

### 7.6 11/25：困难插入连接已解的双端单核HD

中间插入展开后只多出同侧权重源

$$
\mathcal T_\epsilon(z)=\int_{-1}^z\frac{dp}{p}\int_{-1}^z\frac{dq}{q}
(1+p)^\epsilon(1+q)^\epsilon[(1+pq/z)^\epsilon-1].
\tag{G2-15}
$$

对z求导，两个移动上边界给第一项；内部导数消掉pq分母，变回已解的HD：

$$
\boxed{\displaystyle
\theta\mathcal T_\epsilon
=2(1+z)^\epsilon[a_{2\epsilon}(z)-a_\epsilon(z)]
+\epsilon(1+z)^{3\epsilon}z^{-\epsilon/2}H_D(\epsilon,z),
\qquad\mathcal T_\epsilon(-1)=0.}
\tag{G2-16}
$$

这里
$H_D=-(1+z)z^{-\epsilon/2}\int_0^1du\,dv\,u^\epsilon v^\epsilon D_z(u,v)^{\epsilon-1}$，
正是[[1]族已经完成的双端源](../single-opposite/three-loop.md#all-orders)。其三分量常矩阵与Gamma边界给任意阶有限词规则，因而在这里是已解输入。T的零边界由同样的单位区间缩放得到，不增加新的常数；它的任意阶是右端已知系数再补一次字首0，并减去−1端点值。

余下纯一维权重也不需要新源：

$$
\int_{-1}^z\frac{dp}{p}(1+p)^\epsilon[(1-p/z)^\epsilon-1]=-F_\epsilon(z).
\tag{G2-17}
$$

将全部项代回中间插入，并使用E↔F的镜像关系，第三对完整I为

$$
\boxed{\begin{aligned}
I_{6e-11}=I_{6e-25}={}&-\frac12I_{6e-10}-I_{6e-12}
-\frac{Z^3c_1\mathcal V_M(\epsilon,z)^2}{32\epsilon^2c_3}\\
&-\frac{Z^3\rho^{3\epsilon/2}}{32\epsilon^4c_3}
\left[2\mathcal S_\epsilon-\mathcal T_\epsilon(z)-\mathcal T_\epsilon(1/z)
-(E_\epsilon+F_\epsilon)^2\right].
\end{aligned}}
\tag{G2-18}
$$

方括号及前项必须整体展开；原端点差证明实际只含简单极点。新的生成过程也逐项核对了所有虚假的更高负幂抵消。

## 8. 此轮完成范围、核验与下一步

至此，4e-04母图和六个衍生图均有来源明确的全阶生成式：母图使用M/N，衍生图使用已解排序Beta源、三角系统S、明确平方源C以及由已有HD驱动的T。每一个源的边界与任意阶有限取系数规则均已给出。**[2]族按自然全阶生成标准为7/7完成；六个三圈图并未被宣称全阶回到M/N。**

本轮[检查记录](../../../../workflow/checks/same-channel-two-first-generator/README.md)包括42个独立代数/整数控制；源S至ε⁴由上述方程生成后才读取旧目标；六图42个ε⁻⁶…ε⁰位置全部精确恢复。新鲜HyperInt只执行已有显式源的词变换、收敛周期和一维平方积分，发散检查保持开启；没有拟合新函数。84个非整数参数、源导数、镜像与完整I数值控制通过。

原输入、旧结果与旧报告不变；上一版分析主稿保存在[r001](three-loop-r001.md)。全项目固定函数族生成由38/47推进至44/47，6e由30/39到36/39；全阶M/N仍36/47。剩余是[3]的20/21/22。新框架有限阶验收账本仍56/329，formal=not_run，未提交Git。

后续可以从这些自然Beta源研究简洁矩族及相互约化，或转向最后的[3]三图。当前没有证明最小函数族个数、整个高圈的M/N闭包或一般物理复z/i0延拓。
