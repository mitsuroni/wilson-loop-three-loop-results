# [2]族：从4e-04母积分到六个三圈衍生积分

日期：2026-09-29。本页按“已解母图→相邻端点插入→待求生成函数”分析[2]族，作为当前研究入口。**4e-04的全阶M/N生成已完成；六个6e衍生图已有有限项及收敛级数，但符合当前要求的紧凑自然生成与约化尚未完成。** 本页的源定义与降维不改变这个完成边界。

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
| 6e-10、23 | id | 两个对边点之后；镜像为之前 | $(1-s_2^\epsilon)/\epsilon$ | 待完成 |
| 6e-11、25 | id | 两个对边点之间 | $(s_2^\epsilon-s_1^\epsilon)/\epsilon$ | 待完成 |
| 6e-12、24 | swap | 两个对边点之后；镜像为之前 | $(1-s_2^\epsilon)/\epsilon$ | 待完成 |

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

因此**两个未插入的配对都已求解**。尚缺的是带端点幂的取值，例如 $\mathcal Q_{\rm id}(\epsilon;0,\epsilon)$。仅知道λ1=λ2=0这条截面，不能确定横向的λ导数。对ε求导插入的是
$\log D_1+\log D_2-\log z$；对λ1求导插入的是 $\log s_1$。若要在两者之间转换，必须额外推导带边界的恒等式，不能直接把ε平移或求导当作端点插入。

## 3. 三个待求对象及其完整I

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

这些正是旧收敛级数对应的三个物理插入源；此处只是把它们接回母积分，没有把定义当作求值完成。令

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

这条组合已将双核耦合降为一个显式低核原函数平方的一维矩，适合继续寻找自然生成与矩族。**剩余v积分仍需研究；该组合降维并没有分别解出10与12，也没有完成三种源的紧凑全阶约化。**

## 7. 当前成果与下一步应怎样安排

旧结果已完成七图的ε⁻⁶…ε⁰，三圈均只有简单极点；旧收敛级数也给出全ε解析表示。4e-04后来进一步完成了全阶M/N生成。三圈有限项在此前声明的权重六候选空间之外有额外方向；那个有限阶空间增加一个方向足够，不等于全部高阶只需一个新函数族。

现在应把[2]任务组织为：

1. 复用4e-04与其交换配对的零源解，不重新积分母图。
2. 先从(S2-9)的10/23入手，沿4e-04换元保留(S2-11)的权重差，实际求新项及边界。
3. 同时利用(S2-12)研究10与12的和；若该一维平方矩可闭合，就能把12的任务连接到10。
4. 用(S2-10)单独追踪11/25的反对称插入，再检查是否可由前两条源或其参数导数约化。没有证明前不预设三种源的不可约个数。

通过标准仍是：完整I→明确求值的生成源→固定函数族的任意阶有限规则→旧极点/有限项与独立参数检查。只写Q、F或τ的定义不计完成。

本轮[结构核验](../../../../workflow/checks/same-channel-two-descendants/README.md)重新读取七条原输入，完成42个代数/连接/多项式控制和17个数值控制。新一维组合有纸面域分割与分部积分推导；整数控制和数值检查仅为辅助证据。全阶生成计数仍38/47，其中M/N为36/47；finite框架56/329、formal=not_run均不变。未提交Git。

旧[七图解表](solutions.md)、[通用母族结构](unified-structure.md)、[结构系数](structural-solutions.md)、[端点源研究](../../theory/generation/endpoint-sources.md)保留原字节。它们分别记录求值、结构和历史进度；当前生成研究以本页的完成边界为准。
