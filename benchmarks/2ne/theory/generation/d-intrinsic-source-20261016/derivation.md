# 从端点几何平均得到不显含 ρ 的三圈生成表示

2026-10-16。这里继续研究裸几何积分 D=R231=R312，并同时保留 C=R321。目标是把最终展开写成带指标函数的常系数组合，而不是在每阶外面保留 log(ρ) 卷积。原 C/J、物理归一化及已有系数均不改动。`formal=not_run`。

本次得到一个具体候选：先显式求出任意单传播子矩形的二形状参数函数，再由它组装 C、D。D 的所有阶都成为三指标矩的有理常系数组合；每个矩独立满足 crossing 偶性。它解决了当前表示中外乘 ρ 的问题，但尚未把新矩约化成最小的 M/N/P/Q/R 扩充基底。剩余外层积分仍为二维，不能据此声称二维必不可少。

## 1. 为什么从矩形本身开始

一圈 M 母源的关键是把线性幂差按两个端点的几何平均组织。对 D 已完成的有序端点积分，同样的操作留下两个线性原函数和一个互补矩形。因此先研究与 D 无关的单传播子矩形：

\[
\mathcal B_\epsilon(\eta,\kappa)
=-\kappa\int_0^1dx\int_0^1dy\,
[e^\eta x+e^{-\eta}y-\kappa xy]^{\epsilon-1},
\quad 0<\kappa\leq2\cosh\eta.
\tag{R1}
\]

它有两个形状参数，分别描述两条边的相对尺度和双线性项。这里的 B 是矩形函数，不是此前三圈配对类型 B，也不是双 J 函数。原实分支下，孤立矩形在 Re(ε)>−1 收敛。

两次积分可以实际完成。令

\[
u=1-\kappa e^{-\eta},\quad v=1-\kappa e^\eta,
\qquad J_\epsilon(x)=\int_0^x\frac{(1-t)^\epsilon-1}{t}\,dt,
\]

先积 x，再作线性换元，得到

\[
\boxed{\mathcal B_\epsilon(\eta,\kappa)
=-\frac{\kappa^{-\epsilon}}{\epsilon}
[J_\epsilon(u)+J_\epsilon(v)-J_\epsilon(uv)-J_\epsilon(1)].}
\tag{R2}
\]

u、v<1 且 uv≤1，所有项使用同一实分支。ε=0 的可去商给出

\[
\mathcal B_0(\eta,\kappa)
=\operatorname{Li}_2(u)+\operatorname{Li}_2(v)
-\operatorname{Li}_2(uv)-\zeta_2.
\tag{R3}
\]

最重要的低圈连接是

\[
\boxed{\mathcal B_\epsilon\left(-\frac{\log z}{2},
\frac{1+z}{\sqrt z}\right)=\mathcal V_M(\epsilon,z).}
\tag{R4}
\]

此时 κ=2coshη，矩形的另一对角也变为零点，两个形状参数受一个关系约束，恰好回到完整的一圈对边母源。因而 R1 是原 M 几何的一个明确延伸；它不是用 D 的答案反定义的函数。

## 2. 将原 D 的四个端点积分接到这个矩形

沿用

\[
h=(1-t)(1-s)+zts,\quad q=\frac h{\sqrt z},\quad
c_z=\frac{1+z}{\sqrt z},\quad\alpha=\frac1{1+z},
\quad k_\epsilon(t,s)=-c_z q^{\epsilon-1}.
\]

原互补矩形是

\[
G_\epsilon(t,s)=\int_t^1dx\int_0^sdy\,k_\epsilon(x,y).
\]

实际作 x=1−(1−t)ξ、y=sυ，得到

\[
\boxed{G_\epsilon(t,s)=[(1-t)s]^{\epsilon/2}
\mathcal B_\epsilon(\eta(t,s),\kappa(t,s)),}
\]
\[
\eta(t,s)=\frac12\log\frac{1-t}{zs},\qquad
\kappa(t,s)=(1+z)\sqrt{\frac{(1-t)s}{z}}.
\tag{R5}
\]

这里抽出的尺度 [(1−t)s]^(ε/2) 与 z 无关。其余两个线性原函数的端点分别为

\[
A_0=\sqrt z\,t,\qquad q,\qquad B_0=\frac{1-s}{\sqrt z},
\qquad A_0B_0=t(1-s).
\]

原 D 因此是

\[
D=\int_0^1dt\int_0^1ds\,
\frac{(A_0^\epsilon-q^\epsilon)(q^\epsilon-B_0^\epsilon)}
{\epsilon^2(t-\alpha)(s-\alpha)}G_\epsilon(t,s).
\tag{R6}
\]

R5 与 R6 中两处端点几何平均的 z 幂均已相消。这正是应当在展开之前保留的结构。

## 3. 两个幂差产生自然的 cosh 差

定义三个局部对数，它们是积分核的变量，不是外部系数：

\[
L=\frac12\log\frac{zt}{1-s},\qquad
V=\log\frac q{\sqrt{t(1-s)}},\qquad
U=\frac12\log[t(1-t)s(1-s)]+\log q.
\tag{R7}
\]

与一圈 sinh 幂差相似，这里有精确恒等式

\[
(A_0^\epsilon-q^\epsilon)(q^\epsilon-B_0^\epsilon)
=2[t(1-s)]^\epsilon e^{\epsilon V}
[\cosh(\epsilon L)-\cosh(\epsilon V)].
\]

于是完整 D 成为

\[
\boxed{D(\epsilon,z)=2\int_0^1dt\int_0^1ds\,
\frac{e^{\epsilon U}[\cosh(\epsilon L)-\cosh(\epsilon V)]}
{\epsilon^2(t-\alpha)(s-\alpha)}
\mathcal B_\epsilon(\eta(t,s),\kappa(t,s)).}
\tag{R8}
\]

式中的矩形函数已由 R2 显式求值，因此只剩两个外层积分。t=α 时 V=L，s=α 时 V=−L；用

\[
\cosh(\epsilon L)-\cosh(\epsilon V)
=2\sinh\frac{\epsilon(L-V)}2\sinh\frac{\epsilon(L+V)}2
\]

可分别消去两条内部表观极点，包括交点，无需主值。整个 cosh 差从 ε² 开始，也没有真正的 ε 负幂。

## 4. 各阶系数与指标的实际含义

矩形的 raw 参数导数记作

\[
\mathcal B_n=\left.\partial_\epsilon^n\mathcal B_\epsilon\right|_0,
\qquad\mathcal B_\epsilon=\sum_{n\geq0}\frac{\epsilon^n}{n!}\mathcal B_n.
\]

它们不是逐阶待求量。若 j_r(x)=[ε^r]J_ε(x)，则

\[
\mathcal B_n=-n!\sum_{r=1}^{n+1}
\frac{(-\log\kappa)^{n+1-r}}{(n+1-r)!}
[j_r(u)+j_r(v)-j_r(uv)-j_r(1)].
\tag{R9}
\]

每个 j_r 都是已知的一变量有限 GPL 表达。暂以 X 记这组二维矩：

\[
\boxed{\mathfrak X_{m,n,j}(z)=\int_0^1dt\int_0^1ds\,
\frac{U^m(L^{2j}-V^{2j})}{(t-\alpha)(s-\alpha)}
\mathcal B_n(\eta(t,s),\kappa(t,s)),}
\quad m,n\geq0,\ j\geq1.
\tag{R10}
\]

m 数外层尺度对数，n 数矩形传播子指数的导数，j 数两个端点之间的偶次幂差。它们均有原积分操作含义，没有先计算某阶余向量再给它加下标。精确系数公式为

\[
\boxed{D_q=2\sum_{\substack{m,n\geq0,\ j\geq1\\m+n+2j=q+2}}
\frac{\mathfrak X_{m,n,j}}{m!n!(2j)!}.}
\tag{R11}
\]

例如

\[
\begin{aligned}
D_0&=\mathfrak X_{0,0,1},\\
D_1&=\mathfrak X_{1,0,1}+\mathfrak X_{0,1,1},\\
D_2&=\tfrac12\mathfrak X_{2,0,1}+\mathfrak X_{1,1,1}
+\tfrac12\mathfrak X_{0,2,1}+\tfrac1{12}\mathfrak X_{0,0,2}.
\end{aligned}
\tag{R12}
\]

D₃、D₄ 的全部 6 项、9 项组合见[独立推导 K17](independent-kernel.md)。所有外部系数均为有理数，不含 logρ、logz 或 z 的有理函数。X 本身包含运动学依赖，如同 M/N 本身依赖 z。其权重为 m+n+2j+4，R11 每阶统一为 q+6。

每个矩都可独立定义：L^(2j)−V^(2j) 含 (L−V)(L+V)，保持两条内部极点的消去。固定 z>0 时，原完整单传播子核的所有有限阶对数矩都可积，R5 对 ε 求导因此使 B_n 至多增加端点对数。又有 q≥2√[t(1−t)s(1−s)]，U、L、V 也只按端点对数增长；靠近内部两条直线的对数差商有界。

具体地，将 R8 中的三个参数释放为

\[
\mathfrak G(a,b,c;z)=2\int_0^1dt\int_0^1ds\,
\frac{e^{aU}[\cosh(bL)-\cosh(bV)]}{b^2(t-\alpha)(s-\alpha)}
\mathcal B_c(\eta,\kappa),\qquad D(\epsilon,z)=\mathfrak G(\epsilon,\epsilon,\epsilon;z).
\]

令 P=t(1−t)s(1−s)，在 |a|、|b|、|c|<δ<1/8 内，其完整被积函数可用 C(1+|logP|)²P^(−3δ) 控制。矩形因子的界来自 B_c=[(1−t)s]^(−c/2)G_c 与完整单传播子的可积界；两个 sinh 差商控制内部直线。该界在平方域可积，故三参数源在零附近全纯，每个固定指标 X 收敛，逐项展开成立。此论证保留完整差，不允许将两个 cosh 项单独积分后忽略发散。

在 z→1/z、t→1−s、s→1−t 下，U、V、κ 和分母乘积不变，L、η 变号。矩形满足 B_n(−η,κ)=B_n(η,κ)，所以每个 X 自身就是 crossing 偶函数，不需在外面再乘奇对数。

## 5. 同一个矩形也进入 C

原逆序图的有序域给出 C=∫k(t,s)G(t,s)G(s,t)。不引入另一个未求值矩形，直接得到

\[
\boxed{C(\epsilon,z)=-\int_0^1dt\int_0^1ds\,
\frac{c_z}{q}e^{\epsilon U}
\mathcal B_\epsilon(\eta(t,s),\kappa(t,s))
\mathcal B_\epsilon(\eta(s,t),\kappa(s,t)).}
\tag{R13}
\]

这解释了本次矩形函数的复用范围：它在边界给一圈 M，在 C 中成对出现，在 D 中乘两个线性端点差。C、D 因排序方式不同而取同一构件的不同组合。是否能进一步整理出像两圈 N 那样更小的一维源，仍需独立推导。

## 6. ρ 去了哪里，以及目前不应声称什么

R2 的 κ^(−ε) 没有被丢掉。按 R5 的物理参数，

\[
[(1-t)s]^{\epsilon/2}\kappa^{-\epsilon}
=\left[\frac{z}{(1+z)^2}\right]^{\epsilon/2}.
\]

因此把显式 J 重新展开回旧坐标，就会恢复旧 ρ 前因子。新表示的意义是：由独立矩形的几何平均选取归一化，把这一尺度连同内层形状组织成固定函数族，避免在最终系数外面再做 logρ 卷积。它没有消除物理运动学依赖，也没有证明 logρ 在任意函数基底中都会相消。

相比直接定义“ρ^(3ε/2)乘旧 D 母源”，这里额外建立了矩形的独立积分定义、显式求值、原 M 边界、C/D 共同用法及逐指标规则。但 X 仍是候选二维矩族；把 D₀写成 X001 本身不代表新的独立主积分数已经确定。与旧 Euler、双 J 和 P/Q/R 的逐权重约化、跨积分复用范围及一维压缩仍未完成。

下一可执行项是在 R10 的固定定义下，先求低指标矩到既有 Hlog/自然源坐标的精确字典，分离 M/N 及乘积方向，再比较 C 与其他三圈图中的同一矩。这样才能判断哪些指标带来共同的新结构，而不是只为 D 换了一套记号。

实际检查已经完成：40 项 SymPy 恒等式与 5 个负控制，另一套原生 Wolfbook/MCP 的 19 项代数与 4 个负控制均通过。新二维 J 表示对原四维有序数值积分的 7 组比较，以及 D₀ 至 D₄ 的 10 组系数检查通过。涉及的 22 个指标位置另逐一检查 crossing 和系数重建，6 个 M 边界、2 个 C 共用矩形的比较通过。系数数值路线与此前三维原函数路线共享网格，不能称完全独立求积；完整 ε 的二维与四维对照才使用不同的剩余积分域。所有浮点误差由实际精化观察，非严格区间证书。

完整解析步骤见[独立推导](independent-kernel.md)，实际代数、原生 CAS 与数值检查见[核验索引](../../../../../workflow/checks/d-intrinsic-source-20261016/README.md)。本次不更新原物理积分完成数，不宣称新的 canonical Hlog D₄ 或形式证明。
