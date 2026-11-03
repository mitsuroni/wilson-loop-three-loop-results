# 两个纯三对边积分的共享全阶生成结构

2026-10-16。对象是同一 C/I 归一化下的 `6e-20`、`6e-21=6e-22`，取正实 (z>0)。高阶系数提示应保留共同指标的端点矩；本节直接回到原积分，证明两个物理积分共用同一套矩，并把它们组成物理 (z) 的有限微分系统。这里的全阶结论来自积分恒等式和矩递推，有限阶比较是独立检查。

## 1. 从原端点矩到两个物理投影

沿用原来的全区间端点矩，而不是只取正半区间：

\[
U_n=\int_{-1}^{z}p^{n-1}(1+p)^\epsilon dp,
\qquad V_n=\int_{-1}^{z}p^{n-1}(1-p/z)^\epsilon dp.
\]

令 (A=(1+z)^{1+\epsilon})、(b=z^{-\epsilon})，写 (U_n=A u_n,V_n=Abv_n)。它们满足

\[
(n+\epsilon+1)u_{n+1}=-nu_n+z^n,\qquad
(n+\epsilon+1)v_{n+1}=znv_n+(-1)^n,
\quad u_1=v_1=(1+\epsilon)^{-1}.
\tag{P1}
\]

仍取十个二次矩的顺序

\[
\mathbf y_n=(u_n^2,u_nv_n,v_n^2,u_nZ_n,u_nR_n,
v_nZ_n,v_nR_n,Z_n^2,Z_nR_n,R_n^2)^T,
\quad Z_n=z^n,\ R_n=(-1)^n,
\]
\[
h_n=\frac{(-\epsilon)_n}{n!n}(-1/z)^n,
\qquad\mathbf F(x,z)=\sum_{n\ge1}h_n\mathbf y_nx^n.
\tag{P2}
\]

以下未写 x 的 F 均指原分支的 x=1 值。定义两个固定投影

\[
\mathsf D=A^2(F_{u^2}+b^2F_{v^2}),\qquad
\mathsf M=A^2bF_{uv}.
\tag{P3}
\]

已有 W 源中的完整平方是 (\mathcal G(1)=\mathsf D-2\mathsf M)。新的关键连接是：Lambda 不需要另一套独立的二次矩。

为看清这一点，记

\[
B_n(X)=\int_0^1t^{n-1}(1-Xt)^\epsilon dt,\qquad B_n=B_n(1).
\]

把原区间在零处分开，得到

\[
U_n=z^nB_n(-z)-(-1)^nB_n,\qquad
V_n=z^nB_n-(-1)^nB_n(-1/z).
\tag{P4}
\]

对父源

\[
\mathscr F_\epsilon(C;A_1,B_1)=\int_0^1\frac{du}{u}\int_0^1\frac{dv}{v}
(1-A_1u)^\epsilon(1-B_1v)^\epsilon J_\epsilon(Cuv)
\]

逐项代入 P4，三个物理扇区的组合恰好为

\[
\begin{aligned}
\mathscr F_\epsilon(-z;1,1)
+\mathscr F_\epsilon(-1/z;-1/z,-1/z)
-2\mathscr F_\epsilon(1;1,-1/z)
&=\sum_{n\ge1}h_nV_n^2=A^2b^2F_{v^2}.
\end{aligned}
\tag{P5}
\]

换元 (p\mapsto-p/z) 给出 (V_n(1/z)=-(-1)^nz^{-n}U_n(z))。因此 P5 的镜像正好是 (A^2F_{u^2})，两者相加就是 (\mathsf D)。证明保留了两种 n 奇偶性，不使用任何已知 Hlog 答案。

## 2. 一维余项的固定参数曲线

定义

\[
J_\epsilon(x)=\int_0^x\frac{(1-t)^\epsilon-1}{t}\,dt,\qquad
K_\epsilon(x)=\int_0^x\frac{J_\epsilon(t)}t\,dt,
\]
\[
Y_\epsilon(A_1,C)=\int_0^1\frac{(1-A_1u)^\epsilon J_\epsilon(Cu)}u\,du,
\quad
V_\epsilon(C)=\int_0^C\frac{dx}{x}\int_0^x\frac{(1-t)^{2\epsilon}J_\epsilon(t)}t\,dt,
\]
\[
\mathcal B(A_1;B_1,C\mid a,b_1,c)=
\int_0^1\frac{(1-A_1u)^aJ_{b_1}(B_1u)J_c(Cu)}u\,du.
\tag{P6}
\]

这里 (V_\epsilon(C)) 是一个源函数，和端点矩 (V_n) 区分。它也可直接写为一维积分

\[
V_\epsilon(C)=-\int_0^1\frac{\log u}{u}(1-Cu)^{2\epsilon}J_\epsilon(Cu)du.
\]

令 (J_A=J_\epsilon(A_1),J_B=J_\epsilon(B_1),J_D=J_\epsilon(D))。由[完整 Lambda 扇区恒等式](lambda-one-dimensional-insertions.md)及一次分部积分，得到

\[
\mathcal S_\epsilon(A_1,B_1,C;D)=-\mathscr F_\epsilon(C;A_1,B_1)
+\mathcal R_\epsilon(A_1,B_1,C;D),
\tag{P7}
\]
\[
\begin{aligned}
\mathcal R_\epsilon={}&J_AJ_B\left[J_D-\frac{J_A+J_B}{2}\right]\\
&+\mathcal B(A_1;A_1,C\mid\epsilon,\epsilon,\epsilon)
+\mathcal B(B_1;B_1,C\mid\epsilon,\epsilon,\epsilon)\\
&-\frac12\mathcal B(A_1;C,C\mid\epsilon,\epsilon,\epsilon)
-\frac12\mathcal B(B_1;C,C\mid\epsilon,\epsilon,\epsilon)\\
&-\mathcal B(1;A_1,C\mid0,\epsilon,2\epsilon)
-\mathcal B(1;B_1,C\mid0,\epsilon,2\epsilon)\\
&+(J_B-J_D)Y_\epsilon(A_1,C)+(J_A-J_D)Y_\epsilon(B_1,C)
+V_\epsilon(C)+J_DK_{2\epsilon}(C).
\end{aligned}
\tag{P8}
\]

具体说，对旧 Z 源中含 (K_{2\epsilon}) 的一项积分分部，产生

\[
\int_0^1\frac{(1-A_1u)^\epsilon K_{2\epsilon}(Cu)}u du
=L_{2\epsilon}(C)+J_AK_{2\epsilon}(C)
-\mathcal B(1;A_1,C\mid0,\epsilon,2\epsilon).
\]

两个 L 项与原扇区中的 (-2L) 抵消，剩下 P8。这里每个参数曲线已在积分层固定；它比“有限阶可拟合到双 J”更强，但只约化了一维余项。

物理余项是完整镜像之和：

\[
\begin{aligned}
R(z)={}&\mathcal R_\epsilon(1,1,-z;-1/z)
+\mathcal R_\epsilon(-1/z,-1/z,-1/z;-1/z)\\
&-2\mathcal R_\epsilon(1,-1/z,1;-1/z),\qquad
\mathsf R=R(z)+R(1/z).
\end{aligned}
\tag{P9}
\]

## 3. 完整物理积分的全 epsilon 组装

令

\[
\rho=\frac z{(1+z)^2},\quad\sigma=\frac{\rho^{3\epsilon/2}}{\epsilon^3},
\quad P=\frac{[e^{-\gamma_E\epsilon}\Gamma(1-\epsilon)]^3}{\cosh(3\epsilon\log z/2)},
\]
\[
T_0=J_\epsilon(-z)+J_\epsilon(-1/z)-2J_\epsilon(1),\quad
H=-\rho^{\epsilon/2}T_0/\epsilon,
\]
\[
T_J=\int_{-1}^z\frac{[(1+p)^\epsilon-(1-p/z)^\epsilon]J_\epsilon(p/z)}p dp,
\qquad C_W=-2T_0T_J+J_\epsilon(-1/z)T_0^2.
\tag{P10}
\]

从原 W 关系与 P5–P9 有

\[
\mathcal W=-\sigma(\mathsf D-2\mathsf M+C_W),\qquad
\Lambda=\sigma(\mathsf D-\mathsf R).
\tag{P11}
\]

所以最终得到

\[
\boxed{I_{20}=\frac P{16}
\left[\sigma(2\mathsf M-C_W-\mathsf R)-\frac{H^3}{2}\right],}
\tag{P12}
\]
\[
\boxed{I_{21}=I_{22}=\frac P{32}
\left[\frac{H^3}{3}+\sigma\left(\frac{\mathsf D}{2}-2\mathsf M+C_W+\frac{\mathsf R}{2}\right)\right].}
\tag{P13}
\]

I20 的对角平方源在组装时消去，只剩混合矩及一维源；I21 还保留对称对角平方。这是具体的积分层关系，不由系数空间的秩推断。

混合矩还满足更简单的一阶关系：

\[
\boxed{\theta_z\mathsf M=(1+z)^\epsilon
\left\{z^{-\epsilon}[Y_\epsilon(-z,1)-Y_\epsilon(1,-1/z)]
+Y_\epsilon(1,-z)-Y_\epsilon(-1/z,1)\right\}.}
\tag{P14}
\]

其证明是对 P3 求导：F_uv 的对角项与前因子导数相消，留下 (A^2b(F_{uR}+F_{vZ})/(1+z))。由 P4 求和，分别有

\[
AF_{uR}=Y_\epsilon(-z,1)-Y_\epsilon(1,-1/z),\qquad
AbF_{vZ}=Y_\epsilon(1,-z)-Y_\epsilon(-1/z,1).
\]

因此 I20 的耦合生成部分可由已指定边界和一维源的一次 z 积分恢复。这里没有进一步宣称该积分在全 epsilon 下必为同一个双 J 族的有限组合。

## 4. 从辅助变量系统得到物理 z 系统

对原端点积分求导，首先在端点操作合法的参数域内得到

\[
\theta_z u_n=\alpha u_n+cZ_n,\quad
\theta_z v_n=(n+\beta)v_n+cR_n,\quad
\alpha=-\frac{(1+\epsilon)z}{1+z},\quad
\beta=\frac{\epsilon-z}{1+z},\quad c=\frac1{1+z}.
\]

对 (h_n) 也必须求导，(\theta_zh_n=-nh_n)。于是

\[
\theta_z\mathbf F=D_0\mathbf F+D_1\theta_x\mathbf F,\qquad
D_1=\operatorname{diag}(-1,0,1,0,-1,1,0,1,0,-1).
\tag{P15}
\]

D0 的非零行按上述分量顺序为

\[
D_0\mathbf f=(2\alpha f_1+2cf_4,
(\alpha+\beta)f_2+cf_5+cf_6,
2\beta f_3+2cf_7,
\alpha f_4+cf_8,\alpha f_5+cf_9,
\beta f_6+cf_9,\beta f_7+cf_{10},0,0,0)^T.
\]

已有[辅助 x 系统](../pure-three-20261015/moment-closure.md)是

\[
[z\theta_x^2(\theta_x+\epsilon)^2\mathbf1
+x\theta_x(\theta_x-\epsilon)S(\theta_x)]\mathbf F
=\epsilon(1+\epsilon)^2\mathbf y_1x,
\tag{P16}
\]

其中 (S(n)=S_0+nS_1+n^2S_2)，由 P1 的对称平方完全确定。取 (F_k=\theta_x^k\mathbf F|_{x=1})，k=0,1,2,3。P16 的四阶系数矩阵为 (P_4=z\mathbf1+S_2)，它不可逆；零行的零基索引恰好是 `1,3,6,8`。关键在于 D1 的相同四行也是零，故 z 演化不需要这些行的 F4。

其余行可用

\[
P_4F_4=r-P_3F_3-P_2F_2-P_1F_1,
\]
\[
P_3=2z\epsilon\mathbf1+S_1-\epsilon S_2,\quad
P_2=z\epsilon^2\mathbf1+S_0-\epsilon S_1,\quad
P_1=-\epsilon S_0,\quad r=\epsilon(1+\epsilon)^2\mathbf y_1
\tag{P17}
\]

逐行消去 F4。这样得到四组十个 x 导数及一个常数的闭合一阶系统，不需要对奇异 P4 求逆。

再将每组的 uv、v²、vZ 分量乘以 z，其他分量保持，记所得 41 个分量为 Y。系统化为

\[
\boxed{\frac{dY}{dz}=\left[\frac{A_0(\epsilon)}z+\frac{A_{-1}(\epsilon)}{z+1}\right]Y.}
\tag{P18}
\]

[全部稀疏矩阵数据](../../../../../workflow/checks/pure-three-higher-laurent-20261016/physical-z-transport/connection.json)明确给出状态和索引；两个留数矩阵分别有 73、82 个非零元，其 epsilon 多项式次数至多三。构造本身由 [physical_z_system.py](../../../../../workflow/checks/pure-three-higher-laurent-20261016/physical_z_system.py)给出。41 是这套实现的状态数，不是主积分数，也未证明最小。P18 是 Fuchsian 系统，尚未变换成 canonical epsilon form；只有 0、−1 两个有限极点这一事实本身不足以证明所有解均属于某个指定 MPL 代数。

## 5. 边界、延拓与跨阶递推

首先在 x 的共同收敛圆盘使用 P2，再沿原正实 x 路径到达 1。为了同时合法取到 x 的四阶导数并求 z 导数，可先取 (\Re\epsilon>4)。此时在 z=1 的边界由显式收敛级数固定：

\[
(F_k(1,1))_i=\sum_{n\ge1}n^k
\frac{(-\epsilon)_n(-1)^n}{n!n}(\mathbf y_n(1))_i,
\quad k=0,1,2,3.
\tag{P19}
\]

P1 固定全部 y_n，无自由积分常数。P18 的 z=1 是常点，故这些数据确定正实 z 上的解。最后才把整套边界和解作为 epsilon 的亚纯函数延拓到零附近。源的低阶有限 GPL 规则提供独立的局部系数实现。

**不能先把各个辅助 x 导数在 epsilon=0 处逐项置零。** 例如 (F_{ZR}(x,z)=J_\epsilon(x))；先在 (\Re\epsilon>0) 取端点，得到 (\theta_xF_{ZR}(1,z)=-1)，再延拓仍是 −1；先置 epsilon=0 则 J 恒为零。物理 F0 的有限系数与辅助导数边界的取极限次序需要区分。小 epsilon 的直接物理数值检查只使用 F0 的三个平方/混合分量，不借用不合法的高阶 x 导数求和。

写 (A_\sigma(\epsilon)=\sum_{r=0}^3\epsilon^r A_{\sigma,r})，(B_r(z)=A_{0,r}/z+A_{-1,r}/(z+1))。在正确延拓的 Laurent 展开 (Y=\sum_q\epsilon^qY_q) 中，P18 给出固定跨阶关系

\[
Y_q'=B_0(z)Y_q+\sum_{r=1}^3 B_r(z)Y_{q-r}.
\tag{P20}
\]

每阶只由前三阶驱动，并使用 P19 的相应延拓边界；不存在每到一个新权重就增加待拟合函数的步骤。将所得投影代入 P12/P13，保留共同分子和所有前因子，即恢复完整物理阶。分子到 q+3 次才能得到物理 epsilon 的 q 阶。

## 6. 实际检查与结论的范围

| 检查 | 实际结果 |
|---|---|
| 十个对称平方的泛型 z 流、奇异行对齐、整数 epsilon=1…5 的全 z 方程 | 全部精确通过；漏强迫项的控制被拒绝 |
| 离散矩递推与 z 流的全参数兼容性 | 实际 Wolfbook/MCP 差零，漏 h_n 导数的控制被拒绝 |
| epsilon=8.37、12.37，z=0.7、1.4 的 G 及 theta_x²G | 8 项独立原二维核比较通过，最大观测尺度差 5.89e−13 |
| 两种 n 奇偶性的全指标父源投影及完整物理代数组装 | 全部精确通过 |
| P8 与旧完整一维余项，三个物理扇区，源次数 1…7 | 21 项精确 GPL shuffle 差零 |
| epsilon=0.23、0.37，z=1 的完整 I20/I21 | 4 项独立原有序三/四维积分比较，最大观测差 1.02e−12；反转混合项符号被拒绝 |
| P14，epsilon=0.23、0.37，z=0.7、1.4 | 4 项原父源有限差分比较通过，最大观测差 3.15e−9 |

记录：[物理 z 精确检查](../../../../../workflow/checks/pure-three-higher-laurent-20261016/physical-z-transport/exact-verification.json)、[实际兼容性 MCP](../../../../../workflow/checks/pure-three-higher-laurent-20261016/physical-z-transport/flatness-wolfbook.json)、[非整数输运](../../../../../workflow/checks/pure-three-higher-laurent-20261016/physical-z-transport/numeric-verification.json)、[投影及余项](../../../../../workflow/checks/pure-three-higher-laurent-20261016/shared-physical-generator/exact-verification.json)、[完整非整数物理积分](../../../../../workflow/checks/pure-three-higher-laurent-20261016/shared-physical-generator/numeric-verification.json)、[混合矩流 r002](../../../../../workflow/checks/pure-three-higher-laurent-20261016/shared-physical-generator/mixed-flow-verification-r002.json)。r001 混合矩辅助奇偶检查的 V 定义曾抄错，保留其源码、报告和[纠正记录](../../../../../workflow/checks/pure-three-higher-laurent-20261016/shared-physical-generator/mixed-flow-initial-check.txt)；采用的是修正并重新运行的 r002。

这些数值差是实际精化的观测值，不是区间误差证书；`formal=not_run`。较早完成的旧物理零至二阶精确回归、新权重九系数及冻结规则的权重十预测仍由各自历史报告约束。本节没有重新调参或把旧有限阶拟合提升为积分恒等式。

当前得到的结构是：**两个物理积分共用二次端点矩系统；Lambda 的其余部分有固定一维双 J 参数曲线；I20 的混合矩进一步由单 J 源的一阶 z 方程驱动。** 这回答了生成源如何由原积分组织的问题。全 epsilon 的 BJ-only 有限公式、最小主积分个数和一般独立加权有序 Q 的全部端点字仍可继续研究，它们不是本表示已经证明的结论。
