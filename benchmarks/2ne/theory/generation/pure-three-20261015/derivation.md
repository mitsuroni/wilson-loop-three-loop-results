# 三圈纯对边积分：从端点矩得到全阶生成结构

2026-10-15 至 2026-10-16。本轮研究对象是 `6e-20` 与 `6e-21=6e-22`。问题是像 `2e-02`、`4e-04` 那样，从未展开的参数积分导出生成结构，并解释已有 epsilon 平方结果的来源。以下所有新递推都先于、并独立于已知 Hlog 系数；没有通过拟合 M/N 余量来定义新函数。

完整的三核源 `W` 化为带显式初值的十分量四阶微分系统；另一个完整源 `Lambda` 化为一个单重级数，其每个系数由固定的有限递推求出。负参数的解析延拓也有显式收敛级数。两个源均有任意 epsilon 阶的有限 GPL 构造，并保留所有端点差和符号区间。这里的十分量是一个构造出的系统大小，不是十个独立物理主积分。

完整技术推导分别见 [W 的矩闭合](moment-closure.md)、[Lambda 的单重级数](lambda-single-series.md)、[W 的独立源映射及有限词规则](w-source-and-coefficients.md)与 [Lambda 的有限系数及端点规则](w_lambda_review.md)。[实际核验记录](../../../../../workflow/checks/pure-three-generators-20261015/README.md)区分解析推导、精确代数检查、非整数参数数值比较和形式证明状态。

## 1. 保持两个物理积分为最终目标

沿用[已建立的三核物理连接](../natural-beta-source.md)，固定 `z>0`，先在小正实 epsilon 区域推导，再延拓至零附近。定义

\[
D_z(t,s)=(1-t)(1-s)+zts,\qquad
k_\epsilon(t,s)=-(1+z)z^{-\epsilon/2}D_z(t,s)^{\epsilon-1},
\]
\[
F(t,s)=\int_0^t dx\int_0^sdy\,k_\epsilon(x,y),\quad
f(t)=F(t,1),\quad H=F(1,1).
\]

独立三核组合取为

\[
\mathcal W=\int_0^1\!\!\int_0^1k_\epsilon(t,s)f(t)f(s)\,dt\,ds,
\qquad
\mathcal L=\int_0^1\!\!\int_0^1F F_tF_s\,dt\,ds,
\]
\[
\Lambda(\epsilon,z)=\mathcal L(\epsilon,z)+\mathcal L(\epsilon,1/z).
\]

令 `Z=exp(-EulerGamma epsilon) Gamma(1-epsilon)`、`P=Z^3/cosh(3 epsilon log(z)/2)`。始终使用原来的物理归一化：

\[
\boxed{I_{6e-20}=\frac P{16}\left(\mathcal W+\Lambda-\frac{H^3}{2}\right),}
\tag{R1}
\]
\[
\boxed{I_{6e-21}=I_{6e-22}
=\frac P{32}\left(\frac{H^3}{3}-\mathcal W-\frac\Lambda2\right).}
\tag{R2}
\]

`H` 已是一圈已知函数。因此只需解决两个明确的新源，而不是增加任意个逐权重未知量。尤其

\[
I_{6e-20}+4I_{6e-21}=\frac P{16}(H^3/6-\mathcal W),\qquad
I_{6e-20}+2I_{6e-21}=\frac P{32}(\Lambda-H^3/3).
\tag{R3}
\]

本轮没有重新证明原 C/J 到 (R1)–(R2) 的整条链；该连接按上述已有推导复用。新的核验覆盖从这些完整三核源到新表示的步骤。

## 2. 为什么低圈经验仍然有效

在 [4e-04 的直接推导](../4e04-direct-mn.md)中，先做线性幂积分，再使用

\[
p=(1+z)t-1,\qquad q=(1+z)s-1
\]

把相互作用核变成 `1+pq/z`。端点差必须保持完整；反对称项及原函数恒等式随后让二维式降到已知的一维矩源。这是 N 自然出现的步骤。

三圈采用同一换元，但不预设下一步必须再次落回 M/N。设

\[
J_\epsilon(x)=\int_0^x\frac{(1-u)^\epsilon-1}{u}\,du
=-\epsilon x\,{}_3F_2(1,1,1-\epsilon;2,2;x),
\qquad \rho=\frac z{(1+z)^2}.
\]

原矩形原函数直接成为

\[
F(t,s)=-\frac{\rho^{\epsilon/2}}\epsilon
\left[J_\epsilon(-pq/z)-J_\epsilon(p/z)-J_\epsilon(q/z)
+J_\epsilon(-1/z)\right].
\tag{R4}
\]

这同时说明两个原则：新结构由有序积分的端点耦合产生；四个角的减项属于原积分，不能只研究其中一个漂亮的超几何子块。

## 3. W：把二维耦合转成同指标矩的乘积

写

\[
a_\epsilon(p)=(1+p)^\epsilon-(1-p/z)^\epsilon,\quad
T_n=\int_{-1}^{z}p^{n-1}a_\epsilon(p)\,dp.
\]

引入辅助变量 `x`，完整剩余源为

\[
\mathcal G(x)=\int_{-1}^{z}\!\!\int_{-1}^{z}
\frac{a_\epsilon(p)a_\epsilon(q)}{pq}J_\epsilon(-xpq/z)\,dp\,dq
=\sum_{n\ge1}\frac{(-\epsilon)_n}{n!\,n}
\left(-\frac{x}{z}\right)^n T_n^2.
\tag{R5}
\]

最右侧先在 `|x|<min(z,1/z)` 成立。关键是同一个 `n` 同时出现在两条边的矩里；二者不是两个独立求和后相乘的低圈答案。

把 `T_n=U_n-V_n` 拆成两个普通端点矩。对 `p^n(1+p)^(epsilon+1)` 和 `p^n(1-p/z)^(epsilon+1)` 做分部积分，得到

\[
\begin{aligned}
(n+\epsilon+1)U_{n+1}&=-nU_n+(1+z)^{\epsilon+1}z^n,\\
(n+\epsilon+1)V_{n+1}&=znV_n+(1+z)^{\epsilon+1}z^{-\epsilon}(-1)^n.
\end{aligned}
\tag{R6}
\]

初值 `U_1=(1+z)^(epsilon+1)/(1+epsilon)`、`V_1=z^(-epsilon)U_1` 来自初等积分。这里的非齐次端点项是闭合的一部分；遗漏它们会得到错误的生成函数。

归一化 `U_n=A u_n`、`V_n=A z^(-epsilon)v_n`，其中 `A=(1+z)^(epsilon+1)`。四个量 `(u_n,v_n,z^n,(-1)^n)` 构成有限线性递推；它们的对称平方只有十个分量：

\[
\mathbf y_n=(u_n^2,u_nv_n,v_n^2,u_nz^n,u_n(-1)^n,
v_nz^n,v_n(-1)^n,z^{2n},(-z)^n,1)^T.
\]

显式矩阵 `S(n)` 满足 `(n+epsilon+1)^2 y_(n+1)=S(n)y_n`，完整十行见[式 (M8)](moment-closure.md)。令

\[
\mathbf F(x)=\sum_{n\ge1}\frac{(-\epsilon)_n}{n!\,n}
(-1/z)^n\mathbf y_nx^n,\qquad \theta=x\partial_x.
\]

于是得到完整系统

\[
\boxed{\left[z\theta^2(\theta+\epsilon)^2\mathbf1
+x\theta(\theta-\epsilon)S(\theta)\right]\mathbf F
=\epsilon(1+\epsilon)^2\mathbf y_1 x.}
\tag{R7}
\]

`x` 必须置于所有微分算子左侧。零点解析性和 `F(0)=0` 唯一固定解；最后取

\[
\mathcal G=(1+z)^{2\epsilon+2}
(F_1-2z^{-\epsilon}F_2+z^{-2\epsilon}F_3).
\tag{R8}
\]

这个系统不随所需 epsilon 阶数增大。更强的是，epsilon 零阶的对角算子均分解为

\[
\theta^{r_i}(z+d_ix)\theta^{4-r_i}.
\tag{R9}
\]

按 epsilon 阶数递增、同阶分量倒序处理，只需有限次 `dx/x` 积分、除以 `z+d_i x` 和已知低阶运算。由此得到任意阶有限 GPL 规则，辅助变量字母表为 `{0,-z,1,-1/z}`。这不是从有限几个 epsilon 系数猜出的全阶闭合。

对 `z>0`，系统在 `0<x<1` 无奇点。物理值为沿该区间延拓后的 `x→1^-` 极限；不能把初始级数直接代入收敛圆外，也不能额外假定四阶导数在 `x=1` 有限。所有投影、端点差应先组合，再取物理极限。

## 4. 恢复 W 的低圈边界项

令 `j_+=J_epsilon(-z)`、`j_0=J_epsilon(1)`、`j_-=J_epsilon(-1/z)`，则

\[
T_0=j_++j_--2j_0,\qquad H=-\rho^{\epsilon/2}T_0/\epsilon,
\]
\[
T_J=\frac12\left[j_0(j_+-j_-)-j_0^2+j_-^2-K(\epsilon,\epsilon;z)
-K(\epsilon,\epsilon;1/z)\right].
\]

这里的低圈函数已显式求值：

\[
K(a,b;z)=2\mathfrak F(a,b;-z)-2\mathfrak F(a,0;-z)-J_b(1)J_a(-z),
\]
\[
\mathfrak F(a,b;x)=-\frac{ax}{1+b}
{}_4F_3(1,1,1,1-a;2,2,2+b;x),\qquad
J_\epsilon(1)=-\gamma_E-\psi(1+\epsilon).
\]

完整结果因此为

\[
\boxed{\mathcal W=-\frac{\rho^{3\epsilon/2}}{\epsilon^3}
\left[\mathcal G(1)-2T_0T_J+j_-T_0^2\right].}
\tag{R10}
\]

这已经覆盖 (R3) 的第一个物理组合。显式 `epsilon^(-3)` 不代表原图存在三阶极点：完整分子从三阶开始。各辅助分量可从一阶开始，必须按 (R8)、(R10) 完整相加。

## 5. Lambda：求出每一个级数系数

定义完整符号区间源

\[
\mathcal S_\epsilon(A,B,C;d)=\int_0^1\frac{du}{u}\int_0^1\frac{dv}{v}
[J_\epsilon(Cuv)-J_\epsilon(Au)-J_\epsilon(Bv)+J_\epsilon(d)]
\]
\[
\hspace{12mm}\times[(1-Cuv)^\epsilon-(1-Au)^\epsilon]
[(1-Cuv)^\epsilon-(1-Bv)^\epsilon].
\]

保留区间方向后，原源恰为

\[
\mathcal L(\epsilon,z)=-\frac{\rho^{3\epsilon/2}}{\epsilon^3}
\left[\mathcal S(1,1,-z;-1/z)
+\mathcal S(-1/z,-1/z,-1/z;-1/z)
-2\mathcal S(1,-1/z,1;-1/z)\right].
\tag{R11}
\]

继续展开的只有乘积变量 `Cuv`。其余变量的积分已全部求出，得到

\[
\boxed{\mathcal S_\epsilon(A,B,C;d)
=J_AJ_B\left[J_d-\frac{J_A+J_B}{2}\right]
+\sum_{n\ge1}C^n s_n(A,B;d).}
\tag{R12}
\]

`s_n` 的完整七项表达式见[式 (L9)](lambda-single-series.md)。其中只出现二项式系数、它与 `J` 的有限卷积，以及已求值的端点矩

\[
B_n^c(A)=\frac1n{}_2F_1(-c,n;n+1;A),\qquad
Q_n(A)=\frac{J_\epsilon(A)-B_n^\epsilon(A)}n+\frac1{n^2},
\]
\[
H_n(A)=\int_0^1u^{n-1}(1-Au)^\epsilon J_\epsilon(Au)\,du.
\]

最后一个符号不是新的未求积分；它有有限 `3F2` 和式，并满足

\[
\boxed{A(n+\epsilon+1)H_{n+1}
=nH_n+B_n^{2\epsilon+1}-B_n^{\epsilon+1}
-(1-A)^{\epsilon+1}J_\epsilon(A),}
\tag{R13}
\]
\[
H_1=-\frac{(1-A)^{\epsilon+1}J_\epsilon(A)
+J_{\epsilon+1}(A)-J_{2\epsilon+1}(A)}{A(\epsilon+1)}.
\]

`B_n` 也有一阶有理递推。因此 (R12) 每个系数只有有限、显式的计算，不留下逐阶待求二维积分。零指标项已先完整相消再积分，不能单独拆成发散项。

对于物理负参数，进一步取

\[
w=\frac C{C-2},\qquad C=-\frac{2w}{1-w},\qquad
\widetilde s_k=\sum_{n=1}^k(-2)^n\binom{k-1}{n-1}s_n.
\]

于是

\[
\mathcal S=s_0+\sum_{k\ge1}w^k\widetilde s_k.
\tag{R14}
\]

原积分在 `Re C<1` 全纯；这一半平面恰好映到 `|w|<1`。所以任意有限负 `C` 都有真正收敛的表示，不必在 `|C|>1` 处使用发散 Taylor 和。混合区间 `C=1` 使用原级数的 Abel 极限。高阶二项式变换可能有数值消去，这属于计算精度问题，不改变解析等式。

最后把 (R11) 在 `z` 和 `1/z` 的值相加，即得 `Lambda`，再代回原左侧 (R1)–(R2)。

## 6. Lambda 的任意阶系数也只需有限步

全 epsilon 收敛级数和任意 epsilon 阶的有限答案是两种不同要求。为满足后一要求，直接展开原来的三个差因子。采用积分核 `dt/(t-a)` 的 GPL 约定，有

\[
(1-X)^\epsilon=\sum_{m\ge0}\epsilon^mG(\underbrace{1,\ldots,1}_{m};X),\qquad
J_\epsilon(X)=\sum_{m\ge1}\epsilon^mG(0,\underbrace{1,\ldots,1}_{m};X).
\]

每个差因子从一次开始。因此源的第 N 阶只含正整数 `i+j+k=N` 的有限个组合。把同一变量的词 shuffle 后，在最前添零就完成其 `du/u` 积分。内层所得词只有固定字母 `0,1/A` 与移动字母 `1/(Cv)`；Chen 字母变动公式每次降低词长，因此可以在有限步内转成外层变量的 GPL。再 shuffle 和前添零，第二个积分也被消去。

完整算法及初值见 [Lambda 有限规则 (A2)–(A5)](w_lambda_review.md)。泛型外层字母为 `{0,1/B,1/C,A/C}`；代入三个物理区间并加上镜像后，只需 `{0,1,-z,-1/z}`。这是辅助积分变量的字母表，不是已经最小化的物理 z 微分方程字母表。

`A=1` 时，先保留内层常数中的所有 `log(1-A)` 幂，合并趋于相同的字母，再消去对数系数。原有限积分固定这个共同端点极限；不能逐项把发散 GPL 直接置零。源三、四、五阶的两个相关物理区间，已实际检查这些对数系数逐词消去。

这些源阶数对应权重六、七、八。泛型完整系数实际产生 86、342、1136 个 GPL 乘积项，与独立展开原核后的二维求积一致；它们没有读取旧物理系数表。由此，`Lambda` 的单重全阶表示与每阶有限计算规则都已明确。

## 7. 新结构的含义和仍需完成的工作

三圈的新信息有一个直接来源：`W` 中是两个端点矩的同指标乘积，`Lambda` 中是带原函数插入的端点矩。它们在未展开 epsilon 的原积分中已经存在。有限权重的 M/N 非成员检验，可以检验这些源投影后产生了哪些新方向，却不是这些源的定义。

目前不能据此断言“三圈至少需要二维积分”或“每圈恰好增加一个独立主积分”。本轮恰好表明，二维耦合可以用单变量微分系统或单重级数表达；原积分的维数、表示中的变量个数、最小主积分数是不同问题。`W` 的逐阶 GPL 构造也说明，在这个通道出现 M/N 以外方向，并不自动意味着出现椭圆函数。

本轮完成的是两个完整源的生成表示、边界条件与每阶有限规则。尚未求出最小物理 `z` 微分系统，也未证明源不能进一步化成有限个普通超几何函数。新的物理零阶至二阶 Hlog 表没有重新生成；原有结果继续按原来源使用。

下一项可执行工作是在物理参数切片上组合两个源至分子五阶，恢复完整 `I` 的零阶至二阶，与旧结果逐系数比较，再压缩物理 `z` 系统和权重六的新方向。需要分子五阶，是因为 (R10)–(R11) 带有三阶显式除子。已有的泛型系数与端点消去检查不能替代整张物理结果表的回归；旧的 `44/47` 完成盘点本轮不自动改写为 `47/47`。

作为物理连接的独立数值控制，在 `epsilon=1/2,z=1`，新生成式给出 `I20=0.0793498576864`、`I21=I22=0.0454614510299`。与原六维有序匹配积分的八组独立扰乱 Sobol 求积比较，差异分别为 `1.39e-7` 和 `1.01e-7`，小于各自约 `2.86e-7`、`2.31e-7` 的经验标准误差。这是一个物理点的控制，不是一般 z 的证明。完整 `W` 的另三个非整数参数控制达到最大相对差 `9.37e-11`；完整 `L` 的独立原核比较达到绝对差 `6.30e-13`。所有数值误差都只有观测收敛意义，详见核验记录。

本轮新增文件不修改原 C/J、prepared 输入、旧结果对象、Wilson/Lean 或其他研究线。形式证明状态为 `formal=not_run`；CAS 恒等式与数值比较均不冒充原任务的形式证明。
