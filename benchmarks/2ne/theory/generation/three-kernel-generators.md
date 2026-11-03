# 三核行列式与混合组合：从六维排序到二维原函数矩

日期：2026-09-28。按用户要求，参考4e-04的“先处理未展开积分、保留边界、再识别生成权重”路线，研究三核的行列式组合与另一个混合排列组合。

本轮首先得到两个全ε降维恒等式：目标都能由已知低核乘积和两个二维矩表示。二维矩中的矩形单核原函数已显式求成与4e-04相同的三阶超几何原函数。行列式所需的二维矩进一步得到一个由实际边界幂产生的单重对角矩级数；其物理取值的延拓路径明确。本稿不把二维积分表示或需延拓的级数称为有限标准函数闭式。

旧[权重7研究](../weight-seven.md)已完成三核ε¹，而[旧排列稿](../permutation-masters.md)中“q≥1未计算”的段落是较早状态。本轮复用旧证据，不将它重复计作新求值。新路线与检查见[复验](../../../../workflow/checks/three-kernel-generators/README.md)；formal=not_run，当前会话无Wilson MCP工具，不修改Kernel/Lean。

## 1. 先固定两个目标及它们与物理图的连接

沿用

$$
k_\epsilon(t,s)=-(1+z)z^{-\epsilon/2}
[(1-t)(1-s)+zts]^{\epsilon-1},
\qquad k_{ij}=k_\epsilon(t_i,s_j),
$$

以及两边各三个端点递增的域Ω。记

$$
A=R_{123},\quad B=\frac{R_{132}+R_{213}}2,
\quad C=R_{321},\quad D=R_{231}=R_{312}.
$$

两个目标为

$$
\boxed{\begin{aligned}
\Delta_3&=\int_{\Omega_3}\det(k_{ij})\,d^3t\,d^3s
=A-2B-C+2D,\\
X_3&=\int_{\Omega_3}
(2k_{11}k_{22}k_{33}-k_{12}k_{23}k_{31}-k_{13}k_{21}k_{32})\,d^3t\,d^3s
=2A-2D.
\end{aligned}}
\tag{T1}
$$

裸单核为H，裸同序双核为Q。已证明的两条乘积关系给出

$$
\boxed{\begin{aligned}
\Delta_3&=-4A-8B+2HQ-\frac{H^3}{6},\\
X_3&=5A+4B-HQ,\\
A&=\frac{H^3}{36}+\frac{\Delta_3}{6}+\frac{X_3}{3},\\
B&=-\frac{5H^3}{144}+\frac{HQ}{4}-\frac{5\Delta_3}{24}-\frac{X_3}{6}.
\end{aligned}}
\tag{T2}
$$

这只固定结构坐标；两个候选不等于已证明的最小解析函数数。所有母核先取正实z、实小正ε，随后在零附近延拓。角点可积性来自各端点只属于一个核，放宽为立方体后可由独立核乘积控制；有限次对数导数也可支配。这允许下文换序与分部积分。

最后恢复物理I时使用

$$
I_{6e-20}=-\frac{Z^3A}{16\cosh(3\epsilon\ell/2)},\qquad
I_{6e-21}=I_{6e-22}=-\frac{Z^3B}{32\cosh(3\epsilon\ell/2)},
\quad Z=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),\quad\ell=\log z.
\tag{T3}
$$

不将裸核的乘积系数直接套到物理I。

## 2. 真正需要的新工具：任意矩形上的单核原函数

定义已经包含核归一化的矩形积分

$$
F(t,s)=\int_0^t dx\int_0^s dy\,k_\epsilon(x,y),
\qquad f(t)=F(t,1),\qquad H=F(1,1).
$$

对称核保证F对称，且F的混合导数等于原核。相对的右上矩形为

$$
U(t,s)=\int_t^1dx\int_s^1dy\,k_\epsilon(x,y)
=H-f(t)-f(s)+F(t,s).
\tag{T4}
$$

这些原函数不能停留在名字上。令

$$
p=(1+z)t-1,\quad q=(1+z)s-1,\quad
\rho=\frac{z}{(1+z)^2},
$$
$$
J_\epsilon(x)=\int_0^x\frac{(1-u)^\epsilon-1}{u}\,du
=-\epsilon x\,{}_3F_2(1,1,1-\epsilon;2,2;x).
$$

先积一个端点，再用上述线性变量变换，得到

$$
\boxed{\displaystyle
F(t,s)=-\frac{\rho^{\epsilon/2}}\epsilon
\left[J_\epsilon(-pq/z)-J_\epsilon(p/z)-J_\epsilon(q/z)+J_\epsilon(-1/z)\right].}
\tag{T5}
$$

四项来自矩形的四个角，缺一项就会破坏下边界。直接微分可核对：

$$
\partial_p\partial_qJ_\epsilon(-pq/z)
=\frac\epsilon z(1+pq/z)^{\epsilon-1},
$$

其余三项混合导数为零。乘回变量变换因子，正好得到原核；t或s为零时四角项成对消去。t与s均为1时，T5退回已知一圈H。

所有J自变量不超过1；负实自变量沿原实分支取值，达到1时用可积端点极限。ε为零是可去极限，此时J除以ε趋于负二重对数。35位原矩形二重求积及整数ε的混合导数/四角检查均支持此式。

**当前进度：** 矩形单核已实际求值；后面的矩形符号都可由T5替换为明确标准函数。

## 3. 利用中间传播子，把同序三核降到二维

固定同序图的中间端点为t、s。较小的两个端点形成左下矩形F，较大的两个端点形成右上矩形U，两组变量独立。因此

$$
\boxed{\displaystyle
A=\int_0^1dt\int_0^1ds\,k_\epsilon(t,s)F(t,s)U(t,s).}
\tag{T6}
$$

这里六维中的四个变量已经积完。没有将有序域扔掉：两个矩形正是中间端点给出的有序上下限。

再考虑与A相邻的两种配对。固定最大的一对端点，下面两条传播子的两个配对之和等于矩形单核平方的一半；固定最小的一对端点则同理。因此

$$
\begin{aligned}
R_{123}+R_{213}&=\frac12\int k_\epsilon F^2,\\
R_{123}+R_{132}&=\frac12\int k_\epsilon U^2,\\
\boxed{\displaystyle A+B}&=\frac14\int k_\epsilon(F^2+U^2).
\end{aligned}
\tag{T7}
$$

本稿未标范围的二重积分均为单位方形，测度为dt ds。

**当前进度：** 两个需要的三核偶组合都已降为二维的标准原函数矩。此处仍有两个真正需要处理的二维积分，尚不能宣布完成全ε求值。

## 4. 行列式进一步简化：只剩边界原函数的加权单核

由T4，两矩形之差只含边界原函数：

$$
F-U=f(t)+f(s)-H.
$$

定义边界权重矩

$$
\mathcal W=\int k_\epsilon(t,s)f(t)f(s)\,dt\,ds.
$$

注意这里的花体W不是M矩中的对数权重。用f的导数等于核的边缘积分，得到

$$
\int k_\epsilon f(t)=\frac{H^2}{2},\qquad
\int k_\epsilon f(t)^2=\frac{H^3}{3}.
$$

展开矩形差的平方，同时使用T6/T7：

$$
\int k_\epsilon(F-U)^2
=2A+4B
=2\mathcal W-\frac{H^3}{3}.
$$

所以新的全ε桥是

$$
\boxed{\displaystyle A+2B=\mathcal W-\frac{H^3}{6}.}
\tag{T8}
$$

代回两个原目标，而不是只给一个辅助等式：

$$
\boxed{\begin{aligned}
\Delta_3&=2HQ+\frac{H^3}{2}-4\mathcal W,\\
X_3&=3\int k_\epsilon FU+2\mathcal W-HQ-\frac{H^3}{3}.
\end{aligned}}
\tag{T9}
$$

**当前进度：** 行列式只需要一个以边界原函数为权重的单核矩；混合组合另需矩形耦合矩。它们具有共同的已求值原函数来源，但不能因此把两者认作同一个函数。

这也给出检验“新三核函数是否真正需要”的具体对象：先研究这两种权重能否由此前端点源及其边界关系产生，再判定剩余部分。

## 5. 对行列式的加权矩，再做一次带边界分部积分

因为核是F的混合导数，两次分部积分给出

$$
\begin{aligned}
\mathcal W
&=\frac{H^3}{2}-\frac{H^3}{2}
+\int F(t,s)f'(t)f'(s)\,dt\,ds\\
&=\int F(t,s)f'(t)f'(s)\,dt\,ds.
\end{aligned}
\tag{T10}
$$

两个边界常数实际存在，是相消后才离开最终式。现在F只出现一次，其余权重均为初等端点幂差。

沿T5的p坐标，令

$$
a_\epsilon(p)=(1+p)^\epsilon-(1-p/z)^\epsilon,
\quad -1<p<z.
$$

边界原函数的微分为

$$
df=-\frac{\rho^{\epsilon/2}}\epsilon\frac{a_\epsilon(p)}p\,dp.
$$

p为零时差商连续。将T5完整代入T10：

$$
\boxed{\displaystyle
\mathcal W=-\frac{\rho^{3\epsilon/2}}{\epsilon^3}
\left[\mathcal G(1)-2T_0T_J+J_\epsilon(-1/z)T_0^2\right],}
\tag{T11}
$$

其中

$$
\begin{aligned}
\mathcal G(x)&=\int_{-1}^{z}\!\int_{-1}^{z}
\frac{a_\epsilon(p)a_\epsilon(q)}{pq}
J_\epsilon(-xpq/z)\,dp\,dq,\\
T_0&=\int_{-1}^{z}\frac{a_\epsilon(p)}p\,dp,
\qquad
T_J=\int_{-1}^{z}\frac{a_\epsilon(p)}pJ_\epsilon(p/z)\,dp.
\end{aligned}
$$

T0是已知单核，TJ也是已经在4e-04中处理的有序原函数组合。具体令a、b、c分别为J在负z、1、负z倒数处的值，令K为[有序Beta稿B2/B7](ordered-beta-generator.md)中已经求值的交换子，则

$$
T_0=a+c-2b,\qquad
T_J=\frac{b(a-c)-b^2+c^2-K(\epsilon,\epsilon;z)-K(\epsilon,\epsilon;1/z)}2.
\tag{T12}
$$

因此行列式的新内容已经集中到一个明确的乘积自变量生成函数G；其他项都有既有全ε求值。

## 6. G的单重级数由真实边界矩产生

引入辅助参数x用于展开，物理取值最后为1。对J的乘积自变量作幂级数，然后独立积掉p、q：

$$
\boxed{\displaystyle
\mathcal G(x)=\sum_{n\ge1}
\frac{(-\epsilon)_n}{n!\,n}
\left(-\frac{x}{z}\right)^n T_n^2,
\qquad
T_n=\int_{-1}^{z}p^{n-1}a_\epsilon(p)\,dp.}
\tag{T13}
$$

每个系数是同一个边界矩的平方；这不是从M/N或Hlog拟合得到的指标。将区间在零处分开，只需普通Beta和超几何矩：

$$
\boxed{\begin{aligned}
T_n={}&z^n\left[\frac{{}_2F_1(-\epsilon,n;n+1;-z)}n-\mathrm B(n,1+\epsilon)\right]\\
&+(-1)^n\left[\frac{{}_2F_1(-\epsilon,n;n+1;-1/z)}n-\mathrm B(n,1+\epsilon)\right].
\end{aligned}}
\tag{T14}
$$

这使行列式剩余对象成为所有有限项均已求值的单重级数。必须保留其收敛边界：T13最初在

$$
|x|<\min(z,1/z)
$$

绝对收敛。一般z不能把x直接取1后声称仍在该收敛圆内。物理值由积分定义指定沿实区间从零到1的解析延拓与端点极限：此路径上J的最大自变量为x，不跨越大于1的割线。先对小正ε建立，再对ε作局部延拓。

z为1时可以直接在物理点使用收敛级数；奇数n的Tn消失，偶数n由T14给出。三组小x级数、九个独立高精度边界矩及一个物理z等于1的加权矩已数值核对。没有据此宣称一般z下已有快速高精度求和算法。

**当前进度：** 行列式的剩余部分有来自端点幂的单重对角矩生成式；混合组合还保留T6的矩形耦合积分。未证明二者都能压成有限个三阶或四阶超几何函数。

## 7. 旧低阶结果现在约束什么

记A、B和两个目标的下标为ε阶数。已有零阶关系为

$$
A_0+B_0=\frac{H_0^3}{36},\qquad
A_0=2E-\frac{H_0^3}{9}
-\frac12M_{0,1}M_{2,1}+\frac12M_{1,1}^2.
$$

于是

$$
(\Delta_3)_0=4A_0-\frac{H_0^3}{18},\qquad
(X_3)_0=A_0-\frac{H_0^3}{18}.
\tag{T15}
$$

零阶两组合相对于M/N只带同一个已知E方向。权重7的旧计算又给出A的一阶与B的一阶之和属于M/N。由T2直接推出

$$
\boxed{\displaystyle
\Delta_3-4X_3=6HQ-\frac{H^3}{6}-24(A+B).}
\tag{T16}
$$

故这个组合的零阶和一阶都可回到M/N。不能由这两阶猜测它全ε恒属于M/N。已有一阶的每个目标分别可回到声明的双核两源矩空间，这也不等于所有单核或双核源的最小性结论。

## 8. 不依赖低阶拟合的整数参数控制

这里仅检查裸母核，不在Gamma极点处求物理I。ε为1时核为常数、矩阵秩为1，三阶行列式及其一阶参数导数均为零；ε为2时核是双线性多项式，矩阵秩至多2，行列式仍为零。因此在母核解析域内，行列式可提出下面的可去因子：

$$
\Delta_3=(\epsilon-1)^2(\epsilon-2)\mathcal R_3(\epsilon,z).
$$

这并未求出余下函数，但约束了任何候选通式。混合组合在ε为2时却为

$$
X_3(2,z)=-\frac{7(1+z)^4}{2400z^2},
$$

说明不能把它当成行列式的同一个积分。ε为3时，核的二次多项式分解又给出

$$
\Delta_3(3,z)=-\frac{(1+z)^3}{16200z^{3/2}}.
$$

这来自两个有序Vandermonde矩，每个等于1/180。这些控制与原六维积分的多项式求积分别核对；不是用整数样点外推通式。

## 9. 当前边界与下一次具体求解目标

本稿得到未展开积分的解析降维、矩形原函数求值、行列式的对角矩级数，以及原六维数值与多项式控制。本轮尝试三核ε²的原三维/四维路线在900秒限时内未完成第一项；新二维加权矩ε²路线在300秒内也未输出结果。因此没有新的ε²解析系数。两次CAS尝试及其实际完成范围记录在[检查目录](../../../../workflow/checks/three-kernel-generators/README.md)，不把运行中或超时的目标当作已算系数。

下一步优先处理T13的对角矩生成函数及T6的矩形耦合：尝试参数微分、Mellin矩递推和保留边界的分部积分，检查是否可用已有带端点源消去。需要固定阶判别时，优先使用本稿二维表示，比较T16的下一阶是否仍落在明确的低核候选空间；不要只扩大旧高维Hlog积分或把拟合失败说成新函数族已经存在。

原47图结果、此前35/329新框架覆盖和全部历史证据保持。本轮全ε结构研究不重复增加旧有限项覆盖；尚无有限特殊函数闭式、全阶M/N/τ/E闭包或形式证明。
