# 从排序积分导出的生成族：4e-04的全ε标准函数求值

日期：2026-09-27。本轮按用户纠正后的目标推进：从未展开积分产生权重与生成函数，不预设M/N答案。上一轮[行列式稿](determinant-sector.md)的M/N及偶幂源关系属于有限阶后验识别，不是自然基底来源的证明。

本轮的关键步骤是：**对反对称有序积分求z导数并分部积分，消去嵌套原函数；所得普通Beta矩能全部求和为标准₄F₃。** 因此既得到一个由积分结构导出的双参数生成族，也得到4e-04的有限标准函数全ε表达。生成过程没有读取M/N候选或旧系数。

建议先看§2的(B4)：这是新增解析桥；再看§4的完整Q/I及§5的任意阶递推。自然生成族已建立，其分析最小性、所有系数的M/N闭包及向三核/装饰核的推广仍未证明。

## 1. 排序实际留下的两种幂权重

复用此前已从原四维积分得到的[一维反对称原函数](all-epsilon-comparison.md#5-哪部分仍未闭合怎样接回生成与约化主线)及[全部边界连接(D5)](determinant-sector.md#3-精确接到上一轮的一维反对称原函数)。以下不读取其有限阶答案。

记

$$
\begin{aligned}
J_a(x)&=\int_0^x\frac{(1-t)^a-1}{t}\,dt
=-a x\,{}_3F_2\!\left(\begin{matrix}1,1,1-a\\2,2\end{matrix};x\right),\\
T_a(x)&=\int_0^x\frac{J_a(t)}t\,dt .
\end{aligned}
\tag{B1}
$$

物理区域z>0。先取实a,b>0，并在需要时限制到小邻域；所有底数取正实对数，内部零点连续延拓。随后沿参数解析延拓到零附近。这里a,b是两个边界幂的独立记账参数，不将它们误称为原图两个独立维数调节器。

把原反对称积分的q区间在0处分开，正半区令q=zv，负半区令q=−v，自然出现

$$
\begin{gathered}
A_a(v;z)=(1+zv)^a-1,\qquad B_b(v)=(1-v)^b-1,\\
K(a,b;z)=\int_{0<u<v<1}
\frac{A_a(u;z)B_b(v)-B_b(u)A_a(v;z)}{uv}\,du\,dv .
\tag{B2}
\end{gathered}
$$

这是由实际边界幂得到的有序交换子，不是为了匹配某个已知基底而选的权重。由两段换元，上一轮的剩余项满足

$$
\mathcal B_\epsilon=K(\epsilon,\epsilon;z)+K(\epsilon,\epsilon;1/z).
\tag{B3}
$$

更一般地，取不同幂a,b时为K(a,b;z)+K(b,a;1/z)。负半区的反向积分与交换两个权重各给一个负号，故最终是加号。

仅在(B2)给K命名还不是求值。下一节将其积分实际完成。

## 2. 新增解析桥：z微分后消去嵌套原函数

令f(v)=J_a(−zv)、g(v)=J_b(v)，则K=∫(fg'−gf')dv，并有
θ_z f=A_a、θ_z A_a=v∂_v A_a，其中θ_z=z∂_z。因此

$$
\begin{aligned}
\theta_zK
&=\int_0^1\frac{A_a B_b}{v}\,dv
-\int_0^1g(v)\,\partial_v A_a(v;z)\,dv\\
&=2P(a,b;z)-J_b(1)[(1+z)^a-1],\\
P(a,b;z)&=\int_0^1
\frac{[(1+zv)^a-1][(1-v)^b-1]}v\,dv .
\end{aligned}
\tag{B4}
$$

这里两个P相加来自一次分部积分；v=0端点gA=O(v²)，v=1端点必须保留。小参数邻域内B_b的端点幂可积，z在正实紧区间上微分有共同控制，故微分与分部积分合法。这一步没有用有限阶展开。

先取|z|<1，将第一个幂展开，再用普通Beta积分：

$$
\begin{aligned}
P
&=\sum_{n\ge1}\frac{(-a)_n(-z)^n}{n!}
\left[\mathrm B(n,1+b)-\frac1n\right]\\
&=\sum_{n\ge1}\frac{(-a)_n}{(1+b)_n}\frac{(-z)^n}{n}
-J_a(-z).
\end{aligned}
\tag{B5}
$$

所以由排序本身导出的共同生成函数是

$$
\boxed{\displaystyle
\mathfrak F(a,b;x)=
\sum_{n\ge1}\frac{(-a)_n}{(1+b)_n}\frac{x^n}{n^2}
=-\frac{a x}{1+b}\,
{}_4F_3\!\left(
\begin{matrix}1,1,1,1-a\\2,2,2+b\end{matrix};x
\right).}
\tag{B6}
$$

Pochhammer比由Beta矩产生，两个1/n分别来自原dv/v和随后对z作一次原函数。这解释了该函数的参数与指标，不是从答案反推函数形式。标准₄F₃采用[NIST定义及主分支](https://dlmf.nist.gov/16.2)；超出|x|<1时沿负实轴延拓，不穿过[1,∞)割线。

注意T_a(x)=𝔉(a,0;x)，θ_xT_a=J_a。由(B4)/(B5)沿z积分，并用K(a,b;0)=0，得到

$$
\boxed{\displaystyle
K(a,b;z)=2\mathfrak F(a,b;-z)
-2\mathfrak F(a,0;-z)-J_b(1)J_a(-z).}
\tag{B7}
$$

没有未定积分常数，也没有剩余积分。等式先在小正z成立；两边在z>0延拓的定义一致，因而覆盖整个正实区域。

## 3. 边界常数由Beta/Gamma展开确定

令ψ、ψ₁分别为digamma和trigamma，γ_E为Euler常数。由Beta(s,1+a)在s=0的展开，得到

$$
\begin{aligned}
J_a(1)&=-\gamma_E-\psi(1+a),\\
T_a(1)&=\frac{\psi_1(1+a)-\zeta_2-J_a(1)^2}{2},\\
\mathfrak F(a,a;1)&=T_a(1)+\frac12J_a(1)^2
=\frac{\psi_1(1+a)-\zeta_2}{2}.
\end{aligned}
\tag{B8}
$$

最后一式也可由两权重相同时的交换子为零推出：在(B2)令z→−1、a=b，两项逐点相同，再代入(B7)。这里仅用正底数区间的连续极限，不跨物理分支。上述式子先在收敛参数域成立，再在a=0的小邻域解析延拓；[polygamma约定](https://dlmf.nist.gov/5.15)与标准Gamma导数一致。

## 4. 完整4e-04：有限个标准函数，不留参数积分

记ρ=z/(1+z)²、ℓ=log z，并定义同一函数的对称组合

$$
\begin{gathered}
\Psi(a,b;z)=\mathfrak F(a,b;-z)+\mathfrak F(a,b;-1/z)
-2\mathfrak F(a,b;1),\\
H(\epsilon,z)=-\frac{\rho^{\epsilon/2}}{\epsilon}
[J_\epsilon(-z)+J_\epsilon(-1/z)-2J_\epsilon(1)].
\tag{B9}
\end{gathered}
$$

H就是原一圈正规化母核。Ψ仅为crossing组合的缩写，不是又引入一个待求函数。将(B7)/(B8)代回行列式稿(D5)，所有T_ε(−z)、T_ε(−1/z)项按原系数相消，结果为

$$
\boxed{\begin{aligned}
Q(\epsilon,z)
&=\frac12H(\epsilon,z)^2+
\frac{\rho^\epsilon}{\epsilon^2}
[\,2\Psi(\epsilon,\epsilon;z)-\Psi(2\epsilon,0;z)\,],\\
S(\epsilon,z)
&=-\frac{\rho^\epsilon}{\epsilon^2}
[\,2\Psi(\epsilon,\epsilon;z)-\Psi(2\epsilon,0;z)\,],\\
\Delta(\epsilon,z)
&=\frac12H(\epsilon,z)^2+
\frac{2\rho^\epsilon}{\epsilon^2}
[\,2\Psi(\epsilon,\epsilon;z)-\Psi(2\epsilon,0;z)\,].
\end{aligned}}
\tag{B10}
$$

于是项目唯一归一化下的完整I为

$$
\boxed{\displaystyle
I_{4e-04}(\epsilon,z)
=-\frac{e^{-2\gamma_E\epsilon}\Gamma(1-\epsilon)^2}
{8\cosh(\epsilon\ell)}
\left\{\frac{H(\epsilon,z)^2}{2}
+\frac{\rho^\epsilon}{\epsilon^2}
[2\Psi(\epsilon,\epsilon;z)-\Psi(2\epsilon,0;z)]\right\}.}
\tag{B11}
$$

公式先覆盖z>0、实0<ε<1；随后给出固定z下ε=0附近的局部亚纯表达。Q、S、Δ在ε=0的ε⁻²均为可去形式，H中的ε⁻¹亦可去。I在足够小且避开Gamma/cosh极点的邻域取值；ε=0用解析极限，不直接代入0/0。物理复z/i0不在本轮范围。

这是相对于此前分块多重级数与一维剩余积分的实际求值推进：所有出现的J、𝔉及其端点值均已写成有限个₃F₂、₄F₃、polygamma及初等/Gamma函数。式(B10)自动保留Q+S=H²/2；物理伴随I的颜色问题仍按旧约定分开处理。

## 5. 从生成函数直接导出的任意阶权重与递推

这里可以讨论“自然生成族”，而不借助M/N拟合。定义

$$
\mathfrak f_{m,n}(z)=[a^m b^n]\mathfrak F(a,b;-z),\qquad
V(a,b;z)=\theta_z^2\mathfrak F(a,b;-z).
$$

由(B6)的相邻Pochhammer系数比直接求得

$$
\begin{aligned}
(1+z)\theta_z V+(b-az)V&=az,\\
V(a,b;0)&=0,\qquad \mathfrak F(a,b;0)=0,\\
dV&=-bV\,\frac{dz}{z}
+(a+b)V\,\frac{dz}{1+z}
+a\,\frac{dz}{1+z}.
\end{aligned}
\tag{B12}
$$

还要求原点解析，排除积分时的log z齐次项。这是由实际已求值生成函数推出的带边界方程，不是另外假设的微分闭包。与一般[超几何微分方程](https://dlmf.nist.gov/16.8)一致，本例重复参数使其降为对V的一阶方程。

令𝓘₀f=∫₀ᶻf(t)dt/t、𝓘₋₁f=∫₀ᶻf(t)dt/(1+t)。写V=Σa^m b^n v_{m,n}，取负指标与v₀,₀为零，则

$$
\begin{aligned}
v_{m,n}
&=-\mathcal I_0v_{m,n-1}
+\mathcal I_{-1}(v_{m-1,n}+v_{m,n-1})
+\delta_{m1}\delta_{n0}\log(1+z),\\
\mathfrak f_{m,n}&=\mathcal I_0^2v_{m,n}.
\end{aligned}
\tag{B13}
$$

所有被𝓘₀积分的对象在原点消失；此处没有新增正则化常数。这个递推对任意m,n有效，自然权重为m+n+2，只产生运动学字母0、−1。这是生成函数层的字母说明，不是从有限阶输出猜字母表。

对物理需要的对角a=b=ε，令v_r=[ε^r]V(ε,ε;z)，进一步简化为

$$
\begin{gathered}
v_1=\log(1+z),\qquad
v_r=(2\mathcal I_{-1}-\mathcal I_0)v_{r-1}\quad(r\ge2),\\
[\epsilon^r]\mathfrak F(\epsilon,\epsilon;-z)
=\mathcal I_0^2(2\mathcal I_{-1}-\mathcal I_0)^{r-1}\log(1+z).
\end{gathered}
\tag{B14}
$$

例如前两个权重是

$$
\begin{aligned}
[\epsilon]\mathfrak F(\epsilon,\epsilon;-z)&=H_{0,0,-1}(z),\\
[\epsilon^2]\mathfrak F(\epsilon,\epsilon;-z)
&=2H_{0,0,-1,-1}(z)-H_{0,0,0,-1}(z).
\end{aligned}
$$

这里的权重和相对系数由(B4)/(B5)的有序积分导出。另一支𝔉(a,0;−z)则给H_{0,0,−1,…,−1}。因此区别于一圈的部分，明确表现为Beta分母参数b引入的−𝓘₀及另一份𝓘₋₁；它不是事后选几个M/N来匹配。

这给出一个可计算的自然生成族及统一递推，仍不宣称每个𝔣_{m,n}分析独立。若要称“最小基底”，还要讨论其关系、乘积及系数域；M/N可在这之后作为另一套坐标比较。

## 6. 独立生成与实际检查

[计算目录](../../../../workflow/checks/2ne-ordered-beta/README.md)将生成与回归分开：

1. generate.py只按(B13)/(B14)产生至总源次数6的Hlog词；不读旧Q、M/N或拟合结果。
2. 新鲜HyperInt会话只做shuffle/fibration和z→1/z变换，从(B10)/(B11)重建Q、Δ、I的ε⁰…ε⁴；没有重新积分旧二维核，也没有线性拟合。
3. verify.py最后才读取旧结果。**15个完整Hlog系数全部精确相同**，并核对伪极点系数为零。只使用0、−1字母。
4. 六组整数源参数的原有序双积分与(B7)精确一致，包含z符号依赖；通用Pochhammer递推、超几何匹配及边界项负控制通过。
5. 非整数独立a,b的12个有序积分点、13个完整Q点（包括一个复ε点）、端点polygamma恒等式及crossing作数值核对。其容差和限制见实际报告，数值不代替解析推导。

本轮保留所有旧结果、catalog及此前理论证据的原字节；新的全ε表达另存[机器说明](ordered-beta-generator.json)。formal=not_run，未运行Lean/MCP形式整链，也未声称首次发现文献中的同一公式。

## 7. 本轮完成与下一步

已完成的解析链是：
原4e-04 → 一维反对称原函数 → 两种边界幂的有序积分 → z微分/分部积分 → 普通Beta矩 → 标准₄F₃ → 自然系数递推。

因此，4e-04的有限标准函数全ε求值已在上述域完成；“从自然生成族进一步约化到全阶M/N”仍是另一个问题。上一轮的偶源拟合不再作为本轮推导前提。

下一步应检验带一个端点插入及交换配对时(B4)是否仍给出同样的Beta源，记录真正多出的边界/源参数，再连接已有6e目标阶数。裸交换配对S已由(B10)解决。全族最小性、3g和自动扩算ε⁵均不在本次完成范围；原329位置的新生成覆盖仍为此前21，不能重复计数4e-04回归。
