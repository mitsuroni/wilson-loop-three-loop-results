# 36/39：独立端点指数的 Gauss 解与归一化指标矩

2026-10-16。接续 [共同矩形来源桥](../six-e-common-blocks-20261016/source-bridge.md) 与 [单对边全阶表示](../../../families/single-opposite/three-loop.md)。目标是保留传播子及两个端点的独立指数，解释原物理切片背后的生成结构。本次没有改变原 C/J 或发布结果，也没有新求一个物理 Laurent 阶；formal=not_run。

本次完成三件事：独立参数的非齐次 Gauss 方程和全部边界常数；无需剩余求积的三参数闭式；提前合并物理镜像的固定三指标矩。后者已经数值连接到完整物理 I 的权重四至八，但尚未约化为最小的 M/N 扩充基底。

## 1. 从原端点差保留三个独立参数

固定 z>0，记

\[
D_z(t,s)=(1-t)(1-s)+zts,\quad
q_z=D_z/\sqrt z,\quad
k_\lambda(t,s)=-(1+z)z^{-\lambda/2}D_z^{\lambda-1},
\]
\[
u_a(t)=\frac{1-t^a}{a},\qquad u_0(t)=-\log t.
\]

已有原积分来源桥给出

\[
E_\lambda(a,b;z)=\int_0^1dt\,ds\,k_\lambda(t,s)u_a(t)u_b(s).
\tag{1}
\]

它对应独立端点响应；36/39 的物理切片为 λ=a=b=ε。为实际求值，先处理不带归一化的幂矩

\[
F_\lambda(a,b;z)=\int_0^1dt\,ds\,t^a s^b D_z(t,s)^{\lambda-1}.
\tag{2}
\]

于是

\[
\boxed{
E_\lambda(a,b;z)=
-\frac{(1+z)z^{-\lambda/2}}{ab}
[F_\lambda(a,b)-F_\lambda(a,0)-F_\lambda(0,b)+F_\lambda(0,0)].
}
\tag{3}
\]

这里 F 与前次建议中的带 k 的 T 相差明确的公共因子，不能混用。四项必须共同取 a=0 或 b=0 的极限。式(3)的混合有限差分正是两个端点减除的来源。

初始可取 0<Re λ<1、Re a,Re b>0。固定正 z 后，式(2)在参数零点附近全纯：端点要求 Re a,Re b>−1，以及 Re(a+λ),Re(b+λ)>−1。下文超几何式中孤立的 1/λ、1/λ² 是表示的可去奇点，应先组合再取极限。

## 2. 任意端点指数仍给出同一个简单驱动项

用 Mellin–Barnes 展开 D 的两项并分别执行两个 Beta 积分：

\[
F_\lambda(a,b;z)
=\frac{1}{2\pi i}\int_{\cal C}z^w
\frac{\Gamma(-w)\Gamma(1-\lambda+w)
\Gamma(1+a+w)\Gamma(1+b+w)\Gamma(\lambda-w)^2}
{\Gamma(1-\lambda)\Gamma(1+a+\lambda)\Gamma(1+b+\lambda)}
\,dw.
\tag{4}
\]

初始正实参数下可选 λ−1<Re w<0 的向上直线。Gamma 因子的移位比是

\[
\frac{\mathcal M(w+1)}{\mathcal M(w)}
=-\frac{(w+1+a)(w+1+b)}{(w+1)(w+1-\lambda)}.
\tag{5}
\]

记 θ=z∂z。移动一格围道时，算子中的 w(w−λ) 消去 w=0 的简单极点，并把 w=λ 的双极点降为简单极点。Γ(−λ)/Γ(1−λ)=−1/λ 与剩余 λ 因子相抵；按右闭围道方向保留该残数，得到

\[
\boxed{
[\theta(\theta-\lambda)+z(\theta+1+a)(\theta+1+b)]F_\lambda(a,b;z)
=z^\lambda.
}
\tag{6}
\]

右侧与 a、b 无关。这是释放端点指数后仍保持简单结构的关键。令 a=b=λ=ε，准确恢复旧双端源的方程 S2。

作为不依赖围道移位的控制，λ=m=1,2,3,4 时原核是多项式：

\[
F_m(a,b;z)=\sum_{j=0}^{m-1}{m-1\choose j}z^j
B(1+a+j,m-j)B(1+b+j,m-j).
\tag{7}
\]

把这些原积分直接代入式(6)，右侧均为 z^m。非整数参数另以原二维核核对，见最后的检查范围。

### 2.1 完整端点差满足一个由边界驱动的三角系统

令 L_{a,b} 表示式(6)左侧算子（λ 固定），并记

\[
F_0=F_\lambda(0,0),\quad
U_a=\frac{F_0-F_\lambda(a,0)}a,\quad
U_b=\frac{F_0-F_\lambda(0,b)}b,\quad
\mathcal R_{a,b}=\frac{F_\lambda(a,b)-F_\lambda(a,0)-F_\lambda(0,b)+F_0}{ab}.
\]

把四项各自的式(6)相减，直接得到

\[
\boxed{
\begin{aligned}
L_{0,0}F_0&=z^\lambda,\\
L_{a,0}U_a&=z(\theta+1)F_0,\\
L_{0,b}U_b&=z(\theta+1)F_0,\\
L_{a,b}\mathcal R_{a,b}
&=z[(\theta+1+a)U_a+(\theta+1+b)U_b-F_0].
\end{aligned}}
\tag{7a}
\]

这里 U_a、U_b 是单端点响应，与下一节的 Gauss 原函数 U 不同。最底层准确接回既有 M 生成函数：

\[
F_0=-\frac{z^{\lambda/2}}{1+z}
\mathcal V_M(\lambda/2,\lambda/2;z).
\tag{7b}
\]

式(7a)显示无插入、单端点和双端点之间的具体关系；共同源 z^λ 在完整双端差中相消，留下单端响应及必要的 −F₀ 边界。这比单纯定义新矩更强：各层的驱动和低层极限均已固定。原多项式积分的四组检查也分别验证了后两层关系。

## 3. 求出两个齐次解与全部边界常数

作 Euler 规范变换

\[
\nu=1+a+b+\lambda,\qquad Y=(1+z)^\nu F_\lambda.
\]

式(6)成为

\[
[\theta(\theta-\lambda)+z(\theta-a-\lambda)(\theta-b-\lambda)]Y
=z^\lambda(1+z)^\nu.
\tag{8}
\]

两个齐次解为

\[
\Phi={}_2F_1(-a-\lambda,-b-\lambda;1-\lambda;-z),
\]
\[
\Xi=z^\lambda{}_2F_1(-a,-b;1+\lambda;-z).
\tag{9}
\]

θ-Wronskian 为

\[
\Phi\,\theta\Xi-\Xi\,\theta\Phi
=\lambda z^\lambda(1+z)^{a+b+\lambda}.
\tag{10}
\]

将源项除以最高阶系数(1+z)和 Wronskian，恰好只剩 1/λ。因此常数变易不产生新的复杂卷积核，只需两个原函数：

\[
U=\log z+Q,\qquad \theta U=\Phi,\qquad \theta V=\Xi,
\]
\[
\boxed{
Q=-\frac{(a+\lambda)(b+\lambda)z}{1-\lambda}
{}_4F_3\!\left(
\begin{matrix}1,1,1-a-\lambda,1-b-\lambda\\2,2,2-\lambda\end{matrix};-z\right),
}
\tag{11}
\]
\[
\boxed{
V=\frac{z^\lambda}{\lambda}
{}_3F_2\!\left(
\begin{matrix}-a,-b,\lambda\\1+\lambda,1+\lambda\end{matrix};-z\right).
}
\tag{12}
\]

式(11)来自 ∫₀ᶻ(Φ−1)du/u，式(12)先在 Re λ>0 取原函数。先在 |z|<1 逐项验证，再沿负实超几何自变量延拓到所有正 z。

不能把齐次解系数当成自由拟合常数。式(4)的 w=0 与 w=λ 残数分别给出

\[
A_0=B(1+a,\lambda)B(1+b,\lambda),
\]
\[
d_{\lambda,a,b}=\gamma_E-\psi(1-\lambda)
+\psi(1+a+\lambda)+\psi(1+b+\lambda),
\]
\[
F_\lambda
=A_0+\frac{z^\lambda}{\lambda}
\left(\log z-\frac1\lambda+d_{\lambda,a,b}\right)
+O(z)+O(z^{1+\lambda}\log z).
\tag{13}
\]

特别是双极点产生的 log z 及常数 −1/λ 都不能删除。匹配两个 Frobenius 分支后得到完整闭式：

\[
\boxed{
F_\lambda(a,b;z)
=\frac{1}{(1+z)^{1+a+b+\lambda}}
\left\{
A_0\Phi+
\frac{\Xi[\log z+d_{\lambda,a,b}+Q]-\Phi V}{\lambda}
\right\}.
}
\tag{14}
\]

初始非共振参数域中，方程及两个分支系数固定了解；由解析延拓得到参数零点附近的联合极限。式(14)含标准特殊函数，无剩余参数积分。与式(3)结合，独立三参数 E 已求值，而不只是重新命名原二维积分。

## 4. 自然的解析构件：非等上参数 Gauss 与一次 Mellin

式(9)说明释放 a、b 后，应允许 Gauss 的两个上参数彼此独立。原来相等上参数的 P 源仍是其切片，不应提前锁住两个参数。

为明确其矩来源，记

\[
\mathscr P(A,B,C;z)
=\int_0^1\frac{dt}{t}\,
t^B(1-t)^{C-B}\frac{(1+zt)^{-A}-1}{A},
\]
\[
g(B,C)=\frac{\Gamma(1+C)}{\Gamma(1+B)\Gamma(1+C-B)}.
\tag{15}
\]

先在 Euler 积分的收敛域减去常数项，再延拓到零点附近，有

\[
{}_2F_1(A,B;1+C;-z)=1+AB\,g(B,C)\mathscr P(A,B,C;z).
\tag{16}
\]

这个源的指标是独立的对数插入：

\[
\left.\partial_A^i\partial_B^j\partial_C^k\mathscr P\right|_0
=\frac{(-1)^{i+1}}{i+1}
\int_0^1\frac{dt}{t}
\log^{i+1}(1+zt)
\log^j\!\frac{t}{1-t}\,
\log^k(1-t).
\tag{17}
\]

基值为 Li₂(−z)。与旧符号的精确关系为
V_P(a,b;z)=a·𝒫(a,a,b;z)；这里的 𝒫 是 Gauss/Euler 族的参数释放，不能与双 J 的 P/Q/R 切片混淆。

它的一次 Mellin 扩展同样已求值：

\[
\widehat{\mathscr P}(A,B,C,D;z)
=\int_0^1dv\,v^{D-1}\mathscr P(A,B,C;zv)
=\frac{{}_3F_2(A,B,D;1+C,1+D;-z)-1}
{AB D\,g(B,C)}.
\tag{18}
\]

零参数取可去极限。式(11)是式(18)在 D=0 的一类切片，式(12)使用 D=λ。因而整个独立端点解由这一 Euler 族、一次 Mellin 扩展及有限乘积组成。这里只证明了该表示足够，没有证明它们在任意函数系数域下最小。

标准公式约定参考 [DLMF Euler 积分](https://dlmf.nist.gov/15.6.E1) 与 [广义超几何积分](https://dlmf.nist.gov/16.5.E2)；式(15)的 g 已包含普通 ₂F₁ 相对 DLMF 正则化函数的 Gamma 因子。

## 5. 在生成层面合并镜像，得到固定三指标矩

现在对两个端点分别分配 z^(a/2)、z^(b/2)，定义对称延拓

\[
\mathcal E(\lambda,a,b;z)=\frac12
\{z^{(a+b)/2}E_\lambda(a,b;z)
+z^{-(a+b)/2}E_\lambda(a,b;1/z)\}.
\tag{19}
\]

这是一项明确选定的自由参数归一化；它在物理对角线上恰好恢复原 crossing 指数，不宣称所有离开物理切片的归一化都唯一。式(14)与式(3)已给出其显式求值。

利用 kλ(1−t,1−s;1/z)=kλ(t,s)，可在同一原单位方形内写出所有指标。记

\[
\ell=\log z,\qquad
p_m(t;\ell)=
\frac{(\ell/2+\log t)^{m+1}-(\ell/2)^{m+1}}{m+1}.
\]

因为 ∂ₐᵐ[z^(a/2)uₐ(t)]|₀=−p_m(t;ℓ)，有

\[
\mathcal E(\lambda,a,b;z)
=\sum_{p,m,n\ge0}\frac{\lambda^p a^m b^n}{p!m!n!}
\mathcal E_{p,m,n}(z),
\tag{20}
\]
\[
\boxed{
\mathcal E_{p,m,n}
=\frac12\int_0^1dt\,ds\,
k_0(t,s)\log^p q_z(t,s)
\{p_m(t;\ell)p_n(s;\ell)
+p_m(1-t;-\ell)p_n(1-s;-\ell)\}.
}
\tag{21}
\]

每个指标都有固定意义：p 是传播子指数导数；m、n 是两个端点响应的导数。交换端点给出 E_{p,m,n}=E_{p,n,m}；整个矩独立满足 z↔1/z 偶性。自然权重为 p+m+n+4。运动学位于矩定义内部，不需要逐阶在矩外展开 logρ。

式(21)保留两个积分变量；式(14)已给出它的生成函数求值。这里没有把“存在二维矩表示”解释成“二维不可约”。

## 6. 完整物理 Laurent 系数的统一规则

令 Z=e^(−γ_E ε)Γ(1−ε)、c₃=cosh(3εℓ/2)。原来源桥的归一化准确变成

\[
\boxed{
I_{36}=I_{39}
=-\frac{Z^3}{32\epsilon^2c_3}\mathcal E(\epsilon,\epsilon,\epsilon;z).
}
\tag{22}
\]

记

\[
S_q=\sum_{p+m+n=q}\frac{\mathcal E_{p,m,n}}{p!m!n!}.
\]

前四个组合明确为

\[
\begin{aligned}
S_0&=\mathcal E_{0,0,0},\\
S_1&=\mathcal E_{1,0,0}+2\mathcal E_{0,1,0},\\
S_2&=\tfrac12\mathcal E_{2,0,0}
+2\mathcal E_{1,1,0}+\mathcal E_{0,2,0}+\mathcal E_{0,1,1},\\
S_3&=\tfrac16\mathcal E_{3,0,0}
+\mathcal E_{2,1,0}+\mathcal E_{1,2,0}+\mathcal E_{1,1,1}
+\tfrac13\mathcal E_{0,3,0}+\mathcal E_{0,2,1}.
\end{aligned}
\tag{23}
\]

因此最终展开不是一个“取系数”黑箱：

\[
\begin{aligned}
I_{36}=-\frac1{32}\bigg[
&\frac{S_0}{\epsilon^2}+\frac{S_1}{\epsilon}
+S_2+p_2S_0\\
&+\epsilon(S_3+p_2S_1+\zeta_3S_0)\\
&+\epsilon^2(S_4+p_2S_2+\zeta_3S_1+p_4S_0)+O(\epsilon^3)
\bigg].
\end{aligned}
\tag{24}
\]

以 K=2M_{0,1}+6ζ₂=−ℓ² 作纯代数缩写，

\[
p_2=\frac32\zeta_2+\frac98K
=\frac94M_{0,1}+\frac{33}4\zeta_2,
\]
\[
p_4=\frac{57}{40}\zeta_2^2+\frac{27}{16}\zeta_2K+\frac{135}{128}K^2.
\tag{25}
\]

所有外部系数都可写成有理数、ζ 与已有 M 的多项式。K 不是新主函数。任意更高阶由 Z³/c₃ 的固定展开与 S_q 的有限和给出。

原两项极点提供精确低阶连接：

\[
S_0=\frac14M_{2,1}-\frac12\zeta_2M_{0,1},
\]
\[
S_1=\frac12M_{0,1}M_{1,1}-4\zeta_3M_{0,1}
+\frac72\zeta_2M_{1,1}-\frac1{24}M_{1,3}
+\frac14M_{3,1}-\frac34N_{1,2}.
\tag{26}
\]

这些关系由式(22)与既有解析极点匹配得出，本次又用原一维 M/N 定义与式(21)进行数值检查。不是新增的低阶求积或新的形式证明。

## 7. 现在可以明确检验哪几个候选矩

首次超出既定 M/N 空间的是有限项，故应对四个权重六候选

\[
\mathcal E_{2,0,0},\quad
\mathcal E_{1,1,0},\quad
\mathcal E_{0,2,0},\quad
\mathcal E_{0,1,1}
\tag{27}
\]

分别建立 Hlog/MN 字典，再处理它们的物理组合 S₂。不能仅凭“双端”外观就断言只有 E_{0,1,1} 是新方向。

独立参数还暴露了物理对角线看不到的组合。令 a=u+δ、b=u−δ，则

\[
\left.\partial_\delta^2\mathcal E(\lambda,u+\delta,u-\delta)\right|_{\lambda=u=\delta=0}
=2(\mathcal E_{0,2,0}-\mathcal E_{0,1,1}).
\tag{28}
\]

物理 S₂ 只取两者的和，而不同图的端点操作可能取到差。这解释了释放参数对于跨图比较的价值；尚未断言该差属于或超出 M/N。

留出检验可冻结为权重九规则

\[
[\epsilon^3]I_{36}
=-\frac1{32}
[S_5+p_2S_3+\zeta_3S_2+p_4S_1+
(\tfrac35\zeta_5+p_2\zeta_3)S_0].
\tag{29}
\]

本轮未独立求其 canonical Hlog 向量，也未把该规则当作完成一个新物理阶。下一可执行项是式(27)的精确商空间约化及式(28)的跨图作用，再独立核验式(29)。

## 8. 验证范围

[可重放检查与记录](../../../../../workflow/checks/single-opposite-free-endpoints-20261016/README.md)：

- SymPy 54 项精确恒等式、13 项错误式控制；包括一般 Gauss 变换、四个原核多项式、三角端点系统、两种原函数逐级数系数与端点导数。
- 实际 Wolfbook 独立 11 项通过，其中四项为错误源符号控制。
- 六组非整数参数的原二维 F/E 共 12 项比较，最大绝对差 5.96×10⁻¹⁴；另六项非齐次方程数值残差、六项遗漏边界常数控制。
- 六项一维 M/N 对照与三项三阶 Taylor 截断缩放检查。
- 两个 z 点的完整物理 I 从 ε⁻² 到 ε² 共 10 项与已冻结 Hlog 向量比较，最大绝对差 6.95×10⁻¹¹。
- 非等上参数 Euler/Mellin 的八项积分检查及两项小 z 边界缩放。

误差为实际节点精化与浮点比较，不是区间证书。原解析关系按上文完整参数、域和边界推导；CAS 检查不单独充当原任务证明。Wilson 数学 MCP 本轮不可调用，没有修改 Wilson/Lean，也没有形式接纳状态。

本次完成的是“独立指数源的求值及其指标化”，不是最小新主积分基底、跨全部 6e 的闭合证明或二维不可约性证明。
