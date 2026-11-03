# 4e-04在有序对数矩中的全部已知Laurent系数

日期：2026-09-27。本文直接回答“在新生成族中，各系数怎样表示”。沿用[全ε推导](ordered-beta-generator.md)及唯一物理I归一化，列出ε⁰…ε⁴，所有负幂均为零，余项O(ε⁵)。这是已有结果的显式坐标表示，不是重新求积分或增加已验收覆盖。

结果由两部分组成：**已知一圈系数的乘积＋新有序矩的线性组合**。这与全ε式中的H²/2及双源函数两部分一一对应；不需要M/N拟合。

## 1. 本文实际使用的矩

生成族的积分定义是

$$
\boxed{\displaystyle
f_{m,n}(z)=\frac1{m!\,n!}
\int_0^1\frac{dv}{v}\log^n(1-v)
\int_0^v\frac{du}{u}\log^m(1+zu),
\qquad m\ge1,\ n\ge0.}
\tag{L1}
$$

这里f就是前稿的𝔣。阶乘已经包含在定义内，不能按未归一化的M矩习惯再乘一次。展开真实幂源，有

$$
\mathfrak F(a,b;-z)
=\int_{0<u<v<1}
\frac{[(1+zu)^a-1](1-v)^b}{uv}\,du\,dv
=\sum_{m\ge1,n\ge0}a^m b^n f_{m,n}(z).
\tag{L2}
$$

例如f₁,₀=−Li₃(−z)，f₂,₀含½log²(1+zu)，f₁,₁含log(1+zu)log(1−v)。由内层幂级数和普通Beta矩，可直接重建前稿(B6)，不需要从答案反推此积分定义。

以下采用crossing**平均**：

$$
f^+_{m,n}(z)=\frac{f_{m,n}(z)+f_{m,n}(1/z)}2,\qquad
\ell=\log z,\qquad
\lambda=\log\frac{z}{(1+z)^2}.
\tag{L3}
$$

所有ζ_r均为Riemann zeta常数。正实z固定，所有对数取正实主值。

## 2. 五个显式有序矩组合

A₂…A₆只是缩短下面公式的线性组合，**不是新增基底函数**：

$$
A_2=4f^+_{1,1}-4f^+_{2,0}-4\zeta_4.
\tag{L4}
$$

$$
A_3=4(f^+_{1,2}+f^+_{2,1})-12f^+_{3,0}
+16\zeta_2\zeta_3-24\zeta_5.
\tag{L5}
$$

$$
A_4=4(f^+_{1,3}+f^+_{2,2}+f^+_{3,1})-28f^+_{4,0}
+14\zeta_6-16\zeta_3^2.
\tag{L6}
$$

$$
\begin{aligned}
A_5={}&4(f^+_{1,4}+f^+_{2,3}+f^+_{3,2}+f^+_{4,1})
-60f^+_{5,0}\\
&+64\zeta_2\zeta_5+64\zeta_3\zeta_4-180\zeta_7 .
\end{aligned}
\tag{L7}
$$

$$
\begin{aligned}
A_6={}&4(f^+_{1,5}+f^+_{2,4}+f^+_{3,3}
+f^+_{4,2}+f^+_{5,1})-124f^+_{6,0}\\
&+146\zeta_8-128\zeta_3\zeta_5 .
\end{aligned}
\tag{L8}
$$

这些系数来自全ε组合2Ψ(ε,ε)−Ψ(2ε,0)，不是通过线性拟合得到。其一般r≥2规则为

$$
\begin{aligned}
A_r={}&4\sum_{m=1}^{r-1}f^+_{m,r-m}
+(4-2^{r+1})f^+_{r,0}\\
&+(-1)^r\left[
(2^r-2)(r+1)\zeta_{r+2}
-2^r\sum_{i=1}^{r-1}\zeta_{i+1}\zeta_{r-i+1}\right].
\end{aligned}
\tag{L9}
$$

ζ常数来自全ε函数在参数自变量x=1处的已求值边界。求完整I至ε⁴只需r≤6，即新矩总指标m+n≤6。

## 3. 完整物理I的五个系数

令H(ε,z)=Σh_n(z)εⁿ为已知一圈正规化母核。h_n不是新的二圈函数；§4给出它们在同一生成族及其原函数中的明确表示。

把唯一Euler补偿及Gamma/cosh前因子展开：

$$
\begin{gathered}
p_2=\zeta_2-\frac{\ell^2}{2},\qquad
p_3=\frac{2\zeta_3}{3},\\
p_4=\frac7{10}\zeta_2^2-\frac12\zeta_2\ell^2
+\frac5{24}\ell^4,\\
\frac{e^{-2\gamma_E\epsilon}\Gamma(1-\epsilon)^2}
{\cosh(\epsilon\ell)}
=1+p_2\epsilon^2+p_3\epsilon^3+p_4\epsilon^4+O(\epsilon^5).
\end{gathered}
\tag{L10}
$$

记I₄ₑ₋₀₄=C₀+εC₁+ε²C₂+ε³C₃+ε⁴C₄+O(ε⁵)。逐项结果为：

$$
\boxed{C_0=-\frac18\left(\frac{h_0^2}{2}+A_2\right)
=-\frac{(\ell^2+\pi^2)^2}{192}.}
\tag{L11}
$$

$$
\boxed{C_1=-\frac18\left(h_0h_1+A_3+\lambda A_2\right).}
\tag{L12}
$$

$$
\boxed{\begin{aligned}
C_2=-\frac18\bigg[&
h_0h_2+\frac{h_1^2}{2}+\frac{p_2h_0^2}{2}\\
&+A_4+\lambda A_3+
\left(\frac{\lambda^2}{2}+p_2\right)A_2\bigg].
\end{aligned}}
\tag{L13}
$$

$$
\boxed{\begin{aligned}
C_3=-\frac18\bigg[&
h_0h_3+h_1h_2+p_2h_0h_1+\frac{p_3h_0^2}{2}\\
&+A_5+\lambda A_4+
\left(\frac{\lambda^2}{2}+p_2\right)A_3\\
&+\left(\frac{\lambda^3}{6}+p_2\lambda+p_3\right)A_2\bigg].
\end{aligned}}
\tag{L14}
$$

$$
\boxed{\begin{aligned}
C_4=-\frac18\bigg[&
h_0h_4+h_1h_3+\frac{h_2^2}{2}
+p_2\left(h_0h_2+\frac{h_1^2}{2}\right)
+p_3h_0h_1+\frac{p_4h_0^2}{2}\\
&+A_6+\lambda A_5+
\left(\frac{\lambda^2}{2}+p_2\right)A_4\\
&+\left(\frac{\lambda^3}{6}+p_2\lambda+p_3\right)A_3\\
&+\left(\frac{\lambda^4}{24}+\frac{p_2\lambda^2}{2}
+p_3\lambda+p_4\right)A_2\bigg].
\end{aligned}}
\tag{L15}
$$

为了看清一次幂具体包含哪些矩，(L12)完全代入A₂、A₃后就是

$$
\begin{aligned}
C_1={}&-\frac18h_0h_1
-\frac12(f^+_{1,2}+f^+_{2,1})+\frac32f^+_{3,0}
-2\zeta_2\zeta_3+3\zeta_5\\
&-\frac{\lambda}{2}(f^+_{1,1}-f^+_{2,0}-\zeta_4).
\end{aligned}
\tag{L16}
$$

所有更低幂C₋₁、C₋₂、…均严格为零。零阶中A₂=−h₀²/3，所以有序矩组合最终消去，回到(L11)的初等式；此表示没有宣称各矩都独立。

## 4. 一圈乘积也可以直接用同一生成族写出

定义其Euler导数，不引入M符号：

$$
\begin{gathered}
j_r(z)=\theta_z f_{r,0}(z)
=\frac1{r!}\int_0^1\frac{\log^r(1+zt)}t\,dt,\qquad r\ge1,\\
X_r(z)=j_r(z)+j_r(1/z)-2(-1)^r\zeta_{r+1}.
\end{gathered}
\tag{L17}
$$

由原一圈全ε公式直接得到

$$
\begin{aligned}
h_0={}&-X_1=-\frac{\ell^2+\pi^2}{2},\\
h_1={}&-X_2-\frac{\lambda X_1}{2},\\
h_2={}&-X_3-\frac{\lambda X_2}{2}-\frac{\lambda^2X_1}{8},\\
h_3={}&-X_4-\frac{\lambda X_3}{2}-\frac{\lambda^2X_2}{8}
-\frac{\lambda^3X_1}{48},\\
h_4={}&-X_5-\frac{\lambda X_4}{2}-\frac{\lambda^2X_3}{8}
-\frac{\lambda^3X_2}{48}-\frac{\lambda^4X_1}{384}.
\end{aligned}
\tag{L18}
$$

因此(L4)–(L18)构成自包含的解析表示：有序矩f、同一源的Euler导数j、它们的一圈乘积及显式log/ζ常数。它不要求把二圈答案重新拟合到M/N，也不强求只有f的线性组合而删除已知一圈乘积。

## 5. 验证与定位

[检查与完整机器系数](../../../../workflow/checks/4e04-ordered-coefficients/README.md)直接从全ε式卷积：

- 五个A_r的系数、crossing平均的½及端点ζ常数精确核对。
- (L11)–(L15)均反代到完整Hlog，五阶与上一轮独立生成的完整I完全相同。
- A₂=−h₀²/3及h₀初等式通过；crossing错误加倍、漏掉ε²的Gamma/cosh贡献均被检出。

本轮是表达展开与说明，没有新CAS积分或数值拟合；不改变旧r002/catalog、前阶段证据、ε阶数或21/329覆盖。formal=not_run，最小独立基底问题仍开放。
