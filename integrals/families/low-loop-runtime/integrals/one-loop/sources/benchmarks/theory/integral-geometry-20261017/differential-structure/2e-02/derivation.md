# 2e-02：从 J 到 M 与 I

2e-02 是一个对边交换积分。我们要弄清两件事：它的 Laurent 系数为什么能用 \(M\) 表示，循环求和在其中做了什么。先从原始 \(J\) 直接得到 \(M\) 的生成函数，再恢复完整 \(I\)；随后推导微分方程和边界，用它们解释系数递推。这样可以区分“函数怎样随 \(z\) 变化”和“结果怎样压缩成主积分”。

全文取 \(z>0\)，幂函数使用实对数，\(\ell=\log z\)。积分变换和分部积分先在 \(0<\Re\epsilon<1/4\) 建立，再延拓到 \(\epsilon=0\) 的邻域。固定正 \(z\) 时，原二维积分在 \(\Re\epsilon>-1\) 收敛，零点附近的对数展开可逐项积分。

## 1. 原始 J

[原输入](../../../../2ne/workspaces/2e-02/derivation.md)的外到内次序是 \(s_{21}\to s_{41}\)。记 \(u=s_{21},v=s_{41}\)，定义
\[
D_z(u,v)=(1-u)(1-v)+zuv.
\]
于是
\[
\boxed{
J(z,\epsilon)=-(1+z)\int_0^1du\int_0^1dv\,
D_z(u,v)^{\epsilon-1}.
}
\tag{1}
\]

此时分析的是未作循环求和的 \(J\)。分母只在两个角点为零，例如 \(u=0,v=1\) 附近有 \(D_z\asymp (1-v)+zu\)。局部积分形如 \(\int_0^\delta t^{\Re\epsilon}\,dt\)，因此在 \(\Re\epsilon>-1\) 绝对收敛。这也允许我们在共同收敛域先积 \(u\)。

固定 \(v\)，令 \(a=zv,d=1-v\)，则 \(D_z=d+(a-d)u\)，所以
\[
\int_0^1D_z^{\epsilon-1}\,du
=\frac{a^\epsilon-d^\epsilon}{\epsilon(a-d)}.
\]
当 \(a=d\) 或 \(\epsilon=0\) 时取连续极限，右边没有真实的内部极点。

再令 \(\beta=1-v\)、\(r=z/(1+z)\)，有 \(a-d=-(1+z)(\beta-r)\)。式(1)变成
\[
\boxed{
J(z,\epsilon)=\frac1\epsilon\int_0^1
\frac{z^\epsilon(1-\beta)^\epsilon-\beta^\epsilon}
{\beta-r}\,d\beta.
}
\tag{2}
\]
原来的二维分母已经化成一个移动极点 \(\beta=r\)，而分子在此处同时为零。这里始终使用完整分子，不将它拆成两个未经规定的奇异积分。

## 2. 生成 M

式(2)中的两个端点幂，可以组织成一个对称的乘积和一个比值。定义
\[
W=\log[\beta(1-\beta)],\qquad
L=\ell+\log\frac{1-\beta}{\beta}.
\]
直接计算得到
\[
z^\epsilon(1-\beta)^\epsilon-\beta^\epsilon
=2z^{\epsilon/2}e^{\epsilon W/2}
\sinh\!\left(\frac{\epsilon L}{2}\right).
\]
因此将整体尺度幂提出，定义
\[
H_\epsilon(z)=z^{-\epsilon/2}J(z,\epsilon).
\]
它只是原 \(J\) 的归一化写法，还不是完整 \(I\)。代入式(2)，得到
\[
\boxed{
H_\epsilon(z)=\int_0^1
\frac{e^{\epsilon W/2}}{\beta-r}
\frac{\sinh(\epsilon L/2)}{\epsilon/2}\,d\beta.
}
\tag{3}
\]
这就是本例的 \(M\) 生成函数；其结构来自一次实际积分和变量替换，尚未使用循环求和或微分方程。

沿用已有主积分矩的定义
\[
M_{k,j}(z)=\int_0^1\frac{W^kL^j}{\beta-r}\,d\beta,
\qquad k\ge0,\quad j\ge1.
\tag{4}
\]
由于 \(L(r)=0\)，这些矩在 \(\beta=r\) 都取普通连续极限。展开式(3)中的指数和双曲正弦，便有
\[
\boxed{
H_\epsilon(z)=
\sum_{k,j\ge0}
\frac{\epsilon^{k+2j}}
{2^{k+2j}k!(2j+1)!}\,
M_{k,2j+1}(z).
}
\tag{5}
\]

若写 \(H_\epsilon=\sum_{n\ge0}h_n\epsilon^n\)，前三项是
\[
h_0=M_{0,1},\qquad
h_1=\frac12M_{1,1},\qquad
h_2=\frac18M_{2,1}+\frac1{24}M_{0,3}.
\tag{6}
\]
每个固定阶数只有有限个矩，外部系数与 \(z\) 无关。这里的压缩是：原积分中的运动学依赖，被组织进预先定义的 \(M_{k,j}(z)\)。式(5)不宣称这些矩彼此线性独立。

## 3. 循环求和

先从原分母检查反射，而不预先假设 \(J\) 对称：
\[
D_{1/z}(1-u,1-v)=\frac{D_z(u,v)}z.
\]
同时替换两个积分变量，式(1)给出
\[
J(1/z,\epsilon)=z^{-\epsilon}J(z,\epsilon),
\qquad
H_\epsilon(1/z)=H_\epsilon(z).
\tag{7}
\]
所以本例的加权循环组合恰好是
\[
\boxed{
J(z,\epsilon)+z^\epsilon J(1/z,\epsilon)=2J(z,\epsilon).
}
\tag{8}
\]
两项是同一个函数的两次贡献；没有在这一步消去一个额外的端点积分。

现在才引入完整 \(I\)。令
\[
Z(\epsilon)=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon).
\]
按照本例已经固定的[单一 C/I 归一化](../../../../2e-02/current-normalization.md)，循环组合与提取物理尺度后的结果满足
\[
\begin{aligned}
I(z,\epsilon)
&=-\frac{Z(\epsilon)}2\,
\frac{J(z,\epsilon)+z^\epsilon J(1/z,\epsilon)}
{1+z^\epsilon}\\
&=-\frac{Z(\epsilon)}{1+z^\epsilon}J(z,\epsilon).
\end{aligned}
\tag{9}
\]
这里的 \(1+z^\epsilon\) 来自所提取的两种物理尺度；不是为了简化答案而另行引入的因子。利用 \(J=z^{\epsilon/2}H_\epsilon\)，得到
\[
\boxed{
I(z,\epsilon)=
-\frac{Z(\epsilon)}{2\cosh(\epsilon\ell/2)}H_\epsilon(z).
}
\tag{10}
\]

至此已经分清三个对象：\(J\) 是原积分，\(H_\epsilon\) 是剥去整体尺度幂后的同一个单核，\(I\) 则包含原约定的循环组合及归一化。2e-02 的 \(M\) 生成结构出现在式(3)，循环和在式(8)仅使原贡献加倍。

## 4. 微分方程

现在研究这些函数随 \(z\) 的变化。为了直接在原积分内求导，暂记
\[
F_\epsilon(z)=\int_0^1du\,dv\,D_z^{\epsilon-1},
\qquad J=-(1+z)F_\epsilon,\qquad
\theta=z\frac{d}{dz}.
\]
求导时固定 \(\epsilon\)。对于 \(g=D_z^{\epsilon-1}\)，有
\[
\partial_uD_z=(1+z)v-1,\quad
\partial_vD_z=(1+z)u-1,\quad
\theta D_z=zuv.
\]
它们满足
\[
zv\,\partial_vD_z+(1-u)\partial_uD_z
=(1+z)zuv-D_z.
\]
将这个恒等式代入乘积求导公式，得到
\[
\big[(1+z)\theta+z-\epsilon\big]g
=\partial_v(zvg)+\partial_u[(1-u)g].
\tag{11}
\]
这一步把对运动学变量的求导变成积分变量的全导数。

在 \(\Re\epsilon>0\) 时，积分后的两个非零边界为
\[
\int_0^1z(zu)^{\epsilon-1}\,du=\frac{z^\epsilon}{\epsilon},
\qquad
-\int_0^1(1-v)^{\epsilon-1}\,dv=-\frac1\epsilon.
\]
另外两个边界为零；奇角可先截去，再在此共同收敛域取极限。因此
\[
\boxed{
\big[(1+z)\theta+z-\epsilon\big]F_\epsilon
=\frac{z^\epsilon-1}{\epsilon}.
}
\tag{12}
\]
边界项决定了右边的源，不能在求齐次方程时把它丢掉。

由 \(J=-(1+z)F_\epsilon\) 和 \(H_\epsilon=z^{-\epsilon/2}J\)，依次得到
\[
\boxed{
\begin{aligned}
\theta J
&=\frac{\epsilon}{1+z}J-\frac{z^\epsilon-1}{\epsilon},\\
\theta H_\epsilon
&=\epsilon\kappa(z)H_\epsilon
-\frac{2\sinh(\epsilon\ell/2)}{\epsilon},
\qquad
\kappa(z)=\frac{1-z}{2(1+z)}.
\end{aligned}
}
\tag{13}
\]
第二条方程的未知函数是 \(H_\epsilon\)，第一条则直接作用于原 \(J\)。两条方程由明确的归一化连接。

若直接对完整 \(I\) 求导，还必须对式(10)的双曲余弦因子求导，结果为
\[
\theta I
=\epsilon\left[\kappa-\frac12\tanh(\epsilon\ell/2)\right]I
+\frac{Z(\epsilon)}{\epsilon}\tanh(\epsilon\ell/2).
\tag{14}
\]
使用 \(H_\epsilon\) 能让递推系数更简单，但并不改变已选定的原积分。式(13)描述函数的演化；式(3)–(5)才给出它与 \(M\) 的显式连接。

## 5. 边界与递推

一阶方程还需要一个函数值。式(13)中 \(H_\epsilon\) 的齐次解是
\[
C(\epsilon)\left[\frac{z}{(1+z)^2}\right]^{\epsilon/2}.
\]
它也满足反射对称，所以仅要求 \(z\leftrightarrow1/z\) 对称不能固定 \(C(\epsilon)\)。在 \(\Re\epsilon>0\) 时它还在 \(z\to0\) 趋零，因此只指定这个端点的极限也不够。

选择非奇异点 \(z=1\)，从原积分积掉 \(u\)，再令 \(t=2v-1\)，利用被积函数的偶性，得到
\[
H_\epsilon(1)
=-\frac2\epsilon\int_0^1
\frac{[(1+t)/2]^\epsilon-[(1-t)/2]^\epsilon}{t}\,dt.
\tag{15}
\]
这是原积分给出的边界值，而非对微分方程自由选择的常数。

展开 \(H_\epsilon=\sum h_n\epsilon^n\)，式(13)逐阶给出
\[
\boxed{
\theta h_n
=\kappa h_{n-1}
-\mathbf1_{n\ \text{为偶数}}\,
\frac{\ell^{n+1}}{2^n(n+1)!},
\qquad h_{-1}=0.
}
\tag{16}
\]
每阶所需的边界由式(15)直接展开：
\[
h_n(1)=-\frac2{(n+1)!}\int_0^1
\frac{\log^{n+1}[(1+t)/2]-\log^{n+1}[(1-t)/2]}{t}\,dt.
\tag{17}
\]

先完整算出零阶。式(16)给 \(\theta h_0=-\ell\)，所以
\[
h_0(z)=-\frac{\ell^2}{2}+h_0(1).
\]
另一方面，
\[
h_0(1)
=-2\int_0^1\frac{\log[(1+t)/(1-t)]}{t}\,dt
=-4\sum_{m\ge0}\frac1{(2m+1)^2}
=-\frac{\pi^2}{2}.
\]
结合式(6)，得到
\[
\boxed{
M_{0,1}(z)=h_0(z)=-\frac{\ell^2+\pi^2}{2}.
}
\tag{18}
\]
这一次计算同时展示了方程和边界各自的作用：方程固定 \(\ell^2\) 项，原积分的基点值固定 \(\pi^2\) 项。

接下来的两阶满足
\[
\theta h_1=\kappa h_0,\qquad
\theta h_2=\kappa h_1-\frac{\ell^3}{24}.
\]
其常数仍由式(17)固定，而式(6)已经给出同一物理解的 \(M\) 表示。因此我们既知道怎样逐阶积分，也知道这些特定结果怎样组织成主积分。

## 6. 展开 I

式(10)仍显式含有 \(\ell\)，但这本身不意味着超出 \(M\) 空间。由式(18)，对任意整数 \(m\ge0\) 都有
\[
\ell^{2m}=\big[-2M_{0,1}-\pi^2\big]^m.
\tag{19}
\]
因此双曲余弦倒数的每阶系数都能化成 \(M_{0,1}\) 与常数的多项式。这里使用的是允许 \(M\) 乘积的系数空间；将乘积进一步压成某个有限线性基是另一个约化步骤。

为看清归一化的作用，先写
\[
\begin{aligned}
J(z,\epsilon)
&=M_{0,1}
+\frac{\epsilon}{2}\big[M_{1,1}+\ell M_{0,1}\big]
+O(\epsilon^2),\\
Z(\epsilon)&=1+\frac{\zeta(2)}2\epsilon^2+O(\epsilon^3),\\
\frac1{\cosh(\epsilon\ell/2)}
&=1-\frac{\epsilon^2\ell^2}{8}+O(\epsilon^4).
\end{aligned}
\]
原 \(J\) 的一次项带有整体尺度幂产生的 \(\ell M_{0,1}\)。完整归一化移除了这部分尺度依赖；由式(6)、(10)可得
\[
\begin{aligned}
I(z,\epsilon)={}&-\frac12M_{0,1}
-\frac{\epsilon}{4}M_{1,1}\\
&+\epsilon^2\left[
-\frac1{16}M_{2,1}-\frac1{48}M_{0,3}
-\frac{\zeta(2)}4M_{0,1}
+\frac{\ell^2}{16}M_{0,1}\right]
+O(\epsilon^3).
\end{aligned}
\]
再用 \(\ell^2=-2M_{0,1}-6\zeta(2)\)，得到外部不含 \(z\) 或 \(\log z\) 的表示：
\[
\boxed{
\begin{aligned}
I(z,\epsilon)={}&-\frac12M_{0,1}
-\frac{\epsilon}{4}M_{1,1}\\
&-\epsilon^2\left[
\frac1{16}M_{2,1}+\frac1{48}M_{0,3}
+\frac18M_{0,1}^{\,2}
+\frac58\zeta(2)M_{0,1}\right]
+O(\epsilon^3).
\end{aligned}
}
\tag{20}
\]
本例在 \(\epsilon=0\) 没有负幂；其 Laurent 展开实际是 Taylor 展开。更高阶也遵循同一路径：式(5)生成 \(H\) 的矩，\(\log Z=\sum_{n\ge2}\zeta(n)\epsilon^n/n\) 生成常数，式(19)吸收剩余偶次对数。在每个固定阶数，这些操作只产生有限个 \(M\) 及其乘积。

因此，2e-02 已经解释了自身的 \(M\) 闭合：它来自原单核的生成函数恒等式及完整归一化。微分方程给出同一物理解的演化与边界，并不单靠系统的分量数证明这种闭合。此例的循环和仅为原贡献加倍；要研究不同端点结构如何在循环和中合并，应接着看 [4e-02/05](../../single-endpoint-differential/derivation.md)。

---

本文整理已有解析结果，沿用[综合微分推导](../derivation.md)及[超几何—生成函数字典](../../two-loop-hypergeometric-dictionary/derivation.md)，没有改变原输入、归一化或结果选择器。本次公式与文档核对见[整理记录](../../../../../workflow/checks/low-loop-differential-20261018/twoe02-writing-20261019/README.md)；既有局部形式证据不等于原积分全链的形式证明。
