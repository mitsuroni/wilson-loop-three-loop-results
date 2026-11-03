# 从 I₆ₑ₋₂₀ 出发：三核积分中的自然双Beta生成函数

日期：2026-09-28。本页讨论[3]族的 **6e-20、6e-21、6e-22**，以完整物理积分 $I_{6e-20}$ 为推导主线。另两图的完整 $I$ 相等，并与它共享下文的两个积分量。

**上一版求出的5F4只是三图共同出现的一个子块，不是其中某张图的完整通解。** 原稿从辅助量直接开始，没有把它接回图编号，因而无法看清进度。本版从原图出发，主线各步左边始终保留 $I_{6e-20}$；辅助求值集中到附录。

记所有未写出的参数为 $(\epsilon,z)$，先取 $z>0$、实小正 $\epsilon$，再在 $\epsilon=0$ 附近延拓。完整物理前因子为

$$
P(\epsilon,z)=\frac{Z(\epsilon)^3}{\cosh(3\epsilon\log z/2)},
\qquad Z(\epsilon)=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon).
$$

原稿原字节保存在[r001](natural-beta-source-r001.md)；本次[连接与渲染检查](../../../../workflow/checks/three-kernel-i-chain/README.md)不重算既有Laurent系数。原[5F4求值检查](../../../../workflow/checks/natural-beta-source/README.md)继续对应r001，formal=not_run。

## 1. 先说清楚：是哪张图，原积分是什么

使用已由[原图到母核的连接](mother-map.md)固定的核与排序域：

$$
\begin{gathered}
k_\epsilon(t,s)=-(1+z)z^{-\epsilon/2}
\big[(1-t)(1-s)+zts\big]^{\epsilon-1},\\
\Omega_3=\{0<t_1<t_2<t_3<1,\quad0<s_1<s_2<s_3<1\}.
\end{gathered}
$$

6e-20对应两边端点同序配对：

$$
\boxed{\displaystyle
I_{6e-20}=-\frac P{16}
\int_{\Omega_3}k_\epsilon(t_1,s_1)k_\epsilon(t_2,s_2)
k_\epsilon(t_3,s_3)\,d^3t\,d^3s.}
\tag{I1}
$$

前稿把这个六维裸积分记作 $A$；这里直接保留图编号。另一种独立物理结果是

$$
I_{6e-21}=I_{6e-22}
=-\frac P{64}\int_{\Omega_3}
\big[k_{11}k_{23}k_{32}+k_{12}k_{21}k_{33}\big]\,d^3t\,d^3s,
\qquad k_{ij}=k_\epsilon(t_i,s_j).
$$

它等于前稿的 $-PB/32$，其中 $B$ 已含两个排列的平均，不能再次漏乘或多乘一个2。上述连接沿用[三核全ε稿(T3)](three-kernel-generators.md)。

## 2. I₆ₑ₋₂₀：从六维积到二维，再分出两个待处理量

固定中间传播子的端点 $t=t_2,s=s_2$。它左下、右上的两对变量分别形成矩形，定义

$$
\begin{aligned}
F(t,s)&=\int_0^t dx\int_0^sdy\,k_\epsilon(x,y),&
f(t)&=F(t,1),& H&=F(1,1),\\
U(t,s)&=\int_t^1dx\int_s^1dy\,k_\epsilon(x,y)
=H-f(t)-f(s)+F(t,s).
\end{aligned}
$$

把这四个变量实际积掉，完整图变为

$$
\boxed{\displaystyle
I_{6e-20}=-\frac P{16}
\int_0^1dt\int_0^1ds\,k_\epsilon(t,s)F(t,s)U(t,s).}
\tag{I2}
$$

**进度：** 原图六维降为二维。$F,U$ 是原图排序上下限产生的矩形原函数。

利用矩形边界与分部积分，定义两个需要继续处理的量：

$$
\begin{aligned}
\mathcal W(\epsilon,z)
&=\int_0^1dt\int_0^1ds\,k_\epsilon(t,s)f(t)f(s)
=\int_0^1dt\int_0^1ds\,F(t,s)f'(t)f'(s),\\
\mathcal L(\epsilon,z)
&=\int_0^1dt\int_0^1ds\,F(t,s)\partial_tF(t,s)\partial_sF(t,s),\\
\Lambda(\epsilon,z)&=\mathcal L(\epsilon,z)+\mathcal L(\epsilon,1/z).
\end{aligned}
$$

这里 $\Lambda$ 是全ε函数，两项相加不除以2；[系数稿](three-kernel-eps2.md)中的 $\Lambda_q=[\epsilon^q]\Lambda$ 才是各阶系数。由该稿(E1)的边界恒等式，(I2)成为

$$
\boxed{\displaystyle
I_{6e-20}=\frac P{16}\left(\mathcal W+\Lambda-\frac{H^3}{2}\right).}
\tag{I3}
$$

同时，另两图是

$$
\boxed{\displaystyle
I_{6e-21}=I_{6e-22}
=\frac P{32}\left(\frac{H^3}{3}-\mathcal W-\frac{\Lambda}{2}\right).}
$$

这就说明了辅助量与原图的关系：**三图共用 $\mathcal W,\Lambda$；上一版双Beta计算只推进了 $\mathcal W$ 的一部分，尚未处理完 $\mathcal W$，也没有处理完 $\Lambda$。**

还可以直接从物理图提取这两个量：

$$
\begin{aligned}
I_{6e-20}+4I_{6e-21}
&=\frac P{16}\left(\frac{H^3}{6}-\mathcal W\right),\\
I_{6e-20}+2I_{6e-21}
&=\frac P{32}\left(\Lambda-\frac{H^3}{3}\right).
\end{aligned}
$$

第一行消去 $\Lambda$，第二行消去 $\mathcal W$。因此上一版的研究可以定位为：**求解组合 $I_{6e-20}+4I_{6e-21}$ 中的一个子块**，所得结果随后按(I3)回到各图。

## 3. I₆ₑ₋₂₀：把 W 换成实际的端点幂积分

令

$$
\begin{gathered}
p=(1+z)t-1,\qquad q=(1+z)s-1,\qquad
\rho=\frac{z}{(1+z)^2},\\
a_\epsilon(p)=(1+p)^\epsilon-(1-p/z)^\epsilon,\\
J_\epsilon(y)=\int_0^y\frac{(1-v)^\epsilon-1}{v}\,dv,\\
\mathcal J(p,q)=J_\epsilon(-pq/z)-J_\epsilon(p/z)
-J_\epsilon(q/z)+J_\epsilon(-1/z).
\end{gathered}
$$

[已求值的矩形原函数](three-kernel-generators.md)给出

$$
F=-\frac{\rho^{\epsilon/2}}{\epsilon}\mathcal J,
\qquad
df=-\frac{\rho^{\epsilon/2}}{\epsilon}\frac{a_\epsilon(p)}p\,dp.
$$

把它们代入(I3)中 $\mathcal W=\int F\,df(t)\,df(s)$ 的位置，得到

$$
\boxed{\displaystyle
I_{6e-20}=\frac P{16}\left[
\Lambda-\frac{H^3}{2}
-\frac{\rho^{3\epsilon/2}}{\epsilon^3}
\int_{-1}^{z}dp\int_{-1}^{z}dq\,
\frac{a_\epsilon(p)a_\epsilon(q)}{pq}\mathcal J(p,q)
\right].}
\tag{I4}
$$

**进度：** 等号左边仍是完整6e-20；右边只展开了其中的 $\mathcal W$，$\Lambda$ 完整保留。前因子、四个角的边界项都未丢掉。

将 $\mathcal J$ 的后三项分别积分，只产生两个一维边界量：

$$
T_0=\int_{-1}^{z}\frac{a_\epsilon(p)}p\,dp,\qquad
T_J=\int_{-1}^{z}\frac{a_\epsilon(p)}pJ_\epsilon(p/z)\,dp.
$$

它们的全ε求值已在[原稿(T12)](three-kernel-generators.md)接到一圈及4e-04原函数。因此

$$
\boxed{\displaystyle
I_{6e-20}=\frac P{16}\left[
\Lambda-\frac{H^3}{2}
-\frac{\rho^{3\epsilon/2}}{\epsilon^3}
\left\{
\int_{-1}^{z}dp\int_{-1}^{z}dq\,
\frac{a_\epsilon(p)a_\epsilon(q)}{pq}J_\epsilon(-pq/z)
-2T_0T_J+J_\epsilon(-1/z)T_0^2
\right\}\right].}
\tag{I5}
$$

**进度：** 原四角表达只剩第一项仍然耦合两个变量。本次5F4计算发生在这个二重积分内部。

## 4. I₆ₑ₋₂₀：准确标出已算的5F4和没算完的部分

只看(I5)二重积分的 $p>0,q>0$ 区间，取 $p=zu,q=zv$。这一块是

$$
\int_0^1\frac{du}{u}\int_0^1\frac{dv}{v}
\big[(1+zu)^\epsilon-(1-u)^\epsilon\big]
\big[(1+zv)^\epsilon-(1-v)^\epsilon\big]J_\epsilon(-zuv).
$$

乘开两个方括号有四项。上一版只完成了其中的“负项乘负项”：

$$
\int_0^1\frac{du}{u}\int_0^1\frac{dv}{v}
(1-u)^\epsilon(1-v)^\epsilon J_\epsilon(-zuv)
=\mathfrak B(\epsilon,\epsilon,\epsilon;-z).
$$

把(I5)耦合积分的其余部分暂记作 $\mathcal G_{\rm rest}$。为避免用名字掩盖未知量，它的完整定义是

$$
\begin{aligned}
\mathcal G_{\rm rest}={}&
\iint_{[-1,z]^2\setminus(0,z)^2}
\frac{a_\epsilon(p)a_\epsilon(q)}{pq}
J_\epsilon(-pq/z)\,dp\,dq\\
&+\int_0^1\frac{du}{u}\int_0^1\frac{dv}{v}
\Big[(1+zu)^\epsilon(1+zv)^\epsilon
-(1+zu)^\epsilon(1-v)^\epsilon\\
&\hspace{43mm}-(1-u)^\epsilon(1+zv)^\epsilon\Big]
J_\epsilon(-zuv).
\end{aligned}
$$

第一行是另外三个符号区间；后两行是正正区间中尚未求完的三个项。此处使用原 $p,q$ 测度，未省略负区间换元的方向因子。

附录的双Beta求值代回原图，现在得到的完整等式是

$$
\boxed{\begin{aligned}
I_{6e-20}=\frac P{16}\Bigg[
&\Lambda-\frac{H^3}{2}
-\frac{\rho^{3\epsilon/2}}{\epsilon^3}
\Bigg\{
\underbrace{\frac{\epsilon z}{(1+\epsilon)^2}
{}_5F_4\!\left(
\begin{matrix}1,1,1,1,1-\epsilon\\2,2,2+\epsilon,2+\epsilon\end{matrix};-z
\right)}_{\text{本次已求值的双Beta子块}}\\
&\hspace{18mm}+\mathcal G_{\rm rest}
-2T_0T_J+J_\epsilon(-1/z)T_0^2
\Bigg\}\Bigg].
\end{aligned}}
\tag{I6}
$$

**进度：** (I6)没有求完整张6e-20；它准确记录“哪一块已经积完、哪些仍留下”。6e-21/22由(I3)后的共同表示恢复，不能把这个5F4单独当作某一张图。

还有一个必须整体保留的消去：双Beta子块本身从 $\epsilon^1$ 开始，而未拆分的两个端点幂差乘 $J_\epsilon$ 从 $\epsilon^3$ 开始。拆开的项相加后才消去低阶项。因此不能把单个5F4乘上 $\epsilon^{-3}$ 所显现的极点，解释为原图的极点；原三图的负Laurent幂均为零。

## 5. 接下来具体求什么

后续沿(I6)继续求 $\mathcal G_{\rm rest}$ 的伴随端点幂块，并处理 $\Lambda$ 的矩形耦合生成函数。每次求得辅助结果，都代回(I6)或(I3)，保留尚未完成的项；全部完成后，才能给出各图的紧凑全ε通式及自然系数族。

既有三图 $\epsilon^0,\epsilon^1,\epsilon^2$ 的[完整解析解](three-kernel-eps2.md)是另一项已经完成的成果，可作事后回归。这里正在研究的是它们更自然的全ε生成表示。本次改写没有新增系数、改变56/329覆盖或提高形式状态，也没有提交Git。

## 附录A：只求(I6)中的双Beta子块

以下等号左边改为辅助函数 $\mathfrak B$，因为这里只完成(I6)花括号内标出的那一项。读完本附录，应回到(I6)，其余项仍在那里。

将真实出现的三个指数暂时独立化：

$$
\boxed{\displaystyle
\mathfrak B(a,b,c;x)=
\int_0^1\frac{du}{u}\int_0^1\frac{dv}{v}
(1-u)^a(1-v)^bJ_c(xuv),\qquad
J_c(y)=\int_0^y\frac{(1-t)^c-1}{t}\,dt.}
\tag{B1}
$$

物理图在这里取 $a=b=c=\epsilon,\ x=-z$。先在 $\Re a,\Re b>-1$、参数零点邻域、$|x|<1$ 推导，再沿负实轴延拓；$J_c(xuv)=O(uv)$ 保证两个零端点可积。这些参数来自原积分的幂，不从旧M/N系数拟合。

### A.1 两个Beta矩给出5F4

直接展开乘积自变量，交换求和与积分后，两个积分各给一个普通Beta函数：

$$
\begin{aligned}
\mathfrak B
&=\sum_{n\ge1}\frac{(-c)_n x^n}{n!\,n}
\mathrm B(n,1+a)\mathrm B(n,1+b)\\
&=-\frac{cx}{(1+a)(1+b)}
{}_5F_4\!\left(
\begin{matrix}1,1,1,1,1-c\\2,2,2+a,2+b\end{matrix};x\right).
\end{aligned}
\tag{B2}
$$

第一行中每个Beta矩都已经求值。分母中的额外n来自J的原函数；省去它得到另一个函数，不能混用。对负实x的延拓沿原积分的正底数支处理；本轮数值控制在收敛圆内进行，不宣称已验证所有延拓点。

例如c为1时，原J等于负y，两边直接退为

$$
\mathfrak B(a,b,1;x)=-\frac{x}{(1+a)(1+b)}.
$$

这与上参数1减c为零时的超几何终止完全一致。

### A.2 辅助函数的方程和边界

令 $\theta=x\partial_x$，并令 $V=\theta^2\mathfrak B$。若 $V=\sum_{n\ge1}v_nx^n$，相邻级数系数之比直接给出

$$
(n+a)(n+b)v_n=(n-1)(n-1-c)v_{n-1}-c\,\delta_{n,1},\qquad v_0=0.
$$

所以在整个收敛圆内有

$$
\boxed{\displaystyle
\big[(1-x)\theta^2+(a+b+cx)\theta+ab\big]V=-cx,
\qquad V=\theta^2\mathfrak B.}
\tag{B3}
$$

选择零点解析解，并要求

$$
V(0)=0,\qquad
\mathfrak B(0)=0,\qquad
\theta\mathfrak B(0)=0.
$$

参数在零附近时，这些条件及方程固定所需幂级数，没有由旧答案拟合的常数。将参数展开写为 $V=\sum_{i,j,k\ge0}V_{i,j,k}(x)a^ib^jc^k$，源系数递推为

$$
\begin{aligned}
(1-x)\theta^2 V_{i,j,k}
={}&-\theta V_{i-1,j,k}-\theta V_{i,j-1,k}
-x\theta V_{i,j,k-1}\\
&-V_{i-1,j-1,k}
-\delta_{i,0}\delta_{j,0}\delta_{k,1}x,
\end{aligned}
\tag{B4}
$$

负指标为零。右边只含更低源次数；再由零边界积分两次得到B的源系数。这才是可执行的生成规则，不是只给一个参数积分换名字。


在(B2)中取 $a=b=c=\epsilon,\ x=-z$，即得到(I6)标出的正号 $\epsilon z/(1+\epsilon)^2$ 乘5F4。这里求完的是 $\mathfrak B$，回代后的整图状态以(I6)为准。
