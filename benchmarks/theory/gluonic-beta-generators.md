# 纯胶子图的额外函数：共同Beta端点源与扩展生成函数

本稿比较纯2ne、5e-3g、4e-4g和已完成的双TGV图。目标是从原积分解释旧 $\tau^+_{D,4}$、$E$ 及其导数，并寻找能统一生成它们的固定积分族。**先建立共同的未展开源，再计算源矩与旧解比较；不把逐图余项反过来定义为新生成元。**

本轮的一个确定进展是：将旧VB中绑定的Beta指数分离后，已经能分别表示 $\tau^+_{D,4}$ 和 $E$，不再仅得到 $\tau^+_{D,4}+10E$。本稿同时区分有限阶转换、源的全参数恒等式与完整图的全epsilon解；它们不是同一个完成状态。

[逐图出现清单](gluonic-extra-basis-inventory.md)、[计算与核验](../../workflow/checks/gluonic-beta-sources/README.md)。沿用各图现行C/I，正实 $z>0$，formal=not_run；不增加默认47图的Kernel根。

## 1. 先辨认真正相同的对象

统一记

$$
\tau(z)=\tau^+_{D,4}(z),\qquad \Theta=z\frac{d}{dz}.
$$

它不是所有旧 $\tau_{a,q}$ 的统称：许多低阶tau记号已经能够回到M/N。按现行统一系数，真正涉及本轮额外方向的图为：纯6e的11图，5e-3g的15图，4e-4g的05/08，以及已完成双TGV的10、19。具体有理系数见清单。

这里有三个有用的区别。5e-3g-11/12/36/37只需要 $\tau+10E$；纯6e及其他TGV图使用不同的线性组合；4e-4g-05/08需要 $(z-1/z)\Theta\tau/128$，这要求允许运动学导数及有理 $z$ 系数。5e-3g-27/28/29另外含 $\Omega_6$，它不能在只转换tau/E时被省略。

因此，共同理论至少应保留**独立端点权重、有序配对和运动学微分**。仅按对边条数分类，或把同一有限权重符号出现视为全阶恒等式，都不足以确定生成结构。

## 2. 一个从原积分出现的共同Euler模块

对 $w>0$，定义

$$
\boxed{\mathcal F_{\alpha,\beta,\gamma}(w)
=\int_0^1u^\alpha(1-u)^\beta(1-u+uw)^{\gamma-1}du.}
\tag{G1}
$$

初始条件可取 $\Re\alpha,\Re\beta>-1$，之后按需要解析延拓。三个指数分别记录两端Beta密度及线性分母的幂，均由原积分产生。Euler积分实际给出

$$
\mathcal F_{\alpha,\beta,\gamma}(w)
=B(1+\alpha,1+\beta)
{}_2F_1(1-\gamma,1+\alpha;2+\alpha+\beta;1-w).
\tag{G2}
$$

这是标准[Gauss Euler积分](https://dlmf.nist.gov/15.6.E1)，这里使用通常的 ${}_2F_1$，不是已除以Gamma下参数的Olver规范。

这个模块在 $w=1$ 取解析支，并具有以下Gauss方程和边界：

$$
\begin{aligned}
&w(1-w)\mathcal F''
+[1-\beta-\gamma-(3+\alpha-\gamma)w]\mathcal F'
-(1-\gamma)(1+\alpha)\mathcal F=0,\\
&\mathcal F(1)=B(1+\alpha,1+\beta),\qquad
\mathcal F'(1)=-\frac{(1-\gamma)(1+\alpha)}{2+\alpha+\beta}\mathcal F(1).
\end{aligned}
\tag{G3}
$$

参数导数直接产生 $\log u$、$\log(1-u)$、$\log(1-u+uw)$ 的固定幂。因此这不是只把原图重命名成一个未知函数。

### 2.1 TGV：旧VB只取了一条参数方向

[5e-3g端点推导](../tgv/families/t1-opposite/beta-generators.md)中的权重恰好是

$$
\boxed{\Phi_\epsilon(w)=(1-w)\mathcal F_{-\epsilon,-\epsilon,2\epsilon}(w).}
\tag{G4}
$$

原来的两个端点指数都等于 $-\epsilon$，分母指数增量固定为 $2\epsilon$，镜像权重又固定为 $z^{\pm\epsilon}$。旧VB沿着这些绑定的参数展开，因而没有探测全部独立方向。

T3/T4的边界也出现同类模块。例如T4保留的端点函数是

$$
\begin{aligned}
h_j(x,w;\epsilon)
&=\int_0^1[b(1-b)]^{-\epsilon}
[A_j(1-b)+B_jb]^{-1+2\epsilon}db\\
&=A_j^{-1+2\epsilon}
\mathcal F_{-\epsilon,-\epsilon,2\epsilon}(B_j/A_j),
\end{aligned}
$$

其中 $A_j=Q_j(0)$、$B_j=Q_j(1)$，$Q_j$ 见[T4完整边界推导](../tgv/families/t4-opposite/method.md)。分子中的Beta参数多项式对应指数的整数移位，核的导数则满足

$$
\partial_w\mathcal F_{\alpha,\beta,\gamma}
=(\gamma-1)\mathcal F_{\alpha+1,\beta,\gamma-1}.
$$

这些连接涵盖边界模块，并不消去原有的recoil体项或排序权重。尤其 $B_j/A_j$ 还依赖其他端点，后续积分可以产生更多外部奇点；不能由共同Euler模块推断所有TGV图都只有0、−1字母，原Omega方向仍须保留。

### 2.2 纯2ne：双端核可先化为相同模块的谱积分

考虑已存在的带权单核

$$
\mathscr D_{\eta;A,B;C,D}(z)
=\int_0^1du\,dv\,
 u^A(1-u)^Bv^C(1-v)^D D_z(u,v)^{\eta-1},
\quad D_z=(1-u)(1-v)+zuv.
\tag{G5}
$$

先在共同绝对收敛域操作；例如各指数实部足够接近0时，两个角点仍有可积径向界。令 $u=x/(1+x)$、$v=y/(1+y)$，再用乘积变量 $p=xy$。核变成 $(1+zp)/[(1+x)(1+y)]$，内层再取 $x=pt/(1-t)$，得到

$$
\boxed{\mathscr D_{\eta;A,B;C,D}(z)
=\int_0^\infty p^A(1+zp)^{\eta-1}
\mathcal F_{\eta+A+D,\,\eta+B+C,\,-\eta-A-B}(p)\,dp.}
\tag{G6}
$$

这条等式同时适用于纯2ne双端幂和双TGV的两端Beta权重。对于6e-36/39所需的

$$
F_D(\epsilon,z)=\int_0^1du\,dv\,(uv)^\epsilon D_z^{\epsilon-1},
$$

密度是 $\rho_\epsilon(p)=\mathcal F_{2\epsilon,2\epsilon,-2\epsilon}(p)$。利用
$\rho_\epsilon(1/p)=p^{1+2\epsilon}\rho_\epsilon(p)$ 折回单位区间：

$$
\boxed{F_D(\epsilon,z)=\int_0^1\rho_\epsilon(t)
\left[t^\epsilon(1+zt)^{\epsilon-1}+(z+t)^{\epsilon-1}\right]dt.}
\tag{G7}
$$

两个原积分变量中有一个已实际积成标准Gauss函数。端点 $t=1$ 的密度有限，$t\to0$ 的对数行为仍可积。旧tau由

$$
\tau=[\epsilon^4]\frac{H_D(\epsilon,z)+H_D(\epsilon,1/z)}2,
\qquad H_D=-(1+z)z^{-\epsilon/2}F_D
\tag{G8}
$$

生成。接到实际完整图时，仍使用[6e-36/39的原归一化](../2ne/families/single-opposite/three-loop.md)：

$$
\begin{aligned}
I_{6e-36}=I_{6e-39}
={}&\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3(1+z)}
{64\epsilon^4\cosh(3\epsilon\ell/2)}\\
&\times\Big[z^{\epsilon/2}(F_0-2F_1+F_D)(\epsilon,z)\\
&\qquad+z^{-1-\epsilon/2}(F_0-2F_1+F_D)(\epsilon,1/z)\Big].
\end{aligned}
\tag{G9}
$$

$F_0,F_1$已经由M源给出；G7给出这里 $F_D$ 的新共同模块表示。G9保持原物理I，不能把其中的谱积分直接当作已经等于下文某个自由VB取值。

### 2.3 E：相同模块在有序双核中相乘

E的来源是有序双核 $\mathcal Q_{\rm id}$ 的镜像奇端点矩，再乘 $\ell=\log z$，详见[原端点源](../2ne/theory/generation/endpoint-sources.md)。关键的端点原函数为

$$
\begin{aligned}
P_{\eta,\lambda}(x;y)
&=\int_0^x t^\lambda D_z(t,y)^{\eta-1}dt\\
&=(1-y)^{\eta-1}x^{1+\lambda}
\mathcal F_{\lambda,0,\eta}
\left(1+\frac{[(1+z)y-1]x}{1-y}\right).
\end{aligned}
\tag{G10}
$$

于是原四维有序核先积掉两个端点后，成为两个已知 $\mathcal F$ 原函数的乘积或端点差。它与单核tau共享同一种Euler模块，但还保留了排序和乘积积分。这解释了共同函数能够出现，也说明为什么不能据此抹掉单核与有序双核的差别。

### 2.4 四胶子与双三胶子：密度及其导数

4e-4g同边合并留下的密度，直接写成

$$
\begin{aligned}
\kappa_\epsilon(w)={}&w(1-w)^{1-\epsilon}
\mathcal F_{-\epsilon,0,-\epsilon}(w)\\
&+(1-w)w^{1-\epsilon}
\mathcal F_{-\epsilon,0,-\epsilon}(1-w).
\end{aligned}
\tag{G11}
$$

它等于原来的 $wB_{1-w}(1-\epsilon,-\epsilon)+(1-w)B_w(1-\epsilon,-\epsilon)$；与TGV不同，此处两端Beta指数并不相等。标准不完全Beta关系见[DLMF 8.17.7](https://dlmf.nist.gov/8.17.E7)。F3有限项中的 $\Theta\tau$ 是外层核微分后的同一方向，因此共同函数空间还须考虑 $\Theta$ 作用。

双TGV-19的原IBP则留下

$$
\int_0^1du\,dv\, f_\epsilon(u)f_\epsilon(v)D_z^{-1+3\epsilon},
\qquad
f_\epsilon(u)=[u(1-u)]^{-\epsilon}-B(1-\epsilon,1-\epsilon).
\tag{G12}
$$

展开这个端点差乘积后，每项都是G5的参数取值。双TGV-10保留有序三端点的有符号重心密度，属于更深的有序组合；本轮只接通它已知的tau/E部分，没有宣称其整个密度已经约成一个G1。

## 3. 从共同模块导出自由Beta生成函数

继续采用 $r_z=z/(1+z)$、$W=\log[\beta(1-\beta)]$、$L=\ell+\log((1-\beta)/\beta)$。将两端对称和不对称的指数分别记为 $p,q$，定义

$$
\boxed{\begin{aligned}
\mathcal V_{\rm B}^{\rm free}(a,b;p,q,c,d;z)
={}&\int_0^1\frac{e^{aW}}{\beta-r_z}\frac{\sinh(bL)}b\\
&\times\Big[z^d\beta\mathcal F_{p+q,p-q,c}(1-\beta)\\
&\qquad+z^{-d}(1-\beta)\mathcal F_{p+q,p-q,c}(\beta)\Big]d\beta.
\end{aligned}}
\tag{G13}
$$

所有源参数先独立、靠近0；$b=0$及 $\beta=r_z$ 都是可去点。它不是六个新函数族，而是同一个固定Euler核上六种明确权重的生成函数。原VB准确嵌入为

$$
\boxed{\mathcal V_{\rm B}(a,b,t;z)
=\mathcal V_{\rm B}^{\rm free}(a,b;-t,0,2t,t;z).}
\tag{G14}
$$

全零Beta参数时，$\mathcal V_{\rm B}^{\rm free}(a,b;0,0,0,0;z)=-\partial_a\mathcal V_M(a,b;z)$。但G13只保留了sinh部分；完整镜像生成还应保留另一个自然部分。直接从减除指数核作镜像投影，定义

$$
\begin{aligned}
\mathscr U_{\rm B}^{\rm free}
={}&\int_0^1\frac{e^{aW}}{\beta-r_z}\Big[
 z^d\beta\mathcal F_{p+q,p-q,c}(1-\beta)(e^{bL}-1)\\
&\qquad-z^{-d}(1-\beta)\mathcal F_{p+q,p-q,c}(\beta)(e^{-bL}-1)\Big]d\beta.
\end{aligned}
\tag{G17}
$$

它分为 $b\mathcal V_{\rm B}^{\rm free}+b^2\mathcal V_{{\rm B},N}^{\rm free}$。后者把G13的 $\sinh(bL)/b$ 改为 $(\cosh(bL)-1)/b^2$，并将两端权重的和改成差。因为零参数的权重差等于 $-\log((1-\beta)/\beta)$，所以

$$
\boxed{\mathscr U_{\rm B}^{\rm free}(a,b;0,0,0,0;z)
=-b\partial_a\mathcal V_M(a,b;z)-b^2\mathcal V_N(a,b;z).}
\tag{G18}
$$

这里使用原有 $\mathcal V_N=\int e^{aW}R(\cosh(bL)-1)/[b^2(\beta-r_z)]\,d\beta$，$R=\log((1-\beta)/\beta)$。因此M型和N型都由同一个完整镜像源自然出现，不需要分别为两个投影猜测函数族。


为了说明新族到底是什么，令

$$
S_u=\log[u(1-u)],\quad A_u=\log\frac{u}{1-u},\quad
D_+=1-u\beta,\quad D_-=1-u+u\beta.
$$

定义原始对数幂矩，不在矩内放阶乘：

$$
\begin{aligned}
\mathfrak B^{\rm free}_{k,j;m,n,s,h}
={}&\ell^h\int_0^1\frac{W^kL^j}{\beta-r_z}d\beta
\int_0^1du\,S_u^m A_u^n\\
&\times\left[\frac{\beta\log^sD_+}{D_+}
+(-1)^{h+j+1}\frac{(1-\beta)\log^sD_-}{D_-}\right].
\end{aligned}
\tag{G15}
$$

$k,m,n,s,h\ge0$、$j\ge1$；G13用奇数 $j$，完整G17还包括偶数 $j$。权重为 $k+j+m+n+s+h+2$。逐个参数展开就有

$$
\mathscr U_{\rm B}^{\rm free}
=\sum_{k,m,n,s,h\ge0}\sum_{j\ge1}
\frac{a^k b^j p^m q^n c^s d^h}
{k!j!m!n!s!h!}
\mathfrak B^{\rm free}_{k,j;m,n,s,h}.
\tag{G16}
$$

这里没有逐图定义或逐阶猜测权重。对应的Euler/Gauss表示、一阶参数导数、端点解析性和固定积分次序给出任意指定指标的计算办法。**但一个图能否写成G13/G17的有限组合，仍需从其完整源推导；有限阶属于这些矩的空间还不能替代这一步。**

## 4. tau与E现在可以分别转换

下面给出实际求出的完整权重六字典。它们来自先独立积分G15，再做精确向量转换；该转换没有参与G13或G15的定义。所有M/N的自变量为z。

$$
\begin{aligned}
\tau^+_{D,4}={}&- 3 \mathfrak B^{\rm free}_{0,1;1,0,2,0} - 2 \mathfrak B^{\rm free}_{0,1;1,0,1,1}\\
&- \frac{11 N_{2,2}}{8} + \frac{M_{1,1}^{2}}{4}\\
&+ \frac{5 M_{2,3}}{64} + \frac{29 M_{0,1}^{3}}{288}\\
&+ \frac{31 M_{4,1}}{128} - 12 M_{1,1} \zeta_{3}\\
&- \frac{459 M_{0,1} \zeta_{2}^{2}}{40} - \frac{23 M_{0,1}^{2} \zeta_{2}}{48}\\
&+ \frac{3 M_{0,1} M_{2,1}}{8} + \frac{39 M_{2,1} \zeta_{2}}{8}
\end{aligned}
$$

$$
\begin{aligned}
E={}&- \mathfrak B^{\rm free}_{0,1;1,0,1,1} - \frac{M_{2,3}}{4}\\
&+ \frac{N_{2,2}}{4} + \frac{3 M_{0,1}^{3}}{20}\\
&- \frac{79 M_{0,1}^{2} \zeta_{2}}{60} - \frac{53 M_{0,1} \zeta_{2}^{2}}{10}\\
&- \frac{3 M_{2,1} \zeta_{2}}{4} - \frac{M_{0,1} M_{2,1}}{4}
\end{aligned}
$$

两个需要的新矩都位于对称子族 $q=0$；旧VB仍然不能把它们作为独立方向自由取用，因为它绑定了 $p,c,d$。这正是“旧VB只给tau+10E”与“自由Beta源分别给tau、E”的区别。

因此清单中的tau/E贡献已经有统一替换规则。4e-4g的两图可对上面的tau恒等式取 $\Theta$，同时保留其有理运动学前因子；不能把微分后的函数无条件当作原来未微分的M/N空间成员。另含 $\Omega_6$ 的三张TGV图仍须保留该项。

## 5. 正epsilon阶检验与真实边界

| 权重／所测阶 | M/N原秩 | 对称Beta扩展秩 | 再加入目标 | 完整指数与镜像扩展秩 | 再加入目标 |
|---|---:|---:|---:|---:|---:|
| 6／ε^0 | 12 | 14 | 14 | 14 | 14 |
| 7／ε^1 | 22 | 28 | 28 | 28 | 28 |
| 8／ε^2 | 36 | 52 | 53 | 52 | 53 |

下表列出完整指数与镜像扩展对六个代表的检验；等价图由已有全ε关系继承。

| 代表 | 有限项 | ε¹ | ε² |
|---|---|---|---|
| 6e-36 | 可表示 | 可表示 | 可表示 |
| 6e-10 | 可表示 | 可表示 | 可表示 |
| 6e-11 | 可表示 | 可表示 | 在本空间中仍有余项 |
| 6e-12 | 可表示 | 可表示 | 在本空间中仍有余项 |
| 6e-20 | 可表示 | 可表示 | 可表示 |
| 6e-21 | 可表示 | 可表示 | 可表示 |

这些新矩都保持字母0、−1。对称子族已覆盖全部六个ε¹及四个ε²；其权重八缺口在11/25与12/24。不对称指数及偶数L幂的cosh配对都已补测，仍未扩大这三个权重的空间。完整镜像源的结果如表，仍未覆盖的代表为 pure-11, pure-12。

这些检验采用全部所列指标、同权重镜像偶多项式、普通zeta及已记录的权重八周期，系数为有理常数；不允许任意函数分母、任意有理z变换或未列出的特殊函数。成员关系通过完整Hlog向量回代；非成员关系用对全部候选为零、对目标非零的精确线性泛函核验。固定权重方向数不等于全阶函数族最小数。

补测前后的非成员证书也确认两个余项相反：若 $\mathcal A_8$ 表示本节的权重八候选空间，则

$$
[\epsilon^2]I_{6e-11}\equiv-[\epsilon^2]I_{6e-12}
\pmod{\mathcal A_8},
\qquad [\epsilon^2]I_{6e-11}\notin\mathcal A_8.
$$

这是已检验权重的一个共同缺失方向，不是两个已证明独立的全阶新函数族。完整回代记录见[共同余项](../../workflow/checks/gluonic-beta-sources/remaining-direction.json)。

## 6. 下一步应求什么

现在应把“共同Euler模块”“固定矩族”和“完整I的参数取值”三层分开推进。G4、G6、G10、G11、G12已经在未展开积分层建立连接；G13/G17提供自然的扩展候选；权重六字典和高阶检验评估了候选的覆盖范围。

下一项应优先追踪仍未闭合的有序双核端点原函数乘积，尝试在完整镜像G17之外加入由原排序强制产生的Euler模块卷积。先验证6e-11/25与12/24的共享源，再讨论三核，不因某个有限阶可转换就宣布全部全epsilon生成已经完成。4g的进一步高阶还需在其原IBP和密度上展开，不能仅由有限项的导数字典预测。

本轮保留纯2ne的44/47算法生成、36/47既有M/N紧凑和56/329新框架有限阶验收口径，不用跨族字典修改这些计数。所有旧r001、源数据、失败记录及前次VB检查保持。新增的是共同源连接、自由Beta矩及其实际空间检验；未改Kernel/Lean，未提交Git。
