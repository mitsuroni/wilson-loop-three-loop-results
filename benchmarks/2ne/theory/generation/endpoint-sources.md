# 从4e-04推广到端点插入：6e-10/11/12与τ、E的共同来源

日期：2026-09-27。本稿执行4e-04之后的第一项推广：从原有序权重建立带源双核，以6e-10为控制，新生成6e-11/12；同时保留6e-36及其单核源作为τ的控制。不是把旧有限项重新命名为生成函数。

本轮完成：带源积分的四维到二维连接、带边界的端点原函数递推、运动学微分与源移位关系、任意阶系数的积分生成规则，以及四图28个目标位置的实际重建。新增覆盖14位置，累计35/329；另42图、294位置仍待新体系验收。τ与E从各自源矩重算后匹配旧定义，所有实际输出的字母仍为0、−1。

本轮没有得到带源双核的有限维微分闭包，也没有证明其任意阶都属于现有M/N/τ/E小表。检查、完整系数、CAS输出及数值控制见[复验](../../../../workflow/checks/2ne-endpoint-sources/README.md)。formal=not_run；当前会话无可调用的Wilson MCP，未修改Kernel/Lean。

## 1. 我们现在究竟要求哪一个积分

沿用统一母族的核、域和唯一物理归一化：

$$
D_z(t,s)=(1-t)(1-s)+zts,\qquad
k_\epsilon(t,s)=-(1+z)z^{-\epsilon/2}D_z(t,s)^{\epsilon-1},
$$
$$
\Omega_2=\{0<t_1<t_2<1,\quad0<s_1<s_2<1\}.
$$

把两个端点幂暂时取为独立参数：

$$
\boxed{\displaystyle
\mathcal Q_\pi(\epsilon;\lambda_1,\lambda_2;z)
=\int_{\Omega_2}k_\epsilon(t_1,s_{\pi(1)})
k_\epsilon(t_2,s_{\pi(2)})s_1^{\lambda_1}s_2^{\lambda_2}
\,d^2t\,d^2s.}
\tag{S1}
$$

这里下标1、2分别是较小、较大的右端点；配对只能是id或swap。零源时分别退回4e-04已经求值的裸核Q、S。参数先独立，只有还原具体图时才取零或ε。

从[原输入与字槽连接](remaining-bridges.md)得到：

$$
\boxed{\begin{aligned}
B_{10}&=\mathcal Q_{\rm id}(\epsilon;0,0;z)
          -\mathcal Q_{\rm id}(\epsilon;0,\epsilon;z),\\
B_{11}&=\mathcal Q_{\rm id}(\epsilon;0,\epsilon;z)
          -\mathcal Q_{\rm id}(\epsilon;\epsilon,0;z),\\
B_{12}&=\mathcal Q_{\rm swap}(\epsilon;0,0;z)
          -\mathcal Q_{\rm swap}(\epsilon;0,\epsilon;z).
\end{aligned}}
\tag{S2}
$$

因此新增内容很具体：在原双核上分别乘以较大端点的差分、两端点的差分、以及另一配对的较大端点差分。这里尚未乘物理前因子。

所有公式取正实z、正底数实对数。先在正小参数域建立，再在零点的小复邻域展开。核在相对角点的奇性类似两个变量之和的负一次幂；足够小的源幂和有限次对数插入仍有共同可积上界。比如将各源参数与ε的实部绝对值限制在1/8以内，可由完整立方体上的独立核乘积控制原排序域。换序和逐项微分使用这个分析依据，不将prepared的原迭代次序默认当作积测度。

## 2. 先完成两个端点积分，左边始终保留原源积分

令常数因子为

$$
c_\epsilon=-(1+z)z^{-\epsilon/2},\qquad
P_{\epsilon,\lambda}(x;y)
=\int_0^x r^\lambda D_z(r,y)^{\epsilon-1}\,dr.
$$

它已经可以显式求值。取A、B如下，则

$$
\begin{gathered}
A=1-y,\qquad B=(1+z)y-1,\\
P_{\epsilon,\lambda}(x;y)
=\frac{A^{\epsilon-1}x^{1+\lambda}}{1+\lambda}
{}_2F_1\!\left(1-\epsilon,1+\lambda;2+\lambda;-\frac{Bx}{A}\right).
\end{gathered}
\tag{S3}
$$

对id，保留两端点为u、v，即令u为较大左端点、v为较小右端点；分别积掉较小左端点和较大右端点：

$$
\boxed{\displaystyle
\mathcal Q_{\rm id}(\epsilon;\lambda_1,\lambda_2;z)
=c_\epsilon^2\int_0^1du\int_0^1dv\,
v^{\lambda_1}P_{\epsilon,0}(u;v)
\big[P_{\epsilon,\lambda_2}(1;u)-P_{\epsilon,\lambda_2}(v;u)\big].}
\tag{S4}
$$

对swap，保留u为较大左端点、v为较大右端点；分别积掉较小左端点和较小右端点：

$$
\boxed{\displaystyle
\mathcal Q_{\rm swap}(\epsilon;\lambda_1,\lambda_2;z)
=c_\epsilon^2\int_0^1du\int_0^1dv\,
v^{\lambda_2}P_{\epsilon,0}(u;v)P_{\epsilon,\lambda_1}(v;u).}
\tag{S5}
$$

**当前进度：** 同一个源积分已由四维降为二维；两个被消去的积分由明确的超几何原函数取代。两式的v代表不同端点，因此外部源幂不同，不能直接交换。零源时，S4就是4e-04原二维式，S5就是交换配对的对应式。

## 3. 从源原函数得到递推；端点项不能丢

对原函数中的幂乘积求导：

$$
\frac{d}{dx}\big[x^{\lambda+1}(A+Bx)^\epsilon\big]
=x^\lambda(A+Bx)^{\epsilon-1}
\big[(\lambda+1)A+(\lambda+1+\epsilon)Bx\big].
$$

从零积到x，就得到带边界的移位关系：

$$
\boxed{\displaystyle
(\lambda+1)A P_{\epsilon,\lambda}(x;y)
(\lambda+1+\epsilon)B P_{\epsilon,\lambda+1}(x;y)
=x^{\lambda+1}(A+Bx)^\epsilon.}
\tag{S6}
$$

右边是上端点项。下端点项在当前参数域为零。B等于零时使用连续极限，不除以B；直接得到原函数的幂积分。

源展开本身也有明确的内容：

$$
\boxed{\displaystyle
[\epsilon^m\lambda^n]P_{\epsilon,\lambda}(x;y)
=\frac1{m!n!}\int_0^x
\frac{\log^m(A+Br)\log^n r}{A+Br}\,dr.}
\tag{S7}
$$

例如第一个端点源导数，在ε为零时为

$$
\left.\partial_\lambda P_{0,\lambda}(x;y)\right|_{\lambda=0}
=\frac{\log x\log(1+Bx/A)+\operatorname{Li}_2(-Bx/A)}B.
$$

B为零时，它趋于x乘以log x减1再除以A。这解释了端点插入为何会产生二重对数权重；不是先根据旧Hlog答案猜出一个Li₂。

**当前进度：** 插入所需的新原函数已有自然来源、明确边界及源移位关系。将它代入S4/S5之后，剩余二维积分仍须求值；不能把S3的完成误当作整个双核已完成紧凑全ε求值。

## 4. 对运动学变量求导，究竟留下什么

为使crossing后的源仍能表达，在S1的积分权重中允许同时出现

$$
\prod_{i=1}^2s_i^{\lambda_i}(1-s_i)^{\mu_i}.
$$

记对应源为带两组源参数的同一个Q，原S1是所有μ为零的子族。令

$$
\theta=z\partial_z,\qquad
v(s)=s(1-s),\qquad
\eta=-\frac1{1+z}-\frac\epsilon2.
$$

单个核直接满足

$$
\theta k_\epsilon(t,s)
=\partial_s\big[v(s)k_\epsilon(t,s)\big]
+\big[\eta+(1+\epsilon)s\big]k_\epsilon(t,s).
\tag{S8}
$$

把两个核的导数相加，在有序三角形上分部积分。外侧端点的v为零；对角面上两个通量分别带相反法向，且两侧的v相同，因此**相加后**抵消。它们不能在分别处理单个核时各自删掉。导数作用到源权重，给出

$$
v(s)\partial_s\{s^\lambda(1-s)^\mu\}
=\big[\lambda-(\lambda+\mu)s\big]s^\lambda(1-s)^\mu.
$$

于是原带源双核满足

$$
\boxed{\begin{aligned}
\theta\mathcal Q_\pi(\epsilon;\boldsymbol\lambda,\boldsymbol\mu;z)
={}&\big[2\eta-\lambda_1-\lambda_2\big]
\mathcal Q_\pi(\epsilon;\boldsymbol\lambda,\boldsymbol\mu;z)\\
&+\sum_{i=1}^2(1+\epsilon+\lambda_i+\mu_i)
\mathcal Q_\pi(\epsilon;\boldsymbol\lambda+\mathbf e_i,\boldsymbol\mu;z).
\end{aligned}}
\tag{S9}
$$

这是一条真实的源移位微分关系，适用于两种配对。它还说明当前障碍：导数引入了指数加1的源积分。重复操作通常产生更多移位，现有关系没有把它们全部约回有限个已解对象。故不能照搬4e-04的特殊消项，宣称所有嵌套原函数已经消失。

边界条件由S4/S5的原积分定义及S6给出；这里没有将S9作为缺少边界的独立求解器。可复验的整数ε多项式控制、非整数源的数值导数都支持S9。

## 5. 插入情况下，乘积与crossing关系如何改变

同时反射两边端点，源变换为

$$
\boxed{\displaystyle
\mathcal Q_\pi(\epsilon;\lambda_1,\lambda_2;\mu_1,\mu_2;1/z)
=\mathcal Q_\pi(\epsilon;\mu_2,\mu_1;\lambda_2,\lambda_1;z).}
\tag{S10}
$$

因此只保留端点的正幂源，不会自动对crossing封闭。零源的Q、S各自偶对称，带插入者一般没有这个性质。

令单核右端点源为

$$
H_\lambda(\epsilon,z)=\int_0^1\!dt\int_0^1\!ds\,
k_\epsilon(t,s)s^\lambda.
$$

对两个独立单核的两条边分别分割排序域，四个区域给出

$$
\boxed{\displaystyle
H_{\lambda_1}H_{\lambda_2}
=\sum_{\pi\in\{\mathrm{id},\mathrm{swap}\}}
\big[\mathcal Q_\pi(\epsilon;\lambda_1,\lambda_2;z)
+\mathcal Q_\pi(\epsilon;\lambda_2,\lambda_1;z)\big].}
\tag{S11}
$$

零源时它正好退回Q加S等于一半单核平方；源不同时须同时保留两种插入位置。对称关系减少重复工作，不直接证明这些带源对象解析独立。

## 6. τ与E是从哪里读出来的

### 6.1 τ：单核两端同时插入

同一端点原函数给出

$$
\boxed{\displaystyle
H_D(\epsilon,z)
=c_\epsilon\int_0^1s^\epsilon
P_{\epsilon,\epsilon}(1;s)\,ds.}
\tag{S12}
$$

它等于原双端幂权重的二重积分。令

$$
\Phi=\log D_z(t,s)-\frac{\log z}{2},\qquad
U=\log t,\quad V=\log s,
$$

则直接从指数展开得到

$$
\tau^+_{D,q}
=\left[-\frac{1+z}{q!}\int_0^1\!\int_0^1
\frac{(\Phi+U+V)^q}{D_z(t,s)}\,dt\,ds\right]^+.
\tag{S13}
$$

上标加号表示对z与其倒数取平均。q等于4给出统一权重六表保留的τ。本轮实际积分S13的四阶矩，再与旧τ展开核对为零残差。旧Meijer G全ε求值仍有效；本轮增加的是它与双核所用同一端点原函数的明确连接。

### 6.2 E：双核插入的奇部分乘奇对数

S2直接给出

$$
B_{10}=\epsilon b_{10,1}+O(\epsilon^2),\qquad
b_{10,1}=-\int_{\Omega_2}k_0(t_1,s_1)k_0(t_2,s_2)\log s_2.
$$

由S10，奇部分为较大端点与反射后较小端点的对数差，因此

$$
\boxed{\displaystyle
E=-\ell b_{10,1}^-
=\frac{\ell(1+z)^2}{2}\int_{\Omega_2}
\frac{\log(s_2/(1-s_1))}{D_z(t_1,s_1)D_z(t_2,s_2)}\,d^2t\,d^2s,
\qquad\ell=\log z.}
\tag{S14}
$$

这同时给出E的源位置、阶数和crossing性质。它来自未展开插入的首阶，不是有限项减去某个拟合结果后才定义。本轮由新算的首阶奇系数生成E，匹配旧定义；但E的高阶对应物尚未全部约化成一个有限的矩族。

## 7. 任意阶怎样生成，本轮实际算到了哪一阶

设两个核的总对数为

$$
\Phi_\pi=\log D_z(t_1,s_{\pi(1)})+\log D_z(t_2,s_{\pi(2)})-\ell,
\qquad V_i=\log s_i.
$$

将S1沿物理源方向展开，直接有

$$
\begin{aligned}
b_{10,q}&=\frac1{q!}\int_{\Omega_2}k_0k_0
\big[\Phi_{\rm id}^q-(\Phi_{\rm id}+V_2)^q\big],\\
b_{11,q}&=\frac1{q!}\int_{\Omega_2}k_0k_0
\big[(\Phi_{\rm id}+V_2)^q-(\Phi_{\rm id}+V_1)^q\big],\\
b_{12,q}&=\frac1{q!}\int_{\Omega_2}k_0k_0
\big[\Phi_{\rm swap}^q-(\Phi_{\rm swap}+V_2)^q\big].
\end{aligned}
\tag{S15}
$$

各行的两个零阶核使用相应配对。q等于零时全部消失；这证明三图裸组合至少一阶起始。这里是任意阶有效的积分生成规则，不是将每个积分都求值完毕的闭式。

本轮实际计算先解析积掉较小左端点，再展开源差分。以a表示运动学参数，并记

$$
a=\frac1{1+z},\qquad d(u,v)=a(1-u)(1-v)+(1-a)uv,
$$
$$
K_{\rm id}=\frac{\log d(u,v)-\log a-\log(1-v)}{(v-a)d(u,w)},
\qquad
K_{\rm swap}=\frac{\log d(u,w)-\log a-\log(1-w)}{(w-a)d(u,v)}.
$$

积分域为u从零到1、w从零到1、v从零到w。相应首阶核的对数修正为

$$
\begin{aligned}
T_{\rm id}&=\log d(u,w)+\tfrac12\log d(u,v)
+\tfrac12\log(1-v)-\tfrac12\log a-\log(1-a),\\
T_{\rm swap}&=\log d(u,v)+\tfrac12\log d(u,w)
+\tfrac12\log(1-w)-\tfrac12\log a-\log(1-a).
\end{aligned}
$$

两个需要实际积分的系数分别是对应K乘以下多项式：

| 图 | 一阶 | 二阶 |
|---|---|---|
| 6e-10 | $-\log w$ | $-T_{\rm id}\log w-\log^2w/2$ |
| 6e-11 | $\log w-\log v$ | $T_{\rm id}(\log w-\log v)+(\log^2w-\log^2v)/2$ |
| 6e-12 | $-\log w$ | $-T_{\rm swap}\log w-\log^2w/2$ |

这六个积分均在本轮新会话中完成，不读取旧结果作为生成输入。6e-36的二至四阶组合矩及τ的四阶矩也重新积分。消元链的末端候选因子只有z与1加z，实际Hlog结果也只有0、−1字母；相同字母表不意味着相同紧凑M/N函数空间。

## 8. 最后恢复完整物理I，并与旧解比较

三张双核图有相同的显式二阶极点因子与加权crossing参数，但常数不同：

$$
\boxed{\displaystyle
I_g=\frac{\kappa_g Z^3}{\epsilon^2\cosh(3\epsilon\ell/2)}
\frac{e^{\epsilon\ell/2}B_g(\epsilon,z)
+e^{-\epsilon\ell/2}B_g(\epsilon,1/z)}2,
\quad Z=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),}
\tag{S16}
$$

其中10图的κ为负八分之一，11/12图为负十六分之一。因为裸组合从一阶起始，完整I只有简单极点，且

$$
\boxed{\displaystyle
I_{g,-1}=\kappa_g b_{g,1}^+,\qquad
I_{g,0}=\kappa_g\left(b_{g,2}^++\frac\ell2b_{g,1}^-\right).}
\tag{S17}
$$

不能因为最终I是偶函数就删掉第二项；三个漏掉奇部的负控制均不通过。

对6e-36，裸组合为单核乘两个端点差分，从二阶起始。完整卷积为

$$
\begin{aligned}
I_{36,-2}&=-b_{36,2}^+/32,\\
I_{36,-1}&=-(b_{36,3}^++\ell b_{36,2}^-)/32,\\
I_{36,0}&=-\big[b_{36,4}^++\ell b_{36,3}^-
+(3\zeta_2/2-5\ell^2/8)b_{36,2}^+\big]/32.
\end{aligned}
\tag{S18}
$$

四图28个既定位置、其中9个非零，全部与现行结果精确一致。完整表达见[系数表](endpoint-source-coefficients.md)。在既有权重六M/N空间取商，得到

$$
\begin{aligned}
[I_{10,0}]&=\frac1{24}[\tau^+_{D,4}]-\frac1{12}[E],\\
[I_{11,0}]&=\frac1{96}[\tau^+_{D,4}]+\frac7{96}[E],\\
[I_{12,0}]&=-\frac5{96}[\tau^+_{D,4}]+\frac1{96}[E],\\
[I_{36,0}]&=-\frac1{32}[\tau^+_{D,4}].
\end{aligned}
\tag{S19}
$$

方括号只表示模所声明的M/N与周期多项式空间，不能将S19当作完整答案。τ和E的定义来自S13/S14；此处Hlog精确消元只用于验证有限阶表达、确定旧坐标中的系数，不冒充全阶基底来源证明。

## 9. 这一步对主线的推进及剩余工作

现在，裸双核、带插入双核和τ的单核控制已通过同一个端点原函数连接；源的位置、微分、边界、系数预算和物理归一化都显式保留。新鲜计算确认现有τ、E足以补齐本批权重六的M/N缺口，没有出现新函数需求或新字母。

仍未完成的是：把S9的移位层级化为有限闭包，或将S4/S5全部求成类似4e-04的紧凑全ε标准函数组合。普通幂源也不能自动代替一般BA/AB有序字。

下一步按原计划推广已掌握的单核字槽、低核乘积与镜像关系到剩余图，逐项增加既定有限阶验收；源移位的进一步压缩作为明确的研究问题保留。当前不自动提高4e-04阶数、不开展3g或三核全阶求值，也不把本批35/329说成全表完成。
