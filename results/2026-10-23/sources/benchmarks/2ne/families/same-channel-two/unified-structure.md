# [2] 族：从通用 J、I 到双标记有序母积分

日期：2026-09-16。本文回应“Σ 缺少清晰数学结构”的问题，采用[空集与单对边统一笔记](<../../../../GF cal 推导草稿/empty_and_single_opposite_unified_notes.md>) §§2、6–12 的方法：先定义原积分结构及统一生成规则，再计算系数。原笔记只作背景阅读；新结果的必要定义全部写在本文。

适用范围：光锥四边形的纯交换图，恰有两条直接对边传播子，位于同一通道，没有内部顶点、自能或其他对边通道。一般公式允许增加邻边传播子；当前实际对照语料仍为 prepared 0.3 的七个 [2] 根。本文是解析推导与 CAS/数值核验，formal: not_run。

旧 Σ 与 r001/r002 不改写。新的[结构系数表](structural-solutions.md)给出与旧结果逐项相同、但以拓扑母族系数表达的答案。[检查记录](../../../../workflow/checks/same-channel-two-structure/README.md)独立保存。

## 1. 原始通用 J：配对与边序同时属于数据

固定两条对边支撑边为 2、4，选择与当前输入一致的参数方向：

$$D_z(t,s)=(1-t)(1-s)+zts,\qquad z>0.$$

将两条边上的对边端点按参数从小到大编号：

$$0<t_1<t_2<1,\qquad0<s_1<s_2<1.$$

编号记录端点位置，不是传播子标签。另用排列

$$\pi\in S_2=\{\mathrm{id},\mathrm{swap}\}$$

记录传播子配对：

$$
\mathrm{id}: (t_1,s_1),(t_2,s_2),\qquad
\mathrm{swap}: (t_1,s_2),(t_2,s_1).
$$

含 n 条传播子的原积分可写为

$$
J_{\mathscr G}=(1+z)^2
\int_{\Omega_{\mathscr G}}
\prod_{j=1}^{2}D_z(t_j,s_{\pi(j)})^{\varepsilon-1}
\prod_{e\in E_{\rm adj}}
q_e\,[q_e(1-u_e)v_e]^{\varepsilon-1}\,d^{2n}s,
\quad q_e\in\{1,z\}.
$$

Ω 保留全部原边端点排序。写成这一形式前，需逐项核对源多项式及参数映射；不能由图名猜配对，也不能混淆 id/swap 与物理 crossing。当前七图的核与原映射已核对，见每图旧 derivation.md。

令 r_G 是 q_e=z 的邻边数。每条邻边因子为

$$q_e^\varepsilon B_\varepsilon(u_e)A_\varepsilon(v_e),\quad
A_\varepsilon(x)=x^{\varepsilon-1},\ B_\varepsilon(x)=(1-x)^{\varepsilon-1}.$$

两条 D 核的整体号相乘为正；色因子、重数及傅里叶的 (-2)^{-n} 留在 κ_G，不混入这里。

## 2. 六个有序区间：双标记边字

沿用原笔记：

$$
F_w(a,b;\varepsilon)=\int_{a<x_1<\cdots<x_{|w|}<b}\prod_jf_{w_j}(x_j)\,d^{|w|}x,
\qquad\widehat F_w=\varepsilon^{|w|}F_w,\quad\widehat F_\varnothing=1.
$$

两条支撑边各有两个 O，故一般边字是

$$p_0Op_1Op_2,\qquad q_0Oq_1Oq_2.$$

给定对边端点后，它们产生

$$
P_{\mathbf p}(t_1,t_2)=
F_{p_0}(0,t_1)F_{p_1}(t_1,t_2)F_{p_2}(t_2,1),
$$
$$
P_{\mathbf q}(s_1,s_2)=
F_{q_0}(0,s_1)F_{q_1}(s_1,s_2)F_{q_2}(s_2,1).
$$

其余两条边没有 O，给出常数 B_G(ε)=F_{w_1}(0,1)F_{w_3}(0,1)。

这是与 [1] 的关键区别：[1] 每条支撑边分成两个区间；[2] 每条分成三个区间，并额外记录配对排列。尤其不能把中间区间 F_w(t1,t2) 当成两个独立端点权重。

所有邻边共有 2n−4 个端点，故正规化后

$$
J_G=\frac{z^{r_G\varepsilon}}{\varepsilon^{2n-4}}
\widehat B_G(\varepsilon)(1+z)^2
\int_{\Omega_2}D_1^{\varepsilon-1}D_2^{\varepsilon-1}
\widehat P_{\mathbf p}\widehat P_{\mathbf q}\,d^2t\,d^2s.
\tag{2.1}
$$

这里 Ω2={0<t1<t2<1,0<s1<s2<1} 是两个有序二单形的乘积。

更强的拓扑预算是

$$\sum_i|p_i|+\sum_i|q_i|=n-2.$$

原因：每条邻边恰好有一个端点落在支撑边组 {2,4}，另一个落在 {1,3}。因此真正依赖对边端点的字长只有 n−2，而不是 2n−4。

先取 ε>0。原域包含于完整端点立方体；邻边端点幂和每个独立 D 核均可积，故可以用绝对可积性建立此处的 Fubini 分组。prepared 本身不自带这个换序结论。

## 3. 六区间压缩为四个下端点字

同区间 shuffle、路径合成、逆路径与原笔记完全相同。对 0<x<y<1：

$$
\widehat F_q(x,y)=\sum_{q=ab}(-1)^{|a|}
\widehat F_{\operatorname{rev}a}(0,x)\widehat F_b(0,y),
$$
$$
\widehat F_r(y,1)=\sum_{r=cd}(-1)^{|c|}
\widehat F_{\operatorname{rev}c}(0,y)\widehat F_d(0,1).
$$

因此一条双标记边的统一压缩公式为

$$
\boxed{
\widehat F_p(0,x)\widehat F_q(x,y)\widehat F_r(y,1)
=\sum_{\substack{q=ab\\r=cd}}(-1)^{|a|+|c|}\widehat F_d(0,1)
\sum_{\substack{u\in p\shuffle\operatorname{rev}a\\
v\in b\shuffle\operatorname{rev}c}}
\widehat F_u(0,x)\widehat F_v(0,y).
}
\tag{3.1}
$$

shuffle 保留重数；字长不增加。两条边一起压缩，只剩四个带变量的字，其余进入 z 无关的常数系数。

这一步是结构归约，适用于任意字，不是对七张图列一张权重表。检查脚本已对总字长≤3的全部111种三间隔 A/B 插入，在非平凡 ε=2 下精确核对。

## 4. 只定义一次的通用 [2] 母族

定义

$$
\boxed{
\mathcal T^{[2]}_{\pi;\,u_1,u_2;\,v_1,v_2}(\varepsilon,z)
=(1+z)^2z^{-\varepsilon}
\int_{\Omega_2}
\prod_{j=1}^2D_z(t_j,s_{\pi(j)})^{\varepsilon-1}
\prod_{i=1}^2
\widehat F_{u_i}(0,t_i)\widehat F_{v_i}(0,s_i)
\,d^2t\,d^2s.
}
\tag{4.1}
$$

高圈改变的是四个字的允许长度，不改变核或插入规则。设 |u|=|u1|+|u2|+|v1|+|v2|，则

$$
\boxed{
J_G=\frac{z^{(r_G+1)\varepsilon}}{\varepsilon^{2n-4}}
\sum_{|\mathbf u|\le n-2}c^G_{\mathbf u}(\varepsilon)
\mathcal T^{[2]}_{\pi_G;\mathbf u}(\varepsilon,z).
}
\tag{4.2}
$$

c 的来源是完整区间有序字及 shuffle 的整数系数，既不依赖 z，也不从答案拟合。

可以继续沿用固定的非交换系统：

$$
\mathbf U_\varepsilon(x)=\sum_w\widehat F_w(0,x)X_w,\qquad
\mathbf U_\varepsilon'=\mathbf U_\varepsilon\,\varepsilon
\big[x^{\varepsilon-1}X_A+(1-x)^{\varepsilon-1}X_B\big].
$$

取四份张量槽，放在 t1,t2,s1,s2，乘统一双核并在 Ω2 积分，即生成式(4.1)全部系数。也可用
U(0,t1) X_O U(t1,t2) X_O U(t2,1)
直接生成六区间表示。四个槽保存端点位置，X_A/X_B 的非交换性保存 AB/BA。

### 4.1 交换、反射与有限覆盖

两条支撑边交换：

$$
\mathcal T^{[2]}_{\pi;\mathbf u;\mathbf v}
=\mathcal T^{[2]}_{\pi^{-1};\mathbf v;\mathbf u}.
$$

S2 中两个排列均为自身的逆。反射用
(t1,t2,s1,s2)→(1−t2,1−t1,1−s2,1−s1)，配对变为 rπr=π。
一个字反射为交换 A/B 后再逆序；六间隔边字的三段次序同时反转。将右区间再次按(3.1)压缩，得到同字长预算内的母核组合。故全四字母族在 crossing 下封闭；下面的纯 A 子族单独并不因此在 crossing 下封闭。

只扣除支撑边交换、不扣除其他积分恒等式，字长 m 的覆盖上界为

$$
\binom{m+3}{3}2^m+
\mathbf1_{m\ {\rm even}}\left(\frac m2+1\right)2^{m/2}
$$

（已经计入两个排列）。从 m=0 累加到 n−2，n=2,3,4,5,6 时得到 2、10、54、214、786。它们是候选覆盖数，不是物理图数、主积分数或新超越生成元数。

## 5. 从通用 J 连接到唯一物理 I

保持当前规范

$$G(\varepsilon)=e^{-\gamma_E\varepsilon}\Gamma(1-\varepsilon),\qquad
\kappa_G=\frac{\mathrm{sym}_Ga_{1,G}}{2(-2)^n}.$$

prepared C 已含一次 Euler 补偿，不再重复乘入。由

$$I_G=\kappa_GG^n\frac{J_G(z)+z^{n\varepsilon}J_G(1/z)}{1+z^{n\varepsilon}}$$

及(4.2)，令 δ_G=r_G−(n−2)/2，得到

$$
\boxed{
I_G=
\frac{\kappa_GG^n}{\varepsilon^{2n-4}}
\sum_{\mathbf u}c^G_{\mathbf u}(\varepsilon)
\frac{e^{\delta_G\varepsilon\ell}\mathcal T^{[2]}_{\pi;\mathbf u}(z)
+e^{-\delta_G\varepsilon\ell}\mathcal T^{[2]}_{\pi;\mathbf u}(1/z)}
{2\cosh(n\varepsilon\ell/2)},\quad\ell=\ln z.
}
\tag{5.1}
$$

相较 [1]：外部 ε^{-2n+2} 改为 ε^{-2n+4}，中心化 z^{-ε/2} 改为 z^{-ε}，δ 从 r−(n−1)/2 改为 r−(n−2)/2。三个改变都由两个对边核推导，不按圈数猜测。

令 V± 为包含 c_G 的母核组合的偶/奇部分，(5.1)等价于

$$I_G=\frac{\kappa_GG^n}{\varepsilon^{2n-4}}
\frac{\cosh(\delta_G\varepsilon\ell)V^+
+\sinh(\delta_G\varepsilon\ell)V^-}{\cosh(n\varepsilon\ell/2)}.$$

这给出任意目标 Laurent 系数的有限卷积规则。裸写 ε^{-2n+4} 只是统一的极点上界；具体边字和差分可能消去若干阶。当前六张三圈图从母核差分的 ε¹ 开始，实际只有简单极点。

## 6. 与原单对边 T、M、τ 的结构关系

原 T_{a;c} 和新的四字母核满足区域分割恒等式：

$$
\boxed{
T_{a;c}T_{b;d}
=
\mathcal T^{[2]}_{\mathrm{id};a,b;c,d}
+\mathcal T^{[2]}_{\mathrm{id};b,a;d,c}
+\mathcal T^{[2]}_{\mathrm{swap};a,b;d,c}
+\mathcal T^{[2]}_{\mathrm{swap};b,a;c,d}.
}
\tag{6.1}
$$

证明：把两个独立单对边积分的四维区域按 t、s 两组大小关系分成四块；各块重排到 Ω2，权重字随对应端点移动。单对边核的两个负号相乘为正，中心化乘积是 z^{-ε}，因此这里无额外符号或轨道重数。允许换序的可积性同§2。

空字给出

$$
\mathcal T^{[2]}_{\mathrm{id};0}+
\mathcal T^{[2]}_{\mathrm{swap};0}=\frac12H_0^2.
$$

这条和关系只控制组合，不能由“和属于旧函数代数”推出两个加数都属于旧代数。原 M/N 是 T 的系数矩，故它们通过(6.1)控制部分 [2] 组合；余下的配对与排序信息仍需保留。

### 6.1 配对的和与差：乘积部分与行列式核

两种配对还有直接的核层结构。令 Kε(t,s)=D_z(t,s)^{ε−1}，则

$$
K_\varepsilon(t_1,s_1)K_\varepsilon(t_2,s_2)
-K_\varepsilon(t_1,s_2)K_\varepsilon(t_2,s_1)
=\det\!\begin{pmatrix}
K_\varepsilon(t_1,s_1)&K_\varepsilon(t_1,s_2)\\
K_\varepsilon(t_2,s_1)&K_\varepsilon(t_2,s_2)
\end{pmatrix}.
$$

所以在**相同的端点权重**下，id+swap 是两项配对之和，id−swap 是有序二阶行列式核。这个区分直接来自传播子配对，不需要看求解答案。

这里的“配对之差”与 z↔1/z 的奇部分是两种不同操作。物理 crossing 同时反射两边并保持 id/swap，各自的行列式组合还要另外作 crossing 投影。

ε=0 时行列式甚至有显式几何因子：

$$
D_{11}D_{22}-D_{12}D_{21}
=z(t_2-t_1)(s_2-s_1),
$$
$$
\boxed{
\frac1{D_{11}D_{22}}-\frac1{D_{12}D_{21}}
=-\frac{z(t_2-t_1)(s_2-s_1)}
{D_{11}D_{12}D_{21}D_{22}}.
}
$$

它显式保留两边的端点间距。因 t1<t2、s1<s2、z>0，在0<ε<1时相同正权重的id核小于swap核；两者不能直接当作同一个母积分。

下面Q子族可相应分为

$$Q^{\rm sym}_{m,n}=Q_{\rm id;m,n}+Q_{\rm swap;m,n},\qquad
Q^{\rm det}_{m,n}=Q_{\rm id;m,n}-Q_{\rm swap;m,n}.$$

当m=n，区域分割固定 Qsym_m,m=H_m²/2；Qdet仍记录排序耦合。当m≠n，只有
Qsym_m,n+Qsym_n,m=H_mH_n
由乘积关系锁定，不能进一步删去另一种有序加权组合。零阶空字有
Qdet_00(0,z)=−M01²/6，
这是已求值的低阶恒等式，不提升为全ε比例。

这个行列式核及其有序插入，比按图号定义Σ更适合作为研究新函数来源的对象；它不自动证明有限项的最小基底维数。

## 7. 当前七图所需的重复 A 双指标子族

类似原笔记 T_m，定义**与图号无关**的双指标母族：

$$
\boxed{
\mathcal Q_{\pi;m,n}(\varepsilon,z)
=m!n!\mathcal T^{[2]}_{\pi;\varnothing,\varnothing;A^m,A^n}
=(1+z)^2z^{-\varepsilon}\int_{\Omega_2}
D_1^{\varepsilon-1}D_2^{\varepsilon-1}s_1^{m\varepsilon}s_2^{n\varepsilon}\,d^2t\,d^2s.
}
\tag{7.1}
$$

m,n≥0 表示较小、较大端点的重复 A 插入数，**不是圈阶，也不是 Laurent 阶数**。m!n! 抵消两个下区间字的有序阶乘。本批只需 (m,n)=(0,0),(0,1),(1,0)；任意高圈含混合 A/B 时仍使用§4的通用母族。

再令 H_m=m!T_m，则(6.1)给出

$$
\boxed{\mathcal Q_{\mathrm{id};m,n}+\mathcal Q_{\mathrm{id};n,m}
+\mathcal Q_{\mathrm{swap};m,n}+\mathcal Q_{\mathrm{swap};n,m}=H_mH_n.}
\tag{7.2}
$$

特别地 Q_id;00+Q_swap;00=H0²/2。当前在 q=0,1,2 上已用独立新算 Q 与旧 H 做完整 Hlog 恒等检查，不只比较最高 symbol。

### 7.1 不依赖未展开 ε 的系数定义

定义尺度中心化的对数与四维测度：

$$
d\mu_\pi=(1+z)^2\frac{d^2t\,d^2s}{D_1D_2},\qquad
\Phi_\pi=\ln D_1+\ln D_2-\ell,
$$
$$L_1=\ln s_1,\quad L_2=\ln s_2,\quad
\widetilde L_1=\ln(1-s_2),\quad\widetilde L_2=\ln(1-s_1).$$

于是

$$Q_{\pi;m,n}=\int_{\Omega_2}d\mu_\pi\,
e^{\varepsilon(\Phi_\pi+mL_1+nL_2)}.$$

反射同时保持 dμ、Φ，交换端点次序并把 L 换为上面的 tilde L。定义自然的系数函数：

$$
\boxed{\tau^{[2],\pi}_{m,n;q}=[\varepsilon^q]Q_{\pi;m,n}^+,\qquad
\eta^{[2],\pi}_{m,n;q}=\ell[\varepsilon^q]Q_{\pi;m,n}^-.}
\tag{7.3}
$$

显式多项式积分为

$$
\boxed{
\tau^{[2],\pi}_{m,n;q}
=\frac1{2q!}\int_{\Omega_2}d\mu_\pi
\left[(\Phi_\pi+mL_1+nL_2)^q
+(\Phi_\pi+m\widetilde L_1+n\widetilde L_2)^q\right],
}
$$
$$
\boxed{
\eta^{[2],\pi}_{m,n;q}
=\frac{\ell}{2q!}\int_{\Omega_2}d\mu_\pi
\left[(\Phi_\pi+mL_1+nL_2)^q
-(\Phi_\pi+m\widetilde L_1+n\widetilde L_2)^q\right].
}
\tag{7.4}
$$

q 是普通 Taylor 系数，分母 q! 已写在定义中，不再额外乘阶乘。两类函数均 crossing 偶；自然权重分别为 q+4、q+5。因此三圈有限项需要 τ 的 q=2 与 η 的 q=1，均为权重6。

若要分别读出各个矩，可以引入独立源
Qπ(ξ;λ1,λ2)=∫dμ exp(ξΦ+λ1L1+λ2L2)，
原 Q_mn 是 (ξ,λ1,λ2)=(ε,mε,nε) 的限制。对应矩
∫dμ Φ^a L1^b L2^c
由源导数唯一生成；q阶只是 a+b+c=q 的有限组合。这与原 M 的 W、L 双源方式平行，但保留了双核及排序。

这些积分定义解释结构，不能单凭定义宣称算完。此次另已把本批 Q 所需 q≤2 的系数重新积分为 Hlog；完整结果见 [结构基函数](structural-basis-functions.md)。

### 7.2 实际求值的三维形式

先积 t1；令 u=t2，v=s1，w=s2。id 取 s=v、r=w；swap 取 s=w、r=v。记

$$D_a=D_z(u,s),\ D_b=D_z(u,r),\ B=1-s,\ d=(1+z)s-1,$$
$$K_0=(1+z)^2\frac{\ln(D_a/B)}{dD_b},\quad
\Delta=\ln(D_a/B),\quad
K=\ln D_b+\tfrac12(\ln D_a+\ln B)-\ell+m\ln v+n\ln w.$$

则 q=0,1,2 的未对称 Q 系数分别是

$$
\int_{0<v<w<1,\ 0<u<1}K_0
\left\{1,\ K,\ \frac{K^2}{2}+\frac{\Delta^2}{24}\right\}\,du\,dv\,dw.
\tag{7.5}
$$

最后一项 Δ²/24 来自 sinh(x)/x，不能漏掉。d=0 为可去差商点，不作主值解释。

## 8. 七张图的通用母核归约

省略 Q 的 ε,z。定义三个**由邻边位置导出的差分**：

$$
\mathcal A=\mathcal Q_{\mathrm{id};0,0}-\mathcal Q_{\mathrm{id};0,1},\quad
\mathcal B=\mathcal Q_{\mathrm{id};0,1}-\mathcal Q_{\mathrm{id};1,0},\quad
\mathcal C=\mathcal Q_{\mathrm{swap};0,0}-\mathcal Q_{\mathrm{swap};0,1}.
$$

它们的来源依次是 s2之后的 A 端点、s1与s2之间的 A 端点、交换配对下s2之后的 A 端点。此处花体 A/B/C 是全 ε 的插入差分，不是旧解表中的极点系数 A(z)/B(z)/C(z)。

$$J_{4e-04}=z^\varepsilon Q_{\mathrm{id};0,0},$$
$$
J_{6e-10}=\frac{z^{2\varepsilon}}{\varepsilon^2}\mathcal A,\quad
J_{6e-11}=\frac{z^{2\varepsilon}}{\varepsilon^2}\mathcal B,\quad
J_{6e-12}=\frac{z^{2\varepsilon}}{\varepsilon^2}\mathcal C.
$$

另三图由 J23=z^{3ε}J10(1/z)、J25=z^{3ε}J11(1/z)、J24=z^{3ε}J12(1/z) 得到。
两图一组的 κ 相同，故完整 I 相等。各组 κ 依次为 −1/8、−1/16、−1/16；4e-04 为−1/8。

所有 Q_mn 的 q=0 不依赖 m,n；Qπ00 本身 self-mirror。因此差分的 q=0 消失，且所有 q 均有 ηπ00;q=0；三圈前面虽写 ε^{-2}，实际只有 ε^{-1}。

## 9. Σ 的结构来源：母族的有限差分系数

缩写 ρπ_mnq=τ^{[2],π}_{m,n;q}，χπ_mnq=η^{[2],π}_{m,n;q}。它们是上面统一定义的函数，不按图号定义。当前三个有限项恰为

$$
\boxed{\Sigma_A=
\frac{\rho^{id}_{01,2}-\rho^{id}_{00,2}}8
+\frac{\chi^{id}_{01,1}}{16},}
$$
$$
\boxed{\Sigma_B=
\frac{\rho^{id}_{10,2}-\rho^{id}_{01,2}}{16}
+\frac{\chi^{id}_{10,1}-\chi^{id}_{01,1}}{32},}
$$
$$
\boxed{\Sigma_C=
\frac{\rho^{swap}_{01,2}-\rho^{swap}_{00,2}}{16}
+\frac{\chi^{swap}_{01,1}}{32}.}
\tag{9.1}
$$

相应简单极点只需上式第一项把 q=2 改成 q=1，并去掉 χ 项。4e-04 的有限项是 −ρid_00,0/8。全表见 [structural-solutions.md](structural-solutions.md)。

χ 项来自 e^{±εℓ/2} 的一阶展开；G³ 与 sech(3εℓ/2) 没有一次项，所以不会另外改变这些有限项。不能把“系数本身 crossing 偶”误当作可以删掉奇系数乘 ln z。

最低关系：

$$
\rho^{id}_{mn,0}=\frac16M_{0,1}^2,\qquad
\rho^{swap}_{mn,0}=\frac13M_{0,1}^2,\qquad
\chi^\pi_{mn,0}=0.
$$

这将新母族与原 M 从最低层连接起来。高阶是否能回到原 M/N 或旧单对边 τ，仍须对完整系数作精确比较；不能由低阶关系外推全 ε 因子化。

因此 Σ 可以保留为读表方便的别名；理论主表示应是 Q 的系数差分。它们目前在旧单对边权重6空间之外只增加一个方向，是本批物理组合的关系结论，并不是整个 Q 母族所有系数只增加一个方向。

## 10. 可积性、字母结构与完成范围

当前重复 A 子族的每个系数是 dμ 乘有限对数多项式。D 的零点只在两个对角角点，局部 D 与两个非负距离之和可比；两个传播子的变量在完整立方体中独立。对数插入仍可积，限制为 Ω2 不会恶化。固定 m,n 和有限字长时，可取足够小的复 ε 圆盘作一致控制；未声称一个与任意 m,n 无关的统一半径。

对于一般四字母系数，端点递归产生0/1字母。依次积 t1、t2 时增加 D 的线性根；关键 resultant 是

$$\operatorname{Res}_u(D_z(u,v),D_z(u,w))=z(v-w).$$

其余因子为 z、1+z、v、1−v、w、1−w、v−w、(1+z)v−1、(1+z)w−1；后续 v∈(0,w)、w∈(0,1) 仍为线性字母，运动学碰撞位置为0、−1、∞。这一消元顺序解释了本批新系数为何仍是{0,−1} Hlog，而不要求把所有 Hlog 字都列成 ansatz。所有中间表观极点须在完整组合中相消。

以上提供任意有限字长继续求解的结构与线性字母分析；实际重新积分和数值检查只覆盖本批所用 Q、q≤2。没有给出一般高圈最小基底维数、全 ε 有限特殊函数闭式、有限维 canonical z 微分系统或物理复 z/i0 延拓，也没有 Lean 整链证明。

与原统一笔记相同，应分别判断：母族能否覆盖、目标系数是否求值、具体关系是否成立、最小生成元是否确定。此次推进的是前两者的结构连接，并把旧 Σ 的计算结果接回了统一定义。
