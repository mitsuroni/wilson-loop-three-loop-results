# [2]：同通道两条对边传播子的求解方法

本批从 prepared 0.3 的有序原参数根出发。结果为正实 z 的解析计算，formal: not_run。完整 ε 依赖以下给出显式收敛级数；Laurent 到有限项已独立用 HyperInt 求值为 Hlog。级数不是有限个 Gamma/Meijer G 的闭式，也不声称具有快速数值收敛。全文 ln 是正实主值。

## 1. 原积分及归一化

设 D_z(t,s)=(1-t)(1-s)+zts，S={0<t1<t2<1,0<v<w<1}。
id 配对为 D_z(t1,v)D_z(t2,w)，swap 配对为 D_z(t1,w)D_z(t2,v)。

定义中间积分（不是最终求值状态）：

$$
V_\pi[W](\varepsilon,z)=z^{-\varepsilon}(1+z)^2
\int_S D_z(t_1,s_{\pi(1)})^{\varepsilon-1}
D_z(t_2,s_{\pi(2)})^{\varepsilon-1}W(v,w)\,dt_1dt_2\,dv\,dw,
\quad(s_1,s_2)=(v,w).
$$

三种三圈权重：

$$W_A=1-w^\varepsilon,\qquad W_B=w^\varepsilon-v^\varepsilon.$$

从原排序导出：OOA 对应 W_A/ε，OAO 对应 W_B/ε；另一条独立边给 1/ε。BOO、OBO 则为 [1-(1-v)^ε]/ε、[(1-v)^ε-(1-w)^ε]/ε。

| 原编号 | 原边字 | 标准型 | κ=sym·a1/[2(-2)^n] |
|---|---|---|---|
| 4e-04 | OO / OO | id，W=1 | -1/8 |
| 6e-10 | OO / B / OOA | id，A | -1/8 |
| 6e-11 | OO / B / OAO | id，B | -1/16 |
| 6e-12 | OO / B / OOA | swap，A | -1/16 |
| 6e-23 | BOO / A / OO | id，A 的反射 | -1/8 |
| 6e-24 | BOO / A / OO | swap，A 的反射 | -1/16 |
| 6e-25 | OBO / A / OO | id，B 的反射 | -1/16 |

逐条原变量映射、配对及 C 检查存入结果的 reduction。交换两条边并同时反射 t1,t2→1-t2,1-t1 与 v,w→1-w,1-v，保持 S 和 id/swap 类型；D_z(1-t,1-s)=z D_{1/z}(t,s)。因此

$$J_{23}(z)=z^{3\varepsilon}J_{10}(1/z),\quad
J_{24}(z)=z^{3\varepsilon}J_{12}(1/z),\quad
J_{25}(z)=z^{3\varepsilon}J_{11}(1/z).$$

C 已有唯一的 Euler 补偿。令 G=e^{-γ_E ε}Γ(1-ε)，任何一条原 J 的物理组合是

$$I=\kappa G^n\frac{J(z)+z^{n\varepsilon}J(1/z)}{1+z^{n\varepsilon}}.$$

故三对 I 全 ε 相等。不能由这个数学反射跨图更换色因子；此处 κ 已逐图核对。

## 2. 换序及正则性依据

0<ε<1 时 D 核与邻边权重非负。把有序域扩大到两核独立的正方形，∫D_z^{ε-1} 有限；邻边全区间幂积分为 1/ε。故原被积函数绝对可积，Tonelli/Fubini 允许本节具体的换序和分组。此桥是下游论证，不是把 prepared 的根顺序默认改成乘积测度。

ε=0 附近，两独立端点角点是 (t,s)=(1,0),(0,1)，D 与两距离之和可比。取 |ε|≤δ<1/4，|D^{ε-1}|≤C D^{-1-δ}，并有

$$
\left|\frac{1-w^\varepsilon}{\varepsilon}\right|
\le |\ln w|w^{-\delta},\qquad
\left|\frac{w^\varepsilon-v^\varepsilon}{\varepsilon}\right|
\le |\ln v|v^{-\delta}.
$$

在例如 (t,s)=(1,0) 先积 t 后只产生 s^{-δ}，乘额外 s^{-δ}|ln s| 仍可积；其余角点同理。删去排序约束后用两个独立二维核作共同控制，即给出复 ε 小圆盘的一致可积界。有限次 ε 微分新增的对数不改变可积性。

因此 V_id[1] 在 ε=0 正则；三圈 V[W] 从 ε¹ 开始。三圈唯一独立邻边端点仍留下 1/ε，而不是 1/ε²。配上 G 和 cosh 后，负二次及更高极点严格为零。简单极点系数是原域正函数积分乘负 κ，在 z>0 非零。

## 3. 三维对数核：实际用于 HyperInt

先积 t1∈(0,u)，余下 u=t2、0<v<w<1。对 id，取 s=v、r=w；对 swap，取 s=w、r=v。记

$$
D_1=D_z(u,s),\quad D_2=D_z(u,r),\quad B=1-s,\quad d=(1+z)s-1.
$$

精确切片为

$$
\int_0^u D_z(t,s)^{\varepsilon-1}dt
=\frac{D_1^\varepsilon-B^\varepsilon}{\varepsilon d}.
$$

d=0 是可去点，其值为 u B^{ε-1}。差的整体必须保留；不能把单项当成独立收敛的内部极点积分。

设 V[W]=ε U1+ε² U2+O(ε³)，则

$$
K_0=(1+z)^2\frac{\ln(D_1/B)}{dD_2},\qquad
K_{\log}=\ln D_2+\tfrac12(\ln D_1+\ln B)-\ln z,
$$
$$
U_1=\int_0^1du\int_0^1dw\int_0^w dv\,K_0W_1,\qquad
U_2=\int_0^1du\int_0^1dw\int_0^w dv\,K_0(W_2+W_1K_{\log}),
$$

其中 A 型 W1=-ln w，W2=-(ln w)²/2；B 型 W1=ln(w/v)，W2=[(ln w)²-(ln v)²]/2。

计算使用 a=1/(1+z)，D_a=a(1-u)(1-s)+(1-a)us，并严格按 v→u→w 做有限区间积分。新鲜 Maple 会话启用 period 表、发散检查和遇发散中止；全部输出最后转到 z 字母 {0,-1}。没有从旧稿的数值拟合补定积分常数。

4e-04 另从直接二维差商积重新计算 V0，得到

$$V_0(z)=\frac16M_{0,1}(z)^2.$$

旧稿零阶与当前根吻合；旧稿的全阶纯 M 声明不作为推导依据。

## 4. 完整 ε 的显式收敛级数

为避免仅把原积分改名为“完整解”，这里对所有参数积分逐项求值。记

$$
Q_z(t,s)=1-\frac{D_z(t,s)}{1+z}
=\frac{z+t+s}{1+z}-ts,\qquad
Q_z(t,s)^k=\sum_{a,b=0}^k p_k(a,b;z)t^a s^b .
$$

这些系数是有限有理运算：

$$
p_k(a,b;z)=
\sum_{h=\max(0,a+b-k)}^{\min(a,b)}
\frac{k!(-1)^h z^{k-a-b+h}}
{(k-a-b+h)!(a-h)!(b-h)!h!(1+z)^{k-h}}.
$$

给定 a,b,c,d，id 令 β=b、δ=d；swap 令 β=d、δ=b。定义三个已积完的单项式矩：

$$
R_0=\frac1{(a+1)(a+c+2)(\beta+1)(\beta+\delta+2)},
$$
$$
R_A(\varepsilon)=
\frac1{(a+1)(a+c+2)(\beta+1)(\beta+\delta+2)(\beta+\delta+\varepsilon+2)},
$$
$$
R_B(\varepsilon)=
\frac1{(a+1)(a+c+2)(\beta+1)(\beta+\varepsilon+1)(\beta+\delta+\varepsilon+2)}.
$$

R_A、R_B 分别来自 ∫v^β w^δ W_A/ε 与 ∫v^β w^δ W_B/ε；ε=0 已是可去延拓，未保留相消的 1/ε。

令 (x)_k=Γ(x+k)/Γ(x)，取 k=0 时为 1。最终定义的是显式和式：

$$
\boxed{
F_{\pi,X}(\varepsilon,z)=
\left(\frac{(1+z)^2}{z}\right)^\varepsilon
\sum_{k,l\ge0}\frac{(1-\varepsilon)_k(1-\varepsilon)_l}{k!\,l!}
\left[
\sum_{a,b=0}^{k}\sum_{c,d=0}^{l}
p_k(a,b;z)p_l(c,d;z)R_X
\right].
}
$$

X=0,A,B；此处没有剩余积分，也没有待定边界常数。V_id[1]=F_id,0，V_pi[W_X]=ε F_pi,X。这是全 ε 收敛级数表示，不声称为有限特殊函数闭式。

**求和顺序是定义的一部分。** 对每个 (k,l) 先完成有限括号，再对 k,l 求和。不能把有限括号展开成一个任意换序的带符号多重无穷和。

收敛依据：对固定 z>0，0<Q_z≤1，等于 1 仅在零测的角点。广义二项式给 D^{ε-1}=(1+z)^{ε-1}∑(1-ε)_k Q^k/k!。实 0<ε<1 时，含 W_X/ε 的每个积分块非负，单调收敛；小复圆盘上 |(1-ε)_k|≤(1+δ)_k，§2 的 D^{-1-δ} 和权重界支配各块的绝对值。于是**分块双级数绝对且局部一致收敛**，可以逐阶 ε 展开。有限括号内的多项式系数可有符号，不影响这个按积分块证明的结论。

物理 I 因而为

$$
\boxed{I_{4e-04}=-\frac{G^2}{8\cosh(\varepsilon\ell)}F_{\mathrm{id},0}(\varepsilon,z),\quad\ell=\ln z,}
$$

$$
\boxed{
I_{\pi,X}=\frac{\kappa G^3}{2\varepsilon\cosh(3\varepsilon\ell/2)}
\left[e^{\varepsilon\ell/2}F_{\pi,X}(\varepsilon,z)
+e^{-\varepsilon\ell/2}F_{\pi,X}(\varepsilon,1/z)\right].
}
$$

这里 id,A 对应 10/23、id,B 对应 11/25、swap,A 对应 12/24。原始 J 分别为 J4=z^ε F_id,0；J10/11/12=z^{2ε}F/ε；其反射图 J23/25/24=z^ε F(ε,1/z)/ε。

实现 full_series.py 提供每个有限块的精确有理运算。正整数 ε 时 (1-ε)_k 截断，与原顺序 J 作精确比较；非整数有限块和逐级部分和另外数值检查。整数 ε≥1 的检查只针对 J/数学核，不在 Gamma 极点处评价物理 I。

## 5. Laurent 与函数空间

写 Uq±(z)=[Uq(z)±Uq(1/z)]/2。三圈系数由已计算 Hlog 直接得到

$$c_{-1}=\kappa U_1^+,\qquad
c_0=\kappa\left(U_2^++\frac{\ell}{2}U_1^-\right).$$

G³ 与 sech 没有 ε¹ 项，因此不会额外贡献有限项。这一步遗漏 crossing 的奇部分乘 ℓ 会算错 c0。

完整 Hlog 与 M/N 极点见 [solutions.md](solutions.md)、[basis-functions.md](basis-functions.md)。固定权重、crossing 偶、有理 Hlog/ζ坐标的精确线性代数表明：有限项不在此前列出的 M、M+ln z、单对边扩展 τ 多项式空间中。这里只排除被实际测试的有限空间，不是排除任意 M 函数、任意系数域或所有可能的周期恒等式。

旧单对边空间权重6秩13，加本批三个有限项后秩14；选择 c0(6e-10)=Σ_A 为新增生成元已足够。Σ_B、Σ_C 是另外两种有限项的易读记号，均可在旧 τ+Σ_A 中精确展开。它们都附已经求值的 Hlog，不能将符号定义本身当作完成计算。

## 6. 检查与尚未声称的内容

实际运行及哈希见 [检查入口](../../../../workflow/checks/same-channel-two/README.md)。保留原稿规范冲突；不改写本地参考快照或前35个结果。未实施 Kernel/Lean 功能、未注册新根、未执行形式全链。物理复 z/i0 延拓、更高正幂的显式系数和有限特殊函数全 ε 压缩均未完成；局部级数表示本身不覆盖这些额外目标。
