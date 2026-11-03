# 六类代表的归一化与共同矩形构件

2026-10-16。这里分别回答两个问题：旧权重七、八数据能否按现行物理约定联合比较；D 的矩形构件能否从原积分接到纯三对边以外的图。第一项结论是六个代表均无需变号或重新卷积。第二项得到一个完整 ε 的共同构件连接，但不等于证明 D 的三指标 X 已覆盖其他图。`formal=not_run`。

## 1. 六个比较对象使用同一个完整 I 约定

代表为 6e-10、11、12、20、21、36，物理相等图为 10=23、11=25、12=24、21=22、36=39。使用 [2026-09-29 的 targets2](../../../../../workflow/checks/mn-primitive-high-order/targets2.json) 中 `physical[id]["1"]` 与 `physical[id]["2"]`，其权重分别为七、八。

这些向量已经包含单次 Euler 补偿、Gamma³、颜色/轨道权重与 cosh 分母，不是裸几何源系数。[共同符号账本](../../../../../workflow/checks/shared-sign-conventions/corpus-audit.json) 的六个对应行均为 `physical_to_stored_C_ratio=1`；十月的其他族符号迁移不应套到这些纯交换图。

本次重新做了 131 项精确向量相等检查：旧两个高阶版本的共同阶、六个已发布完整 I 的全部低阶、五对相等轨道的低阶，以及 20/21 在十月十六日独立重建的零至二阶。另核对旧目标绑定的 89 次来源哈希，六个错误整体负号均被拒绝。[当前适配器](../../../../../workflow/checks/six-e-common-blocks-20261016/inputs/target-adapter.json) 是恒等映射，没有复制并重写旧结果。39 的独立发布修订为 r001，而 36 为 r002，修订号不能按配对推定。

## 2. 已求值矩形是一项可重复使用的低层构件

沿用 [D 的矩形定义](../d-intrinsic-source-20261016/derivation.md)，记

\[
h_z(t,s)=(1-t)(1-s)+zts,\qquad
k_e(t,s)=-(1+z)z^{-e/2}h_z(t,s)^{e-1},
\]
\[
G_e(t,s)=\int_t^1dx\int_0^sdy\,k_e(x,y),\qquad f_e(s)=G_e(0,s).
\tag{B1}
\]

对正实 z，完整 G 已经由显式四个 J 的组合求值：

\[
G_e(t,s)=[(1-t)s]^{e/2}\mathcal B_e(\eta(t,s),\kappa(t,s)),
\]
\[
\eta(t,s)=\frac12\log\frac{1-t}{zs},\qquad
\kappa(t,s)=(1+z)\sqrt{\frac{(1-t)s}{z}},
\tag{B2}
\]
\[
\mathcal B_e(\eta,\kappa)
=-\frac{\kappa^{-e}}e
\{J_e(1-\kappa e^{-\eta})+J_e(1-\kappa e^\eta)
-J_e[(1-\kappa e^{-\eta})(1-\kappa e^\eta)]-J_e(1)\}.
\tag{B3}
\]

因此 G 的引入不留下一个未经求值的新矩形。边界 t=0 直接给出

\[
\boxed{f_e(v)=v^{e/2}\mathcal B_e\left(-\frac12\log(zv),
(1+z)\sqrt{\frac vz}\right).}
\tag{B4}
\]

设 q=(1+z)v−1，这时 B3 的两个单端点是 −q、−1/z，其乘积为 q/z；又有 v^(e/2)κ^(−e)=ρ^(e/2)。故 B4 与 [旧平方源的 f 闭式 S2-14](../../../families/same-channel-two/three-loop.md) 精确相同。这条连接指出了旧显式 ρ 的几何来源，并将它吸收进已独立定义的矩形尺度。

## 3. 三种双对边插入共用 G 的边界导数

先在小正实 e 的绝对收敛域工作，保留原排序域。由 B1 直接有

\[
-\partial_tG_e(t,s)=\int_0^s k_e(t,y)dy,\quad
\partial_sG_e(t,s)=\int_t^1 k_e(x,s)dx,\quad
f'_e(s)=\int_0^1 k_e(x,s)dx.
\tag{B5}
\]

这些导数也都是初等线性幂差。原 [2] 族两种配对的降维核因此恰好为

\[
\boxed{K_{\rm id}=-G_tG_s,\qquad
K_{\rm sw}=-G_t(f'_e-G_s).}
\tag{B6}
\]

令 u_e(v)=(1−v^e)/e，则三种图的原插入源是

\[
\begin{aligned}
F_A&=\int_0^1dt\int_0^1ds\,u_e(s)(-G_tG_s),\\
F_B&=\int_0^1dt\int_0^1ds\,\frac{s^e-t^e}{e}(-G_tG_s),\\
F_C&=\int_0^1dt\int_0^1ds\,u_e(s)[-G_t(f'_e-G_s)].
\end{aligned}
\tag{B7}
\]

它们分别组装 10/23、11/25、12/24。中间一行仍是两个分别换元后的插入之差，不能将它误作一般双 Mellin 源的原有序变量。

相加第一、三行，先积 t：因为 G(1,s)=0，∫(−G_t)dt=f_e(s)。随后对 s 分部积分，两端的 u_e f² 均消失，得到

\[
\boxed{F_A+F_C=\frac12\int_0^1\frac{dv}{v}
v^{2e}\mathcal B_e\left(-\frac12\log(zv),(1+z)\sqrt{v/z}\right)^2.}
\tag{B8}
\]

这是一项**跨图共享的一维平方块**：来自两个不同配对图的完整源组合，使用与 C、D 同一个 B。它比仅比较有限阶余向量提供了更强的原积分连接。

若记 Z=e^(−γ_E e)Γ(1−e)、c₃=cosh(3e log z/2)，以及 ⟨F⟩ₑ=[z^(e/2)F(z)+z^(−e/2)F(1/z)]/2，则完整物理组合为

\[
\boxed{I_{10}+2I_{12}=-\frac{Z^3}{16e c_3}
\left\langle\int_0^1\frac{dv}{v}v^{2e}\mathcal B_e(\eta_0(v),\kappa_0(v))^2\right\rangle_e.}
\tag{B9}
\]

从 B8 可以直接生成带指标矩，无需先拟合 Hlog。记 B_n=∂ₑⁿBₑ|₀，则

\[
[e^q](F_A+F_C)=\frac12\sum_{m+n+p=q}
\frac{2^m}{m!n!p!}\int_0^1\frac{dv}{v}(\log v)^m
B_n(\eta_0,\kappa_0)B_p(\eta_0,\kappa_0).
\tag{B10}
\]

在完整 I 中，crossing 指数同样可在展开前并入 2log v±log z/2 的内部权重。这里先记录候选矩的自然来源，不另行宣布一个已独立、最小的新函数族。

## 4. 36/39 是同一矩形的双端点 Mellin 插入

为了辨认传播子指数与端点指数的不同作用，先释放两个端点参数 a、b，定义

\[
E_e(a,b;z)=\int_0^1dt\int_0^1ds\,
k_e(t,s)\frac{1-t^a}{a}\frac{1-s^b}{b}.
\tag{B11}
\]

在初始 Re a、Re b>0 的收敛域，利用 k=−∂t∂sG，先对 t 再对 s 分部积分，得到

\[
E_e(a,b;z)=\frac1a\int_0^1s^{b-1}f_e(s)ds
-\int_0^1dt\int_0^1ds\,t^{a-1}s^{b-1}G_e(t,s).
\]

再用 1/a=∫₀¹t^(a−1)dt，合成无表观端点极点的形式：

\[
\boxed{E_e(a,b;z)=\int_0^1dt\int_0^1ds\,
t^{a-1}s^{b-1}[f_e(s)-G_e(t,s)].}
\tag{B12}
\]

必须保留 f 边界。差 f−G 是已知 k 在 [0,t]×[0,s] 上的积分，在 t=0 和 s=0 均消失。这个完整差允许继续到原端点差在零参数附近的解析域；不能单独舍弃发散的边界项。

物理 36/39 的正则几何源就是 a=b=e；以 [原 S1 的 R](../../../families/single-opposite/three-loop.md#regular-source) 记号，E_e(e,e)=−(1+z)z^(−e/2)R。因而

\[
I_{36}=I_{39}=-\frac{Z^3}{64e^2c_3}
\{z^e E_e(e,e;z)+z^{-e}E_e(e,e;1/z)\}.
\tag{B13}
\]

把 B2、B4 代入 B12，得到同一 B 的一项边界值减一项内部值；端点指数的 Taylor 导数产生 log t、log s 插入。它说明 36/39 与 C、D 的确共用矩形低层构件，但它们采用不同的外层积分操作。

## 5. 共同构件与共同主积分空间仍是两层问题

已建立的全参数连接包括：C 使用两个互补矩形的乘积；D 使用两个线性幂差乘一个矩形；[2] 的三种插入使用矩形的边界导数，并在 F_A+F_C 中成为一维边界矩形平方；36/39 使用矩形的双端点 Mellin 差。

由此可以合理地冻结候选积分操作，再比较它们的系数。但**“都由 B 构造”不蕴含“都落在 X 的三指标切片”**。X 固定了 D 的 cosh 差、外层权重和投影；平方、边界导数、独立 Mellin 插入会引入并未包含在该定义中的方向。共同主积分 block 必须在这些操作的精确约化或跨权重字典中识别，不能单凭共享单传播子就宣布全族闭合。

本次 [实际核验](../../../../../workflow/checks/six-e-common-blocks-20261016/inputs/rectangle-bridge-verification.json) 包括 19 项 SymPy 精确恒等式，其中 e=1、2、3 的平方公式分别与原四维有序域直接多项式积分比较；端点公式还使用独立 a、b。六个错误平方系数/遗漏边界控制被拒绝。四组非整数 e、z、a、b 的原核与 B 表示共八项比较通过，112 节点最大观测差 2.07×10⁻¹⁴，并记录 64→112 精化。这些是经验浮点检查，不是区间证书，也不代替正文中全参数的换元和边界推导。
