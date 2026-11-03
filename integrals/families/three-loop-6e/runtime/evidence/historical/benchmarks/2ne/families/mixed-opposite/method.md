# [1,2] / [2,1]：从原始 J 到完整 I

日期：2026-09-17。范围为 prepared 0.3 的 6e-37、6e-38，目标延续前批：完整 ε 解析表示，以及 ε⁻⁶ 至 ε⁰、余项 O(ε)。本文是计算推导；formal=not_run，不扩大 Kernel 的根映射或证明能力。

## 1. 名称与实际拓扑

[1,2] 和 [2,1] 表示两个对边通道各有一条、两条传播子。作为不区分通道名称的族标签，两者相同；实际求解必须继续保留 13/24 通道、端点排序及配对。当前两条根均为 13 一条、24 两条。

定义
$$
D_z(t,s)=(1-t)(1-s)+zts,\qquad
\widetilde D_z(a,b)=z(1-a)(1-b)+ab.
$$

令 $a=s_{11},b=s_{31},t_1=s_{21},t_2=s_{22},s_1=s_{41},s_2=s_{42}$，原域为
$$
[0,1]^2_{a,b}\times\Omega_2,\qquad
\Omega_2=\{0<t_1<t_2<1,\ 0<s_1<s_2<1\}.
$$

两个不同的实际输入是
$$
J_\pi=-(1+z)^3
\int_{[0,1]^2\times\Omega_2}
\widetilde D_z(a,b)^{\varepsilon-1}
D_z(t_1,s_{\pi(1)})^{\varepsilon-1}
D_z(t_2,s_{\pi(2)})^{\varepsilon-1}\,d^6s,
$$
其中 6e-37 对应 $\pi=\mathrm{id}$，6e-38 对应 $\pi=(12)$，记作 swap。分子的负号来自原根，不能丢弃。

本批没有邻边传播子；三条传播子的六个端点互不相同。不同于含邻边端点字的三圈图，本批没有显式 $\varepsilon^{-2}$ 前因子。每图的原外到内顺序另保留在其 derivation.md 中。

## 2. 分组与无极点的依据

固定 $z>0$，每个 $D$ 或 $\widetilde D$ 仅在两个角点为零，且在相应局部非负坐标 $x,y$ 中可与 $x+y$ 上下比较。因此当 $\operatorname{Re}\varepsilon>-1$ 时，
$$
\int_0^\delta\!\int_0^\delta (x+y)^{\operatorname{Re}\varepsilon-1}\,dx\,dy<\infty.
$$

放宽有序域为完整六维立方体，绝对值积分分解为三个有限二维积分；这使用了六个端点没有重复。于是原嵌套积分可用绝对可积性与 Fubini 分组，而不是凭族标签假设因子化。

进一步取 $|\varepsilon|\leq\delta<1$，用 $D^{-1-\delta}$ 加常数界控制整个邻域；任意有限次 ε 导数产生的对数也可由略强的幂界吸收。因此 J 在 ε=0 附近全纯。G 与 crossing 分母在那里也解析且分母非零，I 同样全纯：
$$
c_{-6}=c_{-5}=\cdots=c_{-1}=0.
$$
这里给出数学论证，尚无对应的 Lean 全链收据。

## 3. 复用已有母核

沿用 [单对边定义](../single-opposite/basis.md)：
$$
H_0(\varepsilon,z)=-(1+z)z^{-\varepsilon/2}
\int_{[0,1]^2}D_z(t,s)^{\varepsilon-1}\,dt\,ds.
$$

沿用 [[2] 母族](../same-channel-two/unified-structure.md)，本批只需零插入的两项：
$$
Q_\pi(\varepsilon,z)
=Q_{\pi;0,0}
=(1+z)^2z^{-\varepsilon}
\int_{\Omega_2}
D_z(t_1,s_{\pi(1)})^{\varepsilon-1}
D_z(t_2,s_{\pi(2)})^{\varepsilon-1}\,d^4s.
$$

反射 $(a,b)\mapsto(1-a,1-b)$ 将 $\widetilde D_z$ 变成 $D_z$，积分域及 Jacobian 绝对值不变。因此
$$
\boxed{J_\pi(\varepsilon,z)=z^{3\varepsilon/2}H_0(\varepsilon,z)Q_\pi(\varepsilon,z).}
$$

这一步复用的是已有 H0 与 Q 的精确定义；没有新增 σ 或任意命名的有限项函数。

## 4. Crossing 与唯一 C/I 归一化

两个已有核在本批零插入时分别满足
$$
H_0(\varepsilon,1/z)=H_0(\varepsilon,z),\qquad
Q_\pi(\varepsilon,1/z)=Q_\pi(\varepsilon,z).
$$

Q 的变换具体为 $t_i\mapsto1-t_{3-i}$、$s_i\mapsto1-s_{3-i}$，再重编号；id、swap 各自保留。用 $zD_{1/z}(1-t,1-s)=D_z(t,s)$ 核对比例。这是运动学 crossing，不是交换 id 与 swap。

记
$$
G(\varepsilon)=e^{-\gamma_E\varepsilon}\Gamma(1-\varepsilon),\qquad \ell=\ln z.
$$
prepared C 已含一次 Euler 补偿，不再补乘。通用定义为
$$
I=\kappa G^3\frac{J(z)+z^{3\varepsilon}J(1/z)}{1+z^{3\varepsilon}},
\qquad \kappa=\frac{\mathrm{sym}\,a_1}{2(-2)^3}.
$$

| 积分 | 配对 | $a_1$ | sym | $\kappa$ |
|---|---|---:|---:|---:|
| 6e-37 | id | 1/2 | 2 | −1/16 |
| 6e-38 | swap | 1/4 | 2 | −1/32 |

所以完整结果为
$$
\boxed{I_{6e-37}=-\frac{G^3H_0Q_{\mathrm{id}}}{16\cosh(3\varepsilon\ell/2)},\qquad
I_{6e-38}=-\frac{G^3H_0Q_{\mathrm{swap}}}{32\cosh(3\varepsilon\ell/2)}.}
$$

## 5. 完整 ε 表达如何求值

母核表达保存结构，但交付不止于未求值参数积分。

H0 已有标准特殊函数表示
$$
H_0=-\frac{(1+z)z^{-\varepsilon/2}}
{\Gamma(1-\varepsilon)\Gamma(1+\varepsilon)^2}
G^{3,3}_{3,3}\!\left(z\;\middle|\begin{matrix}
\varepsilon,0,0\\0,\varepsilon,\varepsilon
\end{matrix}\right).
$$

Qid 复用 [已积完的双级数](../same-channel-two/method.md)。这里完整重述零插入式。令
$$
q_z(t,s)=1-\frac{D_z(t,s)}{1+z}
=\frac{z+t+s}{1+z}-ts,\qquad
p_k(a,b;z)=[t^as^b]q_z(t,s)^k.
$$
其显式有限和为
$$
p_k(a,b;z)=
\sum_{h=\max(0,a+b-k)}^{\min(a,b)}
\frac{k!(-1)^h z^{k-a-b+h}}
{(k-a-b+h)!(a-h)!(b-h)!h!(1+z)^{k-h}}.
$$

对 id 取 $(\beta,\delta)=(b,d)$，对 swap 取 $(\beta,\delta)=(d,b)$，则
$$
R_\pi(a,b,c,d)=
\frac{1}{(a+1)(a+c+2)(\beta+1)(\beta+\delta+2)},
$$
$$
\boxed{
Q_\pi=
\left(\frac{(1+z)^2}{z}\right)^\varepsilon
\sum_{k,l\ge0}\frac{(1-\varepsilon)_k(1-\varepsilon)_l}{k!l!}
\left[
\sum_{a,b=0}^{k}\sum_{c,d=0}^{l}
p_k(a,b;z)p_l(c,d;z)R_\pi(a,b,c,d)
\right].
}
$$

每个方括号内先完成有限多项式求和，再求外层双和。$0\le q_z<1$ 在域内部成立；固定正 z、ε=0 的紧邻域内，Pochhammer 系数的绝对值与 $(1+\delta)_k/k!$ 比较，所得支配函数是三个（或两个）独立端点核的可积幂。故按有限块组织的外层双和绝对、局部一致收敛。不能把这一结论误写成任意重排内部带符号单项式后的逐项绝对收敛。

还可直接复用严格的乘积域分割：
$$
\boxed{Q_{\mathrm{id}}+Q_{\mathrm{swap}}=\frac{H_0^2}{2}.}
$$
理由是两个无序 t 与两个无序 s 的乘积域分为四种排序；交换两条核后成对相等，分别给两个 id、两个 swap 区域。因此 6e-38 的全 ε 表达也可只用 H0 和已发布的 Qid 级数，不需计算新的全 ε 母函数。

这是“标准 Meijer G + 显式收敛级数”的解析表示，尚未给 Qid 寻得有限特殊函数闭式，也没有声称双级数是高效数值算法。数值核验使用独立的降维积分。

## 6. 所有极点与有限项

复用并精确反代既有 HyperInt 结果：
$$
H_0(0,z)=M_{0,1}(z)=-\frac{\ell^2+\pi^2}{2},\qquad
Q_{\mathrm{id}}(0,z)=\frac{M_{0,1}^2}{6},\qquad
Q_{\mathrm{swap}}(0,z)=\frac{M_{0,1}^2}{3}.
$$

结合 $G(0)=1$，得到两条相同的有限项：
$$
\boxed{
I_{6e-37}=I_{6e-38}
=-\frac{M_{0,1}(z)^3}{96}+O(\varepsilon)
=\frac{(\ell^2+\pi^2)^3}{768}+O(\varepsilon).
}
$$
这个连等式仅表示所示截断阶相同。纯 M、$\tau_{0,0}^3$、显式 log/π 多项式和完整 Hlog 坐标均保存；这里有限项不需要 τ^[2] 的更高阶、η、Σ 或新 E 方向。不能由此推断所有正幂系数也只含 M01。

## 7. 相同有限项不等于相同全函数

两条完整 I 满足
$$
I_{6e-37}+2I_{6e-38}
=-\frac{G^3H_0^3}{32\cosh(3\varepsilon\ell/2)}.
$$
但一般不相等。一个不依赖数值的检查是：Q 在 ε=1 有限且
$$
Q_{\mathrm{id}}(1,z)=Q_{\mathrm{swap}}(1,z)=\frac{(1+z)^2}{4z}.
$$
若完整 I 恒相等，则消去公共非零因子后应有 $2Q_{\rm id}=Q_{\rm swap}$，与此矛盾。此处只在 ε=1 评估正则的母核/J，不能直接在那里评估含 Gamma 极点的 I。

## 8. 检查与完成边界

[复算入口](../../../../workflow/checks/mixed-opposite/README.md)保存：

- 根多项式、六个端点、原 bounds、配对与 C 的精确检查，拒绝跨通道排序耦合等不适用输入。
- 六次原外到内 J 的整数 ε 精确积分，以及六次终止双级数与它们的相等性。
- 四个非整数参数案例的原六维域随机化 Sobol 对照。
- H0 的 Meijer G/一维求积、Q 的 crossing 与乘积恒等式对照。
- 有限项的 M/τ/log/Hlog 精确等价及独立求积；全 I 减有限项的 O(ε) 缩放。
- 输入、复用母核输出、源码、文档与结果修订的哈希和快照。

以上是计算检查；Sobol 误差与求积差分不属于严格数值区间。当前采用 $z>0$、实 $0<\varepsilon<1$，随后固定正 z 在 ε=0 展开。物理复 z/i0、更高正幂的交付、Lean/MCP 形式整链均未实施。
