# 从固定权重的 τ、E 到两指标 Euler 矩

本页把已经独立求值的自由 Beta 源，整理成便于跨权重追踪的两个矩族。来源是[共同 Euler 模块及 G13—G16](../../../../theory/gluonic-beta-generators.md)，而非用 τ、E 的余项定义新函数。用户此前先从跨图余项发现 N、随后倒推生成源的研究逻辑保持；这里把“发现线索”和“建立原积分定义”分开记录。

结论是：**τ、E 的两个权重六方向，分别对应一个权重六的镜像偶矩，以及一个权重五的镜像奇矩乘上 log z。它们可以归入同一个两参数母源的不同矩与镜像投影。** 这没有证明完整 C/D 只需要这个母源；高阶覆盖由本目录的另一项研究单独检验。

正实分支 z>0，ell=log z；所有矩取普通参数导数，不含阶乘。旧 C/J 和物理归一化不变，formal=not_run。

## 1. 先定义母源，再取指标

记

$$
r_z=\frac z{1+z},\quad
W=\log[\beta(1-\beta)],\quad
R=\log\frac{1-\beta}{\beta},\quad L=\ell+R,
\quad S_u=\log[u(1-u)],
$$

$$
D_+=1-u\beta,\qquad D_-=1-u+u\beta,
\qquad K_z(\beta)=\frac{L}{\beta-r_z}.
$$

原 Euler 模块为

$$
\mathcal F_{p,p,c}(w)
=\int_0^1[u(1-u)]^p(1-u+uw)^{c-1}\,du.
$$

在小参数的共同绝对收敛邻域，例如 |Re p|、|Re c|<1/4，定义一个单侧母源

$$
\boxed{\mathcal P(p,c;z)
=\int_0^1 K_z(\beta)\,\beta\,
\mathcal F_{p,p,c}(1-\beta)\,d\beta.}
\tag{TE1}
$$

它是已知 Euler/Gauss 函数的一维外层积分；若展开 Euler 定义，则是指定先 u 后 beta 的二重积分。参数 p 记录内层两端共同 Beta 权重，c 记录内层线性分母的指数。外层表观极点是可去点：

$$
L(r_z)=0,\qquad
\lim_{\beta\to r_z}K_z(\beta)=-\frac{(1+z)^2}{z}.
\tag{TE2}
$$

因此这里不需要主值处方。端点的对数奇性可积；对上述小参数邻域，角点 u→1、beta→1 的密度由 v^p(v+t)^{c-1} 乘对数控制，故参数微分可以由可积主控实施。

用同一个母源的两个镜像投影定义

$$
\mathcal T(p,c;z)=\mathcal P(p,c;z)+\mathcal P(p,c;1/z),
\quad
\mathcal S(p,c;z)=\mathcal P(p,c;z)-\mathcal P(p,c;1/z).
\tag{TE3}
$$

注意这里不除以 2。取原始幂矩

$$
\boxed{T_{m,s}=\left.\partial_p^m\partial_c^s\mathcal T\right|_{p=c=0},
\qquad S_{m,s}=\left.\partial_p^m\partial_c^s\mathcal S\right|_{p=c=0}.}
\tag{TE4}
$$

等价的固定积分定义为

$$
\begin{aligned}
T_{m,s}&=\int_0^1 K_z(\beta)\,d\beta\int_0^1S_u^m
\left[\frac{\beta\log^sD_+}{D_+}
+\frac{(1-\beta)\log^sD_-}{D_-}\right]du,\\
S_{m,s}&=\int_0^1 K_z(\beta)\,d\beta\int_0^1S_u^m
\left[\frac{\beta\log^sD_+}{D_+}
-\frac{(1-\beta)\log^sD_-}{D_-}\right]du.
\end{aligned}
\tag{TE5}
$$

指标 m,s≥0。增大任一指标都对应在同一个原积分中插入一个明确对数，因此不是按权重给余向量依次命名。Taylor 展开与原始矩的阶乘关系为

$$
\mathcal T=\sum_{m,s\ge0}\frac{p^m c^s}{m!s!}T_{m,s},
\qquad
\mathcal S=\sum_{m,s\ge0}\frac{p^m c^s}{m!s!}S_{m,s}.
\tag{TE6}
$$

它们的均匀权重为 m+s+3：外层 L 提供一个对数，内外两次对数型积分各增加一。与此前六指标矩的精确字典是

$$
T_{m,s}=\mathfrak B^{\rm free}_{0,1;m,0,s,0},
\qquad
\ell S_{m,s}=\mathfrak B^{\rm free}_{0,1;m,0,s,1}.
\tag{TE7}
$$

第二式用来比较旧记号；S 的定义始终是 TE5，不能在 z=1 以除以 ell 的方式定义。

## 2. 镜像奇偶与哪些低指标会退回 M

同时作 z→1/z、beta→1−beta，L 与 beta−r_z 都变号，故 K 不变；D_+、D_- 交换。于是

$$
T_{m,s}(1/z)=T_{m,s}(z),\qquad
S_{m,s}(1/z)=-S_{m,s}(z).
\tag{TE8}
$$

特别是 S_{m,s}(1)=0；ell S_{m,s} 镜像偶。两个族源自同一个 P，奇偶投影并不意味着两个不同的内部 Euler 核。

m=0 时，内层可以立即积掉：

$$
\int_0^1\frac{\beta\log^s(1-u\beta)}{1-u\beta}\,du
=-\frac{\log^{s+1}(1-\beta)}{s+1}.
\tag{TE9}
$$

令原 M 矩为 M_{k,j}=∫_0^1 W^k L^j/(beta−r_z) d beta。用 log(1−beta)=(W+R)/2、log beta=(W−R)/2、R=L−ell，得到所有 s 的闭式

$$
\boxed{
\begin{aligned}
T_{0,s}
&=-\frac1{(s+1)2^s}
\sum_{\substack{0\le v\le s+1\\v\ \text{even}}}
\binom{s+1}{v}
\sum_{h=0}^{v}\binom vh(-\ell)^{v-h}M_{s+1-v,h+1},\\
S_{0,s}
&=-\frac1{(s+1)2^s}
\sum_{\substack{0\le v\le s+1\\v\ \text{odd}}}
\binom{s+1}{v}
\sum_{h=0}^{v}\binom vh(-\ell)^{v-h}M_{s+1-v,h+1}.
\end{aligned}}
\tag{TE10}
$$

例如 T_{0,0}=−M_{1,1}、S_{0,0}=−M_{0,2}+ell M_{0,1}。因此不是本族的每一个矩都新：真正要追踪的是内层 Beta 对数 m>0 打开后留下的方向。这提供了族与旧 M 的明确退化联系。

## 3. τ、E 的完整权重六字典

令 τ=τ^+_{D,4}。由已有直接积分结果独立重放，得到

$$
\boxed{\tau=-3T_{1,2}-2\ell S_{1,1}+R_\tau^{MN},
\qquad E=-\ell S_{1,1}+R_E^{MN}.}
\tag{TE11}
$$

这里没有隐藏拟合别名；旧项完全为

$$
\begin{aligned}
R_\tau^{MN}={}&-\frac{11}{8}N_{2,2}+\frac14M_{1,1}^2
+\frac5{64}M_{2,3}+\frac{29}{288}M_{0,1}^3
+\frac{31}{128}M_{4,1}\\
&-12M_{1,1}\zeta_3-\frac{459}{40}M_{0,1}\zeta_2^2
-\frac{23}{48}M_{0,1}^2\zeta_2
+\frac38M_{0,1}M_{2,1}+\frac{39}{8}M_{2,1}\zeta_2,\\[2pt]
R_E^{MN}={}&-\frac14M_{2,3}+\frac14N_{2,2}
+\frac3{20}M_{0,1}^3-\frac{79}{60}M_{0,1}^2\zeta_2
-\frac{53}{10}M_{0,1}\zeta_2^2\\
&-\frac34M_{2,1}\zeta_2-\frac14M_{0,1}M_{2,1}.
\end{aligned}
\tag{TE12}
$$

也可以单独读出两个种子：

$$
T_{1,2}=-\frac13\left[\tau-2E-(R_\tau^{MN}-2R_E^{MN})\right],
\qquad
\ell S_{1,1}=-E+R_E^{MN}.
\tag{TE13}
$$

这些是从已独立定义且积分的 T/S 到旧 τ/E 的字典。TE13 用于理解和核对，不作为新源定义。

最值得注意的是：**S_{1,1} 本身只有权重五，而且镜像奇；乘一个 ell 后才进入权重六的镜像偶结果。** 因而“物理偶结果首次在权重六超出 M/N”不等于“所有相关新矩都只能从权重六开始”。

## 4. 精确独立性检验及范围

本次另以所有登记到相应权重的原 M/N、ell、普通 zeta 的齐次乘积构成候选空间，允许有理常系数。完全独立的 Fraction/shuffle 消元给出：

| 权重和镜像性 | 原 M/N+ell 候选秩 | 加入的自然矩 | 增广秩 |
|---|---:|---|---:|
| 5，奇 | 5 | S_{1,1} | 6 |
| 6，偶 | 12 | T_{1,2}，ell S_{1,1} | 14 |

三个单目标非成员关系和权重六两个方向的联合独立性均保存了精确对偶泛函。程序逐列检查对所有候选为零、对目标非零。此处“独立”仅在声明的有理常系数、已登记周期空间内；不排除有理 z 系数、运动学微分或其他参数变换后的约化，更不证明全阶最小主函数数。

权重六字典的新重放直接读取旧两次独立 HyperInt 求积的 J/Jb 原输出，使用本次独立实现的 AST 有理数解释器和 shuffle 乘法；没有导入旧拟合器，也没有重新拟合 τ/E。全部旧结果仅只读复用。另对 TE10 的 s=0,1,2,3 共八个矩，与原求积输出逐系数精确比较。

[核验脚本](../../../../../workflow/checks/indexed-extra-sources-20261016/tau-e/verify.py)和[实际结果](../../../../../workflow/checks/indexed-extra-sources-20261016/tau-e/verification.json)记录 79 项通过：完整 τ/E 回代、低高阶求积任务重复、权重与阶乘、镜像和可去点、m=0 切片、非成员对偶和两个错误字典控制。复用 16 个原求积输出，19 个输入核对旧来源账本；本轮新积分运行数为零，formal=not_run。

这一步得到的是**一个可跨权重取矩、具有明确旧族退化及 τ/E 起始字典的自然源**。下一问题是提高 m,s 及必要外层权重后，能否用固定参数曲线连接完整 C/D；该问题不能由本页权重五、六的成功自动解决。
