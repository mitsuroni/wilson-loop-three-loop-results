# [1] 族的结果表示与公共求解依据

2026-09-16。本批范围是当前 prepared 0.3 的 18 个 pure 单对边积分（1 个一圈、2 个二圈、15 个三圈）。目标沿用完整 epsilon 的 I 和全部极点至有限项。这里记录本次求解采用的表示，不把一般高圈闭包猜想当作事实。

## 1. 表示结论

完整 epsilon 层保留有序母积分 T 的来源，并将本批实际需要的对称化组合求值为标准 Meijer G 函数。固定阶数以明确的母核系数函数展示，同时保存完整 Hlog 表达作可检查的后端。旧 M/N 作为经精确系数比较得到的替代表达；不强迫不能归约的项进入纯 M。

旧 tau_(m,q) 只覆盖单端重复 A 的 crossing 偶部分。混合 AB/BA、两端同时插入并不自动属于这个子族。理论中的三圈 10 个覆盖核也不是“10 个最小主函数”的证明。

主要本地依据：[统一理论 §9–17](../../../../references/local/2026-09-16/test0903/empty_and_single_opposite_unified_notes.md)、[gf41 M 求值脚本](../../../../references/local/2026-09-16/test0903/04_GF_work/gf41_M_weight7.mpl)、[gf51 两端插入脚本](../../../../references/local/2026-09-16/test0903/04_GF_work/gf51_two_end_words.mpl)。本批重新执行所需求值，历史脚本/结果只作方法与归一化对照。

## 2. 从当前原参数输入恢复排序

每个相邻传播子逐项重建为 q(1-u)v，q=1 或 z，分子为 -(1+z) 乘所有相邻 q；唯一对边传播子重建为

$$D_z(t,s)=(1-t)(1-s)+zts.$$

每个不同端点只出现一次。原 bounds 是互不相交的有序链；对边端点记为 O，相邻端点为 A(x)=x^(epsilon-1) 或 B(x)=(1-x)^(epsilon-1)。

先取 z>0、epsilon>0。完整立方体上，相邻权重积分有限；对边 D 仅在两个角点为零，局部与两个非负距离之和可比，二维 D^(epsilon-1) 可积。于是原有序域绝对可积，Tonelli/Fubini 允许将原迭代积分按链分组。原外到内顺序仍保存，并以整数 epsilon 的原顺序直接积分作独立连接检查。本次给出了数学换序依据，没有把 prepared 根误称为自带 Fubini 收据。

令 r 为相邻 q=z 的个数，P=2n-2，将每个相邻端点积分乘 epsilon 正规化，得到

$$
J=-(1+z)z^{r\varepsilon}\varepsilon^{-P}
\int_0^1dt\,ds\,D_z(t,s)^{-1+\varepsilon}W(t,s;\varepsilon).
$$

W 来自真实端点排序而非按图号填表。外部未标记链可能贡献
$Q=\Gamma(1+\varepsilon)^2/\Gamma(1+2\varepsilon)$ 等常数。

## 3. 端点权重与 BA 项的处理

记 a=x^epsilon、b=(1-x)^epsilon，B_A(x)=Fhat_BA(0,x)。则下区间权重为

$$
\widehat F_A=a,\quad\widehat F_B=1-b,\quad
\widehat F_{AA}=a^2/2,\quad\widehat F_{BB}=(1-b)^2/2,
$$
$$
\widehat F_{BA}=B_A,\qquad \widehat F_{AB}=a(1-b)-B_A.
$$

上区间由路径合成得到

$$
\widehat F_A(x,1)=1-a,\quad\widehat F_B(x,1)=b,
$$
$$
\widehat F_{AA}(x,1)=(1-a)^2/2,\quad\widehat F_{BB}(x,1)=b^2/2,
$$
$$
\widehat F_{AB}(x,1)=Q-a+B_A,\qquad
\widehat F_{BA}(x,1)=a+b-ab-Q-B_A.
$$

反射整个 BA 有序区，再在 x 处分割，得

$$
B_A(x)+B_A(1-x)=a+b-ab-Q.
$$

因此对称化的 BA 可以精确换成幂权重。程序逐图核对：本批含 BA 的权重只出现在 delta=r-(n-1)/2=0 的情形，对称化后 B_A 完全消失。这是本批的实际简化，不宣称任意高圈 BA 都能删除。

## 4. 六种共同幂权重及标准特殊函数求值

定义统一的已求值函数

$$
\begin{aligned}
\mathcal H_{ab;cd}(\varepsilon,z)
={}&-\frac{(1+z)z^{-\varepsilon/2}}
{\Gamma(1-\varepsilon)\Gamma(1+(a+b+1)\varepsilon)
 \Gamma(1+(c+d+1)\varepsilon)}\\
&\times G_{3,3}^{3,3}\!\left(z\;\middle|\;
\begin{matrix}\varepsilon,-a\varepsilon,-c\varepsilon\\
0,(b+1)\varepsilon,(d+1)\varepsilon\end{matrix}\right).
\end{aligned}
$$

这里的 G_(3,3) 是标准 Meijer G，不能与补偿 Gamma 因子 G(epsilon) 混淆。其定义采用 [DLMF 16.17.1](https://dlmf.nist.gov/16.17#E1)。

推导：对 D^(-1+epsilon) 使用 Mellin–Barnes 分解，积分 t、s 得到 Gamma 因子

$$
\Gamma(-w)\Gamma(1-\varepsilon+w)
\Gamma(1+a\varepsilon+w)\Gamma((b+1)\varepsilon-w)
\Gamma(1+c\varepsilon+w)\Gamma((d+1)\varepsilon-w),
$$

分母即上式三个 Gamma。先取 0<epsilon<1，围道的实部选在 epsilon-1 与 0 之间，左右极点分离。它精确匹配标准 Meijer G，故

$$
\mathcal H_{ab;cd}
=-(1+z)z^{-\varepsilon/2}
\int_0^1dt\,ds\,
D^{-1+\varepsilon}
t^{a\varepsilon}(1-t)^{b\varepsilon}
s^{c\varepsilon}(1-s)^{d\varepsilon}.
$$

右式提供来源和数值对照，左式才是最终全 epsilon 的标准函数求值。没有把未求值积分简单改名来宣称完成。

本批只需要六种：

| 简记 | (a,b;c,d) | 与旧母核的关系 |
|---|---|---|
| H0 | (0,0;0,0) | T_(empty;empty)=T0 |
| H1 | (0,0;1,0) | T_(empty;A)=T1 |
| H2 | (0,0;2,0) | 2 T_(empty;AA)=2 T2 |
| HW | (0,0;1,1) | T_(empty;A)−T_(empty;AB)−T_(empty;BA) |
| HD | (1,0;1,0) | T_(A;A) |
| HX | (1,0;0,1) | T_(empty;A)−T_(A;B) |

镜像满足 H_(ab;cd)(1/z)=H_(ba;dc)(z)；两端交换也对称。H0、HW、HX 因此 crossing 偶。六种足够本批求值，不代表最小独立数或一般 n 圈覆盖数。

## 5. 当前 C/I 与完整结果

记补偿因子 Z(epsilon)=exp(-EulerGamma epsilon) Gamma(1-epsilon)，kappa=sym*a1/[2(-2)^n]。当前 prepared C 已含且仅含一次补偿。令 delta=r-(n-1)/2，则

$$
I=\frac{\kappa Z^n}{\varepsilon^{2n-2}\cosh(n\varepsilon\ell/2)}
\frac{e^{\delta\varepsilon\ell}\mathcal T[W](z)
      +e^{-\delta\varepsilon\ell}\mathcal T[W](1/z)}2,\qquad\ell=\ln z.
$$

T[W] 是归约过程中的线性记号；程序完成 BA 合并并将其替换为上表的显式 Meijer G 函数。每图 JSON 同时保存紧凑六函数式和完全展开定义的 Meijer G 表达式。

I 先在 z>0、0<epsilon<1 建立，固定正 z 后作 epsilon=0 附近的局部亚纯延拓。物理复 z/i0 和本批 Lean 整链另列未执行，不与 CAS 结果混合。

## 6. 系数层：旧 tau 与本次扩展

对六个明确的 H 定义 H_i^±=(H_i(z)±H_i(1/z))/2：

$$
\tau_{0,q}=[\varepsilon^q]H_0,\quad
\tau_{1,q}=[\varepsilon^q]H_1^+,\quad
\tau_{2,q}=[\varepsilon^q]H_2^+.
$$

这与原文的 m![epsilon^q]T_m^+ 完全一致，特别是 H2=2T2，不能漏掉 2!。不额外乘 q!。

额外定义同一求值族的系数

$$
\tau^+_{W,q}=[\varepsilon^q]H_W,\quad
\tau^+_{D,q}=[\varepsilon^q]H_D^+,\quad
\tau^+_{X,q}=[\varepsilon^q]H_X,
$$
$$
\eta_{i,q}=\ell[\varepsilon^q]H_i^-,\qquad i=1,2,D.
$$

eta 的作用是吸收物理 crossing 中的“奇对数乘奇系数”，它们自身均 crossing 偶，不引入 M 分母。偶次 ell 使用
$\ell^2=-2\tau_{0,0}-6\zeta(2)$ 消去。最终外部系数可保持为有理数和 zeta 常数；运动学位于明确基函数内。

每个实际使用的 q≤4 系数都重新积分成 Hlog，不停留在“取某个未求值积分的系数”的定义。JSON 保存其字母 {0,-1}、精确有理系数和 period 幂；基函数附录给出公式及数值检查入口。

## 7. M、N 比较的含义

原定义保持

$$
M_{k,n}=\int_0^1\frac{W^kL^n}{\beta-z/(1+z)}d\beta,\quad
N_{k,n}=\int_0^1\frac{W^kRL^n}{\beta-z/(1+z)}d\beta,\quad
L=\ell+R.
$$

于是 ell*M_(k,n)=M_(k,n+1)−N_(k,n)。允许显式 ell 的 M 多项式表达可通过这个字典改写成 M/N。它不证明所有 N 是纯 M 多项式，也不说明任意新母核都能用 N 表示。

本批比较三个精确定义的空间：固定权重/宇称的全部 M 多项式加 zeta；同一空间再允许 ell（以 ell² 关系约去高次）；上述扩展 tau 系数的多项式。用完整 Hlog/period 坐标做有理线性代数和反代残差检查，不使用数值 PSLQ 猜积分。失败时保存候选列表、秩和不相容见证，结论只针对所声明的坐标空间与常数约定，不宣称任意函数表示或全圈最小闭包。
