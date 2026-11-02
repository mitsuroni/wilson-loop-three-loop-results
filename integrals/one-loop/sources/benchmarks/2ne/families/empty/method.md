# 空集族：从原参数 J 到完整 I

日期：2026-09-16。范围为 prepared 0.3 的 13 个纯交换空集族输入。这里给出纸面解析推导；辅助符号和数值检查见各结果记录，尚未进行本批结果的 Lean/MCP 整链验证。

## 1. 输入与参数约定

读取当前 prepared 清单，按 pure 与空的 direct-opposite channel partition 筛选，再逐条检查传播子、变量、界限和原迭代顺序。族标签本身不作为答案规则。本文用 z=z1>0，J 先取实 epsilon>0；I 先取 0<epsilon<1。固定正 z 后，在 epsilon=0 附近采用解析式的亚纯延拓。物理复运动学的 i0 延拓不在本次计算范围内。

每个输入恰有 2n 个不同端点；每个端点只在一个传播子中出现。程序逐多项式重建

$$
f_j=q_j(1-u_j)v_j,\qquad q_j\in\{1,z\},\qquad
\text{numerator}=\prod_jq_j=z^r.
$$

内部正底数允许拆幂，得到

$$
\text{integrand}
=z^{r\varepsilon}\prod_j(1-u_j)^{\varepsilon-1}v_j^{\varepsilon-1}.
$$

这一步没有使用源文件中仅供追溯的超立方 Jacobian。

## 2. 为什么允许按有序链组织积分

定义 A(x)=x^(epsilon-1)，B(x)=(1-x)^(epsilon-1)。输入的边界图逐项确认为互不相交的链：每个下界是 0，上界是 1 或本链下一个端点。沿每链从小到大读取 A/B，得到字 w。

这里显式补上换序和分组的数学依据：对 epsilon>0，被积函数非负，而且放宽所有排序约束到完整立方体后，

$$
\int_{[0,1]^{2n}}z^{r\varepsilon}\prod_{i=1}^{2n}f_i(s_i)\,d^{2n}s
=z^{r\varepsilon}\varepsilon^{-2n}<\infty.
$$

故原有序域上绝对可积。将约束写为指标函数后，Tonelli/Fubini 连接原指定的迭代积分与这些链的乘积积分。每个原上下界均在导出的链表中保留；链间没有额外约束，重组 Jacobian 为 1。对 Re epsilon>0 的复 epsilon，在任意紧子域也有同类绝对值控制，积分对 epsilon 全纯。

因此

$$
J(z,\varepsilon)=z^{r\varepsilon}K(\varepsilon),\qquad
K=\prod_{\text{chains}}F_w,
$$

$$
F_{c_1\cdots c_m}=
\int_{0<x_1<\cdots<x_m<1}\prod_{j=1}^m f_{c_j}(x_j)\,d^mx,\qquad F_\varnothing=1.
$$

这是本次明确给出的分析桥，不是 prepared 0.3 自带的 Fubini 收据，也没有宣称已形式验证。独立的整数 epsilon 检查另按原始顺序从内向外积分，未使用此分组算法。

## 3. 本批所需的全部字积分

记

$$
G=e^{-\gamma_E\varepsilon}\Gamma(1-\varepsilon),\quad
R_2=\frac{\Gamma(1+\varepsilon)^2}{\Gamma(1+2\varepsilon)},\quad
R_{12}=\frac{\Gamma(1+\varepsilon)\Gamma(1+2\varepsilon)}
{\Gamma(1+3\varepsilon)}.
$$

| 字 | 完整 epsilon 结果 |
|---|---|
| A、B | $1/\varepsilon$ |
| AA、BB | $1/(2\varepsilon^2)$ |
| AAA、BBB | $1/(6\varepsilon^3)$ |
| AB | $R_2/\varepsilon^2$ |
| BA | $(1-R_2)/\varepsilon^2$ |
| AAB、ABB | $R_{12}/(2\varepsilon^3)$ |
| ABA、BAB | $(R_2-R_{12})/\varepsilon^3$ |
| BAA、BBA | $(1-2R_2+R_{12})/(2\varepsilon^3)$ |

推导如下。

同字积分用相同权重的置换对称性或逐次原函数积分，得到
$F_{A^m}=F_{B^m}=1/(m!\varepsilon^m)$。

积分掉最小的 A 后，

$$
F_{AB}=\frac1\varepsilon\int_0^1x^\varepsilon(1-x)^{\varepsilon-1}dx
=\frac{B(1+\varepsilon,\varepsilon)}{\varepsilon}
=\frac{R_2}{\varepsilon^2}.
$$

积分区域分割给 $F_AF_B=F_{AB}+F_{BA}$，从而得 BA。类似地，

$$
F_{AAB}=\frac{B(1+2\varepsilon,\varepsilon)}{2\varepsilon^2}
=\frac{R_{12}}{2\varepsilon^3}.
$$

对 BAB，以中间 A 的位置 x 固定，两侧 B 可以直接积分：

$$
F_{BAB}=\frac1{\varepsilon^2}\int_0^1
x^{\varepsilon-1}[1-(1-x)^\varepsilon](1-x)^\varepsilon dx
=\frac{R_2-R_{12}}{\varepsilon^3}.
$$

对 BBA，先完成两个有序 B：

$$
F_{BBA}=\frac1{2\varepsilon^2}\int_0^1
x^{\varepsilon-1}[1-(1-x)^\varepsilon]^2dx
=\frac{1-2R_2+R_{12}}{2\varepsilon^3}.
$$

变换 x_j→1-x_(m+1-j) 保留有序域且 Jacobian 绝对值为 1，将字变成“逆序后交换 A/B”，于是 ABB=AAB、ABA=BAB、BAA=BBA。所有拆积分均在 epsilon>0 时使用，之后再延拓。超过三字的混合链由脚本拒绝；这组公式不宣称任意圈阶已覆盖。

## 4. 唯一 C、crossing 与 I

设 a1 是完整固定颜色基底的第一坐标，sym 是输入 orbit_multiplicity。正实尺度区内，令

$$
\mathcal N=C_A^{n-1}C_F(\alpha_s/\pi)^n,\qquad
S_{14}=(-2\pi v_{14})^{n\varepsilon},\qquad S_{12}=z^{n\varepsilon}S_{14}.
$$

当前 prepared C 的已补偿表达式整理为

$$
C=\frac{a_1}{(-2)^n}\mathcal N S_{14}G^n.
$$

脚本逐条核对 C.combined 的所有因子及 Gamma/Euler 记录，不能忽略额外符号或重复补偿。沿用当前对称化定义

$$
\frac{\mathrm{sym}}2C[J(z)+z^{n\varepsilon}J(1/z)]
=\mathcal N S_{14}(1+z^{n\varepsilon})I.
$$

取非零共同尺度/耦合，或等价地提取对应微扰系数。由于 K 与 z 无关，

$$
\boxed{
I=\kappa G^nK
\frac{z^{r\varepsilon}+z^{(n-r)\varepsilon}}{1+z^{n\varepsilon}}
=\kappa G^nK
\frac{\cosh[(r-n/2)\varepsilon\ell]}{\cosh(n\varepsilon\ell/2)},
\quad
\kappa=\frac{\mathrm{sym}\,a_1}{2(-2)^n},\quad \ell=\ln z.
}
$$

这里没有假设 J 的 self-mirror 关系；J 和旋转项分别计算。所有实 z>0、0<epsilon<1 下分母为正且 Gamma 有限。I 的 epsilon 亚纯延拓由此确定；对固定 z，局部展开避开其他 Gamma/双曲函数奇点。复 z 的尺度、分支与物理源连接仍需独立工作。

## 5. 展开到有限项

Gamma 的局部对数级数给

$$
\log G=\sum_{k\ge2}\frac{\zeta(k)}k\varepsilon^k,\quad
\log R_2=\sum_{k\ge2}\frac{(-1)^k(2-2^k)\zeta(k)}k\varepsilon^k,
$$

$$
\log R_{12}=\sum_{k\ge2}
\frac{(-1)^k(1+2^k-3^k)\zeta(k)}k\varepsilon^k.
$$

可取共同局部范围 |epsilon|<1/3，并进一步避开给定 z 的双曲分母零点。令指数级数为 exp(sum l_k epsilon^k)=sum a_m epsilon^m，则

$$
a_0=1,\qquad a_m=\frac1m\sum_{k=1}^m k\,l_k a_{m-k}.
$$

脚本在精确数域中卷积这些级数和双曲比值。对 n 圈保守保留正则部分到 2n 阶，乘 epsilon^(-2n) 后得到全部极点及有限项，余项 O(epsilon)。还多算两阶用于余项检查，但标准交付只到零阶。系数中的 ell 就是 log(z)，没有新物理参数；机器表达式中已代回 log(z)。

零系数显式记录。特别地，R2、R12 组合会降低实际极点阶，不能从显式 epsilon^(-2n) 猜最终极点。

## 6. 计算检查的层次

- 输入：逐传播子多项式残差、分子/尺度、所有端点和边界链、完整 C 因子检查。
- 原积分：epsilon=1,2,3，z=2/3 时直接把原始 expanded J 化为多项式，按原顺序精确积分。只检查 J，不在 Gamma 极点上检查 I。
- 字积分：70 位工作精度、三个非整数 epsilon 点的直接一维求积，对照完整 Gamma 表达。
- 原域数值：对原始 bounds 的指标函数作随机化 Sobol 采样，使用归一化 A/B 端点密度；统计误差如实记录，不称为严格区间包围。
- 系数：70 位精度的复 epsilon Cauchy/Fourier 采样，从完整 Gamma 式提取至有限项的系数，与符号递推独立对照。
- 余项：固定三个正 z 点，epsilon 依次减半，检查 (I−Laurent至有限项)/epsilon 向下一系数收敛。余项阶数的分析依据是上述收敛局部级数。

以上是计算与纸面数学证据。旧 Wilson 收据保持原来覆盖范围，本批结果的 formal 字段统一为 not_run。没有为通过测试修改 Kernel、Lean 或默认根输入。

## 7. 保存与复现

每个结果的 rNNN.json 保存输入哈希、完整公式、系数、域、推导和运行记录引用。表达式同时导出 WolframLanguage/InputForm 和 LaTeX，二者来自同一个 SymPy 精确对象；文件仅存数据，不将任意字符串作为代码执行。修订只新增，不覆盖旧文件。

每积分只维护 workspaces/<orbit-id>/derivation.md。公共代码和必要结果在 workflow/checks/empty-family/；临时缓存只放 .wilson。历史推导保留原路径和原文，本轮当前区块与历史内容显式分隔。

执行入口见 [计算说明](../../../../workflow/checks/empty-family/README.md)。
