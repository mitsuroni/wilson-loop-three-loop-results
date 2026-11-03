# 3e-3g-04的正ε阶：ε¹仍在M/N，ε²连接到双端单核

日期：2026-09-29。本轮将当前归一化的完整 $I_{3e-3g-04}$ 计算到ε²，并在已声明的函数空间中检验M/N闭包。ε¹的旧参考算法已隔离重放；ε²由同一已正规化的原积分链实际扩阶得到。旧r001与其至有限项的结果保持原字节。

**结论：ε¹可以用M/N表示；ε²不在通常的权重六M/N多项式空间中，但加入已独立定义的 $\tau^+_{D,4}$ 后精确闭合，其系数为1/8。** 这里没有新运动学字母，结果仍只用0、−1的Hlog；也不需要为本阶另定义一个3g专属函数或加入E。任意有理z系数、自变量变换下的更广约化及全阶生成结构仍未证明。

## 1. 同一个完整I与展开约定

令z>0、$\ell=\log z$。原5维参数积分J、分子P和分母S沿用[当前定义](all-epsilon-status.md)的(TGV9)/(TGV10)。完整图为

$$
I_{3e-3g-04}=
\frac{e^{-2\gamma_E\epsilon}\Gamma(2-2\epsilon)}
{1+z^{2\epsilon}}J_{04}(\epsilon,z).
\tag{H1}
$$

先在积分收敛域建立IBP与换元，再对ε=0亚纯延拓。不把Γ(2−2ε)换成纯交换图的Γ(1−ε)²，Euler补偿恰好一次。

本页使用原有的一维矩：

$$
\begin{gathered}
r=\frac{z}{1+z},\quad
W(\beta)=\log[\beta(1-\beta)],\quad
L(\beta,z)=\ell+\log\frac{1-\beta}{\beta},\\
M_{k,j}=\int_0^1\frac{W^kL^j}{\beta-r}\,d\beta,\qquad
N_{k,j}=M_{k,j+1}-\ell M_{k,j},\quad k\ge0,\ j\ge1.
\end{gathered}
\tag{H2}
$$

内部点β=r可去。$[\epsilon^q]I$表示普通Laurent系数，不额外乘q!。所有旧的负幂和有限项回代保持：

$$
\begin{aligned}
I_{3e-3g-04}={}&
\frac{\zeta_2}{8\epsilon^2}
+\frac{-M_{1,1}+\zeta_3}{8\epsilon}
+\frac{M_{0,1}^2}{48}-\frac{M_{2,1}}8+\frac9{10}\zeta_2^2\\
&+\epsilon[\epsilon^1]I_{3e-3g-04}
+\epsilon^2[\epsilon^2]I_{3e-3g-04}
+O(\epsilon^3).
\end{aligned}
\tag{H3}
$$

## 2. ε¹：需要N型组合，但没有超出M/N

完整一阶为

$$
\boxed{\begin{aligned}
[\epsilon^1]I_{3e-3g-04}={}&
-\frac{5M_{3,1}}{96}-\frac{M_{1,3}}{96}
-\frac{M_{0,1}M_{1,1}}{16}
-\frac{3\zeta_2M_{1,1}}{16}\\
&+\frac{3\zeta_3M_{0,1}}4
+\frac58\zeta_5-\frac{37}{24}\zeta_2\zeta_3
-\frac{N_{1,2}}{16}.
\end{aligned}}
\tag{H4}
$$

若只使用M和ζ、禁止显式log z，在本阶候选空间内会出现秩7到8的缺口。加入 $\ell M_{1,2}=M_{1,3}-N_{1,2}$ 后缺口消失，秩和增广秩都为8。因而这不是发现新的独立函数族。

旧本地参考资料曾在只取奇数第二指标时得到“不闭合”，随后补入ℓM12纠正。本轮没有把这个早期失败继续当作新基底证据；三个积分组成部分已经在当前输入连接下重新计算，最后才与完整M/N空间比较。

## 3. ε²：出现一个已有的双端方向

完整二阶可以写为

$$
\boxed{\begin{aligned}
[\epsilon^2]I_{3e-3g-04}={}&
\frac{223M_{0,1}^3}{11520}
-\frac{13\zeta_2M_{0,1}^2}{1920}
+\frac{169\zeta_2^2M_{0,1}}{320}\\
&-\frac{5M_{0,1}M_{2,1}}{64}
+\frac{M_{1,1}^2}{32}
-\frac{39M_{4,1}}{1024}
+\frac{25M_{2,3}}{1536}
-\frac{N_{2,2}}{32}\\
&-\frac{21\zeta_2M_{2,1}}{64}
+\frac53\zeta_3M_{1,1}
-\frac{19}{24}\zeta_3^2
+\frac{7381}{1680}\zeta_2^3
+\frac18\tau^+_{D,4}(z).
\end{aligned}}
\tag{H5}
$$

这个τ不是从本图余项临时命名的。它来自此前纯交换[1]族的双端单核，明确为

$$
\begin{gathered}
D_z(t,s)=(1-t)(1-s)+zts,\\
X_+=\log\frac{D_z(t,s)\,ts}{\sqrt z},\qquad
X_-=\log\frac{D_z(t,s)(1-t)(1-s)}{\sqrt z},\\
\boxed{\displaystyle
\tau^+_{D,4}(z)=
-\frac{1+z}{48}\int_0^1\!\!\int_0^1
\frac{X_+^4+X_-^4}{D_z(t,s)}\,dt\,ds.}
\end{gathered}
\tag{H6}
$$

等价地，若

$$
H_D(\delta,z)=-(1+z)z^{-\delta/2}
\int_0^1\!\!\int_0^1D_z(t,s)^{\delta-1}t^\delta s^\delta\,dt\,ds,
$$

则τ是其δ⁴系数的z↔1/z偶平均。这里用δ区分源展开参数与本图的ε；两者的相关系数都具有权重六。

表达(H5)还使用了已核验的低指标恒等式

$$
M_{0,5}=\frac43M_{0,1}^3-8\zeta_2M_{0,1}^2,\qquad
N_{0,4}=-\frac4{15}M_{0,1}^3-\frac{32}5\zeta_2M_{0,1}^2
+\frac{96}5\zeta_2^2M_{0,1}.
$$

因此没有将显式ℓM22、ℓM04留作未说明的新方向。完整Hlog表达分别有23和48个非零坐标，见[I¹](../../../workflow/checks/tgv04-high-order/runs/positive-epsilon/I1.txt)、[I²](../../../workflow/checks/tgv04-high-order/runs/positive-epsilon/I2.txt)。

## 4. “不在M/N内”的精确范围

候选使用有理系数、ζ₂/ζ₃/ζ₅及其乘积，完整M族和显式ℓ；N通过(H2)包含，高次ℓ用 $\ell^2=-2M_{0,1}-6\zeta_2$ 消去。枚举相应总权重、总crossing偶的所有乘积，不漏掉奇函数乘奇因子。这里的权重分级及周期坐标沿用已核验的Hlog表示。

| 系数 | 候选空间 | 秩 | 加入该系数后的秩 |
|---|---|---:|---:|
| ε¹ | M及ζ，不允许显式ℓ | 7 | 8 |
| ε¹ | M/N及ζ | 8 | 8 |
| ε² | M/N及ζ | 12 | 13 |
| ε² | 再加log(1+z)的所有允许乘积 | 29 | 30 |
| ε² | M/N再加旧E | 13 | 14 |
| ε² | M/N再加旧τD4 | 13 | 13 |
| ε² | M/N再加τD4与E | 14 | 14 |

最后一行的转换仍不需要E。失败项保留了精确左零空间见证，成功项在完整Hlog/周期坐标中逐项回代为零；并非仅通过数值采样判断。

例如，对权重六M/N候选，令方括号表示完整shuffle展开后的对应坐标系数，线性泛函

$$
\begin{aligned}
\mathcal L(F)={}&
\frac12[\zeta_5G(-1;z)]F
+\frac5{17}[\zeta_2^2G(-1,-1;z)]F\\
&+\frac{10}{17}[\zeta_2^2G(-1,0;z)]F
+[\zeta_5G(0;z)]F
\end{aligned}
$$

在全部该空间候选上为0，但在 $[\epsilon^2]I_{3e-3g-04}$ 上为−1。这给出所声明坐标空间中的不相容证书。G按首字母Cauchy核定义，$G(0;z)=\log z$；本轮所有输出字母仍为0、−1。

没有穷尽一般有理/代数z系数、函数分母、不同自变量上的M/N或其他源变换。因此不能将表中结论扩大为任意解析表示都不可能，也不能据此认定一个全阶最小主积分数。

## 5. 新系数怎样从原积分得到

本轮未展开原发散积分后直接逐项积分。先复用并核对正规化链的原核、原域和IBP恒等式，再对可积的系数核扩阶。

若原分子写成

$$
P=\frac{-z\,\partial_{x_1}S-(1+z)\partial_{x_2}S+\partial_{x_3}S}{4}
+P_{\rm rec},\qquad
P_{\rm rec}=\frac{1+z}{2}\alpha_3
[\alpha_1(1-x_1+zx_3)-z\alpha_2x_3],
$$

则六个端点面与recoil都必须保留。定义边界积分
$F_{i,b}=\int_{\text{其余Wilson变量}}\int_{\Delta_2}d^3\alpha\,
\delta(1-\sum\alpha)(\alpha_1\alpha_2\alpha_3)^{-\epsilon}
(S|_{x_i=b})^{-1+2\epsilon}$。

由一次分部积分以及各面的精确镜像关系，

$$
\begin{aligned}
J_{04}={}&\frac1{4(1-2\epsilon)}
\left[(1+z^{2\epsilon})\mathcal C_{\rm cusp}
-\frac{1+z}{\epsilon}\mathcal D
-z\mathcal H(\epsilon,z)
-z^{2\epsilon-1}\mathcal H(\epsilon,1/z)\right]
+\mathcal O,\\
\mathcal C_{\rm cusp}
={}&\frac{\Gamma(1+\epsilon)\Gamma(1-\epsilon)
-\Gamma(1-\epsilon)^2/\Gamma(1-2\epsilon)}{4\epsilon^4}.
\end{aligned}
\tag{H7}
$$

这里 $\mathcal D=\epsilon(F_{2,0}-F_{2,1})$，$\mathcal H=F_{1,0}$，$\mathcal O$是PRrec的积分；下标i指x₁、x₂、x₃，避免与旧边号21/31/41混淆。H7只组织已有几何分块，不将这些中间符号称为新基底。

| 完整系数 | 对边差商需要 | 困难边界需要 | 正规化recoil需要 |
|---|---|---|---|
| ε¹ | D₂ | H₁ | O₁ |
| ε² | D₃ | H₂ | O₂ |

下标为各块普通Taylor/Laurent系数。每一项在独立的新Maple会话计算，发散检查与遇发散中止均开启。原生dimregPartial再次核对recoil正规化；没有丢弃发散项或在最后用regHlog抹掉极点。

具体地，正规化recoil的分子只有 $A/\epsilon+B$，其余幂产生一个明确对数Lb，因此

$$
O_1=\int(BL_b+AL_b^2/2)\,d\mu_0,\qquad
O_2=\int(BL_b^2/2+AL_b^3/6)\,d\mu_0.
$$

A、B和dμ₀均由同一个原积分图表确定，脚本保存完整表达式。两个困难边界同样由原零阶核乘Lh和Lh²/2计算。最后按(H1)卷积Gamma/crossing因子：求I¹需前因子至ε³，求I²需至ε⁴；没有把裸J系数直接当成完整I系数。

## 6. 对全解析结构的含义

现在已有一个明确的有限阶反例：在上述通常M/N多项式语言中，3e-3g-04的全阶纯M/N闭包不成立，首个已计算的缺口在ε²。

更有用的是这个缺口已经与旧双端单核建立联系。已知
$[\epsilon^0]I_{6e-36}$的M/N外部分为 $-\tau^+_{D,4}/32$，故

$$
[\epsilon^2]I_{3e-3g-04}
+4[\epsilon^0]I_{6e-36}
\quad\text{属于所声明的M/N多项式空间。}
\tag{H8}
$$

这是两个已求系数的精确关系，不是不同圈数完整I的全ε恒等式。它建议下一步从积分层面寻找3g图的高阶权重与双端源HD的连接；尚未证明04任意阶都由HD或现有P/Mellin源覆盖。

## 7. 核验与交付范围

- 6个新系数积分全部完成，最高阶recoil使用约164秒Maple CPU；所有中间输出与日志保留。
- 原核/分子/图表、正规化系数和cusp展开共17项代数控制通过；完整J和I在两阶均通过精确crossing检查。
- 两条紧凑系数及两个低指标恒等式精确回代通过；非成员证书独立检查。
- 45项实积分控制通过：M矩高精度求积、双端τ矩网格收敛、两个对边系数，以及原图表四/五维QMC。I²矩表达与Hlog数值在z=1/2、1、2处差小于7×10⁻¹⁰；QMC仅提供较低精度的独立控制，误差不是严格区间包围。

详细数据和可复算脚本见[本轮检查](../../../workflow/checks/tgv04-high-order/README.md)。ε¹使用了旧本地参考算法，本文不把它说成首次发现；ε²是本轮实际扩阶。四个旧r001、prepared和现有有限项不变，本轮新增正ε系数补充；formal=not_run，未提交Git。
