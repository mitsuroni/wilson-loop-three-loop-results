# 三圈[3]族至ε²：解析系数与M/N之外的两个权重8方向

日期：2026-09-28。对象为6e-20、6e-21、6e-22。本轮已经得到三图完整物理I的ε⁰、ε¹、ε²解析系数，所有负幂为零，余项为O(ε³)。21与22的完整I相同，因此共有六个独立非零表达式。

零阶、一阶与旧结果精确一致；二阶是本轮新求值。结果均已展开为字母0、−1的Hlog及明确周期，不以未求值积分冒充系数。完整表达见[Hlog附表](three-kernel-eps2-hlog.md)、[机器数据](three-kernel-eps2.json)，方法与复验见[检查目录](../../../../workflow/checks/three-kernel-eps2/README.md)。本批不改写旧r001/catalog，补充结果以本页为入口；formal=not_run。

## 1. 先给出结论：二阶的情况与前两阶不同

固定权重、crossing偶、有理多项式系数及既有Hlog/形式周期坐标后：

| 物理系数 | 权重 | 两个独立结果相对所声明M/N空间增加的维数 | 已知表示 |
|---|---:|---:|---|
| ε⁰ | 6 | 1 | 已有E与M/N |
| ε¹ | 7 | 1 | 已有双核两源矩与M/N；两裸核系数之和仍在M/N内 |
| ε² | 8 | 2 | 本轮已求值的两个自然源矩w₂、Λ₂，加M/N及旧函数乘积 |

权重8的精确比较进一步得到：

| 候选空间 | 原秩 | 加入两个三核二阶系数后的秩 |
|---|---:|---:|
| M/N与所需周期 | 36 | 38 |
| 再加旧权重≤6统一基底的权重8乘积 | 40 | 42 |
| 再加现有单核五阶奇系数乘log z | 41 | 43 |
| 再加现有双核三阶两源奇矩乘log z | 45 | 47 |

这里既可用裸A/B系数，也可用完整物理I，两者的联合秩均如表。首层已允许所需的权重8多重ζ常数，不能把缺少一个周期误报成缺少运动学函数。每个失败包含精确不相容见证，而非数值拟合残差。

**本轮可以确认：在已经明确列出的空间之外，需要两个独立的权重8方向。** 尚未枚举单核全部六阶系数、双核四阶偶源矩或任意端点字，因此不能升级为“必须引入两个三核专属函数族”。两个方向目前有来自原积分的自然代表，见下文。

## 2. 两个自然源矩怎样从原积分产生

沿用[全ε研究](three-kernel-generators.md)的矩形原函数F、右上矩形U、边界原函数f及单核H。记同序裸三核为A，另两个配对的偶平均为B。

此前已证明

$$
A=\int k_\epsilon FU,\qquad
\mathcal W=\int k_\epsilon f(t)f(s),\qquad
A+2B=\mathcal W-\frac{H^3}{6}.
$$

本轮再使用三次幂的混合导数：

$$
\partial_t\partial_sF^3
=3F^2k_\epsilon+6F(\partial_tF)(\partial_sF).
$$

定义

$$
\mathcal L(\epsilon,z)=\int_0^1\!\int_0^1
F(\partial_tF)(\partial_sF)\,dt\,ds.
$$

保留四角边界，积分左边等于H的三次方；右上矩形给出crossing伙伴。因此

$$
\boxed{\begin{aligned}
A&=\frac{H^3}{2}-\mathcal W-\mathcal L(\epsilon,z)-\mathcal L(\epsilon,1/z),\\
B&=\mathcal W+\frac{\mathcal L(\epsilon,z)+\mathcal L(\epsilon,1/z)}2-\frac{H^3}{3}.
\end{aligned}}
\tag{E1}
$$

这两个源矩在系数求值和选基以前就已由矩形边界导出。令

$$
w_q=[\epsilon^q]\mathcal W,
\qquad
\Lambda_q=[\epsilon^q]\{\mathcal L(\epsilon,z)+\mathcal L(\epsilon,1/z)\}.
\tag{E2}
$$

Λ使用两项之和，不是除以2的偶投影。它们各阶的实际求值由下面的端点字算法完成，本轮没有只把未知二阶系数换一个名字。

## 3. 为什么新算法能完成上轮超时的计算

取中心变量p、q及ρ为

$$
p=(1+z)t-1,\quad q=(1+z)s-1,
\quad \rho=\frac{z}{(1+z)^2},
\quad -1<p,q<z.
$$

矩形原函数已求值为

$$
F=-\frac{\rho^{\epsilon/2}}\epsilon\,\mathcal J,
\qquad
\mathcal J=J_\epsilon(-pq/z)-J_\epsilon(p/z)-J_\epsilon(q/z)+J_\epsilon(-1/z).
$$

由边界原函数及矩形的两个导数直接得到

$$
\boxed{\begin{aligned}
\mathcal W&=-\frac{\rho^{3\epsilon/2}}{\epsilon^3}
\int_{-1}^z\!\int_{-1}^z\frac{\mathcal J}{pq}
\big[(1+p)^\epsilon-(1-p/z)^\epsilon\big]
\big[(1+q)^\epsilon-(1-q/z)^\epsilon\big],dp,dq,\\
\mathcal L&=-\frac{\rho^{3\epsilon/2}}{\epsilon^3}
\int_{-1}^z\!\int_{-1}^z\frac{\mathcal J}{pq}
\big[(1+pq/z)^\epsilon-(1-p/z)^\epsilon\big]
\big[(1+pq/z)^\epsilon-(1-q/z)^\epsilon\big],dp,dq.
\end{aligned}}
\tag{E3}
$$

每个括号及J都至少一阶，所以分子从ε³起始；二阶结果要求分子完整计算到ε⁵，不能只截取到二阶。

把p、q的区间在零处分开，并缩放到单位方形。两个同号扇区和一个交叉扇区的重数为1、1、−2；交叉扇区的负号来自带方向的dp/p。不得把它们都当正号。

在每个扇区中，所需展开均为明确的迭代积分字：

$$
[\epsilon^m]J_\epsilon(x)
=H\big(x;0,\underbrace{1,\ldots,1}_{m}\big),
\qquad
\frac{\log^m(1-x/a)}{m!}
=H\big(x;\underbrace{a,\ldots,a}_{m}\big).
$$

第一积分使用shuffle展开乘积，然后按dx/x在字首补零；第二变量作fibration后也使用同样操作：

$$
\int_0^1\frac{H(x;w)}x\,dx=H(1;0,w).
\tag{E4}
$$

这条规则只在相应端点极限存在时使用。实现逐条检查：系数不依赖当前积分变量、表达式对该变量的Hlog线性、无未处理常数余项、无不收敛的纯零字。通过这些检查后才取端点；没有关掉发散检查或舍弃字母。

本轮Maple CPU约91秒完成两个源的全部三至五阶分子。通用积分器的四个已完成低阶组合与新端点算法精确一致；其较慢的后续计算在新路线完成后主动停止，保留原日志。最初两个脚本命名冲突也保留为失败会话。

## 4. 完整物理I的三个系数

写H的三次方为

$$
H^3=\sum_{q\ge0}s_q\epsilon^q.
$$

由一圈M生成函数，

$$
\begin{aligned}
s_0&=M_{0,1}^3,\\
s_1&=\frac32M_{0,1}^2M_{1,1},\\
s_2&=\frac38M_{0,1}^2M_{2,1}
-\frac18M_{0,1}^4+\frac34M_{0,1}M_{1,1}^2.
\end{aligned}
\tag{E5}
$$

物理前因子必须保留到同阶：

$$
\frac{Z^3}{\cosh(3\epsilon\ell/2)}
=1+p_2\epsilon^2+O(\epsilon^3),\qquad
p_2=\frac32\zeta_2-\frac98\ell^2
=\frac{33}{4}\zeta_2+\frac94M_{0,1}.
\tag{E6}
$$

于是用已经求值的源矩，可以同时写出两图的前三阶：

$$
\boxed{\begin{aligned}
I_{20}^{(q)}&=\frac{w_q+\Lambda_q-s_q/2}{16}
+\delta_{q,2}p_2 I_{20}^{(0)},\\
I_{21}^{(q)}=I_{22}^{(q)}&=\frac{s_q/3-w_q-\Lambda_q/2}{32}
+\delta_{q,2}p_2 I_{21}^{(0)},\qquad q=0,1,2.
\end{aligned}}
\tag{E7}
$$

零阶的源矩还可写成

$$
w_0=-2E+\frac{M_{0,1}^3}{3}
+\frac12M_{0,1}M_{2,1}-\frac12M_{1,1}^2,
\qquad \Lambda_0=\frac5{18}M_{0,1}^3.
$$

一阶逐项匹配[旧权重7解表](../weight-seven/solutions.md)。二阶的w₂、Λ₂已在[Hlog附表](three-kernel-eps2-hlog.md)完整展开；完整I20与I21的二阶分别有152、153个非零Hlog/周期坐标。这里保留源矩记号是为了显示共同结构，具体函数值并未留待计算。

## 5. 初步基底判断：这次不能只保留一个新方向

令低函数空间包含M/N、所需周期以及旧τ/E的所有已列权重8乘积。物理前因子乘零阶结果属于此空间，因此在该商空间中

$$
\boxed{\begin{aligned}
[I_{20}^{(2)}]&=\frac{[w_2]+[\Lambda_2]}{16},\\
[I_{21}^{(2)}]=[I_{22}^{(2)}]&=-\frac{[w_2]}{32}-\frac{[\Lambda_2]}{64}.
\end{aligned}}
\tag{E8}
$$

两个源方向独立，由精确秩40到42及不相容见证确认；再加入目前可复用的部分高阶低核源，仍留下二维。特别是

$$
[I_{20}^{(2)}+2I_{21}^{(2)}]=\frac{[\Lambda_2]}{32}\ne0.
\tag{E9}
$$

所以此前裸A加B的零阶、一阶仍在M/N中的现象，在二阶不再成立；即使加入旧τ/E乘积也不能恢复这个结论。

后续最值得判断的是：完整单核六阶与双核四阶偶源矩能否覆盖w₂、Λ₂。若能，就继续用已有源族组织；若不能，才扩大具体候选范围。当前结果支持“两项已求值的权重8代表”，不支持未经进一步分析就宣布“两种三核专属的全阶不可约族”。

## 6. 验证与交付范围

24个独立参数系数求积比较通过，另有遗漏区间二阶修正的4个负控制。检查使用原三维/四维端点积分与独立Hlog微分方程级数求值，不依赖新的字首补零算法。零阶三图及一阶两个独立对象与历史解析结果全部精确一致；缺失物理前因子也有非零反控制。

新鲜源生成增加原有限项表的21个位置，累计56/329；ε¹、ε²专项不重复增加这个固定计数。旧结果修订、catalog、历史证据保持。当前仅验证正实z与ε零点邻域；物理复分支、全阶最小性及Lean整链未完成。未提交Git。
