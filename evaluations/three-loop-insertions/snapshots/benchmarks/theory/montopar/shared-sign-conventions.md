# 一圈至三圈的共同符号约定：从作用量到C、J、I

2026-10-06。本文是当前物理符号的共同推导入口。它解决的问题是：相同QCD和Wilson算符下，每张图的整体号如何在原始表达式中确定，而不由图名、答案或RG残差选择。旧冻结源仍保留其数学定义；本页给出它们与统一物理源的明确连接。不是对163个积分重新求值，也不是严格类光分布去调节的形式证明。

本轮从实际源核对一圈2图、两圈12图、三圈149图，共163图。共同规则602项精确检查、坐标时间积分3点控制，逐源589项检查通过。Γ₁／腿Π₁的16个旧固定源需要在颜色收缩阶段转换；还新发现FGV六个非零旧源的星积分相位错误。详见[逐源账本与检查](../../../workflow/checks/shared-sign-conventions/README.md)。当前生产及不可变结果没有在本轮改写，后文直接定义物理C，不把历史转换倍率作为新图的定义。

## 1. 先固定一个作用量和一个Wilson算符

采用Minkowski度规(+---)、D=4−2ε、真空归一化的e^{iS}期望值、Hermitian基本表示生成元。固定

$$[t^a,t^b]=if^{abc}t^c,\qquad \operatorname{Tr}(t^at^b)=T_F\delta^{ab},\qquad T_F=\tfrac12.$$

物理起点为

$$\begin{aligned}
\mathcal L={}&-\frac14F_{\mu\nu}^aF^{a\mu\nu}
-\frac12(\partial_\mu A^{a\mu})^2
+\bar q,i\!\not\!Dq+(\partial_\mu\bar c^a)(D^\mu c)^a,\\
D_\mu q={}&(\partial_\mu-ig_0A_\mu)q,\qquad
F_{\mu\nu}^a=\partial_\mu A_\nu^a-\partial_\nu A_\mu^a
+g_0f^{abc}A_\mu^bA_\nu^c.
\end{aligned}$$

这些定义与[Grozin §4.1，式(4.1)–(4.8)](https://arxiv.org/pdf/hep-ph/0508242#page=30)一致。Feynman规范对应Grozin的a₀=1、DOT的ξ=0。所有动量公式使用

$$A(x)=\int\frac{d^Dp}{(2\pi)^D}e^{-ipx}A(p),\qquad \partial_\mu\longleftrightarrow-ip_\mu.$$

沿参数τ递增的轮廓，定义输运矩阵

$$\frac{dU(\tau,0)}{d\tau}=ig_0\dot X^\mu(\tau)A_\mu(X(\tau))U(\tau,0),\qquad U(0,0)=\mathbf1,$$
$$W_B=\frac1{N_c}\langle\operatorname{Tr}U(1,0)\rangle.$$

这个微分方程明确了P的含义，而不只写一个可能有不同约定的P字符。迭代一次后，两个插入是A(τ₂)A(τ₁)，τ₂>τ₁；一般有序域0<τ₁<⋯<τ_E<1中，Wilson词为

$$t^{a_E}\cdots t^{a_1}.$$

这里把槽位1…E定义为轮廓上的早到晚。若改用晚到早编号，也可以写t¹⋯tᴱ，但位置、动量、中心顶点和排序域必须同时重标。Wilson路径指数的文献背景见[Korchemsky–Marchesini，式(2.1)、(2.9)](https://arxiv.org/pdf/hep-ph/9210281)；晚点在左的具体选择由本页输运方程直接导出。右作用输运方程不是同一固定标签下的替代规则。

## 2. Feynman规则从这个作用量导出

自由传播子已经包括Feynman规则中的i：

$$D_{0,\mu\nu}^{ab}(p)=\frac{-i\delta^{ab}g_{\mu\nu}}{p^2+i0},\quad
S_0(p)=\frac{i\not p}{p^2+i0},\quad
G_{0,\rm gh}^{ab}(p)=\frac{i\delta^{ab}}{p^2+i0}.$$

Wilson发射为ig₀u^μtᵃ，夸克顶点为ig₀γ^μtᵃ。每个闭合夸克圈或鬼圈各计一个−1；胶子圈不使用这条Grassmann规则。其原始图与说明见[Grozin式(3.22)、§4.3](https://arxiv.org/pdf/hep-ph/0508242#page=16)。

三胶子顶点中的符号可以直接检查。三场项为

$$\mathcal L_3=-g_0f^{abc}(\partial_\mu A_\nu^a)A^{b\mu}A^{c\nu}.$$

乘i并将偏导变为−ip，未对称项系数是−g₀fᵃᵇᶜpᵃ_μ。对三个带标记场求泛函导数，六种指派合并成

$$\mathcal V_3^{abc;\mu\nu\rho}(p_1,p_2,p_3)=g_0f^{abc}K^{\mu\nu\rho}(p_1,p_2,p_3),$$
$$K^{\mu\nu\rho}=g^{\mu\nu}(p_1-p_2)^\rho
+g^{\nu\rho}(p_2-p_3)^\mu+g^{\rho\mu}(p_3-p_1)^\nu,$$

其中p₁+p₂+p₃=0，全流入。交换两条完整腿时，f与K各反号，乘积不变。

这解释了文献中表面不同的写法。Grozin将颜色写成if、Lorentz写成ig₀V_G，而式(4.10)的V_G=−K，所以(if)(ig₀)V_G=g₀fK。[Grozin式(4.10)、图18](https://arxiv.org/pdf/hep-ph/0508242#page=32)。DOT式(2.1)–(2.2)定义的截肢顶点对象是Γ_DOT=−ig₀fK；本项目的Feynman插入因此为

$$\boxed{\mathcal V=i\Gamma_{\rm DOT}.}$$

这个字典由同一树级作用量匹配固定，按相同场与动量定义用于圈修正；它不由某个求值结果选择。[DOT定义](https://arxiv.org/pdf/hep-ph/9605348v3#page=4)。

四场项是

$$\mathcal L_4=-\frac{g_0^2}{4}f^{eab}f^{ecd}A_\mu^aA_\nu^bA^{c\mu}A^{d\nu}.$$

其24种场指派给出

$$\begin{aligned}
\mathcal V_4=-ig_0^2\Big[&f^{abe}f^{cde}(g^{\mu\rho}g^{\nu\sigma}-g^{\mu\sigma}g^{\nu\rho})\\
&+f^{ace}f^{bde}(g^{\mu\nu}g^{\rho\sigma}-g^{\mu\sigma}g^{\nu\rho})\\
&+f^{ade}f^{bce}(g^{\mu\nu}g^{\rho\sigma}-g^{\mu\rho}g^{\nu\sigma})\Big].
\end{aligned}$$

没有额外4!。本轮独立由L₄的24项核对这三个张量，未预设顶点号；与[Weinzierl讲义第114页](https://particlephysics.uni-mainz.de/weinzierl/download/script_qft1.pdf#page=114)的完整四顶点相符。鬼顶点由(∂\bar c)D c导出；其动量取对应反鬼腿，颜色顺序与鬼流一并固定。

## 3. 动量与坐标的Wick因子为何相反

对于p²+i0分母，动量时间积分沿p⁰=+ip_E⁰旋转，每个独立动量测度给+i，每个整数次p²分母转换另给(−1)。Grozin式(2.1)–(2.2)给出此方向。[原始推导](https://arxiv.org/pdf/hep-ph/0508242#page=4)。

坐标分母是−y²+i0，时间轮廓却沿y⁰=−iy_E⁰旋转，每个独立内部坐标给−i。一个足以固定号的简单例子是

$$\int_{-\infty}^{\infty}\frac{dt}{R^2-t^2+i0}=-\frac{i\pi}{R},\qquad R>0.$$

正根位于上半平面，其残数为−1/(2R)，闭合上半围道即得到−iπ/R。不能把动量的+i复制到坐标星积分。三点有限iη数值积分也独立支持这一方向。

于是含m条腿的坐标标量星核，令β_j>0、B=Σβ_j、Q=−Σ_{i<j}α_iα_j(X_i−X_j)²，采用共同边界值，有

$$\int d^Dy\prod_{j=1}^m[-(X_j-y)^2+i0]^{-\beta_j}
=-i\pi^{D/2}\frac{\Gamma(B-D/2)}{\prod_j\Gamma(\beta_j)}
\int d^m\alpha\,\delta(1-\sum\alpha_j)
\prod_j\alpha_j^{\beta_j-1}(Q+i0)^{D/2-B}.$$

这里的−i与m无关。两中心给(−i)²；不是按m或“是否有夸克”决定符号。参数化、换序和去调节先在合适解析域定义，再作共同亚纯延拓；本页未证明整个严格类光图有共同绝对收敛条带。

裸胶子坐标传播子为

$$\int\frac{d^Dp}{(2\pi)^D}\frac{-i e^{-ipX}}{p^2+i0}
=-\frac{\Gamma(1-\epsilon)}{4\pi^{2-\epsilon}}(-X^2+i0)^{-1+\epsilon}.$$

对三腿星核，三条坐标传播子的(−1)³乘中心−i给+i；用两个独立动量计算，同样是(−i)³·i²·(−1)³=+i。对四腿两条路线均给−i。这是对同一对象的两种计算，而不是分别采用两个相位惯例。

动量微分也由指数固定：e^{-iℓX}时ℓ=i∂_X；中心入射p=−ℓ时p=−i∂_X。改成p变量后指数为e^{+ipX}，仍给p=−i∂_X。内部坐标积分、顶点腿置换、导数方向和链式法则的负号必须分别保存。

## 4. 为什么3e-3g的现行结果保持，而Γ₁旧源需要转换

颜色恒等式始终是

$$\frac{f^{abc}}{N_c}\operatorname{Tr}(t^at^bt^c)=+\frac i2C_AC_F.$$

一般SU(N)的证明为

$$f^{abc}t^at^bt^c=\frac12f^{abc}[t^a,t^b]t^c
=\frac i2C_A\sum_ct^ct^c=\frac i2C_AC_F\mathbf1.$$

这是群论恒等式，SU(2)矩阵控制只是附加检验。

在本页早到晚槽位下，P给反序词，所以同一中心顶点标记下的实际收缩为

$$\frac{f^{abc}}{N_c}\operatorname{Tr}(t^ct^bt^a)=-\frac i2C_AC_F.$$

选择晚到早的腿编号可继续显示+i/2，同时重标中心腿、位置与动量即可；完整贡献一样。Π的颜色δ不会改变此恒等式，内部夸克圈负号也不属于这个外部色迹。

树图的无量纲Q用独立位置定义。记T_K为K(∂)作用于Q，本轮四个树父图及44个mixed父核均得到T_K=2A₀²P。标量星核相位+i、K(−ℓ)=−iK(∂)，再由∂Q^{-1+2ε}=−(1−2ε)Q^{-2+2ε}∂Q，实际物理密度为

$$\left(-\frac i2C_AC_F\right)i^3g_0^4
\frac{i\Gamma(1-2\epsilon)}{64\pi^{4-2\epsilon}}
(-i)\,[-(1-2\epsilon)]T_K Q^{-2+2\epsilon}.$$

因此直接提取的正确树C为

$$\boxed{C_{\rm tree}=\frac{C_AC_F}{4}
\left(\frac{\alpha_s}{\pi}\right)^2\Gamma(2-2\epsilon)(\pi A_0)^{2\epsilon}.}$$

它与当前可信3e-3g一致。旧文稿把色迹写成同序+i/2，又把标量星FT写成−i，两个号差相互补偿，恰好保留上述C×J。旧颜色恒等式本身正确；旧动量源、P和星FT的整链写法没有统一。不能只把其中一个字段改号，再把可信最终I翻转。

mixed使用同一完整Wilson词、一中心FT及裸旁观传播子，每个裸线额外保留r/A₀=−1/2。其完整颜色反序与星FT变化同样相消，现行修复后的44图C×J全部保持。这一结论从共同作用量与实际P产生，不把“mixed不翻号”作为定义。

Γ₁及腿Π₁的当前固定源已经使用正确的Wick链，但仍保留早到晚色迹。因此从初始源把正确P词收缩后，物理C直接成为下文的负值；J及原域保持。这里没有再给Π或夸克圈加负号。Γ₁与树图差异来自旧树文稿／实现的成对号差，不能归结为“所有非纯胶子多一个负号”。

## 5. Π算符、闭圈符号与Dyson连接

把完整1PI插入定义成含i对象

$$\Pi_L^{ab,\mu\nu}=-i\delta^{ab}(k^2g^{\mu\nu}-k^\mu k^\nu)\mathcal P_L(k^2),\qquad
D=D_0+D_0\Pi D_0+D_0\Pi D_0\Pi D_0+\cdots.$$

这里\(\mathcal P_L\)已经包含内部夸克、鬼和胶子的相位／符号。如果R=δ(D^{-1})，则R=−Π并有δD=−D₀RD₀；不能把R与Π混名而同时保留这个外负号。

取横向投影P_T，D₀=−iP_T/k²，Π=−ik²P_T\mathcal P。直接相乘给D₀ΠD₀=+iP_T\mathcal P/k²=−D₀\mathcal P。横向完整分母展开是1/(1+\mathcal P₁+\mathcal P₂+⋯)，单Π₁／单Π₂各一次负\(\mathcal P\)，两个Π₁串联给\(\mathcal P_1^2\)。它按插入次数计算，不按内部圈数或图名计算。

独立夸克泡仍从未积原图的闭圈−1开始，得到

$$\Pi_q^{\mu\nu}=i\mathcal A_q(k^2g^{\mu\nu}-k^\mu k^\nu)(-k^2-i0)^{-\epsilon},$$
$$\mathcal A_q=-\frac{8g_0^2T_Fn_f}{(4\pi)^{2-\epsilon}}
\Gamma(\epsilon)B(2-\epsilon,2-\epsilon).$$

这一负号已经进入\(\mathcal B_1=[(10-6ε)C_A-8(1-ε)T_Fn_f]/(3-2ε)\)。再给完整Π₁图补一个闭费米圈负号会重复计数，并错误翻转没有夸克的部分。胶子／鬼内圈和Π₂的标量仍采用已核对全D文献输入；本页固定定义与连接，不冒充对其全部原始内图的新积分计算。[独立泡与外圈推导](../../self-energy/pi1-chain/source-to-I-review.md)、[Π₂来源连接](../../self-energy/pi2/source-to-I-review.md)、[Ruijl等式(2.1)、(2.4)](https://arxiv.org/pdf/1703.08532v2#page=4)。

## 6. Γ₁内部的负号到底放在哪里

DOT全D算符接入项目时，定义\(\mathcal V_1=i\Gamma_{\rm DOT}^{(1)}\)，与本页树级字典一致。恢复显式内圈采用

$$T_c=\frac1i\int\frac{d^Dk}{(2\pi)^D}
\left[\frac{N_c}{d_1d_2d_3}+\sum_r\frac{M_{c,r}}{d_ad_b}\right].$$

既有全D约化已逐母积分核对这个字典与DOT完整张量；圈测度、η和耦合不能再乘一次。本次重新检查其相位来源：

| 原始内部图 | 原规则到canonical字典 |
|---|---|
| 胶子三角 | 三条−i传播子给+i，内部f收缩给C_A/2；写成(1/i)N要求N_tri=−G/2 |
| 鬼三角 | 鬼流、两个取向、闭鬼−1与各顶点顺序同时固定；既有源重建给−(G_a+G_b)/2 |
| 剑鱼泡 | 四顶点−i与两条−i传播子相乘为+i；颜色／对称因子／Lorentz合并后M_r=−9B_r/4 |
| 夸克三角 | 三个ig顶点、三个i传播子和闭费米−1相乘为+1；两个内部颜色方向结合反向Dirac链给iT_Ff，故N_q=−Tr_D(六个slash) |

表中的负分子与外部Wilson色迹反序是不同操作。内部完整字典一旦固定，不能再为q或ghost添加额外闭圈号。[全D分子与DOT连接](../../self-energy/gamma1/parametric-theory.md#23-为什么这确实是所引用的全-d-顶点)、[胶子／鬼／剑鱼重建](../../../workflow/checks/gamma1-gluon-numerator/README.md)。

固定旧源的三角有三条独立动量、六个二次分母和三次分子，i³(−1)⁶(−i)³=+1；泡图为i³(−1)⁵(−i)=+1。色迹i/2与内圈1/i相消，旧源C为正。正确P色迹为−i/2，相消结果直接是−1/2，因而两个物理C均为负。这个号由原作用量、P与原内圈字典导出，未读取RG或任何最终Laurent值。

## 7. FGV的新反例：此处没有成对补偿

四槽迹与两个f收缩给

$$\frac{\operatorname{Tr}(t^at^bt^ct^d)}{N_c},(f^{abe}f^{cde},f^{ace}f^{bde},f^{ade}f^{bce})
=C_A^2C_F(-\tfrac14,0,\tfrac14).$$

反转四槽词后，这三个实颜色系数保持。完整Lorentz组合为C_A²C_F r²ℓ，ℓ=(−L₁+L₃)/(4r²)，恰为当前J中的lorentzZ。四个Wilson顶点、四顶点、四个坐标传播子和一个中心积分的相位为

$$i^4\cdot(-i)\cdot(-1)^4\cdot(-i)=-1.$$

因此从原图提取的物理C是

$$\boxed{C_{\rm FGV}=-\frac{C_A^2C_F}{16}
\left(\frac{\alpha_s}{\pi}\right)^3\Gamma(2-3\epsilon)(\pi A_0)^{3\epsilon}.}$$

当前源账本写star_Wick_rotation=I、product=1，C恰为上述式的相反数。02、05、06、07、08、09的ℓ非零，所以这是六图共同的整体号失配；01、03、04仍为零。九图Lorentz多项式、原域和α测度未因符号审核改变。

独立Euclidean控制也支持该判断：e^{-S_E}的单四顶点展开给−g²Q，四个Wilson插入给i⁴=1；SU(2)取方向(e₁,e₁,e₂,e₂)，完整颜色×Lorentz为−3/4，物理密度因此为+3/4。旧相位给−3/4。它不使用Minkowski星积分相位、不使用RG和已知I。

当前冻结数学结果仍证明／计算旧根；本轮不手改旧表或上游。后续正确生产必须从本节−i星母式生成负C，再做同样的prepared与I提取；不以逐图附加−1作为永久生成规则。

## 8. 各族直接得到的C与历史对象的关系

令A₀=−2v₁₄>0、λ_b=α_s/π，先不补Euler；裸g₀的μ^{nε}按原字典保留。下表使用μ=1简写尺度，实际记录中的μ^{2nε}始终保留。

| 图类 | 统一规则直接提取的物理raw C（其余标量／颜色说明见各主稿） | 与当前存储C×J的关系 |
|---|---|---|
| 纯交换，n=1,2,3 | \(\mathcal C_\rho\lambda_b^n\Gamma(1-ε)^n(\pi A_0)^{nε}(-2)^{-n}\) | 保持，47图 |
| 树TGV | 上述正C_tree | 保持，4图；统一重推旧文稿的成对相位 |
| mixed，当前n_pure=1 | \(+C_F\lambda_b^3\Gamma(1-ε)\Gamma(2-2ε)(\pi A_0)^{3ε}/16\) | 保持，44图 |
| FGV | 上述负C_FGV | 六非零旧根相反；三零图保持零 |
| 双TGV | \(\mathcal C_\rho\lambda_b^3\Gamma(2-3ε)(\pi A_0)^{3ε}/16\)，两个导数号保留在P₁/P₂ | 保持，26图，包含9零源 |
| 单Π₁、串联、Π₂ | 同一含i定义、正Dyson链和完整横向FT所导出的[现行raw C](../../self-energy/conventions.md#2-三个分支的-c-与唯一-i) | 现行r002保持，17图；旧版导入按插入数转换 |
| TGV腿Π₁ | \(-C_AC_F\mathcal B_1\lambda_b^3\Gamma(1-ε)^2\Gamma(2-3ε)(\pi A_0)^{3ε}/[64ε(1+ε)\Gamma(2-2ε)]\) | 当前12个固定早到晚源相反 |
| 中央Γ₁ | \(-C_FC_Ac\lambda_b^3\Gamma(1-3ε)(\pi A_0)^{3ε}/128\)，c=C_A或T_Fn_f | 当前4个固定早到晚源相反 |

这张关系表用于阅读历史数据；新物理图的定义来自前七节。每个原J的变量、依赖界、积分顺序、i0和内圈定义都保持；对应整个I的号关系由积分线性和相同规范化导出。它不自动转移旧根的形式证明状态。

## 9. 从raw贡献到唯一I，不再引入另一套符号

在固定原域下，物理贡献首先是\(\mathcal A_\rho=C_{\rho,P}J_\rho\)。连通web颜色、C4、镜像／crossing和尺度提取分别执行，不能重复计数。纯交换／Π₁加裸交换的ordinary trace与web颜色不同，例如CF(CF−CA/2)→−CA CF/2；这个−1/2来自连通颜色，不是闭费米圈号。

补Euler恰为e^{-nγ_Eε}一次，并验证相应Γ乘积的线性γ_E项。普通MS耦合反项随后代入真实裸耦合展开，其−β₀/ε等符号来自Z_α；轮廓Z减除来自log W_R=log W_B−log Z_C。两者都不能用来补未确定的原图符号。

若C已写成a₁\(\mathcal N_nS_{14}R_n\)，则

$$I_\rho=\frac{m_\rho a_1}{2}R_n
\frac{J_\rho(z)+z^{nε}J_\rho(1/z)}{1+z^{nε}},$$
$$\frac{m_\rho}{2}C[J(z)+z^{nε}J(1/z)]=\mathcal N_n(S_{14}+S_{12})I_\rho.$$

这条恒等式是旧、新对象比较的唯一提取规则；C、J之间移动一个常数允许，但完整乘积和原域必须相同。采用web根的解表不再二次扣除W₁W₂。RG通过只核验所列结构，不能替代上述符号推导，也不能证明全部单图极点或有限项。

## 10. 验收与仍开放的工作

2026-10-06原生迁移完成：共同规则已贯穿实际qgraf-pipeline→Montopar→原生prepared→I。163源、所有raw/IBP C及17个自能值通过；814/94回归与8个安装路径CLI重放通过，发布22个物理I及Γ₁03的z表示。319个系数与31图217个组装位置相符。新I的assembly_sign=1；FGV、腿Π₁、Γ₁号在生成C前由本页规则产生，旧结果仅通过同J/域/分支和线性提供积分值，不作为符号定义。默认47根与40个绑定历史文件保持，formal=not_run，未提交Git。[验收与版本选择](../../../workflow/checks/shared-sign-production/README.md)。下两段记录先前理论阶段的验收范围。

共同约定和163源的符号连接已从明确作用量、P、FT和原C/J给出。Γ₁的−1来源现在可以与可信树父核在同一规则中比较；没有使用“非纯胶子统一变号”。FGV新增反例保存在实际检查中，既有来源审计、数学I以及计算完成计数未被重写。

本轮理论验收不包含外部生产回写、新数学根登记、严格复运动学全域证明或重新求解全部内圈。原生生产整合要把真实Wilson词和星FT在源端固定，再输出上述C；已有冻结对象应版本化转换，不能修改原字节冒充旧证据。本页的规则审核与生产迁移是两个可分别验收的工作。

[原始文献及PDF哈希](../../../workflow/checks/shared-sign-conventions/references.json) · [602项共同规则检查](../../../workflow/checks/shared-sign-conventions/rules-checks.json) · [163源实际C/J连接](../../../workflow/checks/shared-sign-conventions/corpus-audit.json)。
