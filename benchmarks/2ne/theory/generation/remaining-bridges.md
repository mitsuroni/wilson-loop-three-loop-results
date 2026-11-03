# Step 3 完成：后16图的全ε定义桥

日期：2026-09-25。本文将[统一母族G1–G6](mothers.md)连接到同通道双核、双通道、混合核数及三核的现行结果。连同首批空集与单对边31图，**47图现在都已从原有序端点生成规则连接到现行完整ε的I**。连接包括参数域、配对、字槽、尺度和颜色；不是只比较有限项。

记上横线为z→1/z，ℓ=log z、Z=exp(−γ_E ε)Γ(1−ε)。所有等式先在z>0、0<Re ε<1的公共收敛域建立，再沿用各族已经给出的ε=0附近延拓。本文的精确代数和有限积分检查不是Lean证明；formal=not_run。逐图完整I见[当前生成表](mother-map.md)，机器数据采用[mother-map-r002.json](mother-map-r002.json)。

## 1. 双核Q是字母族的重复A子族

沿用原[双标记理论§7](../../families/same-channel-two/unified-structure.md#7-当前七图所需的重复-a-双指标子族)，统一定义直接给

$$
Q_{\pi;m,n}
=m!n!\mathcal T^{[2]}_{\pi;\varnothing,\varnothing;A^m,A^n}
=(1+z)^2z^{-\epsilon}\int_{\Omega_2}
\prod_{i=1}^2D_z(t_i,s_{\pi(i)})^{\epsilon-1}s_1^{m\epsilon}s_2^{n\epsilon}
\,d^2t\,d^2s.
\tag{C1}
$$

Q中两个下标分别属于较小、较大的s端点。本批只需00、01、10；它们不是ε阶数。阶乘来自Fhat_(A^m)=x^(mε)/m!，不能在任意长度推广时漏掉。

左边t上的插入先用G4交换两条有序边，π同时变成π⁻¹。S₂中两种排列都等于自身的逆。B插入则用Fhat_B(0,x)=1−(1−x)^ε及同时反射：

$$
\mathcal T^{[2]}_{\pi;\varnothing,\varnothing;B,\varnothing}
=Q_{\pi;0,0}-\overline Q_{\rho\pi\rho;0,1},\qquad
\mathcal T^{[2]}_{\pi;\varnothing,\varnothing;\varnothing,B}
=Q_{\pi;0,0}-\overline Q_{\rho\pi\rho;1,0},\quad\rho=(12).
\tag{C2}
$$

反射会交换较小、较大端点，不能只交换A/B而保留插入位置。对S₂有ρπρ=π；零插入Q各自crossing偶，含插入的Q不作此假设。上述规则只依赖字槽，不依赖图号。

## 2. 旧F级数与Q差分：必须保留ε

旧结果中的F_base、F_id_after、F_id_middle、F_swap_after是[完整块级数](../../families/same-channel-two/full-series.json)的函数名。按其权重定义，精确字典为

$$
\begin{aligned}
F_{\rm base}&=Q_{\rm id;0,0},\\
F_{\rm id,after}&=\frac{Q_{\rm id;0,0}-Q_{\rm id;0,1}}{\epsilon},&
F_{\rm swap,after}&=\frac{Q_{\rm swap;0,0}-Q_{\rm swap;0,1}}{\epsilon},\\
F_{\rm id,middle}&=\frac{Q_{\rm id;0,1}-Q_{\rm id;1,0}}{\epsilon}.
\tag{C3}
\end{aligned}
$$

证明可以在每个有限求和块完成。若两个核的单项式为t₁^a s_(π(1))^b t₂^c s_(π(2))^d，令(β,δ)=(b,d)或(d,b)，则Q的单项式矩是

$$
\mathcal R_{m,n}
=\frac1{(a+1)(a+c+2)(\beta+1+m\epsilon)
(\beta+\delta+2+(m+n)\epsilon)}.
\tag{C4}
$$

精确相减得到

$$
\begin{aligned}
\frac{\mathcal R_{0,0}-\mathcal R_{0,1}}{\epsilon}
&=\frac1{(a+1)(a+c+2)(\beta+1)(\beta+\delta+2)(\beta+\delta+\epsilon+2)},\\
\frac{\mathcal R_{0,1}-\mathcal R_{1,0}}{\epsilon}
&=\frac1{(a+1)(a+c+2)(\beta+1)(\beta+\epsilon+1)(\beta+\delta+\epsilon+2)}.
\tag{C5}
\end{aligned}
$$

这正是旧F_A、F_B的两个分母。公共尺度、Pochhammer权重及求和指标也完全相同；检查程序读取实际保存的四条F级数逐式核对，而不是只按函数名建立别名。分块绝对局部一致收敛的依据沿用[双核方法§2、4](../../families/same-channel-two/method.md)。

C3中的ε=0按可去极限理解：权重分别为(1−w^ε)/ε与(w^ε−v^ε)/ε，其极限为−log w与log(w/v)。本步骤没有重新计算这些极限积分的Laurent系数。统一式中的ε⁻²乘Q差分，正好还原旧表达式中的ε⁻¹乘F。

## 3. 后16图实际产生的B组合

继续采用统一式

$$
J_G=z^{(a+r/2)\epsilon}\epsilon^{-2(n-r)}B_G,\qquad
I_G=\frac{\kappa_GZ^n}{\epsilon^{2(n-r)}\cosh(n\epsilon\ell/2)}
\frac{e^{\delta_G\epsilon\ell}B_G+e^{-\delta_G\epsilon\ell}\overline B_G}{2},
\quad\delta_G=a+\frac{r-n}{2}.
\tag{C6}
$$

下面的B由四边间隔字的压缩生成；a、r、κ逐条来自原输入，完整值见生成表。Q_id、Q_swap均简写零插入00。

| 图 | 从统一规则得到的B_G | 与旧求值的连接 |
|---|---|---|
| 4e-04 | $Q_{\rm id}$ | F_base |
| 6e-10 | $Q_{\rm id}-Q_{\rm id;0,1}$ | ε F_id,after |
| 6e-11 | $Q_{\rm id;0,1}-Q_{\rm id;1,0}$ | ε F_id,middle |
| 6e-12 | $Q_{\rm swap}-Q_{\rm swap;0,1}$ | ε F_swap,after |
| 6e-23 | $Q_{\rm id}-\overline Q_{\rm id;0,1}$ | ε F_id,after(1/z) |
| 6e-24 | $Q_{\rm swap}-\overline Q_{\rm swap;0,1}$ | ε F_swap,after(1/z) |
| 6e-25 | $\overline Q_{\rm id;0,1}-\overline Q_{\rm id;1,0}$ | ε F_id,middle(1/z) |
| 4e-06 | $H_0^2$ | 两个零插入单核 |
| 6e-32、6e-34 | $H_1(H_0-H_1)$ | 两通道各自端点权重 |
| 6e-33 | $(H_0-H_1)^2$ | 两通道的相同差分 |
| 6e-37 | $H_0Q_{\rm id}$ | 单核乘零插入双核 |
| 6e-38 | $H_0Q_{\rm swap}$ | 单核乘另一配对双核 |
| 6e-20 | $R_{123}$ | 三核零插入 |
| 6e-21 | $R_{132}$ | 三核零插入 |
| 6e-22 | $R_{213}$ | 三核零插入 |

双核的三对镜像图有δ=±1/2，B同时互为crossing。代入C6后10=23、11=25、12=24的完整I相等；这一步还核对了每对相同的κ。双通道32与34的B、外部尺度、κ均相同，所以原J与完整I都相同。来源推导分别见[双核](../../families/same-channel-two/method.md)、[双通道](../../families/double-opposite/method.md)。

## 4. 混合通道的Q_swap定义也要连接

旧混合结果把Q_swap定义为H₀²/2−Q_id。与统一母族连接需要区域分割：两个独立单核的t、s分别有两种排序，四个区域按交换核标签成对相等，故

$$Q_{\rm id}+Q_{\rm swap}=\frac12H_0^2.\tag{C7}$$

分组合法性由共同收敛域中的绝对可积性保证。两个单核负号相乘为正，尺度z^(−ε)也恰好一致。这个全ε关系使[混合族现行定义](../../families/mixed-opposite/method.md)的Q_swap与G1的swap配对核成为同一个函数。

因此37、38的完整I分别为−Z³H₀Q_id/[16 cosh(3εℓ/2)]与−Z³H₀Q_swap/[32 cosh(3εℓ/2)]。二者已有有限项相同，但不能用Q_swap=2Q_id作为全ε关系：在ε=1，两个正则母核都等于(1+z)²/(4z)，足以反驳该全ε比例。此处没有在Γ(1−ε)极点处计算I。

## 5. 三核R：符号与crossing排列

原[三标记定义](../../families/same-channel-three/method.md)正是统一母族的零插入槽：

$$
R_\pi=\mathcal T^{[3]}_{\pi;\varnothing,\varnothing,\varnothing;\varnothing,\varnothing,\varnothing},\qquad
\overline R_\pi=R_{\rho\pi\rho},\quad\rho(i)=4-i.
\tag{C8}
$$

三个核的公共因子是−(1+z)³z^(−3ε/2)，负号没有吸收到κ中。对123保持；132与213互换。两者各自等于逆排列并不意味着各自crossing偶。

零插入统一块级数按三个t和三个s的有序单项式矩求积，恰好给出旧R的三重级数。检查除读取现行R定义外，还对S₃全部六种配对的有限块，用“先展开整个六变量多项式”的独立路线比较旧实现。由C6得到

$$
I_{20}=-\frac{Z^3R_{123}}{16\cosh(3\epsilon\ell/2)},\qquad
I_{21}=I_{22}=-\frac{Z^3(R_{132}+R_{213})}{64\cosh(3\epsilon\ell/2)}.
\tag{C9}
$$

21、22的J互为crossing，完整I相同；推导中没有把两个原J强行合并。

## 6. 完成的证据与下一步

本轮从原多项式及bounds重新核对47根，统一字压缩的全部生成项与首批记录一致；47条现行完整I符号残差均为0，首31条生成表达式保持原样。后16条另核对完整J符号连接，并在两个非平凡整数ε/有理z点上将命名母积分求值与原次序J比较，共32例。

12项单元检查包括：四个旧F级数定义；72个双核有限块（保留符号ε）；48个三核有限块；32个双核单字母槽的反射/交换；零插入乘积，以及漏ε、错crossing、错配对、有限阶比例冒充全ε关系等负控制。实际输出见[本批检查](../../../../workflow/checks/2ne-mothers-complete/README.md)。

Step 3到此完成的是**全47图的统一生成及与已有全ε求值的连接**。没有新增物理结果或修改旧结果，也未声称有限特殊函数闭式、最小主积分数、一般复z或形式整链。下一步Step 4登记全ε反射、域分割和乘积关系，并明确它们与固定ε阶关系的区别；Step 5、6再做系数生成和329位置的新链验收。
