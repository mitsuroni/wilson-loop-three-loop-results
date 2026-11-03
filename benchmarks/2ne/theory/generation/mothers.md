# Step 3：由有序端点生成同一个母积分系统

日期：2026-09-25。本文接续 [有序字理论](ordered-words.md) 的 W11、W12。完成范围是：47个当前纯交换根的统一母族生成表与全部现行完整ε结果的精确连接。本文先建立共同定义和首31图求值桥；后16图的Q/R/F连接见[补充推导](remaining-bridges.md)。本文给出分析推导与辅助检查，formal=not_run。

核心区分是：**生成母族按核数、配对和端点字定义；最终基底由这些母族在指定阶数的关系约化得到。** 因而M、N、τ、E不参与本步骤的输入分类，也不用于反推生成规则。

## 1. 只定义一次的单通道母族

沿用A/B字及正规化字积分Fhat。令r为一个对边通道的核数，π是从左边由小到大的端点到右边由小到大端点的配对排列，定义

$$
\begin{gathered}
\Omega_r=\{0<t_1<\cdots<t_r<1,\ 0<s_1<\cdots<s_r<1\},\\
\mathcal T^{[r]}_{\pi;\mathbf u;\mathbf v}(\epsilon,z)
=[-(1+z)]^r z^{-r\epsilon/2}
\int_{\Omega_r}\prod_{i=1}^r D_z(t_i,s_{\pi(i)})^{\epsilon-1}
\prod_{i=1}^r\widehat F_{u_i}(0,t_i)\widehat F_{v_i}(0,s_i)
\,d^rt\,d^rs .
\tag{G1}
\end{gathered}
$$

r=0采用空母积分1；所有全边字常数留在系数内。正实z及Re ε>0下，Fhat在端点区间有统一有限界，未加权对边核按W12的独立端点上界可积，故定义成立。连接I时先取公共条带0<Re ε<1；ε=0附近的延拓须由已有求值或后续减除建立。

三个不同的下标各有作用：r记录核数；π记录配对；字槽记录相邻端点的排序。不能用圈阶替代它们，也不能删掉暂未出现在语料中的排列。例如r=3有六种配对，当前三条图只选中了其中三种。

为了写成真正的生成系统，引入两个非交换形式变量X_A、X_B：

$$
\mathbf U_\epsilon(x)=\sum_w\widehat F_w(0,x)X_w,\qquad
\frac{d\mathbf U_\epsilon}{dx}
=\mathbf U_\epsilon\,\epsilon
\big[x^{\epsilon-1}X_A+(1-x)^{\epsilon-1}X_B\big].
\tag{G2}
$$

空字系数为1，其余在x=0为0。乘法在右侧添加末字母，与W5一致。给2r个端点各一份独立字槽，用G1的核积分这些生成级数，逐槽取字系数就得到全部母积分。这里按有限字长截断操作，不假设无限非交换级数作为算子收敛。

对于n≤3，单通道留在端点上的总字长≤n−r。因此[1]只需总字长≤2，[2]≤1，[3]为零插入。r=1的原始有序字对共有17个，利用两边交换后可用10个覆盖标签；这只是覆盖上界，不是独立主函数数。

## 2. 两个通道的方向与精确尺度

47图从原多项式重建后的实际方向全部为：24通道使用D_z，13通道使用D_tilde。逐条识别记录见 [输入表](input-map.json)。关键恒等式是

$$
\widetilde D_z=zD_{1/z},\qquad
[-(1+z)]^r\prod_i\widetilde D_z^{\epsilon-1}
=z^{r\epsilon}[-(1+1/z)]^r\prod_iD_{1/z}^{\epsilon-1}.
\tag{G3}
$$

所以含任意相同端点权重的13通道积分正好等于
z^(rε/2)乘同一母族在1/z的取值；24通道是z^(rε/2)乘母族在z的取值。此路线上不反射端点，因而不改字或π。也可先反射端点，但必须同时反射间隔字与排列，不能只把核名换掉。

两边交换则给出一般全ε恒等式

$$
\mathcal T^{[r]}_{\pi;\mathbf u;\mathbf v}
=\mathcal T^{[r]}_{\pi^{-1};\mathbf v;\mathbf u}.
\tag{G4}
$$

同时反射两条边时，位置次序倒转，π变为ρπρ，ρ(i)=r+1−i；下端点字会先变成上区间字，再由W8压回。只有零插入时可以直接写
T_pi(ε,1/z)=T_(ρπρ)(ε,z)。r=2的id/swap分别保持，r=3的132与213交换。

## 3. 全部47图的公共J及I

把四条边的间隔权重用W11压缩。每项都是整数shuffle系数、全边常数和两个通道各自的端点字。两通道分别占据边{1,3}与{2,4}，所以压缩后的每一项在剩余变量上因子化。邻边传播子连接的两端早已由W2拆开；没有跨通道不等式。原域到分组积分的合法性来自W12的绝对可积论证。

$$
\begin{aligned}
B_G(\epsilon,z)
&=\sum_\lambda c_{G,\lambda}(\epsilon)
\mathcal T^{[r_{13}]}_{\pi_{13};\mathbf u_{13};\mathbf v_{13}}(\epsilon,1/z)
\mathcal T^{[r_{24}]}_{\pi_{24};\mathbf u_{24};\mathbf v_{24}}(\epsilon,z),\\
J_G(\epsilon,z)&=z^{(a+r/2)\epsilon}\epsilon^{-2(n-r)}B_G(\epsilon,z),\qquad r=r_{13}+r_{24}.
\tag{G5}
\end{aligned}
$$

a仅计q=z的邻边数。c是全边字常数的有限整系数组合，既不依赖z，也不从已有答案拟合。本批所需常数长度≤3，使用[空集Gamma表](../../families/empty/method.md#3-本批所需的全部字积分)，记

$$R_2=\frac{\Gamma(1+\epsilon)^2}{\Gamma(1+2\epsilon)},\qquad
R_{12}=\frac{\Gamma(1+\epsilon)\Gamma(1+2\epsilon)}{\Gamma(1+3\epsilon)}.$$

唯一物理归一化W13随之变成

$$
I_G=\frac{\kappa_G Z^n}{\epsilon^{2(n-r)}\cosh(n\epsilon\ell/2)}
\frac{e^{\delta_G\epsilon\ell}B_G(\epsilon,z)+e^{-\delta_G\epsilon\ell}B_G(\epsilon,1/z)}2,
\quad\delta_G=a+\frac{r-n}{2},\quad\ell=\log z.
\tag{G6}
$$

其中Z=exp(−γ_E ε)Γ(1−ε)、κ=sym·a₁/[2(−2)^n]。Euler只计一次。δ与r都不能误用历史单对边文档中也命名为r的“q=z邻边数”。[逐图生成表](mother-map.md)列出每条的b=a+r/2、显式极点预算、κ及B；[JSON](mother-map-r002.json)另保留每个字槽、排列和当前结果路径。

| 原族 | 同一系统中的组成 | 本轮生成 | 本轮与现行完整I逐式比较 |
|---|---|---:|---:|
| [] | 全边常数，r=0 | 13 | 13 |
| [1] | T^[1]，总字长≤2 | 18 | 18 |
| [2] | T^[2]，π∈S₂，总字长≤1 | 7 | 7 |
| [1,1] | 两个T^[1]的有限乘积和 | 4 | 4 |
| [1,2] | 零插入T^[1]与T^[2]乘积 | 2 | 2 |
| [3] | 零插入T^[3] | 3 | 3 |

乘积结构是本批纯交换几何的结论。存在内部3g顶点时，Feynman参数可能产生跨端点的相关权重，不能直接套用G5。

## 4. 空集与单对边的31条完整求值桥

空集的B仅为四条全边常数之积。G6立即给Gamma乘双曲函数比；13条均与现行结果精确相等。

对单对边，简写a_x=x^ε、b_x=(1−x)^ε、K_x=Fhat_BA(0,x)。所有允许的下端点字为

$$
\begin{aligned}
\widehat F_A&=a_x,&\widehat F_B&=1-b_x,&
\widehat F_{AA}&=a_x^2/2,&\widehat F_{BB}&=(1-b_x)^2/2,\\
\widehat F_{AB}&=a_x(1-b_x)-K_x,&\widehat F_{BA}&=K_x,&
K_x+K_{1-x}&=a_x+b_x-a_xb_x-R_2.
\tag{G7}
\end{aligned}
$$

最后一式由全BA区间的分割与几何反射得到，参见[现行推导§3](../../families/single-opposite/basis.md#3-端点权重与-ba-项的处理)。crossing将a与b交换，K换为a+b−ab−R₂−K。本批需K的图在G6中δ=0，两个K精确消去；单独的J一般仍含K。检查器拒绝在crossing投影前把K丢掉。

剩下幂权重，采用已经求值的同一个函数H_(ab;cd)：它等于G1的r=1核对权重t^(aε)(1−t)^(bε)s^(cε)(1−s)^(dε)的积分泛函。其标准Meijer G表达、围道与域见[现行求值式](../../families/single-opposite/basis.md#4-六种共同幂权重及标准特殊函数求值)。实际只需以下固定字典，不增加新的自由函数：

| 幂权重 | 现行求值名 | 统一母族中的来源 |
|---|---|---|
| 1 | H0 | T^(1)_(∅;∅) |
| s^ε | H1 | T^(1)_(∅;A) |
| s^(2ε) | H2 | 2T^(1)_(∅;AA) |
| s^ε(1−s)^ε | HW | T^(1)_(∅;A)−T^(1)_(∅;AB)−T^(1)_(∅;BA) |
| t^εs^ε | HD | T^(1)_(A;A) |
| t^ε(1−s)^ε | HX | T^(1)_(∅;A)−T^(1)_(A;B) |

两端交换及同时交换a↔b、c↔d给其余取值。H0、HW、HX为crossing偶；H1、H2、HD要保留镜像。18条完整I均从压缩生成的权重得到，并与现行修订作符号比较，未调用低阶M/N关系来弥补全ε差异。

例如6e-36和6e-39分别给B=H0−2H1+HD与其镜像，δ分别为1和−1，因此投影后的I相同。6e-04与6e-05虽然J相同，κ却相差2，完整I也相差2。两个检查阻止把“相同母族”误用为“同一个物理积分”。

## 5. 从生成规则到已积完的统一级数

G1说明函数从何处产生；另可给全部47图一个共同的显式级数求值规则。此处直接在原有序端点上写出，以避免把母积分的命名当作求值。

令b_k(ε)=(1−ε)_k/k!。对每条对边核d_e定义q_e=1−d_e/(1+z)：

$$
q_D(x,y)=\frac{z+x+y}{1+z}-xy,\qquad
q_{\widetilde D}(x,y)=\frac{1+z(x+y)}{1+z}-xy.
\tag{G8}
$$

在正实z的内部域有0≤q_e<1。展开每个对边核的(1−q_e)^(ε−1)，以及每个B端点的(1−x)^(ε−1)；A端点保留x^(ε−1)。给对边核分配指标k_e，给B端点分配指标m_v。定义有限多项式

$$
P_{\mathbf k}(\mathbf x)=\prod_{e\in E_O}q_e(\mathbf x)^{k_e}
=\sum_{\boldsymbol\nu}p_{\mathbf k,\boldsymbol\nu}(z)\mathbf x^{\boldsymbol\nu},\qquad
\alpha_v=\nu_v+(\epsilon-1)\mathbf1_{A}(v)+m_v\mathbf1_B(v).
$$

一条按小到大排列为v₁…v_L的边，其单项式积分已经显式积完：

$$
\int_{0<x_1<\cdots<x_L<1}\prod_jx_j^{\alpha_j}\,d^Lx
=\prod_{j=1}^L\frac1{\sum_{h=1}^j(\alpha_h+1)}.
\tag{G9}
$$

原bounds正好是这四条有序链。由原分子及核尺度得到

$$
\boxed{
J_G=(-1)^r(1+z)^{r\epsilon}z^{a\epsilon}
\sum_{\mathbf k,\mathbf m\ge0}
\left(\prod_e b_{k_e}(\epsilon)\right)
\left(\prod_{v\in B}b_{m_v}(\epsilon)\right)
\left[
\sum_{\boldsymbol\nu}p_{\mathbf k,\boldsymbol\nu}(z)
\prod_{i=1}^4\prod_{j=1}^{L_i}
\frac1{\sum_{h=1}^j(\alpha_{v_{ih}}+1)}
\right].
}
\tag{G10}
$$

每个方括号都是有限有理运算，不再含参数积分。该式无需新增一个按图号定义的函数。r=0时没有对边指标；r=2、3的零插入情形包含已有双重、三重块级数。

收敛先在0<Re ε<1建立。对条带内的ε紧集，取0<σ<min Re ε；Pochhammer系数的模可由常数乘(1−σ)_k/k!支配。对方括号，先还原为完整非负q多项式幂的积分，再取模，而不是对展开后每个带符号单项式分别求和。B端点及对边核的绝对和于是受(1−x)^(σ−1)和d_e^(σ−1)控制；A端点同样受x^(σ−1)控制。放宽排序域后，独立端点上界按各传播子分解且可积。因此按方括号分块的级数绝对、局部一致收敛。z也可限制在任意正实紧区间。

G10先给共同收敛域上的精确解析表示；向ε=0的亚纯延拓沿用具体求值或减除论证，不能直接在未减除级数内展开发散端点。正整数ε时各二项式级数终止，可作为额外的精确实现检查；这不同于在Gamma极点上检查I。本轮没有声称这一高重级数适合高效数值计算，也没有给非整数截断的认证尾界。

## 6. 对“统一基底来源”的实际推进

目前可以从原图唯一确定：核数、通道、配对、字长预算、生成母积分及固定归一化。空集常数、单端插入和双端插入都是同一规则的不同系数槽。τ是后续取ε系数的记号；M/N是某些系数可采用的矩坐标。增加允许字长会扩大候选覆盖空间，不需要修改母族定义。

这仍不证明M/N在任意阶闭合，也不证明出现新字槽就产生新独立函数。端点A/B字母与运动学Hlog字母属于不同变量，不能从前者推定后者只有0、−1。关于字母表、独立维数与最小性的结论，应在后续关系、微分和系数约化中分别验证。

后16图的Q/R及乘积表达现已[逐式接到G5/G6](remaining-bridges.md)，Step 3全47图连接完成。下一步Step 4整理全ε关系；Step 5、6才重建固定阶基底与全部329个系数位置。[首批检查](../../../../workflow/checks/2ne-mothers/README.md)与[完成检查](../../../../workflow/checks/2ne-mothers-complete/README.md)分别保留；现行结果、catalog、历史证据及Kernel均保持。
