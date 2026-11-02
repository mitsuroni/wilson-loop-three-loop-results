# 两圈组装前审计：12 个 I 的物理权重

日期：2026-09-24。完成用户所指的“剩余工作1”：来源、符号、颜色、轨道重数、尺度和展开预算。**现有 6+4+2 个积分足以继续，但 Pi1 不能按旧候选式的正号加入物理指数。** 来源中最外层负号与当前算符的插入约定重复；下面从传播子展开推导转换，保留原 I、C/J 和所有结果字节。

本次交付是[组装字典](assembly-bridge-r001.json)，不是已完成的两圈有限结果。实际生产 API 重放与辅助精确检查见[复核记录](../../../workflow/checks/rg-two-loop-audit/README.md)。没有改上游代码、注册 Kernel 根或新增形式收据。

后续执行更新：本页来源审计之后，[两圈完整组装与RG/cusp](derivation.md)已完成，且完整有限函数与作者WRen.m解析一致。本页的“未完成组装”描述保留来源审计当时的阶段；最新结果见[result-r001.json](result-r001.json)。

当前约定更新：本页审计对象是旧r001。项目已发布[统一自能2.0](../../self-energy/conventions.md)，当前I吸收了这里推导的物理号，直接正号求和；不得对当前r002再乘−1。

## 1. 起点和审计范围

采用与[一圈](../one-loop/derivation.md)一致的归一化迹 Wilson 圈：

$$
W_B=\frac1{N_c}\left\langle\operatorname{Tr}
P\exp\left(ig_0\oint_C A_\mu^a T^a dx^\mu\right)\right\rangle,
\qquad D^{(0)}_{\mu\nu}(k)=\frac{-ig_{\mu\nu}}{k^2+i0}.
$$

真空泡由真空归一化消去。沿轮廓有序积分已经吸收指数展开的阶乘；不同边上的两个端点没有额外的 1/2。四条边均类光、轮廓无自交，取解析不变量起始区 s,t<0，z=s/t>0。度规为 (+−−−)，幂按共同 Feynman 边界值解释。

两圈来源是实际 qgraf-pipeline 的组合生成器、分类器和结构化导出 API；不是重新运行 QGRAF 二进制。冻结的 12 个代表与当前上游输出逐字段一致。独立枚举有序端点和 C4 作用，得到以下筛选表。

| 轨道 | 按轮廓排列的附着边 | 原颜色 / 指数颜色 | Lorentz 或参数核 | m | 旧 I 到物理贡献的乘数 |
|---|---|---|---|---:|---:|
| 4e-01 | 3,3,4,4 | CF(CF−CA/2) → −CFCA/2 | z² f1^(ε−1) f2^(ε−1) | 4 | +1 |
| 4e-02 | 2,3,4,4 | 同上 | −z(1+z) f1^(ε−1) f2^(ε−1) | 4 | +1 |
| 4e-03 | 2,3,3,4 | 同上 | z f1^(ε−1) f2^(ε−1) | 4 | +1 |
| 4e-04 | 2,2,4,4 | 同上 | (1+z)² f1^(ε−1) f2^(ε−1) | 2 | +1 |
| 4e-05 | 2,2,3,4 | 同上 | −(1+z) f1^(ε−1) f2^(ε−1) | 4 | +1 |
| 4e-06 | 1,2,3,4 | 同上 | (1+z)² f1^(ε−1) f2^(ε−1) | 1 | +1 |
| 3e-3g-01 | 3,4,4 | iCFCA/2；相位与 C 合并 | α 单纯形、完整 P/S^(2−2ε) | 4 | +1 |
| 3e-3g-02 | 3,3,4 | 同上 | 同上；与01是不同来源轨道 | 4 | +1 |
| 3e-3g-03 | 2,4,4 | 同上 | 同上；保留同边接触面 | 4 | +1 |
| 3e-3g-04 | 2,3,4 | 同上 | 同上；保留完整余项 | 4 | +1 |
| 2e-Pi1-01 = BA | 3,4 | CF B1(ε)；两味道通道 | 完整 K，见原主稿式(14) | 4 | **−1** |
| 2e-Pi1-02 = BO | 2,4 | 同上 | 同上 | 2 | **−1** |

表中全部是非零、保留于指数的轨道，没有再次丢项。f、P、S、K 的实际多项式已由端点坐标与 Gram 矩阵独立重建，并与保存的 J、实际 montopar 函数比较。这里 v_ij=v_i·v_j；不采用部分外部笔记的 v_ij=2v_i·v_j 记号。

| 未进入这12个代表的对象 | 排除阶段与理由 |
|---|---|
| 裸传播子两端在同一类光直边 | Lorentz 收缩给出 v_i²=0 |
| 三个 TGV 外腿均在同一类光直边 | 三个顶点项均含该边的 v_i²；参数核为零 |
| 同边 Pi1 | 在共同维数解析延拓下没有非零不变量，为无尺度项；不能仅用“一个 v_i²”代替完整张量论证 |
| 无质量四胶子 tadpole | 内部动量积分无尺度，在维数正规化中为零 |
| 两裸交换的非交错配对 | 在 W2 中保留 CF²；与交错图的 CF² 部分一起组成 W1²/2，在 log W 中消去 |
| Pi1 内部胶子、鬼、夸克图 | 已在 B1 中求和，不再另加同一批泡图或闭合费米子负号 |

独立的带边标记配置总数分别为 **19、16、6**，合计41；C4 商后分别为6、4、2。这里既未作镜像商，也未把“解相等”当成删图理由。这是指定两圈、无质量、Feynman 规范图集的完整性检查，不推广为任意高圈筛选证明。

## 2. 为什么这些颜色直接给 w2

由生成元恒等式

$$
T^aT^bT^a=(C_F-C_A/2)T^b
$$

得到交错配对的完整迹 CF²−CFCA/2；另外两种配对均为 CF²。对四个轮廓点，交换核的阿贝尔部分满足

$$
\frac12(W_1)^2
=\frac{g_0^4C_F^2}{8}\int_{C^4}K_{12}K_{34}
=g_0^4C_F^2\int_{\tau_1<\tau_2<\tau_3<\tau_4}
(K_{12}K_{34}+K_{13}K_{24}+K_{14}K_{23}).
$$

第二步按24种点序分区，每种无标号配对出现8次；使用核对称性，并在共同调节下作解析延拓。脚本枚举了这24个排列及其8+8+8分组。因而投影后的 −CFCA/2 已是 w2 的交换颜色；再次减 W1²/2 会引入本不应出现的 CF² 项。TGV 的反对称颜色和单 Pi1 直接是连通颜色贡献。[两圈指数化的文献对应：DHKS §2.2](https://arxiv.org/pdf/0709.2368)。

Pi1 的完整系数为

$$
\mathcal B_1(\varepsilon)=
\frac{(10-6\varepsilon)C_A-8(1-\varepsilon)T_Fn_f}{3-2\varepsilon}.
$$

故两圈物理指数的颜色为 CFCA 和 CFTFnf。该结论限于两圈；不是把旧 a1 规则提升为全阶指数化定理。

## 3. Pi1 绝对符号：从完整传播子推导

用 p1(k²) 记文献的无量纲一圈逆传播子修正，避免与本项目参数积分 J 混名。1996年(C.1)+(C.3)在 Feynman 规范给出

$$
p_1(k^2)=-\mathcal A_1(\varepsilon)(-k^2-i0)^{-\varepsilon},
\quad
\mathcal A_1=
\frac{g_0^2}{2(4\pi)^{2-\varepsilon}}
\frac{\Gamma(\varepsilon)\Gamma^2(1-\varepsilon)}{\Gamma(2-2\varepsilon)}
\mathcal B_1(\varepsilon).
$$

这也保留了夸克与胶子/鬼项的相对号。[DOT 式(2.7)、附录C](https://arxiv.org/pdf/hep-ph/9605348)。令 P_T=g−kk/k²，恢复项目使用的 Feynman 传播子相位，完整横向传播子写成

$$
D_T=\frac{-iP_T}{k^2+i0}\frac1{1+p_1+\cdots},
\qquad \delta D_T=-D_T^{(0)}p_1.
$$

分母 1+Π 的约定可直接核对 [Ruijl 等式(2.1)、(2.4)](https://arxiv.org/pdf/1703.08532#page=4)。这里是逆函数的一阶展开，并未引用 cusp 数值。

当前[Pi1主稿](../../self-energy/pi1/theory.md)式(3)中的算符恰为

$$
\Pi_{\rm current}^{\mu\nu}
=-i(k^2g^{\mu\nu}-k^\mu k^\nu)p_1
=i\mathcal A_1(k^2g^{\mu\nu}-k^\mu k^\nu)(-k^2-i0)^{-\varepsilon}.
$$

在横向子空间，用 P_T²=P_T，直接得到

$$
D^{(0)}\Pi_{\rm current}D^{(0)}
=\left(\frac{-i}{k^2}\right)(-ik^2p_1)
\left(\frac{-i}{k^2}\right)P_T
=\frac{ip_1}{k^2}P_T
=-D_T^{(0)}p_1=\delta D_T.
$$

所以相应 Wilson 图是两个原有 ig0 发射顶点乘这个传播子修正；外面没有新的负号。与之相比，当前来源及主稿式(1)定义的是

$$
\mathcal A_{\rm src}=-(ig_0v_i)(ig_0v_j)
D^{(0)}\Pi_{\rm current}D^{(0)},
\qquad
\boxed{\mathcal A_{\rm Wilson}=-\mathcal A_{\rm src}.}
$$

来源定位有三层：冻结的 `qgraf_pipeline/internal_types.py:get_vertex_info` 对 N=2 设置 `sign=-1`；`emitter.py` 实际把它写到图式最外面；`montopar/self_energy.py` 的 `source_external_sign=-1` 忠实翻译了这个输入。当前 C/J 与旧 I 对原输入的求值并未因此被否定；失配发生在**该输入与物理 Wilson 图的连接**。

如果坚持最外 −1 的图式，插入符号必须改成 −Π_current；如果保留当前已核对的 Π_current，物理图就应去掉外负号。两种源记法等价，但不能同时保留当前两处符号。**本轮不修改两处生产约定，也不重新定义 I；组装时明确乘 −1。** 三圈单插入、Pi2 和串联的来源约定应另行审计，不能把本次两圈桥直接当作它们已修复的证据。

### 3.1 对用户所贴两式的再次核对：显式保留全部 i

用户再次指出最初的“外负号 + 含 i 的 Π”定义。本次另核对 Grozin, *Lectures on QED and QCD* 的式(3.11)、(3.15)–(3.16)、(4.20)–(4.25)：它明确画出含 i 的 1PI 插入，无需从省略相位的传播子约定推测。[原文](https://arxiv.org/pdf/hep-ph/0508242)。

令 a(k²) 为用户 Π_scalar 除以 i 后的整个标量（包括耦合、Gamma、颜色和动量幂）。由该文 G(1,1) 的 Gamma 式(2.23)，胶子/鬼与夸克的系数相加为

$$
a(k^2)=\frac{g_0^2}{(4\pi)^{D/2}}
(-k^2-i0)^{D/2-2}
\frac{\Gamma(2-D/2)\Gamma^2(D/2-1)}{\Gamma(D-2)}
\left[\frac{(3D-2)C_A}{2(D-1)}
-\frac{2(D-2)T_Fn_f}{D-1}\right].
$$

取 TF=1/2、nf=NF，这恰好就是用户公式除以 i 的结果。**标量公式本身没有错。** 作者的含 i 插入因此正是用户的 Π：

$$
\Pi_{\rm user}=ik^2P_Ta(k^2).
$$

下面用 D_F 表示包含 −i 的传播子，区别于教材剥离 −i 后的 D。在横向子空间，完整传播子和一阶修正为

$$
D_{F,T}=\frac{-iP_T}{k^2+i0}\frac1{1-a},
\qquad
\delta D_{F,T}=D_{F,T}^{(0)}a
=D_F^{(0)}\Pi_{\rm user}D_F^{(0)}.
$$

最常见的混淆来自逆矩阵公式。若定义 R=δ(D_F^{-1})，则确实有

$$
\delta D_F=-D_F^{(0)}R D_F^{(0)},
\qquad R=-\Pi_{\rm user}.
$$

因此“负的 D0 Π D0”可以是一种正确记法，但其中 Π 必须是 R，不能同时取用户给出的 +i 标量。两个 ig 发射顶点已在 Wilson 展开中给出其自身的负号；闭合夸克/鬼环的符号也已在 a 中，均没有理由再整体乘 −1。

结论保持：用户两式作为指定数学积分的定义完全确定，既有 montopar 仍在求这个对象；若要求它等于标准 Wilson 圈的该图贡献，就必须转换整体号。更准确的诊断是**混用了 1PI 插入与逆传播子修正的符号约定**，不是简单地说“自能图永远不能有外负号”。本次未改任何 I 或生产源，亦未靠 RG 数值决定该号。

## 4. 完整张量与六边对相消

此处的负号与是否保留导数张量无关。相邻图的 1/(1+ε) 是完整张量的几何结果，仍须保留。对闭合轮廓，先保留共同端点值 H_ij，再作双导数积分的端点求和：

$$
\sum_{i,j=1}^4
(H_{i+1,j+1}-H_{i+1,j}-H_{i,j+1}+H_{ij})=0,
\qquad H_{i+4,j}=H_{ij}=H_{i,j+4}.
$$

恒等式在取光锥极限之前成立。不能先分别把奇异端点写成 0×∞。本轮还用已保存的 BA、BO 全ε解作了连接检查。令 G=e^(−γE ε)Γ(1−ε)，相邻完整值与只含度规项的辅助值之差为

$$
I_{B_A}-I_{B_A,g}
=-\frac{G^2}{32\varepsilon^2(1-2\varepsilon)(1+\varepsilon)},
\qquad I_{B_A,g}=\frac{G^2}{32\varepsilon^3(1-2\varepsilon)}.
$$

BO 中的额外端点项正好是其相反数；剩下的度规辅助值为

$$
I_{B_O,g}=\frac{G^2 H_0(2\varepsilon,z)}
{16\varepsilon(1-2\varepsilon)\cosh(\varepsilon\log z)}.
$$

因此已包含 m=4、2 的两个完整 I 满足

$$
I_{B_A}+I_{B_O}=I_{B_A,g}+I_{B_O,g}.
$$

这里只用辅助分解核验同一个完整 I 的闭合轮廓求和；不新增另一套待求 I，也不逐图删除张量项。

## 5. 一次 Euler、尺度和重数

除去 λ0²、颜色和 π^(2ε)(−t+i0)^(2ε)，三类 prepared 前因子分别是

$$
\begin{array}{c|c}
4e&-G^2/8\\
3e\text{-}3g&(1-2\varepsilon)G_2/4\\
2e\text{-Pi1（当前源）}&G^2/[16\varepsilon(1-2\varepsilon)]
\end{array},\qquad
G_2=e^{-2\gamma_E\varepsilon}\Gamma(1-2\varepsilon).
$$

颜色依次提出 CFCA、CFCA、CF B1；Pi1 的物理乘数仍单列为 −1。TGV 的色迹 iCFCA/2 与 C 中 I³ 相乘为 CFCA/2，导数因子为 2(1−2ε)，与当前正前因子一致。完整 P 的符号由同一 Fourier 导数方向重建；不能只看色迹中的 i 判断正负。

**G2 不等于 G²**，虽然二者 Euler 补偿次数相同。其比值从 ε² 就不同，故组装时各自保留。

沿用一圈字典

$$
\lambda_0 S_\varepsilon=\mu_{\overline{\rm MS}}^{2\varepsilon}
Z_\alpha\lambda,\quad S_\varepsilon=(4\pi)^\varepsilon e^{-\gamma_E\varepsilon},
\quad\mu_x^2=\frac{e^{2\gamma_E}}4\mu_{\overline{\rm MS}}^2,
$$

并令 Xs=μx²(−s+i0)、Xt=μx²(−t+i0)。三类 raw C 都具有同一个 π^(2ε) 尺度幂；恢复原 Euler、作裸耦合替换再与 prepared 比较，均得到相同的 Xs、Xt。这里暂不乘 Zα，普通反项留给下一阶段。

对任意轨道，旧 I 满足

$$
\frac m2 C\,[J(z)+z^{2\varepsilon}J(1/z)]
=\mathcal N_2(S_{14}+S_{12})I.
$$

右边已经包含 m/2、crossing 和颜色投影。最终求和各轨道一次，包括 I01=I02 的两张 TGV；不再乘 m，也不再补 Euler。

## 6. 可交给下一步的组装式

在 λ=αs/π 的展开中，hat w2 表示已恢复 MS-bar 尺度、尚未加入普通反项的 λ² 系数。现有 I 原样代入：

$$
\boxed{
\widehat w_2=C_F(X_s^{2\varepsilon}+X_t^{2\varepsilon})
\left\{C_A\left[\sum_{r=1}^6 I_{4e-r}
+\sum_{r=1}^4 I_{3e-3g-r}\right]
-\mathcal B_1(\varepsilon)(I_{B_A}+I_{B_O})\right\}.}
$$

与旧计划候选式的差别是 Pi1 前面的物理转换号，不是改变其 I。该式为本次审计得到的连接公式；本轮尚未计算它的完整 Laurent 总和或有限部分。

为衡量原号错误的后果，只作了最高极点诊断：交换图与 TGV 的 ε⁻⁴ 相消；在共同尺度 L=0 的 ε⁻³ 检查中，旧正号候选加上标准耦合反项后，相对 RG 必需的 3CFβ0/8 留下残差

$$
\Delta c_{-3}=C_F\left(\frac5{12}C_A-\frac13T_Fn_f\right).
$$

使用由 Dyson 展开独立得到的 −1 桥后，这一残差为零。这只是辅助诊断，**不是用已知 RG 结果决定符号，也不代表低阶极点和有限 RG 已匹配**。

展开预算已足够：12 个 I 都有有限项；B1 必须卷积到 ε³；一圈组装保存至 ε²，覆盖反项所需的 ε¹。接下来执行完整和、加入 −β0 hat w1/(4ε)，再做严格极点减除与双路线 RG/cusp 提取。现有来源修复、物理组装和形式验证三种状态分别记录。
