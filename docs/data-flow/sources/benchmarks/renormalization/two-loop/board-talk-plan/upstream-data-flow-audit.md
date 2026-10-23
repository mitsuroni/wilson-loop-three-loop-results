# 两圈前置数据流：图、动量表达式、C/J 与 I

2026-10-17，WILSON-5.6。先核对项目事实，不改 PDF。旧板书 v1 把四段上游解释压得过短；相关资料实际存在，排版交付不等于已经讲清这条推导链。

本轮实际重放冻结的组合生成/分类/导出 API，与旧冻结的十二条结构化定义完全一致。独立枚举、当前 prepared 连接、自能版本绑定与三类前因子关系共 110 项检查通过。证据见 [audit.json](../../../../workflow/checks/two-loop-upstream-flow-20261017/audit.json)；不是新增多维积分或形式证明。

## 1. 先区分每一层的对象

```text
Wilson observable + g0^4 order + finite oriented contour
 -> allowed attachments/pairings
 -> zero/scaleless removal + exponent-reducible product separation
 -> C4 orbit representatives, multiplicities, ordered endpoints
 -> complete momentum-space contribution A_r
 -> Fourier / common-center / bubble integration
 -> A_r = C_raw,r J_r on the original ordered parameter domain
 -> explicit color/sign conversion + one Euler preparation
 -> A_prep,r = C_r J_r
 -> cyclic sum and scale extraction DEFINE the normalized I_r
 -> EVALUATE J_r to obtain exact or Laurent expressions for I_r
 -> full color/scale sum -> bare w2 -> renormalization -> F2
```

三件事必须分开：

- 旧图式把原始贡献标成 `I_(...)`，不是下游归一化的 `I_r`。讲稿宜把前者统一记为 A_r。
- C_raw、prepared C、轨道归一化 I 属于不同层。加入 Euler 补偿后应明确写 A_prep，不能无说明地改变原振幅。
- J→I 的归一化公式只定义对象，不是积分求解。得到生成函数与 Laurent 系数还需要真正计算 J。

## 2. 为什么是 6+4+2

### 2.1 从耦合阶数开始

每个有限 Wilson 边插入给一个 g0。在 g0^4 阶：

| Wilson 插入数 | 内部结构 | 图族 |
|---:|---|---|
| 4 | 两条裸胶子交换 | 4e |
| 3 | 树级三胶子顶点，再给 g0 | 3e-3g |
| 2 | 一圈两点插入，再给 g0^2 | 2e-Pi1 |

真空泡由真空归一化移去；无质量四胶子 tadpole 无尺度。Pi1 内的胶子、鬼、夸克泡已合成一个 D 维算符，不额外增加轮廓轨道。

实际代码是 `qgraf_pipeline.generator.generate_diagrams` → `classifier.classify_orbits` → `definition_exporter`。该组合生成路径替代 QGRAF 枚举；本轮与旧两圈审计都没有运行 QGRAF 二进制。不能讲成“QGRAF 直接输出十二张图”。

### 2.2 先说明为什么算 log W

四个有序轮廓点的交错配对有

\[
c_\times=\frac1{N_c}\sum_{a,b}\operatorname{Tr}(T^bT^aT^bT^a)
=C_F^2-\frac12 C_FC_A.
\]

另外两种非交错配对只有 C_F^2。全部阿贝尔部分组成 W1^2/2，在 log W 中消去，因此交换类保留指数颜色 -C_F C_A/2。这项两圈论证不能用一个没有解释的高圈 2PI 标签代替；选好指数颜色后不再扣 W1^2/2。

### 2.3 独立带边枚举，然后才做 C4 商

- 4e：四个点的弱递增边号有 35 种；交错配对为 (1,3)、(2,4)，排除 e1=e3 或 e2=e4 的同边裸线零项，余 19。
- 3e-3g：三个点的弱递增边号有 20 种；删去三个端点全在同一类光边的 4 个零项，余 16。
- Pi1：两条不同边的选择有 C(4,2)=6 种；同边完整张量贡献为无尺度项。

只取循环旋转 C4 商。镜像只记录关系，不再合并：

| ID | 沿轮廓的附着边 | m |
|---|---|---:|
| 4e-01 | 3,3,4,4 | 4 |
| 4e-02 | 2,3,4,4 | 4 |
| 4e-03 | 2,3,3,4 | 4 |
| 4e-04 | 2,2,4,4 | 2 |
| 4e-05 | 2,2,3,4 | 4 |
| 4e-06 | 1,2,3,4 | 1 |
| 3e-3g-01 | 3,4,4 | 4 |
| 3e-3g-02 | 3,3,4 | 4 |
| 3e-3g-03 | 2,4,4 | 4 |
| 3e-3g-04 | 2,3,4 | 4 |
| 2e-Pi1-01 | 3,4 | 4 |
| 2e-Pi1-02 | 2,4 | 2 |

带边配置为 19+16+6=41，代表为 6+4+2=12。35/20 不是全部 QCD 图数。源码先规范化部分边分布，本轮 generator 的 19 个 synthetic 对象也不是上述 41 个带边配置。

解相等不删轨道：4e-02/05、3e-3g-01/02 均各自保留。内部的 `tree` 标签只指子图为树级，不表示整个 Wilson 贡献是一圈。

依据：[冻结定义](../../../../workflow/checks/rg-two-loop-audit/sources/upstream/orbit_definitions.json)、[旧独立枚举](../../../../workflow/checks/rg-two-loop-audit/check.py)、[本轮重放](../../../../workflow/checks/two-loop-upstream-flow-20261017/audit.py)。

## 3. 费曼规则怎样给出动量空间起点

固定 X_i(s)=x_i+s v_i、v_i^2=0、Fourier 相位 exp(-ikX)、测度 d^Dk/(2pi)^D、D=4-2eps。先列规则，再写三类图：

1. Wilson 插入：ig0 v_i^mu T^a ds。
2. 胶子传播子：-i g_mu,nu delta_ab/(k^2+i0)，连同端点差相位。
3. 三胶子顶点：全流入约定的 g0 f^{abc} K_mu,nu,rho，K 为三项循环动量差。
4. 自能插入：D0 Pi1 D0，Pi1 已含 1PI 相位及内部物种系数。
5. 归一化反序颜色迹与原有序域；既已排序，不再乘 1/E!。

本项目保留有限边 ds 与相位，不能直接拿无限直线的 1/(v.k) 替代。沿轮廓早到晚编号，较晚的生成元在迹的左边。

交换类交错图的完整起点为

\[
\mathcal A_r^E=(ig_0)^4c_\times
\int_{\mathcal D_r}d^4s\,(u_1\!\cdot u_3)(u_2\!\cdot u_4)
\int_k\int_p\frac{(-i)^2e^{-ik(X_1-X_3)-ip(X_2-X_4)}}
{(k^2+i0)(p^2+i0)}.
\]

顶点图让 k_a 从中心流向 X_a，流入顶点的动量为 -k_a：

\[
\begin{aligned}
\mathcal A_r^V={}&(ig_0)^3g_0\int_{\mathcal D_r}d^3s\,
\frac{f^{abc}}{N_c}\operatorname{Tr}(T^cT^bT^a)\\
&\times\mathsf R_r\left[
\int_{k_1,k_2,k_3}(2\pi)^D\delta^{(D)}(\sum k_a)
\frac{(-i)^3u_1^\mu u_2^\nu u_3^\rho
K_{\mu\nu\rho}(-k_1,-k_2,-k_3)e^{-i\sum k_aX_a}}
{\prod_a(k_a^2+i0)}\right].
\end{aligned}
\]

R_r 表示先对独立位置求导，再限制到路径。颜色收缩为 -i C_A C_F/2，动量分子在同一指数上变成 -i mathscr K。

自能的起点为

\[
\mathcal A^\Pi_{ij}=\int d\xi d\eta\,
\frac{\operatorname{Tr}[(ig_0v_i^\mu T^a)(ig_0v_j^\rho T^b)]}{N_c}
\int_k e^{-ikX}[D_0\Pi_1D_0]^{ab}_{\mu\rho}(k).
\]

当前物理插入外部系数为 +1。旧冻结源额外写 -1；现行 r002 的 C、I 已在导入边界翻号，J 保持，不能再补负号。采用入口是 `benchmarks/self-energy/current-r002.json`。

完整三类链见 [momentum-to-parameters.tex](../report/momentum-to-parameters.tex)，顶点的逐步规则与颜色推导见 [momentum.tex](../../../tgv/3e-3g/lecture-notes/momentum.tex)。旧 `integrals.md` 的 I 标签、TGV 相位印刷括号及旧自能源号需要解释；新讲稿不能只照抄其渲染文字。

## 4. 从动量空间到 C_raw,J：三类不是同一个操作

设 A0=-t=-2v14>0。共同单传播子变换：

\[
\Delta_0(X)=-\frac{\Gamma(1-\epsilon)}{4\pi^{2-\epsilon}}
(-X^2+i0)^{-1+\epsilon}.
\]

### 4.1 两次 Fourier 变换：交换类

sigma_ab=-(X_a-X_b)^2+i0，f1=sigma13/A0，f2=sigma24/A0，
n=4(u1.u3)(u2.u4)/A0^2。保持同一振幅：

\[
\mathcal A_r^E=
\underbrace{\frac{g_0^4c_\times\Gamma(1-\epsilon)^2A_0^{2\epsilon}}
{64\pi^{4-2\epsilon}}}_{C_{\rm raw}^E}
\underbrace{\int_{\mathcal D_r}d^4s\,n_r f_1^{\epsilon-1}f_2^{\epsilon-1}}_{J_r}.
\]

不加 simplex，参数就是原路径位置。例如 4e-01：
n=z^2，f1=z(1-s31)s41，f2=z(1-s32)s42，界限为 s31<s32、s41<s42。

### 4.2 公共中心、单纯形、高斯积分与导数：顶点类

依次用守恒 delta 引入 Y、变换三条传播子、用 alpha 单纯形合并分母、完成 Y 的平方并积分、作用 mathscr K。定义

\[
Q_\star=-\sum_{a<b}\alpha_a\alpha_b(X_a-X_b)^2,\quad
S_r=\frac{\mathsf R_rQ_\star}{A_0},\quad
P_r=\frac{\mathsf R_r(\mathscr KQ_\star)}{2A_0^2}.
\]

由此

\[
C_{\rm raw}^V=\frac{g_0^4C_AC_F\Gamma(2-2\epsilon)A_0^{2\epsilon}}
{64\pi^{4-2\epsilon}},\quad
J_r=\int_{\mathcal D_r}d^3s\int d\mu_\epsilon(\alpha)
\frac{P_r}{S_r^{2-2\epsilon}},
\]

dmu=d^3alpha delta(1-sum alpha)(alpha1 alpha2 alpha3)^(-eps)。
例如 3e-3g-01 的 S=az(1-t)(bx+cy)，P=az^2(1-t)(b-c)/4，0<x<y<1。
Gamma(2-2eps)、分子和分母幂来自同一微分，不能套交换类 Gamma 因子。

### 4.3 泡积分、横向张量与再 Fourier：自能类

内部泡给 Gamma(eps) Gamma(1-eps)^2/Gamma(2-2eps)，全颜色系数为

\[
\mathcal B_1=\frac{10-6\epsilon}{3-2\epsilon}C_A
-\frac{8(1-\epsilon)}{3-2\epsilon}T_Fn_f.
\]

当前物理号给出

\[
C^\Pi_{\rm raw}=-\frac{g_0^4C_F\mathcal B_1\Gamma(1-\epsilon)^2A_0^{2\epsilon}}
{256\pi^{4-2\epsilon}\epsilon(1-2\epsilon)},\quad
J_{ij}=\int d\xi d\eta\left[d f^{2\epsilon-1}
+\frac{\partial_\xi\partial_\eta f^{2\epsilon}}{4\epsilon(1+\epsilon)}\right].
\]

等价显式核为 [(1+2eps)d f+(1-2eps)u w] f^(2eps-2)/[2(1+eps)]。
f=-X^2/A0、d=-2vi.vj/A0、u=-2vi.X/A0、w=-2vj.X/A0，幂保留共同 i0。
混合导数需积出四角点，不能直接删“纵向项”。

### 4.4 实际文件和准备边界

| 族 | 原核与原因子 | 下游 prepared |
|---|---|---|
| 4e | `benchmarks/2ne/sources/montopar-ordered/4e/montopar-initial.json` | `benchmarks/2ne/prepared/4e.json` |
| 3e-3g | `benchmarks/tgv/3e-3g/sources/canonical.json` 的 original_C、J | `benchmarks/tgv/3e-3g/prepared/3e-3g.json` |
| Pi1 | current-r002 选择的 sources/2e-Pi1-r002.json | 同一选择器绑定的 prepared/2e-Pi1-r002.json |

4e prepared C 含固定基底颜色投影和一次 exp(-2 gamma_E eps)，完整颜色坐标仍保存。
有序 s 与历史 unit_cube 带各自界限/Jacobian，不能混拼。
IBP 的 sum C_t J_t 是同一图的分解；combined_C 已含全局 C，不能再乘一次。

## 5. C,J 到 I：定义与求值分开

交换图先定义投影后的 A_E,exp（c_cross -> -CFCA/2）。随后

\[
C_r=e^{-2\gamma_E\epsilon}C_{{\rm raw},r},\qquad
\mathcal A_r^{\rm prep}=C_rJ_r.
\]

### 5.1 先把旋转的图逐个写出来

这里的“旋转 90 度”指沿四边形循环移动一条边的附着位置：
\(1\to2\to3\to4\to1\)，连同动量、端点参数及原有序积分域一起重标记。
它不要求一般的 null 四边形具有欧氏正方形的刚性旋转对称。
取

\[
s=2v_1\!\cdot v_2=2v_3\!\cdot v_4,\qquad
t=2v_1\!\cdot v_4=2v_2\!\cdot v_3,\qquad z=s/t.
\]

循环移动一次使 \(s\leftrightarrow t\)，所以 \(z\to1/z\)；移动两次后运动学变量恢复。
这并不表示两次旋转后的**带边图配置**一定与原图相同。

为解释尺度权重，先保留一般圈数 \(n\)，本节实际应用为 \(n=2\)。
在 \(-s,-t>0\) 的共同起始区域定义

\[
S_t=[\pi(-t)]^{n\epsilon},\qquad
S_s=[\pi(-s)]^{n\epsilon}=z^{n\epsilon}S_t.
\]

将不含这两个运动学尺度的公共因子记为 \(K_r\)，代表图写成
\(\mathcal A_r=K_rS_tJ_r(\epsilon,z)\)。
旋转一次后的图自然以 \(s\) 为参考尺度，因此

\[
\mathcal A_{Rr}=K_rS_sJ_r(\epsilon,1/z)
=K_rS_t\,z^{n\epsilon}J_r(\epsilon,1/z).
\]

所以 \(z^{n\epsilon}\) 来自**旋转后参考尺度的转换**；仅代入 \(z\to1/z\) 会漏掉它。
后续解析延拓需沿用共同支路，不能把此起始区域的幂律身份当作任意复数幂的无条件乘法。
记 \(a_r=J_r(\epsilon,z)\)、\(b_r=z^{n\epsilon}J_r(\epsilon,1/z)\)，则

| 循环移动 | 图像 | 运动学变量 | 除以公共因子 \(K_rS_t\) 后的贡献 |
|---|---|---|---|
| \(0^\circ\) | \(r\) | \((s,t),z\) | \(a_r\) |
| \(90^\circ\) | \(Rr\) | \((t,s),1/z\) | \(b_r\) |
| \(180^\circ\) | \(R^2r\) | \((s,t),z\) | \(a_r\) |
| \(270^\circ\) | \(R^3r\) | \((t,s),1/z\) | \(b_r\) |

\(360^\circ\) 回到原图，不是第五个贡献。表中四个位置先按群作用列全，再去除相同的图配置；
不同图的积分表达可以相同，不能据此删图。

### 5.2 为什么系数是 \(m_r/2\)

\(m_r\) 是该代表图的**不同循环像的个数**。若保持图配置不变的旋转有
\(|H_r|\) 个，则 \(m_r|H_r|=4\)。
本项目的识别包含允许的内部哑标签重命名，不包含镜像商。
四次群作用会将每个不同图重复列出 \(|H_r|\) 次，所以应求和

\[
\begin{aligned}
\sum_{\text{distinct cyclic images}}\mathcal A
&=\frac{1}{|H_r|}\sum_{k=0}^{3}\mathcal A_{R^kr}\\
&=\frac{m_r}{4}K_rS_t(2a_r+2b_r)\\
&=\boxed{\frac{m_r}{2}C_r
\left[J_r(\epsilon,z)+z^{n\epsilon}J_r(\epsilon,1/z)\right]},
\qquad C_r=K_rS_t.
\end{aligned}
\]

此处 \(m_r\) 也等于首次回到同一图所需的四分之一周旋转次数，具体为：

| \(m_r\) | 不同的图像与首次返回 | 去重后的和，除以 \(C_r\) | 当前两圈例子 |
|---|---|---|---|
| 4 | \(r,Rr,R^2r,R^3r\)；\(R^4r=r\) | \(a_r+b_r+a_r+b_r=2(a_r+b_r)\) | 4e-01 |
| 2 | \(r,Rr\)；\(R^2r=r\) | \(a_r+b_r\) | 4e-04、2e-Pi1-02 |
| 1 | 只有 \(r\)；\(Rr=r\) | \(a_r=\tfrac12(a_r+b_r)\)，因 \(a_r=b_r\) | 4e-06 |

\(m_r=1\) 时，同一图的旋转等价性要求
\(J_r(\epsilon,z)=z^{n\epsilon}J_r(\epsilon,1/z)\)；
前面的 \(1/2\) 正好消除同一图的重复，并未只计半张图。
其他轨道也可能在特殊运动学或求值后出现相等函数，\(m_r\) 仍按图配置计算。

**写这个式子的目的：**只需计算每个循环轨道的一张代表图，就能恢复该轨道的全部贡献；
式中的两项保留了 \(s\) 与 \(t\) 两种参考尺度，\(m_r/2\) 保证每张不同的图只计一次。
因此应先展示旋转图和上表，再把逐图和缩写成这个公式。

英文板书口稿可接在四张旋转图后：

> A quarter-turn exchanges the two kinematic scales. The rotated graph therefore
> contributes \(z^{n\epsilon}J_r(\epsilon,1/z)\) when we express it using the original
> reference scale. We now sum the distinct cyclic images. A four-member orbit gives
> two copies of each term, whereas a two-member orbit gives one copy of each.

### 5.3 再提出共同尺度，定义唯一的 I

回到两圈，令 lambda_b=g0^2/(4pi^2)、S14=(pi A0)^(2eps)、S12=z^(2eps)S14；
它们分别是上述 \(S_t,S_s\) 的 \(n=2\) 版本，\(A_0=-t\)。
取 h_r=C_A（交换/顶点）或 B1（自能），则

\[
\boxed{\frac{m_r}{2}C_r[J_r(\epsilon,z)+z^{2\epsilon}J_r(\epsilon,1/z)]
=\lambda_b^2C_Fh_r(S_{14}+S_{12})I_r.}
\]

左端已经是不同旋转图的和；右端再提出颜色、耦合及对称的外尺度。
其中

\[
S_{14}+S_{12}=S_{14}(1+z^{2\epsilon})
\]

解释了 \(I_r\) 定义里为什么会出现分母 \(1+z^{2\epsilon}\)：
它来自共同尺度的提取，不是额外的图平均，也不是再乘一次旋转重数。
将两项以这个尺度归一化，还使留下的形状组合在 \(z\leftrightarrow1/z\) 下显式不变。
Z=exp(-gamma_E eps) Gamma(1-eps)。记 B[J]=[J(z)+z^(2eps)J(1/z)]/(1+z^(2eps))：

| 族 | I 的定义 |
|---|---|
| Exchange | -m Z^2 B[J]/16 |
| Vertex | m exp(-2 gamma_E eps) Gamma(2-2eps) B[J]/8，本族 m=4 |
| Pi1 | -m Z^2 B[J]/[32 eps(1-2eps)] |

所以得到 I 后，不再乘 C、Euler 或轨道倍数。

### 5.4 归一化之后，还要实际求积分

**实际求值**才给具体的生成解/系数。4e-01 可以从头走完：

\[
J_{01}=z^{2\epsilon}\int_0^1ds_{32}\int_0^1ds_{42}
\int_0^{s_{32}}ds_{31}\int_0^{s_{42}}ds_{41}
[(1-s_{31})(1-s_{32})s_{41}s_{42}]^{\epsilon-1}
=\frac{z^{2\epsilon}}{4\epsilon^4}.
\]
\[
I_{01}=-\frac{Z^2}{4}
\frac{(z^{2\epsilon}+1)/(4\epsilon^4)}{1+z^{2\epsilon}}
=-\frac{Z^2}{16\epsilon^4}.
\]

其他图依次需要单核一维化、Beta/simplex 积分、路径分部积分与接触面、对数矩母源，或系数级 Hlog 积分。3e-3g-04 的已知 eps^2 系数与紧凑全 epsilon 表示的问题分开。

旧结果 JSON 不一定单文件列齐正 eps 阶。例：4e-01/r002 本体到有限项，精确生成式与后续研究支持高阶展开，由报告构建器汇总。TGV04 正阶还消费 `workflow/checks/tgv04-high-order/runs/positive-epsilon/I1.txt`、`I2.txt`。报告有 eps^2 数据不等于每个旧 r001/r002 都存了这些字段。

## 6. 下一稿需要怎样改变

先推得十二个轨道，再列费曼规则与三类完整动量起点，随后每类保持同一 A 的连续等式链得到 C_raw,J。
prepared 与 I 另作明确接口，并把 4e-01 的 J 求值和 I 归一化接在同一算例里。
之后才进入其他积分、全表、物理组装与重整化。v1 后半段可以复用，上游需按这条依赖重新写；本轮没有重写旧 PDF。
