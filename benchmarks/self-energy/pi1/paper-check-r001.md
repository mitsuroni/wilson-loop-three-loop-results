# 相邻 Π₁ 图与论文式的精确比较

日期：2026-09-23。对象：已发布的 `2e-Pi1-01`，原 s 域、完整颜色、raw 归一化。结论：**完整核相对于度规核的 1/(1+ε) 正确；当前 r001 的绝对符号与用户给出的论文式相反。** 因此本次论文控制不能记为通过。原始 JSON、生产代码及 r001 理论文档保持原字节，本报告记录新增审计结论。

## 1. 比较对象和术语

使用用户本轮明确给出的公式，记作 P：

$$
P=\frac18\left(\frac{\alpha_s}{\pi}\right)^2 C_F
(-2\pi v+i0)^{4-D}
\frac{\Gamma^2(D/2-1)[(3D-2)C_A-2(D-2)N_f]}
{(D-4)^3(D-1)(D-3)}.
$$

本报告的直接比较取 v=v34=v₃·v₄、D=4−2ε、TF=1/2、Nf=nf，同一相邻代表，无 C4 求和、无 Euler 补偿。先在 −2v14>0、z=v34/v14>0 的起始区域比较，再按同一个 Feynman 处方共同延拓；没有把复杂幂作不受条件约束的拆分。

完整自能张量 `(ℓ²g−ℓℓ)Π_scalar` 本身是横向的。本报告把 `ℓ²g` 称作度规项，把 `−ℓℓ` 所产生的坐标导数称作导数项（或纵向动量项）。这里讨论的是保留后者产生的单图修正，不能说删去了完整的横向算符。

通过 INSPIRE `open_search` 核实的文献是 I.A. Korchemskaya、G.P. Korchemsky，*On lightlike Wilson loops*，Phys. Lett. B 287 (1992) 169–175，[DOI](https://doi.org/10.1016/0370-2693(92)91895-G)。本地保存的[原文 OCR](/C:/Users/Lenovo/Documents/learning/gbrain-master/brain-host/sources/korchemsky-1992-ll-wl/full.md:121)在第3节将同样的 D 依赖结构列作 W_b，原式尺度记为 L；用户本轮给出的 `2πv` 用作本次直接对照的尺度。在线出版社全文返回403，未视觉核对原 PDF；也没有在所读原文中找到作者明确声称“逐图抛弃某张量项”的句子。论文度规/约化代表的解释与原作者的明确表述应区分。

## 2. 从实际生产 J 得到精确积分

本次不是重新抄写理论里的 C/J：检查器直接加载外部 montopar 当前源码生成结果，核对实际 workspace JSON 与 Wilson 冻结包一致，再读取 `tensor_decomposition` 的两项与完整 `integrand`。

令 x=1−s31、t=s41。相邻线有

$$
f=zxt,\qquad d=z,\qquad uw=zf.
$$

因此在开域 x,t,z>0，

$$
\mathcal S_g=z^{2\varepsilon}x^{-1+2\varepsilon}t^{-1+2\varepsilon},
\qquad \mathcal S_{\partial}=-\frac{\varepsilon}{1+\varepsilon}\mathcal S_g,
\qquad \mathcal S_{\rm full}=\frac{\mathcal S_g}{1+\varepsilon}.
$$

Re ε>0 时两个端点积分均绝对收敛，每一个等于 1/(2ε)。原 s31 外、s41 内的次序保持；本例是实际执行两次一维积分，不需要假设一个未知的换序桥。得到

$$
J_g=\frac{z^{2\varepsilon}}{4\varepsilon^2},\qquad
J_{\partial}=-\frac{z^{2\varepsilon}}{4\varepsilon(1+\varepsilon)},\qquad
J_{\rm full}=\frac{z^{2\varepsilon}}{4\varepsilon^2(1+\varepsilon)}.
$$

这些式子按亚纯函数继续到 ε=0 邻域。修正不只在积分后成立：本相邻例的核在积分前已具有相同比值。当前两点核可由此确认；不能据此把所有对边或四端点图都替换为相同比值。

## 3. 乘上当前 C：存在精确的整体负号差

统一记

$$
B(\varepsilon)=\frac{(10-6\varepsilon)C_A-4(1-\varepsilon)n_f}{3-2\varepsilon},
$$

$$
K=\frac1{64}\left(\frac{\alpha_s}{\pi}\right)^2 C_F B(\varepsilon)
\frac{\Gamma^2(1-\varepsilon)}{\varepsilon^3(1-2\varepsilon)}
(-2\pi v_{34}+i0)^{2\varepsilon}.
$$

当前生产 C 乘上上节的实际 J，精确给出

$$
(CJ)_{r001,g}=K,\qquad (CJ)_{r001,\partial}=-\frac{\varepsilon}{1+\varepsilon}K,
\qquad (CJ)_{r001,\rm full}=\frac{K}{1+\varepsilon}.
$$

但是用户给出的论文式满足

$$
(D-4)^3=(-2\varepsilon)^3=-8\varepsilon^3,\qquad D-3=1-2\varepsilon,
\qquad P=-K.
$$

故本次实际验收比值为

$$
\boxed{\frac{(CJ)_{r001,g}}P=-1,\qquad
\frac{(CJ)_{r001,\rm full}}P=-\frac1{1+\varepsilon}.}
$$

Γ函数、圈阶、味数、尺度幂和数值系数均匹配；差异是一个整体负号。它不会被 C4 重数、相位向量反向、Euler 补偿或六边对导数相消修复。

## 4. 符号来源与应有的约定转换

当前生产实现忠实采用此前用户／QGRAF给出的起点：最外层负号，Wilson 顶点 ig，两侧传播子 −i，Π_scalar=+iAq^(−ε)。度规项再用 ℓ²=−q 和 Fourier 相位 i。其相位链为

$$
(-1)_{\rm outer}\,i^2\,(-i)^2\,i\,(-1)_{\ell^2=-q}\,i=-1.
$$

再用 v14/(−2v14)=−1/2，给出当前 C₂ 的正号。上游 `qgraf-pipeline/src/qgraf_pipeline/internal_types.py:get_vertex_info` 对 N=2 固定返回 `sign=-1`；当前 montopar 的 C 账本则同时记录 `source_external_sign=-1` 和 `Pi1_scalar_phase=I`。这就是两个定义接合处需要统一的符号。

若目标就是与 P 同定义的物理 Wilson 图贡献，其完整结果应为

$$
\boxed{(CJ)_{\rm paper,full}=\frac{P}{1+\varepsilon}
=\frac{2P}{6-D}=-\,(CJ)_{r001,\rm full}.}
$$

对应的必要转换是 **C_paper=−C_r001，J 不变**。在动量输入层，可以取消当前独立的最外负号，或在保留它时重新定义插入符号 Π；必须同步写清两种符号代表的是哪个对象，不能同时改两处而把负号又抵消。若保留现有含外负号的 I 定义，则应明确它与论文 W 贡献相差负号，不能直接说两者数值相同。

本轮任务是核验。这里没有静默更改用户此前的动量起点、十三图 C、上游图生成器或不可变 r001。因为它们共享同一单插入约定，后续若统一到论文 W，应整体审阅13图的公共 C 和上游 Π 约定，不能只对01加一个特例负号。

## 5. 对上轮验证结论的限定

上轮测试确认的是实现与所选动量起点之间的代数连接，不能代替与论文的绝对符号对照。此前讨论总结的第2.2及5.1–5.5节已明确保留这项缺口；本次以用户给出的明确尺度和实际发布数据，得到确定的 `−1` 比值，不能将文献验收标为通过。

完整核/度规核=1/(1+ε) 及二点六边对导数相消仍成立；本次发现的是公共符号约定问题，不要求撤销这些相对恒等式，也不对全部十三图授予文献验证状态。

本次得到相邻控制图的精确 raw C×J，未展开 Laurent 系数或求其余图。[检查器及19项复核](../../../workflow/checks/self-energy-pi1-paper-check/README.md)。状态：`paper_comparison=failed_overall_sign`、`formal=not_run`。
