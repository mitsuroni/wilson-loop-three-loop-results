2026-10-06共同规则已融入正文：§1保留顶点标记、使用真实晚到早Wilson迹，§5直接生成负C，§10据该C提取物理I；内部canonical分子、两种J、域和IBP保持。当前原生包为 `montopar-r002.json`，I为current选择的r002，原r001数据保留。

# Gamma1 Montopar production contract

2026-10-03. Schema **1.8**, dispatch `gamma1`, family `3e-Gamma1`.
The full theoretical derivation is reproduced below, with links to its fixed
Wilson-workspace evidence. Source, parameterization, and evaluation are separate:
this implementation generates C/J; it does not evaluate I.

```powershell
montopar --source ../qgraf-pipeline/outputs/3loop/orbit_definitions.json --type 3e-Gamma1 --gamma1-ibp --package-initial
```

`--gamma1-ibp` is an alias of `--tgv-ibp`. A single source representative is
written for each of the four diagrams. No hypercube or canonical-leg relabeling
changes the input contour order. `gamma1_source.py` validates the complete
momentum contract and pins the embedded DOT operator hash; malformed Gamma
sources cannot enter the pure/tree branches.

## Output and consumption

- `C.combined` is the common factor including CF CA. `C.channels.g/q.combined`
  adds CA or TF nf respectively. Every scalar term stores its full `combined_C`.
- `J_raw.terms` contains seven channel/block integrals: N0 and L0 in each
  colour, plus three gluonic bubble blocks. The two triangle simplex measures
  give seven integration variables; a bubble has five. Zero algebraic blocks
  remain defined in raw; they are not new physical diagrams.
- `J_sector.terms` gives the equivalent fifteen integrals with explicit
  u,v >= 1 triangle sectors. All t substitutions and Jacobians are included.
  The scalar expressions contain no Max. The three partitions carry no new
  graph multiplicity.
- `J_raw_s` preserves the complete ordered raw structure. `J_raw.operators`
  preserves the underlying differential-operator structure, while each term's
  integrand has already applied all ambient derivatives. Numerators are exact,
  self-contained factored polynomials; no G/T/W or other undefined abbreviation
  remains in them.
- Optional `J_ibp.terms` contains the face/recoil integrals and
  `J_ibp.sector_terms` provides their no-Max sector versions. Each term has its
  own channel, C, variables, bounds, order, face provenance and common regulator
  group. 344/334/244 each span three faces; 234 spans six faces and a bulk
  remainder. Contacts retain both normal contributions. Derivatives are applied
  before restricting coordinates to a face. The original s order is inherited.
- Output validation reconstructs all data from the embedded source before
  packaging; dropping a colour, changing a kernel, order, branch, or coefficient
  is rejected. Type packages use schema 1.8. Other families retain their schemas.

The amplitude is `sum(term.combined_C * integral(term.integrand))`.
The legacy Stage 3 consumer of a single `J_raw.integrand` is insufficient for
this two-colour multi-measure object. A downstream solver must explicitly
iterate the terms and preserve their measures; no single-integrand shortcut
is exported.

## Normalization and checks

Raw output includes neither exp(-3 EulerGamma eps), C4 multiplicity, colour a1
projection, crossing assembly, external-leg Pi1 nor a counterterm. The explicit
coupling dictionary is g0=mu^eps g_b and alpha_s=g_b^2/(4 Pi), with mu counted
once. The DOT loop normalization is already represented by the Gram-free
integral and is not multiplied a second time.

The runtime uses only Python standard-library exact Fraction arithmetic. CAS
and numerical libraries appear only in independent checks. Source validation,
component-level ambient derivatives, finite-I0 IBP reconstruction, C4 geometry,
explicit sector Jacobians, no-site-packages execution, real CLI generation, and
package tampering are checked separately. No Lean/MCP proof status is implied.

## Full theory (r002 mathematical specification)

The following r002 text predates this production adapter. Its mathematical
formulae are unchanged; its statements that production is pending describe
that earlier checkpoint, not the current schema 1.8 capability.

# Γ₁：从完整动量源到 C、J 的参数化理论

2026-10-03，r002。本文完成中央一圈三胶子顶点 Γ₁ 的参数化推导：从[四条固定动量源](D:/wilson-agent/benchmarks/self-energy/gamma1/sources/qgraf-r001.json)出发，给出两个颜色通道各自的完整 $C,J$，并保留可直接生成有序边界项的偏导结构。对象是 **3e-Gamma1**，不含外腿 Π₁，不含顶点反项。本文不计算 $I$ 的 Laurent 系数。

算符以 [DOT 全 D 理论输入](D:/wilson-agent/benchmarks/self-energy/gamma1/theory.md)为依据：Davydychev–Osland–Tarasov，[hep-ph/9605348v3](https://arxiv.org/abs/hep-ph/9605348)，式(2.2)、(2.4)、(3.4)–(3.9)、(3.17)–(3.22)，结合[1999勘误](https://doi.org/10.1103/PhysRevD.59.109901)。本次以独立张量约化，精确证明下文无 Gram 分母的分子重构该算符；不再仅以 UV 留数或有限数值点认定等价。

四图的**理论 C/J 已齐备**，[固定数据](D:/wilson-agent/benchmarks/self-energy/gamma1/sources/theory-cj-r001.json)含展开后的微分多项式、原域和 IBP 通量算子；这不是生产 Montopar 的导出包。实际 Montopar 专用模块仍待接入，当前拒绝未实现 Γ₁ 的保护保留。$I$ 仍为 0/4，`formal=not_run`。[检查和适用范围](D:/wilson-agent/workflow/checks/gamma1-parametric/README.md)。旧设计另存于[历史 r001](D:/wilson-agent/benchmarks/self-energy/gamma1/history/parametric-design-r001.md)。

## 1. 沿用树顶点的几何，重新计算圈修正的权重

3e-3g 已经提供了需要保留的骨架：三条 Wilson 外腿连接一个中心，外部动量变为位置偏导，而沿本腿切向的偏导变成路径全导数。Γ₁ 改变了中心的张量和参数权重，但没有改变三个外接点的几何。我们的任务是把这部分新增内容算完整，再接回同一个星核。

采用 $(+---)$ 度规、$D=4-2\epsilon$、Feynman 规范（DOT 的 $\xi=0$）、无质量夸克和 $\operatorname{Tr}_{D}(1)=4$。$v_{ij}=v_i\cdot v_j$ 没有额外的 2。中心所有动量流入，$p_i=-\ell_i$，$p_1+p_2+p_3=0$；因此外部相位是 $e^{+i\sum p_iX_i}$。必须区分 Wilson 边向量的 $v_i^2=0$ 与中心动量的 $p_i^2$：后者不能设零。

取
$$
A_0=-2v_{14},\qquad z=\frac{v_{12}}{v_{14}},\qquad
\xi_j=\frac{X_j}{\sqrt{A_0}},\qquad n_j=\frac{v_{e_j}}{\sqrt{A_0}}.
$$
先在 $A_0>0,z>0$ 的解析区域提取尺度。所有坐标偏导在三个独立的 $\xi_j$ 上进行，$n_j,A_0$ 固定；求导后才代入轮廓或接触面。特别是，不能先令同边距离为零再做法向微分。

源中 Wilson 顶点与三传播子的显式相位相消，颜色迹给出
$$
i^3(-i)^3=1,\qquad
\frac{1}{N_c}f^{abc}\operatorname{Tr}(t^ct^bt^a)=-\frac{i}{2}C_FC_A.
$$
内部圈有 $C_A$ 与 $T_Fn_f$ 两部分。它们具有不同的张量权重，最后得到 $C_gJ_g+C_qJ_q$，不是同一个 $J$ 乘两个自能系数。

## 2. 用无 Gram 分母的分子重构完整 Γ₁

直接将 DOT 的约化形状函数 Fourier 变换，会遇到 Gram 分母和差商。它们来自张量约化，未必适合下一步积分。更合适的做法是保留圈积分的多项式分子，并用 DOT 检验其完整性。

以下分子恒等式对任意三个切向量成立；进入最终 $J$ 时统一使用第 1 节的无量纲 $n_j$。记
$$
d_1=(k-p_2)^2+i0,\quad d_2=(k+p_1)^2+i0,\quad d_3=k^2+i0.
$$
去掉 $g_0^2$ 和内部颜色因子的顶点采用明确字典
$$
\mathcal T_c(p)=\frac1i\int\frac{d^Dk}{(2\pi)^D}
\left[\frac{N_c(k,p)}{d_1d_2d_3}
+\sum_{r=1}^3\frac{M_{c,r}(p)}{[k^2+i0][(k+p_r)^2+i0]}\right],
\qquad c=g,q.
$$
这里所有图权、内部颜色比和闭圈负号都已放入 $N,M$，不得再额外乘一次。该字典与源的 $\mathcal V_1=i\Gamma^{(1)}_{\rm DOT}$ 相连。

### 2.1 胶子、鬼圈和剑鱼图

用一个向量空间上的矩阵写胶子分子，可以保留 $D$ 维迹而避免长串指标。定义
$$
\mathsf M(n,p,k)=n\otimes(k+2p)+(-2n\cdot k-n\cdot p)\mathbf1
+(k-p)\otimes n,
$$
其中 $(a\otimes b)^\mu{}_{\nu}=a^\mu b_\nu$、$\operatorname{tr}\mathbf1=D$。三个树顶点的收缩为
$$
G=\operatorname{tr}\!\left[
\mathsf M(n_1,p_1,k)
\mathsf M(n_2,p_2,k-p_2)
\mathsf M(n_3,-p_1-p_2,k+p_1)\right].
$$
两个鬼圈方向给
$$
G_a=[n_1\cdot(k+p_1)][n_3\cdot(k-p_2)](D:/wilson-agent/benchmarks/self-energy/gamma1/n_2/cdot k),
\qquad
G_b=(n_1\cdot k)[n_2\cdot(k-p_2)][n_3\cdot(k+p_1)].
$$
于是完整三角分子是
$$
\boxed{N_g=-\frac12(G+G_a+G_b).}
$$
四胶子顶点与三胶子顶点构成的三个剑鱼图，可在积分后等价表示为不依赖 $k$ 的一次分子：
$$
B_1=(n_1\cdot n_2)(n_3\cdot p_1)-(n_1\cdot n_3)(n_2\cdot p_1),
$$
$$
B_2=(n_2\cdot n_3)(n_1\cdot p_2)-(n_1\cdot n_2)(n_3\cdot p_2),
$$
$$
B_3=(n_1\cdot n_3)(n_2\cdot p_3)-(n_2\cdot n_3)(n_1\cdot p_3),
\qquad \boxed{M_{g,r}=-\frac94B_r.}
$$
胶子三角的内部颜色比为 $C_A/2$；剑鱼色收缩给 $3C_A/2$，再结合 Lorentz 收缩和 $1/2$ 对称因子得到上述系数。[独立图权与分子推导](D:/wilson-agent/workflow/checks/gamma1-gluon-numerator/README.md)给出详细收缩及控制。这里的分解允许已经积分为零的项重新组织；真正与源连接的对象是完整积分恒等式。

### 2.2 夸克圈

两种取向合并并提取 $T_Fn_f$ 后，可以取
$$
\boxed{N_q=-\operatorname{Tr}_{D}
[\not n_1(\not k+\not p_1)\not n_3(\not k-\not p_2)\not n_2\not k],
\qquad M_{q,r}=0.}
$$
Dirac 迹按 $\{\gamma^\mu,\gamma^\nu\}=2g^{\mu\nu}$ 和 $\operatorname{Tr}_{D}(1)=4$ 计算。这里的整体负号是 canonical 字典的一部分；旧预研文件保存的裸 trace 恰为它的相反数。

上述矩阵迹与 Dirac 迹都只是有限多项式的紧凑记法。[理论数据](D:/wilson-agent/benchmarks/self-energy/gamma1/sources/theory-cj-r001.json)已实际展开 $N_c(\bar k,p)$ 和 $\Box_kN_c(\bar k,p)$，因此下一阶段无需重新计算这些迹。

### 2.3 为什么这确实是所引用的全 D 顶点

这一点用解析约化验证。设 $x=p_1^2,y=p_2^2,a=p_1\cdot p_2$，将 $k$ 分成 $p_1,p_2$ 张成的平面部分与正交部分。秩不超过三时只需要
$$
\langle k_\perp^\mu\rangle=0,\qquad
\langle k_\perp^\mu k_\perp^\nu\rangle
=\frac{k_\perp^2}{D-2}g_\perp^{\mu\nu}.
$$
再用
$$
k\cdot p_1=\frac{d_2-d_3-x}{2},\qquad
k\cdot p_2=\frac{y+d_3-d_1}{2}
$$
约掉分母。剩余三角和泡积分分别为 $\varphi,\kappa_1,\kappa_2,\kappa_3$；单传播子的无质量多项式积分在维数正规化中为零。泡积分若带缺失分母的幂，只需零、一、二次矩，均可由 Lorentz 分解精确求出。

按此方法，胶子和夸克两个通道各四个母积分系数，与 DOT 全张量收缩的对应系数在一般 $D,x,y,a$ 下逐项相等。这检验包含 A、B、C、F、H 的全部循环项，未删去横向部分。中间约化在非零 Gram 的一般区域证明，随后按完整解析函数延拓；最终参数表示不含 Gram 分母。八条精确系数等式见[夸克](D:/wilson-agent/workflow/checks/gamma1-parametric/quark-exact.json)及[胶子](D:/wilson-agent/workflow/checks/gamma1-parametric/gluon-exact.json)。

## 3. 先积内部圈，把三角化成平移后的星核

三个外腿 Schwinger 参数为 $a_i$，内部为 $b_i$，$B=b_1+b_2+b_3$。在欧氏起始区域，二次型矩阵为
$$
\mathsf A=\begin{pmatrix}
a_1+a_3+b_2&a_3&b_2\\
a_3&a_2+a_3+b_1&-b_1\\
b_2&-b_1&B
\end{pmatrix}.
$$
完成内部平方，令 $t_i=b_i/B$，得到
$$
\bar k=-t_2p_1+t_1p_2,\qquad
A_1=a_1+Bt_2t_3,\quad A_2=a_2+Bt_1t_3,\quad A_3=a_3+Bt_1t_2.
$$
行列式和位置多项式分别是
$$
U=B(A_1A_2+A_1A_3+A_2A_3),
$$
$$
F=B(A_3\rho_{12}+A_2\rho_{13}+A_1\rho_{23}),\qquad
\rho_{ij}=(X_i-X_j)_E^2.
$$
标量高斯为 $(4\pi)^{-3D/2}U^{-D/2}e^{-F/(4U)}$。三次分子只产生两个 Wick 矩：
$$
N_{c,0}(t,p)=N_c(\bar k,p),\qquad
L_{c,0}(t,p)=(\Box_kN_c)(\bar k,p).
$$
全部 $D$ 依赖都保留在迹和 $\Box_k$ 中。先作内部 Minkowski 积分时，这两个矩的相对号为
$$
\mathcal T_c^\triangle=(4\pi)^{-D/2}\int_{\Delta_t}d^2t
\left[-\Gamma(1+\epsilon)\frac{N_{c,0}}{\Delta_t^{1+\epsilon}}
+\frac{\Gamma(\epsilon)}4\frac{L_{c,0}}{\Delta_t^\epsilon}\right],
$$
$$
\Delta_t=-t_2t_3p_1^2-t_1t_3p_2^2-t_1t_2p_3^2-i0.
$$
这也解释了为什么不能把全顶点看成树顶点乘一个标量泡系数。

## 4. 参数域和两种等价的圈权重

平移 $a_i\to A_i$ 后，原条件 $a_i>0$ 给出真正的域
$$
0<B<B_{\max},\qquad
B_{\max}=\min\left(\frac{A_1}{t_2t_3},\frac{A_2}{t_1t_3},\frac{A_3}{t_1t_2}\right).
$$
取 $A_i=1/(\sigma\alpha_i)$、$\sum\alpha_i=1$，记
$$
P=\alpha_1\alpha_2\alpha_3,\qquad
m(\alpha,t)=\max(\alpha_1t_2t_3,\alpha_2t_1t_3,\alpha_3t_1t_2).
$$
则 $B_{\max}=1/(\sigma m)$，而
$$
d^3A\,(A_1A_2+A_1A_3+A_2A_3)^{-D/2}
=d\sigma\,d^2\alpha\,\sigma^{-2\epsilon}P^{-\epsilon}.
$$
内部 $B$ 积分只有
$$
\frac{B_{\max}^{1+\epsilon}}{1+\epsilon},\qquad
\frac{B_{\max}^{\epsilon}}{\epsilon}.
$$
这给出有限双 simplex 的表示。若希望去掉 max，按哪个限制最小分三个扇区。在第 $k$ 扇区，把另两个指标记为 $i,j$，令
$$
t_k=\frac{\alpha_k}{Z},\quad t_i=\frac{\alpha_i u}{Z},\quad
t_j=\frac{\alpha_jv}{Z},\quad Z=\alpha_k+\alpha_i u+\alpha_jv,\quad u,v\ge1.
$$
Jacobian 为 $P/Z^3$，$B_{\max}=Z^2/(\sigma Puv)$。两种表示完全等价；三个扇区只划分一次原 $t$ 域，不是三个额外物理图。

## 5. 尺度积分、Fourier 相位与完整前因子

解析延拓后定义独立坐标上的星核
$$
q_\alpha(\xi)=-\sum_{i<j}\alpha_i\alpha_j(\xi_i-\xi_j)^2+i0.
$$
在有量纲变量中 $Q=A_0q_\alpha$。三次项的尺度母函数必须减去一个常数：
$$
G_0(Q)=\int_0^\infty d\sigma\,\sigma^{-1-3\epsilon}
(e^{-\sigma Q/4}-1)=\Gamma(-3\epsilon)(Q/4)^{3\epsilon},
$$
$$
G_1(Q)=\int_0^\infty d\sigma\,\sigma^{-3\epsilon}e^{-\sigma Q/4}
=\Gamma(1-3\epsilon)(Q/4)^{-1+3\epsilon}.
$$
对固定 $Q>0$，共同收敛条带为 $0<\Re\epsilon<1/3$；三次微分消去减掉的常数。这不等于宣称全部轮廓端点在这一条带绝对收敛。完整积分以共同调节及亚纯延拓定义，不能逐项丢弃发散边界。

由U′=igAU，颜色迹是晚到早而顶点标签保持。源相位需要单独核对：颜色的 $-i/2$ 与 canonical 内部圈的 $1/i$ 给 $-1/2$。三角有三个独立动量、六个分母，Wick 因子 $i^3(-1)^6=-i$；三次 Fourier 微分因子 $(-i)^3=i$，相乘为正。泡图有五个分母与一次微分，$i^3(-1)^5\times(-i)=+1$。内部 Laplacian 项的外微分次数较低，综合它的 Wick 号后也取最终公式中的同号。详见[独立相位与归一化核验](D:/wilson-agent/workflow/checks/gamma1-normalization/README.md)。

利用 $\Gamma(-3\epsilon)=-\Gamma(1-3\epsilon)/(3\epsilon)$，统一提取
$$
\boxed{C_g=-\frac{g_0^6C_FC_A^2\Gamma(1-3\epsilon)A_0^{3\epsilon}}
{8192\pi^{6-3\epsilon}},\qquad
C_q=-\frac{g_0^6C_FC_AT_Fn_f\Gamma(1-3\epsilon)A_0^{3\epsilon}}
{8192\pi^{6-3\epsilon}}.}
$$
其中 $8192=2\cdot4^6$：$1/2$ 来自色迹，三圈高斯的 $4^{-6+3\epsilon}$ 与尺度的 $4^{-3\epsilon}$ 相乘。Laplacian 自带 $1/4$ 恰好抵消 $G_1$ 多出的 4；泡项没有该 $1/4$，所以其 $J$ 留有系数 4。

这里已还原并重新积分 DOT 的内部圈，**不再额外乘** $\eta/(4\pi)^{D/2}$。否则会把该圈计入两次。

## 6. 最终的两个 J

用 $\mathcal N_c(t,\partial)=N_{c,0}(t,p)|_{p\to\partial_\xi}$、$\mathcal L_c(t,\partial)=L_{c,0}(t,p)|_{p\to\partial_\xi}$、$\mathcal M_{c,r}=M_{c,r}|_{p\to\partial_\xi}$。这里不再附加 $i$，它已经计入上一节的相位链。偏导只作用于 $q_\alpha$，不作用于参数权重。

每个代表图 $o$ 的原始贡献是
$$
\boxed{\mathcal A_o=C_gJ_{o,g}+C_qJ_{o,q}.}
$$
有限双 simplex 形式为
$$
\boxed{\begin{aligned}
J_{o,c}^{\triangle}={}&\int_{\mathcal D_o}d^3s
\int_{\Delta_\alpha}d^2\alpha\int_{\Delta_t}d^2t\,P^{-\epsilon}
\left[-\frac{m^{-1-\epsilon}}{3\epsilon(1+\epsilon)}
\mathcal N_c\,q_\alpha^{3\epsilon}
+\frac{m^{-\epsilon}}\epsilon\mathcal L_c\,q_\alpha^{-1+3\epsilon}\right],\\
J_{o,c}^{b}={}&\frac{4B(1-\epsilon,1-\epsilon)}\epsilon
\int_{\mathcal D_o}d^3s\int_{\Delta_\alpha}d^2\alpha\,P^{-\epsilon}
\sum_{r=1}^3\alpha_r^{-\epsilon}\mathcal M_{c,r}\,q_\alpha^{-1+3\epsilon},\\
J_{o,c}={}&J_{o,c}^{\triangle}+J_{o,c}^{b}.
\end{aligned}}
$$
因此夸克没有泡项；胶子泡项可以直接写成
$$
J_{o,g}^{b}=-\frac{9B(1-\epsilon,1-\epsilon)}\epsilon
\int_{\mathcal D_o}d^3s\int_{\Delta_\alpha}d^2\alpha\,P^{-\epsilon}
\sum_r\alpha_r^{-\epsilon}B_r(\partial)q_\alpha^{-1+3\epsilon}.
$$
三角是七维参数积分，泡项是五维；公共尺度和泡图内部比例参数已经积掉。若使用矩形扇区，三角的等价式为
$$
\boxed{\begin{aligned}
J_{o,c}^{\triangle}={}&\int_{\mathcal D_o}d^3s\int_{\Delta_\alpha}d^2\alpha
\sum_{k=1}^3\int_1^\infty du\int_1^\infty dv\,
\frac{P^{-2\epsilon}Z^{-1+2\epsilon}}{(uv)^{1+\epsilon}}\\
&\times\left[-\frac{\mathcal N_c(t,\partial)q_\alpha^{3\epsilon}}{3\epsilon(1+\epsilon)}
+\frac{Puv}{\epsilon Z^2}\mathcal L_c(t,\partial)q_\alpha^{-1+3\epsilon}\right].
\end{aligned}}
$$
每个扇区分别代入自己的 $Z,t$。这就是从完整 Γ₁ 到参数空间的 $C,J$：所有分子、域、权重、颜色和前因子均已定义，没有未确定形状函数。该选择优先保留 IBP 结构，不声称变量数或主积分数已经最小。

## 7. 四个原有序域和共同几何

令 $x,y,w$ 对应源中三个外接槽位。下表的积分顺序严格为外到内。

| 图 | 几何 | 外到内顺序与域 | $\rho_{12}/A_0$ | $\rho_{13}/A_0$ | $\rho_{23}/A_0$ |
|---|---|---|---|---|---|
| 01 | 344 | $x\in[0,1],w\in[0,1],y\in[0,w]$ | $z(1-x)y$ | $z(1-x)w$ | $0$ |
| 02 | 334 | $y\in[0,1],x\in[0,y],w\in[0,1]$ | $0$ | $z(1-x)w$ | $z(1-y)w$ |
| 03 | 244 | $x\in[0,1],w\in[0,1],y\in[0,w]$ | $D_z(x,y)$ | $D_z(x,w)$ | $0$ |
| 04 | 234 | $x\in[0,1],y\in[0,1],w\in[0,1]$ | $(1-x)y$ | $D_z(x,w)$ | $z(1-y)w$ |

其中 $D_z(s,t)=(1-s)(1-t)+zst$，表内 $\rho$ 表示从欧氏区域延拓的 $-(X_i-X_j)^2$。只有完成父空间微分后，才使用表中多项式计算函数值。具体端点始终由
$$X^{(e)}(s)=\sum_{j<e}v_j+s v_e$$
确定；$v_{12}=v_{34}=zv_{14}$、$v_{23}=v_{14}$、$v_{13}=v_{24}=-(1+z)v_{14}$。

## 8. 从最终 J 直接生成边界项和余项

先在多项式算子层面做分解。令
$$
\delta_i=n_i\cdot\partial_i=\partial_{s_i},\qquad
\mathcal O=\sum_i\delta_i\mathcal O_i+\mathcal R_{\mathcal O},
$$
其中 $\mathcal O$ 依次取 $\mathcal N,\mathcal L,\mathcal M_r$。这些是常系数算子；$\mathcal O_i$ 可以继续含二阶环境导数。动量守恒的替换在这里来自母核的平移恒等式 $(\partial_1+\partial_2+\partial_3)q_\alpha=0$，同样适用于其导数；它不是任意独立坐标函数上的额外约束。沿指定顺序逐个把 $\delta_i$ 置零，取两次多项式之差除以 $\delta_i$，即可得到可重现的分解，完全不需要除以几何函数。

对 344/244，$n_2=n_3$ 且 $n_2^2=0$。每个非零单项对三条切向量各线性，所以必含重复方向与某动量的点积；这些点积全属于 $(\delta_2,\delta_3)$：
$$
n_2\cdot p_1=n_3\cdot p_1=-\delta_2-\delta_3,\qquad
n_2\cdot p_2=n_3\cdot p_2=\delta_2.
$$
因此可只按 $\delta_2,\delta_3$ 分解，$\mathcal O_1=0$、$\mathcal R_{\mathcal O}=0$。同理，334 只需 $\delta_1,\delta_2$，可取 $\mathcal O_3=0$、余项为零。这一论证适用于完整 Γ₁，包括 F、H，并非把树级零项规则套到一圈。

把第 6 节中的 $\mathcal N,\mathcal L,\mathcal M$ 分别换成它们的 $\mathcal O_i$，保留所有参数权重，得到通量 $\mathcal F_i$；换成各自余项则得到体核 $\mathcal R$。于是
$$
J=\int_{\mathcal D_o}d^3s\left(\sum_i\partial_{s_i}\mathcal F_i+\mathcal R\right).
$$
344/244 的完整式只有三个面：
$$
\boxed{J=\int_0^1dx\left[
\int_0^1dy\,\mathcal F_3(x,y,1)
-\int_0^1dw\,\mathcal F_2(x,0,w)
+\int_0^1dw\,(\mathcal F_2-\mathcal F_3)(x,w,w)\right].}
$$
334 也只有三个面：
$$
\boxed{\begin{aligned}
J={}&\int_0^1dx\int_0^1dw\,\mathcal F_2(x,1,w)
-\int_0^1dy\int_0^1dw\,\mathcal F_1(0,y,w)\\
&+\int_0^1dy\int_0^1dw\,(\mathcal F_1-\mathcal F_2)(y,y,w).
\end{aligned}}
$$
参数积分隐含于通量中，不改变其权重。334 保留 $w$ 最内层，边界测度继承父域的原迭代顺序。接触面上的两个法向贡献必须相减后完整保留。234 采用 $\delta_1,\delta_2,\delta_3$ 的顺序分解，保留六个面及非零体余项。所有面上的环境导数先完成，再代入面约束。[固定数据](D:/wilson-agent/benchmarks/self-energy/gamma1/sources/theory-cj-r001.json)保存了四图两个颜色逐项的通量和余项算子。

## 9. 若要展开成普通有理参数核

保留微分不妨碍数值或符号实现。记
$$
R_i=\xi_i-\sum_j\alpha_j\xi_j,\qquad H_{ij}=\alpha_i\delta_{ij}-\alpha_i\alpha_j.
$$
则
$$
\partial_{i\mu}q=-2\alpha_iR_{i\mu},\qquad
\partial_{i\mu}\partial_{j\nu}q=-2H_{ij}g_{\mu\nu},\qquad
\partial^3q=0.
$$
对于任意常方向导数 $D_a,D_b,D_c$，记 $q_a=D_aq$、$q_{ab}=D_aD_bq$，有
$$
D_aD_bD_cq^\beta=
\beta(\beta-1)(\beta-2)q_aq_bq_cq^{\beta-3}
+\beta(\beta-1)(q_{ab}q_c+q_{ac}q_b+q_{bc}q_a)q^{\beta-2}.
$$
二阶度规收缩须使用 $g^\mu{}_{\mu}=D$。因此直接核只含 $q^{-3+3\epsilon}$ 与 $q^{-2+3\epsilon}$ 两种幂；IBP 通量只需二阶或零阶导数。它们可以按同一组公式求值，无须再做动量积分。

## 10. raw、prepared 与下一步的 I

$C_g,C_q$ 是单个有序代表的 raw 前因子：不含 C4、不含 Euler 准备、不做颜色投影。若 $\lambda_b=g_0^2/(4\pi^2)$，公共尺度为
$$
C_c=-\lambda_b^3\,\mathrm{color}_c(\pi A_0)^{3\epsilon}
\frac{\Gamma(1-3\epsilon)}{128}.
$$
prepared 只额外加入一次 $e^{-3\gamma_E\epsilon}$，并保存原 $C$。原 C4 轨道重数为 4；完成几何旋转映射后，同一归一化规则下的每个颜色分量为
$$
I_{o,c}=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-3\epsilon)}{64}
\frac{J_{o,c}(z)+z^{3\epsilon}J_{o,c}(1/z)}{1+z^{3\epsilon}}.
$$
这给出以后求 $I$ 的目标，不是本轮已经求出的解析解。实际生产接入时仍须逐源重放 C4 的端点顺序与尺度交换；不能把此代数组装公式当成已经发布的 prepared 输出。一般复运动学中的幂始终沿源的同一解析路径定义。

当前工作流位置是：**qgraf 完整源 4/4 → 完整参数化理论 C/J 4/4 → Montopar 专用实现待接入 → I 0/4**。下一可执行任务是据本文接入专用模块，输出原有序域 raw 核及可选三面／六面加余项的 IBP 核，完成导出回归后从 344 的两颜色开始求 $I$。默认 47 根、旧积分结果及形式证明状态不因本文自动改变。
