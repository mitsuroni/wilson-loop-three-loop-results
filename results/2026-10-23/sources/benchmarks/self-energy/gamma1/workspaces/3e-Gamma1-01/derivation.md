# 3e-Gamma1-01：344 首图的两颜色 I

## 当前物理积分的入口（2026-10-06）

从共同作用量出发，U′=igAU确定真实Wilson词；FT/Dyson规则在提取C前生效。当前[原生C/J](../../sources/montopar-r002.json)给出下式，J就是本文实际积分的原核，域、分支和IBP全保留：

$$I_{P,c}=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-3\epsilon)}{64}\frac{J_c(z)+z^{3\epsilon}J_c(1/z)}{1+z^{3\epsilon}}.$$

用下文已经求出的J代入此定义，得到[当前物理I](../../results/3e-Gamma1-01/r002.json)；组装系数为+1，Euler只补一次。以下带日期的求值过程与r001公式属于旧根历史，计算J的步骤可按同域积分线性复用，旧形式状态不继承。

2026-10-03，用户在暂停汇报及完整流程说明后明确恢复求解。起点为正式 [Montopar raw 包](../../sources/montopar-r001.json)的第01轨道，SHA256 `3b92e23125cb6b654ad63d5c19c17fead9cfeb3bff56739dd9f07c6eda87e65a`。本工作区仅保存推导；脚本及计算记录位于 [gamma1-01 检查目录](../../../../../workflow/checks/gamma1-01/)。

**2026-10-03完成：两颜色全部极点至有限项已解析求出。** 夸克最高五阶，胶子最高六阶；最后HH母项化成普通一维polylog积分并实际积完。正确MB展开的42个系数、HH的14个符号连接及4个Mellin对、原27项80位实卷积、完整12项组装控制通过。正式结果见 [01/r001](../../results/3e-Gamma1-01/r001.json)，02由完整反射得到同值；本分支有限阶2/4，全ε闭式仍0/4。`formal=not_run`。

## 0. 从初始积分读起：与 3e-3g-01 的完整流程对照

本节在2026-10-03暂停期间按用户要求补写；随后用户明确恢复，现已完成首图有限阶求值，结果见§2–4。先区分两个容易混用的量：qgraf 给出的原单图贡献记为 $\mathcal A_{01}$，它含耦合、颜色、量纲尺度及全部积分；本文最终的 $I_{01,c}$ 是在固定规范下提取这些因子并组装旋转后得到的无量纲函数。下面每次换表示，都要说明积分变量和前因子去了哪里。

两张01图的三个外接点完全一样，都是344：一个点在边3，两个有序点在边4。令
$$
x=s_{31},\quad y=s_{41},\quad w=s_{42},\qquad
\int_{\mathcal D_{01}}d^3s:=\int_0^1dx\int_0^1dw\int_0^wdy.
$$
若四边形边向量为 $v_i$，起点取原点，则
$$
X_1=v_1+v_2+xv_3,\qquad
X_2=v_1+v_2+v_3+yv_4,\qquad
X_3=v_1+v_2+v_3+wv_4.
$$
$X_1,X_2,X_3$ 是三条胶子腿的位置，不能把它们的下标理解成四边形的边号。树级01的中心是树三胶子顶点；Γ₁-01把同一个中心替换成裸的一圈1PI三胶子顶点。外腿自能 Π₁和顶点反项都不包含在本图中。

两者的共同动量空间起点，可以用中心顶点 $\mathcal V$ 写成一个完整公式：
$$
\begin{aligned}
\mathcal A_{01}[\mathcal V]={}&\int_{\mathcal D_{01}}d^3s\,
\frac1{N_c}\operatorname{Tr}\left[
(ig_0v_3^{\mu_1}t^{a_1})(ig_0v_4^{\mu_2}t^{a_2})
(ig_0v_4^{\mu_3}t^{a_3})\right]\\
&\times\int\prod_{j=1}^3\frac{d^D\ell_j}{(2\pi)^D}\,
(2\pi)^D\delta^{(D)}(\ell_1+\ell_2+\ell_3)
e^{-i(\ell_1X_1+\ell_2X_2+\ell_3X_3)}\\
&\times\prod_{j=1}^3
\frac{-ig_{\mu_j\nu_j}\delta^{a_jb_j}}{\ell_j^2+i0},
\mathcal V^{b_1b_2b_3,\nu_1\nu_2\nu_3}
(-\ell_1,-\ell_2,-\ell_3).
\end{aligned}
$$
三条传播子各有自己的动量；δ消去一个之后，剩下两个独立动量积分。取 $\mathcal V=\mathcal V_0$ 时，树顶点为
$$
\mathcal V_0^{abc,\mu\nu\rho}(p_1,p_2,p_3)
=g_0f^{abc}\left[g^{\mu\nu}(p_1-p_2)^\rho
+g^{\nu\rho}(p_2-p_3)^\mu
+g^{\rho\mu}(p_3-p_1)^\nu\right],
\quad p_i=-\ell_i.
$$
于是原单图是 $g_0^4$ 阶、两个独立动量积分。取 $\mathcal V=\mathcal V_1$ 时，中心内部另有一个圈动量 $k$，所以总共三个独立动量积分、$g_0^6$ 阶。原源使用已经含内部圈的有效顶点；为了参数化，我们把它等价重写成下面的原始圈积分，而不是在已经积分的顶点外再额外乘一次圈测度。

先把树级01真正算到底。令 $A_0=-2v_{14}$、$z=v_{12}/v_{14}>0$，并写 $a=\alpha_1,b=\alpha_2,c=1-a-b$。三传播子的 Fourier/星型参数化给出同一个几何二次式
$$
q=-\sum_{i<j}\alpha_i\alpha_j(\xi_i-\xi_j)^2+i0,
\quad\xi_i=X_i/\sqrt{A_0},\qquad
q\big|_{344}=z(1-x)a(by+cw).
$$
在推导中先对独立的 $\xi_i$ 求导，再代入344轮廓。由于 $v_4^2=0$，两个同边点的平方距离为零；另外两个距离给出上式。这里的 α 是用来结合三条外传播子的参数，满足 $a,b,c\ge0$、$a+b+c=1$，并不是三条独立的单位区间。

实际树级输入的完整参数积分为
$$
\boxed{J_{\rm tree}(z,\epsilon)=
\int_{\mathcal D_{01}}d^3s
\int_0^1da\int_0^{1-a}db\,
(abc)^{-\epsilon}
\frac{z^2(1-x)a(b-c)}4\,
q^{-2+2\epsilon}.}
$$
树顶点的分子只含一次动量，Fourier 变换把它转成一次坐标导数。对这个具体几何，可以直接看见它就是路径全导数：
$$
\frac{z^2(1-x)a(b-c)}4q^{-2+2\epsilon}
=\frac{z(\partial_w-\partial_y)}{4(1-2\epsilon)}
q^{-1+2\epsilon}.
$$
因为 $\partial_yq=z(1-x)ab$、$\partial_wq=z(1-x)ac$，此式只需链式法则。参数权重 $(abc)^{-\epsilon}$ 不依赖路径变量，所以可以沿原有序域作 IBP。若定义
$$\Phi(x,y,w)=\int_{\Delta_\alpha}d^2\alpha\,(abc)^{-\epsilon}q^{-1+2\epsilon},$$
则三面形式为
$$
J_{\rm tree}=\frac{z}{4(1-2\epsilon)}\int_0^1dx
\left[\int_0^1dy\,\Phi(x,y,1)
+\int_0^1dw\,\Phi(x,0,w)
-2\int_0^1dw\,\Phi(x,w,w)\right].
$$
最后一项来自同边两点重合的接触面；系数−2由两个法向通量共同产生。这张图的路径 IBP 没有体余项。实际来源和算子见[树级规范输入](../../../../tgv/3e-3g/sources/canonical.json)及[prepared](../../../../tgv/3e-3g/prepared/3e-3g.json)；旧理论文§5b/5c的局部 Fourier 符号问题不作为本节重新取号的依据。

同一个积分也可以沿尺度直接降维，这正是树级01的完整求值路线。设
$$
r=1-x,\quad y=wh,\quad
a=A,\quad b=(1-A)\beta,\quad c=(1-A)(1-\beta).
$$
两个 Jacobian 分别为 $w$ 和 $1-A$。完整核变为
$$
\frac{J_{\rm tree}}{z^{2\epsilon}}=
\frac14\int_0^1dr\,r^{-1+2\epsilon}
\int_0^1dw\,w^{-1+2\epsilon}
\int_0^1dA\,A^{-1+\epsilon}
\int_0^1d\beta\int_0^1dh\,
[\beta(1-\beta)]^{-\epsilon}(2\beta-1)
(1-\beta+\beta h)^{-2+2\epsilon}.
$$
值得观察的是所有 $(1-A)$ 幂恰好抵消。因此 $r,w,A$ 三个积分分别是 $1/(2\epsilon),1/(2\epsilon),1/\epsilon$。再积 $h$：
$$
\int_0^1dh\,(1-\beta+\beta h)^{-2+2\epsilon}
=\frac{1-(1-\beta)^{-1+2\epsilon}}{\beta(-1+2\epsilon)}.
$$
最终只剩一个明确的一维积分：
$$
\frac{J_{\rm tree}}{z^{2\epsilon}}=
\frac1{16\epsilon^3}\int_0^1d\beta\,
\beta^{-1-\epsilon}(1-\beta)^{-\epsilon}(2\beta-1)
\frac{1-(1-\beta)^{-1+2\epsilon}}{-1+2\epsilon}.
$$
端点差整体保持，然后使用 Beta 恒等式及其解析延拓，得到
$$
J_{\rm tree}=\frac{z^{2\epsilon}}{16(1-2\epsilon)\epsilon^4}
\left[2\Gamma(1+\epsilon)\Gamma(1-\epsilon)
-B(1-\epsilon,1-\epsilon)\right].
$$
按已经固定的 Euler/C4 规范化，最后为
$$
\boxed{I_{\rm tree}=
\frac{e^{-2\gamma_E\epsilon}\Gamma(1-2\epsilon)}{32\epsilon^4}
\left[2\Gamma(1+\epsilon)\Gamma(1-\epsilon)
-B(1-\epsilon,1-\epsilon)\right].}
$$
到此树级01的所有积分都已消去，得到全 ε 闭式；展开它才产生极点和有限项。[树级现行方法](../../../../tgv/3e-3g/method.md)给出其独立控制和唯一 I 的定义。

现在把相同初始公式中的树顶点换成 Γ₁。采用参数化主理论的去色、去耦合及相位字典，中心与三个切向量收缩后的圈函数为
$$
\boxed{\mathcal T_c(p)=\frac1i\int\frac{d^Dk}{(2\pi)^D}
\left[\frac{N_c(k,p)}{d_1d_2d_3}
+\sum_{j=1}^3\frac{M_{c,j}(p)}{[k^2+i0][(k+p_j)^2+i0]}\right],}
$$
$$d_1=(k-p_2)^2+i0,\quad d_2=(k+p_1)^2+i0,\quad d_3=k^2+i0.$$
内部颜色分别是 $C_A$ 和 $T_Fn_f$；结合外部色迹产生 $C_FC_A^2$、$C_FC_AT_Fn_f$。$c=g,q$ 标记两个函数，不意味着两者共享同一个标量顶点。这里 $n_1=v_3/\sqrt{A_0}$、$n_2=n_3=v_4/\sqrt{A_0}$。夸克分子完全由
$$
N_q=-\operatorname{Tr}_D[
\not n_1(\not k+\not p_1)\not n_3(\not k-\not p_2)\not n_2\not k],
\quad\operatorname{Tr}_D1=4,\qquad M_{q,j}=0
$$
定义。胶子/鬼圈的紧凑、无 Gram 分母表达是
$$
\mathsf M(n,p,k)=n\otimes(k+2p)+(-2n\cdot k-n\cdot p)\mathbf1+(k-p)\otimes n,
$$
$$
G=\operatorname{tr}[\mathsf M(n_1,p_1,k)\mathsf M(n_2,p_2,k-p_2)
\mathsf M(n_3,p_3,k+p_1)],\qquad p_3=-p_1-p_2,
$$
$$
G_a=[n_1\cdot(k+p_1)][n_3\cdot(k-p_2)](n_2\cdot k),\qquad
G_b=(n_1\cdot k)[n_2\cdot(k-p_2)][n_3\cdot(k+p_1)],
\quad N_g=-\tfrac12(G+G_a+G_b).
$$
胶子另有三个剑鱼项 $M_{g,j}=-9B_j/4$，其中
$$
\begin{aligned}
B_1&=(n_1\cdot n_2)(n_3\cdot p_1)-(n_1\cdot n_3)(n_2\cdot p_1),\\
B_2&=(n_2\cdot n_3)(n_1\cdot p_2)-(n_1\cdot n_2)(n_3\cdot p_2),\\
B_3&=(n_1\cdot n_3)(n_2\cdot p_3)-(n_2\cdot n_3)(n_1\cdot p_3).
\end{aligned}
$$
三角部分由三条外传播子和三个 $d_j$ 构成，总计六个分母；剑鱼部分有三条外传播子和两个内部传播子，共五个分母。这就是 Γ₁ 比树级多出的具体积分内容。该多项式圈表示已逐母积分系数还原源中固定的 DOT 全 D 顶点；不需要把论文形状函数里的 Gram 分母带入求值。

为了积掉内部 $k$，给其三个分母引入内部 Schwinger 参数 $b_j=B t_j$，$t_j\ge0,\sum t_j=1$。配方后圈动量中心为
$$\bar k=-t_2p_1+t_1p_2.$$
三次分子有两个非零高斯矩：
$$N_{c,0}=N_c(\bar k,p),\qquad L_{c,0}=(\Box_kN_c)(\bar k,p).$$
前者来自平移后的分子，后者来自圈动量的二阶协方差；只代入 $k=\bar k$ 会漏掉后者。它们仍是外部动量的多项式，分别成为三次和一次坐标偏导算子。

内部圈也改变了外部 Schwinger 域。若外参数为 $a_j^{\rm ext}>0$，平移后
$$
A_1=a_1^{\rm ext}+Bt_2t_3,\quad
A_2=a_2^{\rm ext}+Bt_1t_3,\quad
A_3=a_3^{\rm ext}+Bt_1t_2.
$$
所以内部尺度并非独立地积到无穷大，而是
$$
0<B<\min\left(\frac{A_1}{t_2t_3},\frac{A_2}{t_1t_3},\frac{A_3}{t_1t_2}\right).
$$
取 $A_j=1/(\sigma\alpha_j)$ 后，这个上界为 $1/(\sigma m)$，其中
$$m=\max(at_2t_3,bt_1t_3,ct_1t_2).$$
因此 m 的来源是原来各条外 Schwinger 参数必须为正的限制。积掉 B 和公共尺度 σ 后，仍保留两个 α 比例和两个 t 比例，不能把 t 当成又一套外部路径参数。

把相位、颜色、耦合、外尺度及公共 Gamma 因子提取为
$$
C_g=\frac{g_0^6C_FC_A^2\Gamma(1-3\epsilon)A_0^{3\epsilon}}{8192\pi^{6-3\epsilon}},
\qquad
C_q=\frac{g_0^6C_FC_AT_Fn_f\Gamma(1-3\epsilon)A_0^{3\epsilon}}{8192\pi^{6-3\epsilon}},
$$
于是 $\mathcal A_{01}=C_gJ_g+C_qJ_q$。这里 raw C 不含后续的 Euler 因子或 C4 旋转重数。全部动量积分已经消去后，真正待求的参数积分是
$$
\boxed{\begin{aligned}
J_c^\triangle={}&\int_{\mathcal D_{01}}d^3s
\int_{\Delta_\alpha}d^2\alpha\int_{\Delta_t}d^2t\,(abc)^{-\epsilon}
\left[-\frac{m^{-1-\epsilon}}{3\epsilon(1+\epsilon)}
\mathcal N_c q^{3\epsilon}
+\frac{m^{-\epsilon}}\epsilon\mathcal L_c q^{-1+3\epsilon}\right],\\
J_c^b={}&\frac{4B(1-\epsilon,1-\epsilon)}\epsilon
\int_{\mathcal D_{01}}d^3s\int_{\Delta_\alpha}d^2\alpha\,(abc)^{-\epsilon}
\sum_{j=1}^3\alpha_j^{-\epsilon}\mathcal M_{c,j}q^{-1+3\epsilon},
\qquad J_c=J_c^\triangle+J_c^b.
\end{aligned}}
$$
$\mathcal N_c=N_{c,0}|_{p\to\partial_\xi}$、$\mathcal L_c=L_{c,0}|_{p\to\partial_\xi}$、$\mathcal M_{c,j}=M_{c,j}|_{p\to\partial_\xi}$。Fourier 所需的 i 已放进 C/J 的相位字典，不在这些算子上再乘一次。偏导只作用于 q 的幂；q 在微分前仍取独立坐标定义。$\Delta_\alpha$ 的显式积分为 $\int_0^1da\int_0^{1-a}db$，$\Delta_t$ 同样有两个独立变量。三角项有7个参数变量，剑鱼项有5个。

树级路径全导数结构在 Γ₁-01 仍然存在。因 $n_2=n_3$ 且 $n_2^2=0$，上述算子可以完整分解到 $\partial_y,\partial_w$。因此
$$
J_c=\int_{\mathcal D_{01}}d^3s\,(\partial_y\mathcal F_{c,2}+\partial_w\mathcal F_{c,3})
$$
也给出三个面：
$$
J_c=\int_0^1dx\left[\int_0^1dy\,\mathcal F_{c,3}(x,y,1)
-\int_0^1dw\,\mathcal F_{c,2}(x,0,w)
+\int_0^1dw\,(\mathcal F_{c,2}-\mathcal F_{c,3})(x,w,w)\right].
$$
通量内部仍含 α/t 积分与二阶或零阶坐标导数。这些面不是已经积完的数值结果。01在这一步没有路径体余项；具体算子及方向约定见[参数化理论§8](../../parametric-theory.md)。

本次实际求值使用与三面 IBP 等价的直接核，先沿 $r=1-x,y=wh$ 积掉路径参数。链式法则把三次微分写成 $q^{-3+3\epsilon}$ 和 $q^{-2+3\epsilon}$ 项；已核对各分子的齐次性使 r、w 的测度均变为 $r^{-1+3\epsilon},w^{-1+3\epsilon}$，各给 $1/(3\epsilon)$。h 的分子至多一次，而分母是 $a(bh+c)$ 的幂，因此 h 也能解析积分。剩下的四维积分才是当前主要的求值对象。

夸克通道在这一步可以写成一个没有未定系数的普通积分。令 $t_1=t,t_2=u,t_3=v=1-t-u$，则
$$
\begin{aligned}
m&=\max(auv,btv,ctu),\qquad F=cu-bv,\\
K&=abv(t+u)+acu(t+v)-bct^2,
\end{aligned}
$$
$$
\boxed{\begin{aligned}
j_q(\epsilon)={}&\frac{16(\epsilon-1)}{9\epsilon^3}
\int_0^1da\int_0^{1-a}db\int_0^1dt\int_0^{1-t}du\,
a^{2\epsilon-1}b^{-\epsilon-1}c^{-\epsilon}(cu-bv)m^{-\epsilon}\\
&\times\left[(b+c)^{3\epsilon-1}-c^{3\epsilon-1}\right]
\left(1-\frac{\epsilon K}{(1+\epsilon)m}\right).
\end{aligned}}
$$
这个公式可直接与前面的树级一维积分比较：路径积分同样产生了端点差，但是 m 和 K 把 α 与 t 耦合在一起。树级能把 A 积成 $1/\epsilon$ 的完全分离在这里尚未实现。把 m 的三个分支分域可以去掉 Max，却不会自动消去这四个剩余变量。

当前夸克求值路线是：三个 Max 扇区→共同调节下的端点差→Euler 分部积分正规化→ε展开→逐阶 HyperInt 解析积分→乘唯一规范化因子。Euler 正规化在展开前把边界发散提取为显式 ε 极点，不能先逐项展开原始发散核。全部极点和有限项已完成，细节与结果在下文§2。与树级不同，我们目前先求了有限阶 Laurent 系数，没有得到全 ε Gamma 闭式。

胶子求值另外采用内部传播子抵消。利用 $2k\cdot p_1=d_2-d_3-p_1^2$、$2k\cdot p_2=p_2^2+d_3-d_1$ 的多项式关系，在完整分子中组织出
$$N_g+\frac{D-2}{4}N_q=d_1H_1+d_2H_2+d_3H_3+R.$$
这里的抵消恒等式在传播子二次多项式层面使用，完整积分仍保留共同的 $+i0$ 边界值处方。除以 $d_1d_2d_3$ 后，前三项变成只有两个内部传播子的泡图；R 对 k 只有一次，故 $\Box_kR=0$。原剑鱼与这些泡图组成 $I_p$，剩余三角分子产生 $I_R$，从而
$$\boxed{I_g=-\frac{1-\epsilon}{2}I_q+I_p+I_R.}$$
**这里的 R 是内部圈分子抵消后的余项，并不是前述路径 IBP 的体余项。** $I_p$ 已完整积成 Gamma/Beta 组合，$I_R$ 仍保留耦合参数积分。分子与路径密度的完整多项式见[pinch-result.json](../../../../../workflow/checks/gamma1-01/bubble/decomposition/pinch-result.json)和[remainder-path.json](../../../../../workflow/checks/gamma1-01/bubble/decomposition/remainder-path.json)。

最后，路径尺度已证明 $J_c(z)=z^{3\epsilon}j_c$。按正式四重旋转及一次 Euler 因子得到
$$\boxed{I_{01,c}=\frac{e^{-3\gamma_E\epsilon}\Gamma(1-3\epsilon)}{64}\,j_c.}$$
这一步只是已求 j 的规范化及卷积，不能补出尚未求得的积分系数。暂停时夸克通道和 $I_p$ 已完成；$I_R$ 的最高两阶有独立解析确认，但单极点与有限项的旧 MB 候选因展开深度错误已撤回。原共同端点差拆开时出现零 sdd、最初漏看六阶极点，以及 MB 展开深度错误，是三个不同问题；准确范围与恢复任务见下文§4。

| 连续步骤 | 树级 3e-3g-01 | 3e-Gamma1-01 |
|---|---|---|
| 原动量积分 | 两个独立外动量，树中心 | 同样的两个外动量，加中心内部 k |
| 参数空间起点 | 三个路径＋两个 α 比例，共5维 | 三角三路径＋两 α＋两 t，共7维 |
| 外接几何 | 344、原序 x→w→y | 完全相同 |
| 星核 q | $z(1-x)a(by+cw)$ | 完全相同 |
| 偏导与路径 IBP | 一次导数；三个面，无体余项 | 三次/一次导数；三个面，无体余项 |
| 路径、尺度降维后 | 一维 β 端点差 | 三角剩四维 α/t；泡项进一步闭合 |
| 当前如何求完 | Beta/Gamma全 ε 闭式 | 两颜色有限阶均完成；胶子末项由Mellin卷积的一维polylog积分闭合 |

## 1. 固定源、顺序与唯一 I

三端点为 $(3,1),(4,1),(4,2)$。用 $x=s_{31},y=s_{41},w=s_{42}$，原有序域为
$$
\int_0^1dx\int_0^1dw\int_0^wdy.
$$
令 $a=\alpha_1,b=\alpha_2,c=1-a-b$，$z=v_{12}/v_{14}>0$，有
$$q=z(1-x)a(by+cw).$$
两颜色分别是 $C_FC_A^2$ 与 $C_FC_AT_Fn_f$。raw 的公共 Γ/耦合/尺度因子和是否含 C4、Euler 始终采用[参数化理论](../../parametric-theory.md)，不再次插入圈归一化。

取 $r=1-x,y=w h$，Jacobian 为 $w$。正式核中每个 $q^{-n+3\epsilon}$ 项的分子都精确满足
$$K_n=z^n r^{n-1}w^{n-2}\widehat K_n(a,b,t,h).$$
因此两个路径尺度分别给 $1/(3\epsilon)$，再积 $h$ 的线性分母即可。所有块均得到
$$J_{01,c}(z)=z^{3\epsilon}j_c(\epsilon).$$
九个原核齐次性恒等式和端点密度保存在 [path-reduced.json](../../../../../workflow/checks/gamma1-01/path-reduced.json)。这里先按受控调节推导，再采用完整源的共同解析延拓。

C4旋转的两个尺度分别贡献两次 $z^{3\epsilon}$ 和两次1，所以统一规范化为
$$
\boxed{I_{01,c}(\epsilon)=\frac{e^{-3\gamma_E\epsilon}\Gamma(1-3\epsilon)}{64}\,j_c(\epsilon).}
$$
因此本图的每个颜色 I 与 z 无关。原源4重计权保持，没有合并不同轨道的物理权重。

## 2. 夸克核的四参数形式

记 $t_1=t,t_2=u,t_3=v=1-t-u$，并定义
$$
m=\max(auv,btv,ctu),\qquad F=cu-bv,
$$
$$K=abv(t+u)+acu(t+v)-bct^2.$$
把正式 N0、L0 两块合并后，完整夸克积分为
$$
\boxed{\begin{aligned}
j_q(\epsilon)={}&\frac{16(\epsilon-1)}{9\epsilon^3}
\int_{\Delta_\alpha}da\,db\int_{\Delta_t}dt\,du\,
a^{2\epsilon-1}b^{-\epsilon-1}c^{-\epsilon}
F\,m^{-\epsilon}\\
&\times\big[(b+c)^{3\epsilon-1}-c^{3\epsilon-1}\big]
\left(1-\frac{\epsilon K}{(1+\epsilon)m}\right).
\end{aligned}}
$$
修正项前面的 ε 必须保留。原 N0 的 $q^{-3+3\epsilon}$ 系数在此几何恒零，但 N0 本身不为零，也不能删去圈动量协方差项。

为了进行解析展开，按最大的 $\alpha_i/t_i$ 分成三域，取
$$\alpha_i=\frac{t_iv_i}{S},\qquad S=\sum_i t_iv_i,\qquad\max_i v_i=1.$$
每域其余两个 $v_i$ 在 $[0,1]$ 内，Jacobian 为 $t_1t_2t_3/S^3$，且 $m=t_1t_2t_3/S$。再取 $t_1=\tau,t_2=(1-\tau)y$，四个单位区间映射为正半轴。完整有理因子及 ε 幂基底已作12项精确换元核验。

端点差不能在正 ε 下拆成两个普通收敛积分。采用共同辅助解析调节后，分别正规化六个图表，再对完整和取延拓。六图表的 ε=0 发散均已由显式 Euler 分部积分移除；没有关闭 HyperInt 的发散检测。

正规化后可写
$$
j_q=\frac{16(\epsilon-1)}{9\epsilon^5}
\sum_{n\ge0}T_n\epsilon^n.
$$
已实际解析求得
$$
\begin{aligned}
T_0&=\frac1{24},\qquad T_1=-\frac1{18},\\
T_2&=-\frac{143}{216}+\frac5{24}\zeta_2,\\
T_3&=-\frac{4109}{1296}+\frac{19}{18}\zeta_2+\frac74\zeta_3,\\
T_4&=-\frac{93089}{7776}+\frac{881}{216}\zeta_2+\frac23\zeta_3
+\frac{71}{80}\zeta_2^2,\\
T_5&=-\frac{1923857}{46656}+\frac{18455}{1296}\zeta_2+\frac{17}{9}\zeta_3
-\frac{23}{60}\zeta_2^2-\frac{21}{4}\zeta_2\zeta_3+\frac{125}{4}\zeta_5.
\end{aligned}
$$
$T_5$ 已经实际积分完成。先对六个块积分 u、v，合并7743个中间词项，再在 r、h 的共同基中压为897项，最后积完两变量。没有对单块先作可能改变分支抵消的 RegHlog 转换。包含分项解析积分的记录见 `quark_uv_periods.m`、`quark_finite_fibration_period.m` 和 [quark-complete.json](../../../../../workflow/checks/gamma1-01/quark-complete.json)。前四个极点另由保留辅助调节及留数的 Mellin–Barnes 路线独立求得；单极点和有限项另有独立复轮廓数值检查。

对应的已知规范化系数为
$$
\begin{aligned}
I_q={}&-\frac1{864\epsilon^5}+\frac7{2592\epsilon^4}
+\frac1{\epsilon^3}\left(\frac{131}{7776}-\frac{19\zeta_2}{1728}\right)\\
&+\frac1{\epsilon^2}\left(\frac{3251}{46656}-\frac{59\zeta_2}{5184}-\frac{17\zeta_3}{288}\right)\\
&+\frac1\epsilon\left(\frac{68435}{279936}-\frac{127\zeta_2}{15552}
+\frac{47\zeta_3}{864}-\frac{827\zeta_2^2}{11520}\right)+I_{q,0}+O(\epsilon).
\end{aligned}
$$
有限项为
$$
\boxed{I_{q,0}=\frac{1365323}{1679616}+\frac{2921\zeta_2}{93312}
+\frac{305\zeta_3}{2592}-\frac{739\zeta_2^2}{34560}
-\frac{11\zeta_2\zeta_3}{64}-\frac{1331\zeta_5}{1440}.}
$$
其数值为 $-0.350329707054862310\ldots$。独立 MB 的两组轮廓对 raw $j_{q,0}$ 的误差分别约 $5.2\times10^{-14}$ 和 $3.4\times10^{-14}$；见 [独立交叉检查](../../../../../workflow/checks/gamma1-01/triangle/quark-complete-crosscheck.json)。该通道已完成所需有限阶，但没有宣称全 ε 标准函数闭式。

## 3. 胶子通道的精确简化

对未积圈动量分子先使用传播子恒等式，得
$$\boxed{I_g=-\frac{1-\epsilon}{2}I_q+I_p+I_R.}$$
$I_p$ 包括三个内部传播子抵消后的泡积分和原剑鱼；$I_R$ 的分子对内部圈动量仅为一次，故 $L_R=0$。此等式由完整分子回代得到，不是拟合已知系数。

$I_p$ 已全部积分为 Γ 函数组合并展开到有限项，精确表达式、六个系数及独立Euler控制见 [pinch-result.json](../../../../../workflow/checks/gamma1-01/bubble/decomposition/pinch-result.json)。其中
$$I_p=-\frac1{576\epsilon^5}+0\,\epsilon^{-4}+O(\epsilon^{-3}).$$
原剑鱼已经包含在 $I_p$ 中，不可再加一次 [bubble/result.json](../../../../../workflow/checks/gamma1-01/bubble/result.json) 的独立子块。

$I_R$ 的三个路径变量也已积完。拆上下端点后出现不随 ε 调节的中间发散，朴素 ε-only 正规化试验未予采用。有效路线保留共同 $b^\delta$；Mellin–Barnes 中先继续 $\delta\to0$ 并核对两条一维留数的相消，再继续 ε。已经独立求得
$$
\begin{aligned}
I_R={}&\frac1{4608\epsilon^6}+\frac1{1152\epsilon^5}
+\frac1{\epsilon^4}\left(-\frac5{1152}+\frac{13\zeta_2}{3072}\right)\\
&+\frac1{\epsilon^3}\left(-\frac{31}{1152}+\frac{\zeta_2}{384}+\frac{121\zeta_3}{4608}\right)
+O(\epsilon^{-2}).
\end{aligned}
$$
上式展示最初确认的四阶。恢复后，修正版r003已把全部原母核展开到正确深度；42个系数与未展开原核独立Cauchy提取相符，且没有浮点污染。三个原一维块与零维留数已经积完；二维单极点利用共同核在 $x,y,1-x-y$ 上的置换对称，把每个单坐标ψ移到外变量，内圈只需有理cot卷积，外圈用Barnes导数与MZV矩积完。其直接二维轮廓数值与解析式误差约 $1.1\times10^{-14}$。

目前完整R和胶子通道的极点系数如下。记 $Z_n=\zeta(n)$：

| ε次幂 | $[\epsilon^k]I_R$ | $[\epsilon^k]I_g$ |
|---|---|---|
| −6 | $1/4608$ | $1/4608$ |
| −5 | $1/1152$ | $-1/3456$ |
| −4 | $13Z_2/3072-5/1152$ | $13Z_2/3072-65/10368$ |
| −3 | $Z_2/384+121Z_3/4608-31/1152$ | $11Z_2/6912+121Z_3/4608-271/7776$ |
| −2 | $2659Z_2^2/61440-17Z_2/2304+11Z_3/2304-31/288$ | $2659Z_2^2/61440-95Z_2/20736+7Z_3/192-6943/46656$ |
| −1 | $113Z_2^2/15360+253Z_2Z_3/3072-97Z_2/2304-Z_3/288+867Z_5/2560-2587/6912$ | $1537Z_2^2/46080+253Z_2Z_3/3072-619Z_2/15552-19Z_3/3456+867Z_5/2560-154267/279936$ |

六阶项经原参数端点的统一可积界及有效RHW解析积分确认；原三扇区又独立重现 $j_R[-4]=-5/18+5Z_2/24$。较高阶已完成证据只有在完整核与r003精确相等后才复用。前五阶及新单极点的连接、唯一Gamma卷积与组装见 [corrected-front-poles-r003.json](../../../../../workflow/checks/gamma1-01/corrected-front-poles-r003.json)，单极点组装另有[五项精确复核](../../../../../workflow/checks/gamma1-01/triangle/R-single-assembly-check-r003.json)。没有再加一次原剑鱼，也保留了 $(1-\epsilon)I_q$ 的移位贡献。

有限项也已闭合。其余二维单坐标项解析积完后，置换对称及两次轮廓IBP将8个混合块化为6个额外单坐标块与一个HH母核。这里的Mellin变量不是原路径x、y：若记为X、Y、T=1−X−Y，共同核为 $\prod_{s=X,Y,T}\Gamma(s)\Gamma(1-s)$。将分母和ψ标记分配到三个单变量Mellin变换，使用
$$
\int\frac{dX}{2\pi i}\frac{dY}{2\pi i}
F_1(X)F_2(Y)F_3(1-X-Y)
=\int_0^\infty f_1(\rho)f_2(\rho)f_3(\rho)\,d\rho,
$$
再取新的卷积实变量 $u=\rho/(1+\rho)$，得到普通绝对收敛的一维积分
$$
\mathcal H=\int_0^1\frac{du}{(1-u)^2}
G_0(u)\sum_{i,j}n_{ij}H_i(u)H_j(u).
$$
这里 $u$ 与§2的内部参数 $t_2$ 不同；$n_{ij}$ 是原母核分子的27个精确有理系数，不是待定拟合参数。令 $L=\log u$、$D_u=-u(1-u)\partial_u$，完整反演为
$$
G_0=1+\operatorname{Li}_2(1-u)+\frac{L^2}{2}+\frac{L}{1-u},
$$
$$
H_0=1-\zeta_3+\operatorname{Li}_3(u)-L\operatorname{Li}_2(u)
-\frac{L^2\log(1-u)}2+\frac{L^3}6+\frac{L+L^2/2}{1-u},\qquad H_i=D_u^iH_0.
$$
反演的导数恒等式、边界条件和测度逐项核对，u=0最坏为可积的log⁷，u=1为 $O((1-u)^3)$。27项合并成 $G_0(AH_0+B)/(1-u)^2$ 后，HyperInt单变量实际积分约36.5CPU秒，得到
$$
\boxed{\mathcal H=\frac{56}3+\frac{28\pi^2}{27}+\frac{112\zeta_3}9
-\frac{14\pi^4}{15}-80\zeta_5-\frac{\zeta_3^2}4+\frac{193\pi^6}{8505}.}
$$
[完整Mellin连接和端点论证](../../../../../workflow/checks/gamma1-01/bubble/decomposition/hh-convolution-proof.md)给出27项和紧凑A/B核。独立80位实积分采用原27项而非紧凑核，与闭式差约 $4.23\times10^{-72}$；另一条独立75位控制与闭式差约 $2.99\times10^{-72}$。完整二维有限项与直接复轮廓积分相符，误差约 $2.3\times10^{-13}$。

加上原零维与一维块后，raw余项有限系数为
$$
j_{R,0}=\frac{33809Z_2^3}{5040}+\frac{13Z_2^2}{60}+\frac{3Z_2Z_3}4
+\frac{184Z_2}9-\frac{115Z_3^2}{18}+\frac{62Z_3}{27}
+\frac{101Z_5}{12}-\frac{50665}{648}.
$$
按完整 $G_3=e^{-3\gamma_E\epsilon}\Gamma(1-3\epsilon)$ 卷积，特别保留
$$[\epsilon^6]G_3=\frac{44469Z_2^3}{560}+\frac{81Z_3^2}2,$$
最终胶子有限项为
$$
\boxed{I_{g,0}=\frac{150911Z_2^3}{516096}+\frac{457Z_2^2}{27648}
+\frac{Z_2Z_3}{12}-\frac{18325Z_2}{93312}+\frac{1177Z_3^2}{9216}
-\frac{347Z_3}{5184}+\frac{1907Z_5}{2880}-\frac{3201427}{1679616}.}
$$
数值为 $0.0725695427120798869\ldots$，夸克有限项为 $-0.350329707054862310\ldots$。所有极点加上述有限项给出两颜色完整Laurent到ε⁰，余项 $O(\epsilon)$；[全部组装](../../../../../workflow/checks/gamma1-01/complete-r003.json)与[独立12项控制](../../../../../workflow/checks/gamma1-01/triangle/R-finite-assembly-check-r003.json)均保留，没有宣称全ε闭式。

## 4. 反射及当前验收边界

01与02的原核满足同时反射
$$
(x,y,w)_{01}=(1-w',1-y',1-x'),\quad
(\alpha_1,\alpha_2,\alpha_3)\mapsto(\alpha'_3,\alpha'_2,\alpha'_1),\quad
(t_1,t_2,t_3)\mapsto(t'_3,t'_2,t'_1).
$$
在共同有限 I0 及置换对称的 α/t 截断下，核和域的等价允许原序变换；按相同延拓取极限后，两个颜色均有 $I_{02,c}=I_{01,c}$。点态控制见[独立核验](../../../../../workflow/checks/gamma1-01/independent/README.md)。两个源轨道仍各自计权。

01两颜色已经完成并发布[不可变r001](../../results/3e-Gamma1-01/r001.json)。上述完整核/原域反射将结果严格连接到[02推导](../3e-Gamma1-02/derivation.md)与[02/r001](../../results/3e-Gamma1-02/r001.json)；两源轨道仍分别计权，没有合并成一个物理图。本分支极点至有限项2/4，下一03的244对边几何；全ε闭式0/4。生产与冻结raw未改，未提交Git，`formal=not_run`。

### 暂停时的纠错与独立核验状态

旧 MB 单极点与有限项候选已撤回。生成器把有理分母展开后隐藏了显式 ε 负幂，本机 `MBexpand` 因而低估 Gamma 展开所需深度。直接从未展开原核用 Cauchy 系数作核对时，二维 ε⁻² 系数相符，而 ε⁻¹ 和 ε⁰ 明确不符；见 [展开审计](../../../../../workflow/checks/gamma1-01/triangle/R-old-expansion-audit.json) 与 [原因及修复方案](../../../../../workflow/checks/gamma1-01/bubble/decomposition/mb-expansion-audit.md)。旧候选的数值吻合只核验了错误截断核，不能算原积分证据。恢复后的精确log-Gamma展开r003已通过[独立42系数检查](../../../../../workflow/checks/gamma1-01/triangle/Cauchy-check-log-gamma-expanded-r003.json)，原一维末两阶亦已实际解析并经55位原核轮廓核对；旧失败仍保留，不迁移其通过标签。

不拆路径端点的五维 RHW 路线已经实际解析验证 `j_R[-6]=1/72` 和 `j_R[-5]=1/18`。RHW3 单图表的 ε⁻³ 贡献也已算得 `23 ζ3/36−1/18−2 ζ2/9`，仅为一个图表，不能当作完整 R 系数。两个较大的 RHW1/RHW2 图表的 ε⁻⁴ 作业在用户要求暂停时中止；已有中间文件和日志保留，未完成作业不计通过。
