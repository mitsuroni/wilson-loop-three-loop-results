2026-10-06共同规则已融入正文：保留中心腿/动量标记，U′=igAU确定晚到早Wilson迹；§2的G₀因此为负，§6和§13直接给出负C/R及唯一物理I。A₁内部闭圈、J、域与IBP不另变号。原r001输入保留，当前原生C/J与I为r002。

# 2026-10-01 TGV腿Pi1生产参数化

schema 1.7和`tgv_pi1`专用分派已实现以下原直接理论。原r001的“尚未生产”和“IBP未完成”文字保留为历史状态；当前IBP的完整推导及适用边界见[tgv-pi1-ibp.md](tgv-pi1-ibp.md)。

运行：`D:/py/python.exe -m montopar.cli --source ../qgraf-pipeline/outputs/3loop/orbit_definitions.json --type "3e-3g+Pi1" --tgv-ibp --package-initial`。

主输出是原s排序域和两个独立α，完整两个分母项，共五维。raw不乘Euler、不投影颜色、不乘C4重数、不组装crossing。J_ibp只由显式开关添加；每行保存自己的C、核、变量、域、面来源和共同调节组。原源中downstream_status是历史快照，实际当前状态以输出translation_status=translated为准。

本次没有求I。I的定义沿用原理论§13.3；四次旋转的q/L/E/F尺度及标记腿已在测试中核对。全部非慢速测试541项通过，单一外部CAS慢测试未运行。

---

# TGV腿Π₁：Montopar通用参数化理论

日期：2026-09-25。实现规格r001。以用户网页推导为基础，经本地复核及接口补全；原稿保存在[输入快照](../../../workflow/checks/tgv-pi1-parametric/web-input.md)，原下载文件不变。

**本稿为本分支直接参数化的实现依据；IBP输出理论尚未完成。** 源算符与12图由[动量源理论](theory.md)、[冻结源包](sources/qgraf-r001.json)固定；本稿不修改其含义。已得到可实施的raw C/J定义与输入验收规则；生产Montopar尚未实现，不能把本文或检查用几何JSON当作已发布prepared输入。与中央Γ₁分支无关。

本地验证见[检查报告](../../../workflow/checks/tgv-pi1-parametric/README.md)。下文先给完整推导，§13集中规定实现接口、分支与唯一I；程序优先输出§6.1的两个分母项。

本文采用任务输入的度规、Fourier 约定、裸自能、颜色迹及 Wilson 生成元顺序。研究对象是单个有序代表图。本文的 A 始终是[动量源理论](theory.md) §1 的原始贡献；不乘 C4 重数、镜像重数、自能对称因子或额外负号，也不引入反项。以下结果是按给定输入推导的通用参数化，不是 Laurent 系数计算。

## 1. 记号

$$
d=\frac D2=2-\varepsilon,\quad
\beta=1-\varepsilon,\quad h=1-2\varepsilon,\quad
\nu=1-3\varepsilon,\quad m=\nu+1=2-3\varepsilon.
$$

对标记腿 r，以 (r,s,t) 表示 (1,2,3)、(2,3,1)、(3,1,2) 中对应的循环排列。这里 s,t 是本节的腿标签，不是 Wilson 轮廓参数。

所有坐标微分先在独立的 X1、X2、X3 上完成；切向量 u1、u2、u3 在微分时固定。


### 1.1 独立可读的算符、几何与约定

本文采用(+---)度规、$D_0=-ig/(k^2+i0)$、$e^{-ikX}$及正号Dyson插入。令$q=-k^2-i0$，则

$$
\begin{aligned}
\Pi_1^{ab,\mu\nu}(k)&=i\delta^{ab}\mathcal A_1(\varepsilon)
(k^2g^{\mu\nu}-k^\mu k^\nu)q^{-\varepsilon},\\
\mathcal A_1(\varepsilon)&=\frac{g_0^2}{2(4\pi)^{2-\varepsilon}}
\frac{\Gamma(\varepsilon)\Gamma(1-\varepsilon)^2}{\Gamma(2-2\varepsilon)}
\mathcal B_1(\varepsilon),\\
\mathcal P_1(k^2)&=-\mathcal A_1(\varepsilon)q^{-\varepsilon}.
\end{aligned}
$$

最后一行连接[2.0共同约定](../conventions.md)，不是再翻号；$\mathcal B_1$见§6.2，包含胶子、鬼和无质量夸克的总一圈系数。裸$g_0$有质量维ε。

$X^{(e)}(s)=\sum_{k<e}v_k+s v_e$，$u_i=v_{e_i}$。$v_i^2=0$、$\sum_i v_i=0$；$v_{14}=v_{23}=v_{\rm ref}$、$v_{12}=v_{34}=zv_{\rm ref}$、$v_{13}=v_{24}=-(1+z)v_{\rm ref}$。$v_{ij}$是点积，没有额外2。全篇统一用$A_0=-2v_{\rm ref}$；先在$A_0>0,z>0$的解析区域提取尺度。

所有非整数幂由同一个Feynman解析边界值定义。在一般非重合坐标及可交换积分/微分的收敛域建立恒等式，再作维数与必要辅助指数的亚纯延拓，最后施加类光多边形。不能把单独发散项当作普通Lebesgue积分交换次序；也不据此删除端点分布。本实现不做轮廓IBP，因此不需要凭空增加或删除轮廓接触项。

## 2. 颜色、相位及完整动量收缩

由反对称性及生成元对易关系，

$$
f^{ABC}t^At^B=\frac12f^{ABC}[t^A,t^B]=\frac{iC_A}{2}t^C,
\qquad
\frac1{N_c}f^{ABC}\operatorname{Tr}(t^Ct^Bt^A)=-\frac i2C_AC_F.
$$

因此三个 Wilson 顶点与有向颜色迹的乘积为

$$
(ig_0)^3\left(-\frac i2C_AC_F\right)=-\frac{g_0^3}{2}C_AC_F.
$$

令

$$
\mathcal G_0:=-\frac{g_0^4 C_AC_F}{2}\mathcal A_1(\varepsilon),
$$

并保留完整动量测度

$$
d\Lambda_X=\prod_{j=1}^3\frac{d^D\ell_j}{(2\pi)^D}
(2\pi)^D\delta^{(D)}(\ell_1+\ell_2+\ell_3)
e^{-i\sum_j\ell_j\cdot X_j}.
$$

定义

$$
\mathcal L_\ell=
\sum_{(i,j,k)=(1,2,3),(2,3,1),(3,1,2)}
(u_i\cdot u_j)\,u_k\cdot(\ell_i-\ell_j),
$$

$$
\mathcal W_r=
-(u_s\cdot\ell_r)(u_t\cdot\ell_s)
+(u_t\cdot\ell_r)(u_s\cdot\ell_t)
+(u_s\cdot u_t)\ell_r\cdot(\ell_s-\ell_t).
$$

其中 W_r 是将 Γkin(ℓ) 的第 r 个切向量换成 ℓr 后的完整收缩。于是 Γkin(−ℓ) 与自能分子的收缩为

$$
\mathcal N_r^{\rm mom}=-\ell_r^2\mathcal L_\ell
+(u_r\cdot\ell_r)\mathcal W_r.
$$

原图严格变为

$$
\boxed{
\mathcal A_{a,b,c}^{[r]}=
i\mathcal G_0\int_{\mathcal D_s}d^3s\int d\Lambda_X\,
\frac{(-\ell_r^2-i0)^{-\varepsilon}\mathcal N_r^{\rm mom}}
{(\ell_r^2+i0)^2\prod_{j\ne r}(\ell_j^2+i0)}.
}
$$

外面的 i 来自完整插入腿与另外两条裸腿的乘积。自能中的颜色及内部积分已经全部包含在 A1 中。

利用动量守恒，亦可检查

$$
\mathcal W_r=(u_s\cdot u_t)(\ell_t^2-\ell_s^2)
+(u_s\cdot\ell_s)(u_t\cdot\ell_s)
-(u_s\cdot\ell_t)(u_t\cdot\ell_t).
$$

这不是零。本文不通过轮廓 IBP 删除该部分。

## 3. 完整坐标空间传播子

记 Δ(X)=−X²+i0，并定义

$$
c_0=\frac{\Gamma(1-\varepsilon)}{4\pi^{2-\varepsilon}},\quad
c_A=\frac{\Gamma(1-2\varepsilon)}
{4^{1+\varepsilon}\pi^{2-\varepsilon}\Gamma(1+\varepsilon)},\quad
c_B=\frac{\Gamma(-2\varepsilon)}
{4^{2+\varepsilon}\pi^{2-\varepsilon}\Gamma(2+\varepsilon)}.
$$

由输入的 Fourier 母式直接得到

$$
D_{0,\mu\nu}^{AB}(X)=-\delta^{AB}c_0g_{\mu\nu}\Delta(X)^{-\beta}.
$$

设 q=−k²−i0。在 Feynman 边界值的解析延拓意义下，

$$
\mathcal D_{\mu\nu}^{AB}(k)
=i\delta^{AB}\mathcal A_1
\left[g_{\mu\nu}q^{-1-\varepsilon}+k_\mu k_\nu q^{-2-\varepsilon}\right].
$$

因此

$$
\mathcal D_{\mu\nu}^{AB}(X)
=-\delta^{AB}\mathcal A_1c_Ag_{\mu\nu}\Delta^{-h}
+\delta^{AB}\mathcal A_1c_B\partial_\mu\partial_\nu\Delta^{2\varepsilon}.
$$

所需导数及系数关系为

$$
\partial_\mu\partial_\nu\Delta^{2\varepsilon}
=-4\varepsilon g_{\mu\nu}\Delta^{2\varepsilon-1}
+8\varepsilon(2\varepsilon-1)X_\mu X_\nu\Delta^{2\varepsilon-2},
\qquad \frac{c_B}{c_A}=-\frac1{8\varepsilon(1+\varepsilon)}.
$$

最终完整核是

$$
\boxed{
\mathcal D_{\mu\nu}^{AB}(X)=
-\delta^{AB}\mathcal A_1\frac{c_A}{2(1+\varepsilon)}
\left[(1+2\varepsilon)g_{\mu\nu}\Delta^{-h}
-2hX_\mu X_\nu\Delta^{-h-1}\right].
}
$$

第二项不能由“传播子改幂”获得；第一项的系数本身也包含动量 kk 项 Fourier 后的贡献。

### 3.1 横向性

$$
\partial^\mu(g_{\mu\nu}\Delta^{-h})=2hX_\nu\Delta^{-h-1},
$$

$$
\partial^\mu(X_\mu X_\nu\Delta^{-h-1})
=(D-2h-1)X_\nu\Delta^{-h-1}
=(1+2\varepsilon)X_\nu\Delta^{-h-1}.
$$

故完整核满足 ∂μDμν=0。该关系继承自原动量核的横向性，并按同一解析延拓解释为分布恒等式。

## 4. 中心 y 积分与两类参数核

定义

$$
Z_j=X_j-y,\qquad \Delta_j=-Z_j^2+i0,
$$

$$
T_{\mu\nu}(Z)=\frac1{2(1+\varepsilon)}
\left[(1+2\varepsilon)g_{\mu\nu}\Delta(Z)^{-h}
-2hZ_\mu Z_\nu\Delta(Z)^{-h-1}\right],
$$

$$
\mathscr K_{\mu\nu}^{[r]}(X)=\int d^Dy\,
T_{\mu\nu}(Z_r)\prod_{j\ne r}\Delta_j^{-\beta}.
$$

对于任何张量 Bμν，定义三胶子坐标微分的完整收缩

$$
\begin{aligned}
\mathscr V_r[B]={}&
u_t\cdot(\partial_r-\partial_s)(u_r^\mu u_s^\nu B_{\mu\nu})\\
&+(u_s\cdot u_t)(\partial_s-\partial_t)^\nu(u_r^\mu B_{\mu\nu})\\
&+u_s\cdot(\partial_t-\partial_r)(u_r^\mu u_t^\nu B_{\mu\nu}).
\end{aligned}
$$

于是

$$
\mathcal A^{[r]}=i\mathcal G_0c_0^2c_A
\int_{\mathcal D_s}d^3s\,\mathscr V_r[\mathscr K^{[r]}].
$$

其相位来自 Γkin(−ℓ)=−iΓkin(∂)，以及三条坐标腿的总负号。

### 4.1 标量星型积分

$$
d\mu_\alpha=d\alpha_1d\alpha_2d\alpha_3\,
\delta(1-\alpha_1-\alpha_2-\alpha_3),\qquad \alpha_i\ge0,
$$

$$
\bar X=\sum_i\alpha_iX_i,\quad R_i=X_i-\bar X,\quad
Q=-\sum_{i<j}\alpha_i\alpha_j(X_i-X_j)^2+i0.
$$

因为

$$
\sum_i\alpha_i\Delta_i=-(y-\bar X)^2+Q,
$$

对 Σ=a1+a2+a3 有

$$
\boxed{
\int d^Dy\prod_i\Delta_i^{-a_i}
=-i\pi^d\frac{\Gamma(\Sigma-d)}{\prod_i\Gamma(a_i)}
\int d\mu_\alpha\prod_i\alpha_i^{a_i-1}\,Q^{d-\Sigma}.
}
$$

这里的相位是 −i。坐标核为 −y²+i0，中心时间分量沿 y0=−iyE0 转动；它不能照搬输入中动量 Fourier 母式的 +i。

### 4.2 两类实际需要的权重

令

$$
w_r=\alpha_r^{-2\varepsilon}\prod_{j\ne r}\alpha_j^{-\varepsilon}.
$$

度规项的指数是 ar=h、其余 aj=β，其标量星积分含 Γ(ν)Q^(−ν)，权重为 wr。

张量项的指数是 ar=h+1、其余 aj=β，其标量星积分含 Γ(m)Q^(−m)，权重为 αr wr。两个裸腿的权重不变，插入腿权重与树图不同。

张量矩积分由 Zr=Rr−(y−Xbar) 得到：

$$
\begin{aligned}
&\int d^Dy\,Z_{r\mu}Z_{r\nu}\Delta_r^{-h-1}
\prod_{j\ne r}\Delta_j^{-\beta}\\
&=-\frac{i\pi^d}{\Gamma(h+1)\Gamma(\beta)^2}
\int d\mu_\alpha\,\alpha_rw_r
\left[\Gamma(m)\frac{R_{r\mu}R_{r\nu}}{Q^m}
-\frac{g_{\mu\nu}}2\Gamma(\nu)\frac1{Q^\nu}\right].
\end{aligned}
$$

其中度规矩的 −1/2 可直接由

$$
\int d^D\eta\,\eta_\mu\eta_\nu(Q-\eta^2)^{-\Sigma}
=\frac{g_{\mu\nu}}D
\left[Q\int d^D\eta(Q-\eta^2)^{-\Sigma}
-\int d^D\eta(Q-\eta^2)^{1-\Sigma}\right]
$$

推出。组合两项后，中心核为

$$
\boxed{
\mathscr K_{\mu\nu}^{[r]}
=-\frac{i\pi^d\Gamma(\nu)}
{2(1+\varepsilon)\Gamma(h)\Gamma(\beta)^2}
\int d\mu_\alpha\,w_r
\left[(1+2\varepsilon+\alpha_r)\frac{g_{\mu\nu}}{Q^\nu}
-2\nu\alpha_r\frac{R_{r\mu}R_{r\nu}}{Q^{\nu+1}}\right].
}
$$

这里的 +αr gμν 正是张量矩产生的度规贡献。

### 4.3 二阶导数路线的交叉检查

也可从 cB∂rμ∂rνΔr^(2ε) 出发。其标量星指数为 ar=−2ε、其余 aj=β，权重为 αr^(−1−2ε)∏j≠rαj^(−ε)，Q 幂为 Q^(3ε)。该中间表达按解析延拓定义。

利用

$$
\partial_{r\mu}\partial_{r\nu}Q^{-q}
=2q\alpha_r(1-\alpha_r)g_{\mu\nu}Q^{-q-1}
+4q(q+1)\alpha_r^2R_{r\mu}R_{r\nu}Q^{-q-2}
$$

并取 q=−3ε，得到与上面完全相同的中心核。最终实际实现不必保留这个负指数中间星积分。

## 5. 将剩余坐标导数全部作用完

定义 Pi=αiRi。基础导数为

$$
\partial_{i\mu}Q=-2P_{i\mu},\qquad
\partial_{i\mu}R_{r\nu}=(\delta_{ir}-\alpha_i)g_{\mu\nu}.
$$

引入三个已经没有微分算子的标量：

$$
\mathscr L=\sum_{(i,j,k)\,\mathrm{cyclic}}
(u_i\cdot u_j)u_k\cdot(P_i-P_j),
$$

$$
\begin{aligned}
\mathscr E_r={}&(1-\alpha_r+\alpha_s)(u_r\cdot u_t)(u_s\cdot R_r)\\
&-(1-\alpha_r+\alpha_t)(u_r\cdot u_s)(u_t\cdot R_r)\\
&+D(\alpha_t-\alpha_s)(u_s\cdot u_t)(u_r\cdot R_r),
\end{aligned}
$$

$$
\begin{aligned}
\mathscr F_r=(u_r\cdot R_r)\big[&
-\alpha_s(u_s\cdot R_r)(u_t\cdot R_s)
+\alpha_t(u_t\cdot R_r)(u_s\cdot R_t)\\
&+(u_s\cdot u_t)R_r\cdot(\alpha_sR_s-\alpha_tR_t)\big].
\end{aligned}
$$

D=4−2ε 必须保留在 Er 中；它来自对张量分子的导数收缩。

将上一节中心核中方括号记为 Bμν，则实际微分结果为

$$
\boxed{
\mathscr V_r[B]
=2\nu\left[
\frac{(1+2\varepsilon+\alpha_r)\mathscr L-\alpha_r\mathscr E_r}{Q^m}
-2m\alpha_r\frac{\mathscr F_r}{Q^{m+1}}\right].
}
$$

这个式子对任意独立坐标和固定切向量成立，尚未使用类光或同边条件。

相位相乘为：三个坐标腿给 −1，顶点微分给 −i，中心 Wick 转动给 −i，总乘积是 +1。因此

$$
\boxed{
\begin{aligned}
\mathcal A^{[r]}={}&
\frac{\mathcal G_0\Gamma(m)}
{4^{3+\varepsilon}\pi^{4-2\varepsilon}\Gamma(2+\varepsilon)}
\int_{\mathcal D_s}d^3s\int d\mu_\alpha\,w_r\\
&\times\left[
\frac{(1+2\varepsilon+\alpha_r)\mathscr L-\alpha_r\mathscr E_r}{Q^m}
-2m\alpha_r\frac{\mathscr F_r}{Q^{m+1}}\right].
\end{aligned}
}
$$

## 6. 无量纲化与最终 Csrc/J

先在 A0=−2vref>0 的起始解析区域定义

$$
q_{ij}=-\frac{(X_i-X_j)^2}{A_0},\quad
q_0=\sum_{i<j}\alpha_i\alpha_jq_{ij},\quad
\widehat Q=q_0+i0,
$$

$$
G_{ij}=\frac{u_i\cdot u_j}{A_0},\qquad
T_{ij}=\frac{u_i\cdot R_j}{A_0},\qquad
W_{ij}=\frac{R_i\cdot R_j}{A_0}.
$$

同样定义 L=𝓛/A0²、Er=𝓔r/A0²、Fr=𝓕r/A0³，即

$$
L=\sum_{(i,j,k)\,\mathrm{cyclic}}
G_{ij}(\alpha_iT_{ki}-\alpha_jT_{kj}),
$$

$$
E_r=(1-\alpha_r+\alpha_s)G_{rt}T_{sr}
-(1-\alpha_r+\alpha_t)G_{rs}T_{tr}
+D(\alpha_t-\alpha_s)G_{st}T_{rr},
$$

$$
F_r=T_{rr}\left[-\alpha_sT_{sr}T_{ts}
+\alpha_tT_{tr}T_{st}
+G_{st}(\alpha_sW_{rs}-\alpha_tW_{rt})\right].
$$

这些完全是 z、Wilson 参数及 Feynman 参数的多项式。

### 6.1 最终几何积分

固定本文 J 的归一化为

$$
\boxed{
\begin{aligned}
J_{a,b,c}^{[r]}(\varepsilon,z)
={}&\int_{\mathcal D_s}d^3s\int d\mu_\alpha\,
\alpha_r^{-2\varepsilon}\prod_{j\ne r}\alpha_j^{-\varepsilon}\\
&\times\left[
\frac{(1+2\varepsilon+\alpha_r)L-\alpha_rE_r}{\widehat Q^{2-3\varepsilon}}
-\frac{2(2-3\varepsilon)\alpha_rF_r}{\widehat Q^{3-3\varepsilon}}
\right].
\end{aligned}
}
$$

亦可写成一个共同分母：

$$
\boxed{
J^{[r]}=\int_{\mathcal D_s}d^3s\int d\mu_\alpha\,w_r
\frac{N_r}{\widehat Q^{3-3\varepsilon}},\quad
N_r=\widehat Q[(1+2\varepsilon+\alpha_r)L-\alpha_rE_r]
-2(2-3\varepsilon)\alpha_rF_r.
}
$$

两个分母幂的合并是代数通分，没有作 IBP。它不是只保留原树顶点分子的单核假设。

### 6.2 源系数

$$
\boxed{
\mathcal A_{a,b,c}^{[r]}=C_{\rm src}(\varepsilon,v_{\rm ref})
J_{a,b,c}^{[r]}(\varepsilon,z).
}
$$

直接保留输入 A1 的形式是

$$
C_{\rm src}=
-\frac{g_0^4C_AC_F\mathcal A_1}{2}
\frac{\Gamma(2-3\varepsilon)}
{4^{3+\varepsilon}\pi^{4-2\varepsilon}\Gamma(2+\varepsilon)}
A_{\rm ref}^{3\varepsilon}.
$$

代入输入完整裸自能后，

$$
\boxed{
C_{\rm src}=
-\frac{g_0^6C_AC_F\mathcal B_1(\varepsilon)}
{4096\pi^{6-3\varepsilon}}
\frac{\Gamma(1-\varepsilon)^2\Gamma(2-3\varepsilon)}
{\varepsilon(1+\varepsilon)\Gamma(2-2\varepsilon)}
A_{\rm ref}^{3\varepsilon}.
}
$$

这里

$$
\mathcal B_1(\varepsilon)=
\frac{(10-6\varepsilon)C_A-8(1-\varepsilon)T_Fn_f}{3-2\varepsilon}.
$$

数值及 Gamma 化简逐项为

$$
-\frac1{2}\times\frac1{2(4\pi)^{2-\varepsilon}}
\times\frac1{4^{3+\varepsilon}\pi^{4-2\varepsilon}}
=-\frac1{4096\pi^{6-3\varepsilon}},
\qquad
\frac{\Gamma(\varepsilon)}{\Gamma(2+\varepsilon)}
=\frac1{\varepsilon(1+\varepsilon)}.
$$

若 g0=μ^ε gb、αs,b=gb²/(4π)，则

$$
\boxed{
C_{\rm src}=
-C_AC_F\mathcal B_1(\varepsilon)
\left(\frac{\alpha_{s,b}}\pi\right)^3
\frac{\Gamma(1-\varepsilon)^2\Gamma(2-3\varepsilon)}
{64\varepsilon(1+\varepsilon)\Gamma(2-2\varepsilon)}
(\pi\mu^2A_{\rm ref})^{3\varepsilon}.
}
$$

μ^(6ε) 显式保留在最后的尺度因子内。

上述因子分解先在 A0>0、z>0 建立，然后按原 Feynman 处方整体解析延拓。继续到其他区域时，不能分别对 A0 的幂和 Qhat 的幂重新选择独立主值分支；其支路必须继承未提取尺度时的 Q。全文没有将这些幂换成绝对值。

### 6.3 Euler 因子

Csrc 中没有额外乘入 Euler 因子。Gamma 组合本身的对数含 +3γEε。

若另定义 prepared 量为 Aprepared=e^(−3γEε)Asrc，则 Cprepared=e^(−3γEε)Csrc，J 不变。若要求原始 A 不变而仅重新分配两个因子，则必须同时取 C′=e^(−3γEε)Csrc、J′=e^(+3γEε)J。不能只修改 C 又仍宣称乘积是同一原始 A。

## 7. 仅用标量多项式实现全部十二图

边向量的无量纲 Gram 矩阵为

$$
\mathbb G_{ij}=\frac{v_i\cdot v_j}{A_0}
=-\frac12\begin{pmatrix}
0&z&-(1+z)&1\\
z&0&1&-(1+z)\\
-(1+z)&1&0&z\\
1&-(1+z)&z&0
\end{pmatrix}.
$$

若 Xi=Σk cik vk、ui=v_ei，则 c_i 由 Vi+si vei 直接读出。

$$
\bar c=\sum_i\alpha_i c_i,\quad R_i\leftrightarrow c_i-\bar c,
$$

$$
G_{ij}=\mathbb G_{e_i e_j},\quad
T_{ij}=\sum_k\mathbb G_{e_i k}(c_{jk}-\bar c_k),
$$

$$
q_{ij}=-(c_i-c_j)^T\mathbb G(c_i-c_j).
$$

也可完全免去 Ri·Rj 的显式向量计算：

$$
d_i=\sum_k\alpha_kq_{ik},\qquad
\boxed{W_{ij}=q_0+\frac{q_{ij}-d_i-d_j}{2}.}
$$

这里 q0 不含 i0；W 是几何多项式。

令三个附着参数统一按腿顺序暂记为 s1,s2,s3，并定义

$$
D_z(s,t)=(1-s)(1-t)+zst.
$$

| 源图 | 边序列 | r（依源图顺序） | q12 | q13 | q23 | 原排序约束 |
|---|---|---|---|---|---|---|
| 01–03 | (3,4,4) | 3,2,1 | z(1−s1)s2 | z(1−s1)s3 | 0 | s2≤s3 |
| 04–06 | (3,3,4) | 3,2,1 | 0 | z(1−s1)s3 | z(1−s2)s3 | s1≤s2 |
| 07–09 | (2,4,4) | 3,2,1 | Dz(s1,s2) | Dz(s1,s3) | 0 | s2≤s3 |
| 10–12 | (2,3,4) | 3,2,1 | (1−s1)s2 | Dz(s1,s3) | z(1−s2)s3 | 无同边约束 |

所有参数仍在 [0,1] 内。改变 r 只改变标记权重及 Er、Fr，不改变附着点顺序。

## 8. 01 图完整代入例子

记

$$
s=s_{3,1},\quad t=s_{4,1},\quad u=s_{4,2},\quad x=1-s,
\qquad b_\alpha=\alpha_2t+\alpha_3u.
$$

以 V4 为坐标原点，

$$
X_1=-xv_3,\quad X_2=tv_4,\quad X_3=uv_4,
\qquad (u_1,u_2,u_3)=(v_3,v_4,v_4),\quad r=3.
$$

保持原域

$$
\int_{\mathcal D_s}d^3s=
\int_0^1ds\int_0^1du\int_0^u dt.
$$

此处 x 只是 1−s 的缩写，并未改变积分变量。

$$
q_{12}=zxt,\quad q_{13}=zxu,\quad q_{23}=0,
\qquad q_0=z\alpha_1xb_\alpha.
$$

代入通式得到三个显式分子：

$$
\boxed{L_{01}=\frac{z^2\alpha_1x}{4}(\alpha_3-\alpha_2),}
$$

$$
\boxed{E_{01}=\frac{z^2\alpha_1x}{4}
\left[D(\alpha_2-\alpha_1)-(1-\alpha_3+\alpha_2)\right],}
$$

$$
\boxed{F_{01}=\frac{z^3\alpha_1^2x^2}{8}
\left[(\alpha_1-\alpha_2)b_\alpha+\alpha_2u\right].}
$$

因此

$$
\begin{aligned}
J_{01}={}&\int_0^1ds\int_0^1du\int_0^u dt
\int d\mu_\alpha\,
\alpha_1^{-\varepsilon}\alpha_2^{-\varepsilon}\alpha_3^{-2\varepsilon}\\
&\times\left[
\frac{(1+2\varepsilon+\alpha_3)L_{01}-\alpha_3E_{01}}{(z\alpha_1xb_\alpha+i0)^m}
-\frac{2m\alpha_3F_{01}}{(z\alpha_1xb_\alpha+i0)^{m+1}}
\right].
\end{aligned}
$$

虽然 u2·u3=0 且 q23=0，E01、F01 一般都不为零。这也展示了为什么必须在施加同边条件前完成坐标微分。

## 9. 10 图完整代入例子

记 x=1−s2,1、w=1−s3,1、t=s4,1。将两项 1−s 的变量替换实际用于积分时，Jacobian 绝对值各为 1，反转的上下限各自抵消负号，因此域仍是单位立方体。

以 V4 为原点，

$$
X_1=-xv_2-v_3,\quad X_2=-wv_3,\quad X_3=tv_4,
\qquad (u_1,u_2,u_3)=(v_2,v_3,v_4),\quad r=3.
$$

$$
q_{12}=x(1-w),\quad q_{13}=x(1-t)+z(1-x)t,\quad q_{23}=zwt,
$$

$$
q_0=\alpha_1\alpha_2x(1-w)
+\alpha_1\alpha_3[x(1-t)+z(1-x)t]
+\alpha_2\alpha_3zwt.
$$

令 A=α1x、B=α1+α2w、C=α3t，则

$$
G=\frac12\begin{pmatrix}0&-1&1+z\\-1&0&-z\\1+z&-z&0\end{pmatrix},
$$

$$
T=\frac12\begin{pmatrix}
1-B-(1+z)C&w-B-(1+z)C&-B+(1+z)(t-C)\\
x-A+zC&-A+zC&-A-z(t-C)\\
(1+z)(A-x)-z(B-1)&(1+z)A-z(B-w)&(1+z)A-zB
\end{pmatrix}.
$$

因此，完全显式的标量分子可写成

$$
L_{10}=\frac12\left[-(\alpha_1T_{31}-\alpha_2T_{32})
-z(\alpha_2T_{12}-\alpha_3T_{13})
+(1+z)(\alpha_3T_{23}-\alpha_1T_{21})\right],
$$

$$
E_{10}=-\frac12\left[
z(1-\alpha_3+\alpha_1)T_{13}
+(1+z)(1-\alpha_3+\alpha_2)T_{23}
+D(\alpha_2-\alpha_1)T_{33}\right],
$$

$$
F_{10}=T_{33}\left[-\alpha_1T_{13}T_{21}
+\alpha_2T_{23}T_{12}
-\frac12(\alpha_1W_{31}-\alpha_2W_{32})\right],
$$

其中 W31、W32 由上一节的距离多项式公式给出。这是有限次数多项式的分块写法，不含未作用的坐标导数。

$$
J_{10}=\int_0^1dx\int_0^1dw\int_0^1dt\int d\mu_\alpha\,
\alpha_1^{-\varepsilon}\alpha_2^{-\varepsilon}\alpha_3^{-2\varepsilon}
\left[\frac{(1+2\varepsilon+\alpha_3)L_{10}-\alpha_3E_{10}}{(q_0+i0)^m}
-\frac{2m\alpha_3F_{10}}{(q_0+i0)^{m+1}}\right].
$$

## 10. 排序域与可选超立方映射

本文没有做轮廓 IBP，因此没有删除端点或同边接触面。前面的导数均在独立坐标上作用于传播子核，未作用于排序阶跃函数。

若实现只接受单位立方体，可对同边的一对槽位取

$$
s_{e,2}=b,\quad s_{e,1}=ab,\quad 0\le a,b\le1,
\qquad ds_{e,1}ds_{e,2}=b\,da\,db.
$$

接触面 s_e,1=s_e,2 对应 a=1，不能删去。

对 Feynman 单纯形可取

$$
\alpha_1=\xi,\quad\alpha_2=(1-\xi)\eta,\quad
\alpha_3=(1-\xi)(1-\eta),\qquad
 d\mu_\alpha=(1-\xi)d\xi d\eta.
$$

若有更多同边槽位，映射 s_e,k=∏j=k^me t_j 的 Jacobian 为 ∏j=2^me t_j^(j−1)。本任务最多只有两个同边点。

## 11. 两个颜色通道

定义

$$
b_g(\varepsilon)=\frac{10-6\varepsilon}{3-2\varepsilon},\qquad
b_f(\varepsilon)=-\frac{8(1-\varepsilon)}{3-2\varepsilon}.
$$

把 Csrc 写为 CA CF B1 Kcommon，则

$$
\boxed{\mathcal A^{[r]}=
C_FC_A^2\,[b_gK_{\rm common}J^{[r]}]
+C_FC_AT_Fn_f\,[b_fK_{\rm common}J^{[r]}].}
$$

几何 J 完全相同，但

$$
b_g=\frac{10}{3}+\frac{2}{9}\varepsilon+O(\varepsilon^2),\qquad
b_f=-\frac83+\frac89\varepsilon+O(\varepsilon^2).
$$

若 H=Kcommon J=Σn Hn ε^n，则颜色通道 χ 的第 n 阶系数是 Σk≥0 bχ,k Hn−k，必须保留卷积。

## 12. 检查与实现状态

### 耦合及尺度

A1 提供 g0²，Wilson 顶点及中央树顶点提供 g0⁴，总阶数 g0⁶。

坐标和切向量一同缩放 λ 时，Q→λ²Q、𝓛和𝓔→λ⁴、𝓕→λ⁶，两项均缩放 λ^(6ε)，与 A0^(3ε) 一致。固定 u 而仅平移 Xi 不影响 Ri、Q 或所有分子。

### 置换

循环重标全部腿、位置、颜色和标记 r 时，核协变。奇置换使有向 Lorentz 分子反号；颜色结构和 Wilson 生成元顺序也须同步处理。不能把 r 移到第一腿而继续沿用原颜色顺序的所有符号。

### 动量度规项控制

只保留原自能的 k²gμν 项时，插入腿坐标核为 −A1 cA gμνΔ^(−h)。使用本文同一个 Csrc，控制积分为

$$
\boxed{J_g^{[r]}=2(1+\varepsilon)
\int_{\mathcal D_s}d^3s\int d\mu_\alpha\,w_r
\frac{L}{\widehat Q^m}.}
$$

因而完整张量相对该控制的差为

$$
J^{[r]}-J_g^{[r]}=
-\int_{\mathcal D_s}d^3s\int d\mu_\alpha\,w_r
\left[\frac{(1-\alpha_r)L+\alpha_rE_r}{\widehat Q^m}
+\frac{2m\alpha_rF_r}{\widehat Q^{m+1}}\right].
$$

这一区别还说明“完整坐标核中显式 gμν 那一项”不等于“原动量投影中单独的 k²gμν 项”。

### 本地实际执行的复核

独立构造了12图的Gram几何、距离核与W恒等式，核对01例子的L/E/F；另从坐标分量对中心张量直接微分，在D=3、4、5及各标记腿作9个一般位置精确样点。它们均与§5的显式分子相符。加上Gamma递推、相位、横向系数和两项错误控制，连同10图T矩阵、第二导数路线及64个C4 Gram条目，共179项通过。

这9个样点不是任意坐标恒等式的形式证明；通式的依据仍是§5的乘积求导。可直接看出E的D项来自散度：对R_r^μR_r^ν求∂_i,ν得到(D+1)(δ_ir−α_i)R_r^μ；另两条顶点项合并抵消其中的1，留下D(α_t−α_s)。F来自分母求导后的二次张量收缩。维数D不得提前改成4。

网页原稿声明有随附脚本，本次只收到Markdown，未获得或运行那些脚本；上面的结果均来自本地另写检查。未运行完整积分数值比较、Laurent展开、RG或Lean/MCP证明，formal=not_run。

## 13. Montopar实现契约

### 13.1 分派与输入拒绝条件

新增专用`tgv_pi1`路径，精确匹配`3e-3g+Pi1`。不得复用仅支持两点骨架的旧self_energy分派，也不得落入中央Gamma1分支。源验收须检查：

- 总三圈、g_power=6、三个Wilson端点；内部恰一个树TGV和一个裸Pi1，标记腿r∈{1,2,3}。
- 四段传播子；标记腿两段共用l_r，通过Pi1连通；另两腿各一段。重建连通与指标配对，不仅信任marked_contour_position。
- 中央动量为−l₁、−l₂、−l₃；三相位对应各端点，唯一中心δ携带(2π)^D归一化。错号或丢δ应拒绝。
- sign_convention为2.0、source_external_sign=+1、Dyson系数1；旧源翻号适配不可再次施加。未知约定应拒绝，不能猜测。
- 同边槽位按源序保留，核对bounds及iterated_order_outer_to_inner。拒绝缺排序链、循环依赖、槽位与相位不一致的输入。

01–12与r的对应取§7及源元数据；不从编号推断物理规则。编号表仅用于回归。允许颜色或指标重命名，但必须保留完整连通语义。

### 13.2 主输出：一个C、一个J、两个核项

raw必须保持以下语义；字段名称在实现时接入现有输出模型，本文不提前占用schema版本号。

| 内容 | 要求 |
|---|---|
| 来源 | 源文件/对象哈希、orbit、符号版本、标记腿、算符圈数1、总圈数3 |
| C_src | §6.2的全色全Gamma表达式，含g₀⁶与A₀^(3ε)，无Euler补偿、无C4/镜像重数 |
| 原Wilson变量 | 源变量名、依赖bounds、外到内次序逐字段保留 |
| Feynman域 | α₁、α₂、α₃非负且和为1，权重α_r^(−2ε)其余α^(−ε) |
| 几何 | q₀、L、E_r、F_r为有理系数的z/ε/s/α多项式，D已替为4−2ε |
| 第一核项 | K_r=(1+2ε+α_r)L−α_rE_r，分母(q₀+i0)^(2−3ε) |
| 第二核项 | −2(2−3ε)α_rF_r，分母(q₀+i0)^(3−3ε) |
| 颜色 | 完整CF CA B₁及CF CA²、CF CA TF nf两个通道；B₁只乘一次 |
| 状态 | raw参数化，不是已求I；无旧族输出的隐式变化 |

C在两项之外只乘一次，两项共同定义一个J，不生成两个独立物理图。若消费者要求一个分母，可导出派生形式，但必须保存处方：

$$
N_r(\delta)=(q_0+i\delta)K_r-2m\alpha_rF_r,\qquad
J=\lim_{\delta\downarrow0}\int w_r\frac{N_r(\delta)}{(q_0+i\delta)^{m+1}}.
$$

若只存实多项式$N_r(0)$，还必须另存$i\delta K_r$及其边界值规则，不能默默丢弃后声称与两项式严格相同。主实现选择两项式，可避免这个额外义务。实现不得把Feynman i0当成可约掉的普通浮点小数。

若输出接口要求独立积分变量，先消去α₃=1−α₁−α₂，定义α₁外层0..1，α₂内层0..1−α₁，Jacobian=1；必须对权重和L/E/F同步替换。由三个Wilson变量加两个单纯形变量，共五个独立参数。将α积分置于原Wilson嵌套之内是本表示的明示定义，不重排原Wilson链。§10立方体版本仅是可选派生表示，必须存替换及两个Jacobian，不替换ordered主对象。

### 13.3 raw、prepared与唯一I的关系

定义共同无色系数

$$
R_{\rm src}(\varepsilon)=
-\frac{\Gamma(1-\varepsilon)^2\Gamma(2-3\varepsilon)}
{64\varepsilon(1+\varepsilon)\Gamma(2-2\varepsilon)},\qquad
R(\varepsilon)=e^{-3\gamma_E\varepsilon}R_{\rm src}(\varepsilon).
$$

令$\lambda_b=\alpha_{s,b}/\pi$、$S_{14}=(\pi\mu^2A_0)^{3\varepsilon}$。源代表图满足

$$
\mathcal A_r=\lambda_b^3 C_AC_F\mathcal B_1 S_{14}R_{\rm src}J_r.
$$

prepared代表图为原源乘一次$e^{-3\gamma_E\varepsilon}$，保留original_C。它与原始裸图是显式归一化转换，不是声称改变C而乘积不变；若恢复原图需乘回逆转换。raw阶段不能内置该补偿。

这12个C4轨道各有四个旋转代表。将v₁₂与v₁₄互换的两个旋转单独计入，设$S_{12}=z^{3\varepsilon}S_{14}$，则本分支沿现有项目的对称尺度提取方式定义

$$
\begin{aligned}
\mathcal W_r^{\rm prepared}
&=2\lambda_b^3 C_AC_F\mathcal B_1 S_{14}R
[J_r(\varepsilon,z)+z^{3\varepsilon}J_r(\varepsilon,z^{-1})]\\
&=\lambda_b^3 C_AC_F\mathcal B_1(S_{14}+S_{12})I_r(\varepsilon,z),\\
I_r(\varepsilon,z)
&=2R\frac{J_r(\varepsilon,z)+z^{3\varepsilon}J_r(\varepsilon,z^{-1})}
{1+z^{3\varepsilon}}.
\end{aligned}
$$

这是后续组装的唯一I定义，非raw Montopar额外输出的旋转副本。旋转时必须同步移动边标记和附着槽位；不能用镜像代替C4旋转。此处采用冻结源C4商的每轨道4重数，不把源中canonical diagram计数1误当物理重数。实际组装实现还须重放四个旋转验证该缩放关系。

两个单项式颜色通道的归一化系数为$b_g I_r$与$b_f I_r$，不是两个不同的几何积分。分味道做精确颜色重构后，剩余方向分别为CF CA²和CF CA TF nf，剥离B₁后的首坐标均为1。求有限项时与$b_g,b_f$的全ε展开卷积，不预设本族极点阶数。

### 13.4 编码顺序与验收

1. 新分派、拓扑/源字段验证，非法样例必须拒绝。
2. 复用现有多边形Gram/位置构造器，生成q₀/G/T/W并按§5–7得到L/E/F；不硬编码12图长分子。
3. 输出有序两项J与可乘C账本；保留完整颜色、Gamma、尺度和符号。可用本地[几何夹具](../../../workflow/checks/tgv-pi1-parametric/geometry.json)逐图交叉比对，但夹具不是生产导出。
4. 回归12图、01和10显式例子、三种插入腿；验证D项、张量项、相位以及两颜色。增加丢kk、提前D=4、重复B₁、额外负号、误用树权重、错序和重复Euler的负控制。
5. 对旧pure/TGV/FGV/mixed/double-TGV与Pi1/chain/Pi2执行相关生产回归；旧输出应保持。再做实际CLI导出并审阅C/J。

本地理论审计支持开始编码；尚未取得生产回归或实际CLI结果。无新积分求值、无新Laurent系数，未修改Kernel或现行15个三圈自能I。原网页稿、旧源与旧理论留存。

## 14. IBP形式：下一阶段，尚未推导完成

2026-09-25补充：用户要求对齐现有3e-3g的边界项与余项结构。本稿§5–6为坐标导数全部作用后的直接表示；L、E_r、F_r不是IBP边界项或余项的分类。之前“可开始实现”仅适用于直接版C/J，不代表IBP版就绪。当前优先级改为先补齐IBP理论，再建设完整Montopar输出。

推导起点应回到§4的$\mathscr V_r[\mathscr K^{[r]}]$。对固定切向量和独立端点坐标，

$$
\frac{d}{ds_i}K(X_1,X_2,X_3)=u_i\cdot\partial_{X_i}K.
$$

交叉导数$u_j\cdot\partial_{X_i}$不能直接当作$s_i$导数。需用明确的张量收缩和平移恒等式整理，目标是在保留同一测度及归一化的条件下构造

$$
\mathcal K_r=\sum_i\partial_{s_i}\mathcal B_{r,i}+\mathcal R_r.
$$

此式是待构造的目标，不是本稿已证明的分解。插入腿的二次张量分子也被微分，其产生项必须完整保留；不能直接套树级3e-3g的边界系数，也不预设余项为零。

IBP验收要求：

1. 给出每个全导数的原函数与显式余项，逐点回代重构§6两项核，保持D=4−2ε及i0。
2. 在原排序域积分，列出端点面及同边接触面；同一接触面上的法向符号和贡献必须合并。不能将依赖上限当作三个独立区间处理。
3. 导出边界积分及体积余项的C/J账本，所有项求和回到原始直接C/J；不改变颜色、Euler准备或最终I定义。
4. 在重复附着与三边附着代表上核验，再扩展12图；检查直接形式与IBP形式的等价关系，而非仅检查表达式可生成。

当前179项检查及公式渲染只覆盖直接表示，不是IBP验收。下一步实施本节理论推导；本次仅记录决定，没有计算新的边界项或修改生产程序。
