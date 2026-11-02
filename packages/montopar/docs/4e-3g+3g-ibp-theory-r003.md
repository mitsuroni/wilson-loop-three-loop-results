> Implementation update, 2026-09-20: this full theory is now implemented by `--tgv-ibp`. See [output contract](double-tgv-ibp.md). The final status paragraph below records the earlier document-only milestone; it is preserved as history.

# 4e-3g+3g完整理论：从初始Feynman积分到双层轮廓IBP（r003）

2026-09-20。本版是后续montopar改造的完整理论依据。正文自行定义全部符号，从五传播子动量积分开始，推导raw C/J，再推导IBP C/J；无需先读r001或r002才能补全步骤。r001、r002及其证据保留原字节。

本版的核心问题是：**三胶子顶点中的哪一部分动量因子能变成Wilson参数的全导数，这些全导数在有序积分域上产生什么边界？** 参数化仅负责把共同标量核写清楚，IBP的来源是原始顶点与Wilson相位。

推导依次经过：初始动量式 → 开放顶点收缩 → 顶点的切向/余项分解 → Fourier标量母核 → 五α参数化 → 统一尺度与raw C/J → 两层有序域IBP → 最终可输出的C/J。

目前原始C/J已由montopar输出；本版IBP层尚未进入生产代码。本文不新增Laurent系数或形式证明收据。

## 1. 几何、腿身份、有序域与颜色

### 1.1 几何约定

取D=4−2ε、Lorentz度规(+,-,…,-)，四边形的边向量满足

$$v_e^2=0,\qquad \sum_{e=1}^4v_e=0,\qquad
V_1=0,\quad V_e=\sum_{r<e}v_r.$$

本项目约定vᵢⱼ=vᵢ·vⱼ，不额外含2。四条外胶子腿用i=1,…,4编号，所在轮廓边用eᵢ编号：

$$u_i=v_{e_i},\qquad X_i=V_{e_i}+s_i u_i,\qquad u_{ij}=u_i\cdot u_j.$$

腿(1,2)接左三胶子顶点，腿(3,4)接右顶点，第五条胶子连接两个内部顶点。腿编号、边编号和边内端点次序是三个不同概念。

取参考尺度及无量纲运动学

$$v_{\rm ref}=v_1\cdot v_4\ne0,\qquad A_0=-2v_{\rm ref},\qquad
A_{\rm ref}=-2v_{\rm ref}+i0,\qquad z=\frac{v_1\cdot v_2}{v_{\rm ref}}.$$

从A₀>0、z>0的解析区域开始。四边的归一化Gram矩阵为

$$\frac{v_e\cdot v_f}{v_{\rm ref}}=
\begin{pmatrix}
0&z&-1-z&1\\
z&0&1&-1-z\\
-1-z&1&0&z\\
1&-1-z&z&0
\end{pmatrix}_{ef}.$$

该矩阵用于化简外部切向量内积；其秩不能替代内部D维度规迹。

### 1.2 原有序域

每条边保存从小到大的腿序ρₑ,₁,…,ρₑ,ₘₑ。积分域为

$$\mathcal D_\rho=\prod_e
\{0\le s_{\rho_{e,1}}\le\cdots\le s_{\rho_{e,m_e}}\le1\}.$$

每条链的最大参数先作外层积分；全图外到内次序随输入保存。在共享一条边的不同顶点外腿之间，也要保留排序约束。不能将Dρ拆成“左顶点域×右顶点域”，除非这确实是原图的域结构。

### 1.3 原始颜色与单一规范代表

每条腿的颜色Aᵢ独立。按轮廓次序π排列Wilson生成元，定义

$$\mathcal C_\rho=\frac1{N_c}\operatorname{tr}
(T^{A_{\pi(1)}}T^{A_{\pi(2)}}T^{A_{\pi(3)}}T^{A_{\pi(4)}})
f^{A_1A_2B}f^{A_3A_4B}.$$

采用Tr(TᵃTᵇ)=δᵃᵇ/2、[Tᵃ,Tᵇ]=ifᵃᵇᶜTᶜ，以及

$$\sum_a(T^a)_{ij}(T^a)_{kl}=
\frac12\left(\delta_{il}\delta_{jk}-\frac1{N_c}\delta_{ij}\delta_{kl}\right),
\qquad C_A=N_c,\quad C_F=\frac{N_c^2-1}{2N_c}.$$

对轮廓顺序a<b<c<d，配对(ab)(cd)、(ac)(bd)、(ad)(bc)的规范有向颜色依次为−CA²CF/4、0、+CA²CF/4。顶点内反向会同时改变颜色与Lorentz顶点符号。

颜色零图先筛除。本篇Wρ表示一张规范代表的贡献，不含C4重数、a₁投影或Euler补偿；IBP不会改变它的颜色或重新计算接触面颜色。

## 2. 初始动量空间表达式

外腿动量为ℓ₁,…,ℓ₄，内部动量为q，采用

$$\ell_1+\ell_2+q=0,\qquad \ell_3+\ell_4-q=0.$$

写

$$[d\ell]=\frac{d^D\ell}{(2\pi)^D},\qquad
\Delta(p)=\frac{-i}{p^2+i0}.$$

三胶子运动学顶点约定为

$$\Gamma_{\mu\nu\rho}(k_1,k_2,k_3)=
g_{\mu\nu}(k_1-k_2)_\rho+
g_{\nu\rho}(k_2-k_3)_\mu+
g_{\rho\mu}(k_3-k_1)_\nu.$$

对应的完整Lorentz收缩是

$$\mathcal L=
u_1^\mu u_2^\nu u_3^\alpha u_4^\beta
\Gamma_{\mu\nu\rho}(-\ell_1,-\ell_2,-q)
g^{\rho\sigma}
\Gamma_{\alpha\beta\sigma}(-\ell_3,-\ell_4,q).$$

四个Wilson顶点给(ig)⁴，两个三胶子顶点各给g，故耦合为g⁶，Wilson相位i⁴=1。保持五个标量传播子的−i，初始式为

$$\boxed{\begin{aligned}
\mathcal W_\rho={}&g^6\mathcal C_\rho\int_{\mathcal D_\rho}d^4s
\int\prod_{i=1}^4[d\ell_i]\,[dq]\;
(2\pi)^D\delta^D(\ell_1+\ell_2+q)\\
&\times(2\pi)^D\delta^D(\ell_3+\ell_4-q)
\exp\!\left[-i\sum_i\ell_i\cdot X_i\right]
\left[\prod_{i=1}^4\Delta(\ell_i)\right]\Delta(q)\,\mathcal L.
\end{aligned}}$$

用第一个δ积分掉q，得到q=−(ℓ₁+ℓ₂)；第二个δ仍然留下总体动量约束：

$$\boxed{\begin{aligned}
\mathcal W_\rho={}&g^6\mathcal C_\rho\int_{\mathcal D_\rho}d^4s
\int\prod_{i=1}^4[d\ell_i]\;(2\pi)^D\delta^D(\ell_1+\ell_2+\ell_3+\ell_4)\\
&\times e^{-i\sum_i\ell_i\cdot X_i}
\prod_{i=1}^4\frac{-i}{\ell_i^2+i0}
\frac{-i}{(\ell_1+\ell_2)^2+i0}\,\mathcal L.
\end{aligned}}$$

总体δ不能丢。L中已经收缩四个Wilson切向量，传播子之外不再额外乘四个带自由指标的uᵢ。

## 3. 收缩顶点，并在动量式中找出IBP来源

### 3.1 保留内部指标的两个向量

先只收缩左顶点的两条Wilson外腿。利用q=−(ℓ₁+ℓ₂)，三个顶点项分别贡献

$$u_{12}(\ell_2-\ell_1),\qquad
-u_2(u_1\cdot\ell_1+2u_1\cdot\ell_2),\qquad
u_1(2u_2\cdot\ell_1+u_2\cdot\ell_2).$$

因此

$$\boxed{\mathbf V_L=
u_{12}(\ell_2-\ell_1)
+u_1(2u_2\cdot\ell_1+u_2\cdot\ell_2)
-u_2(u_1\cdot\ell_1+2u_1\cdot\ell_2).}$$

右顶点同样给出

$$\boxed{\mathbf V_R=
u_{34}(\ell_4-\ell_3)
+u_3(2u_4\cdot\ell_3+u_4\cdot\ell_4)
-u_4(u_3\cdot\ell_3+2u_3\cdot\ell_4),\qquad
\mathcal L=\mathbf V_L\cdot\mathbf V_R.}$$

粗体V是动量的一次函数，仍保留内部Lorentz指标。这使两个顶点之间的度规收缩可追踪。

### 3.2 哪些动量因子可以积分成边界

Wilson参数sᵢ只通过相位中的Xᵢ=Vₑᵢ+sᵢuᵢ出现，因此

$$\boxed{(u_i\cdot\ell_i)e^{-i\sum_j\ell_j\cdot X_j}
=i\,\partial_{s_i}e^{-i\sum_j\ell_j\cdot X_j}.}$$

这就是轮廓IBP的来源：**自腿动量与自腿切向量的收缩，可以转成相应路径参数的全导数。** 传播子不显含s；本模块也没有额外pure传播子形成的s依赖旁观因子。

据此，在左顶点选择

$$\mathbf T_L=u_2(u_1\cdot\ell_1)-u_1(u_2\cdot\ell_2),$$

并定义W_L=u₁⊗u₂−u₂⊗u₁，其中(u⊗v)k=u(v·k)。逐项相减得到

$$\boxed{\mathbf R_L=
u_{12}(\ell_2-\ell_1)+2W_L(\ell_1+\ell_2),\qquad
\mathbf V_L=\mathbf T_L+\mathbf R_L.}$$

由动量守恒，后一项也可写成−2W_Lq。右顶点为

$$\mathbf T_R=u_4(u_3\cdot\ell_3)-u_3(u_4\cdot\ell_4),$$

$$\boxed{\mathbf R_R=
u_{34}(\ell_4-\ell_3)+2W_R(\ell_3+\ell_4)
=u_{34}(\ell_4-\ell_3)+2W_Rq,\qquad
W_R=u_3\otimes u_4-u_4\otimes u_3.}$$

左右余项都保留共享内部动量q的信息。这说明两个顶点不能当成互不相关的两份3e-3g标量公式。

T/R分解是明确的一种选择，不是唯一IBP。R不表示正交于轮廓的纯横向分量；它表示这次提取自腿路径导数之后剩下的顶点动量结构。

### 3.3 退化为什么自动消掉一层余项

若左顶点两条外腿在同一类光边上，u₁=u₂=u且u²=0，则u₁₂=0、W_L=0：

$$\boxed{\mathbf R_L=0,\qquad
\mathbf V_L=u\,[u\cdot\ell_1-u\cdot\ell_2].}$$

这一结论发生在参数化之前。右顶点同理。更一般地，非零且共线的类光外腿也有同样退化；一般复运动学中仅有u₁₂=0而无W_L=0，不足以断言余项消失。

## 4. Fourier变换：将动量顶点变成微分算子

### 4.1 共同坐标标量核

在两个动量δ中分别引入内部坐标Y_L、Y_R：

$$(2\pi)^D\delta^D(k)=\int d^DY\,e^{iY\cdot k}.$$

这样五条传播子分别对应位移X₁−Y_L、X₂−Y_L、X₃−Y_R、X₄−Y_R、Y_R−Y_L。采用

$$\int[d\ell]\,\frac{-i}{\ell^2+i0}e^{-i\ell\cdot X}
=-\frac{\Gamma(\nu)}{4\pi^{D/2}}\frac1{(-X^2+i0)^\nu},
\qquad \nu=1-\epsilon.$$

定义

$$\Delta_{1,2}=-(X_{1,2}-Y_L)^2+i0,quad
\Delta_{3,4}=-(X_{3,4}-Y_R)^2+i0,quad
\Delta_5=-(Y_L-Y_R)^2+i0.$$

于是共同坐标标量核为

$$\Phi(X_1,\ldots,X_4)=\int d^DY_Ld^DY_R\prod_{r=1}^5\Delta_r^{-\nu}.$$

### 4.2 两个一次算子与共同负号

令∂ᵢ表示对独立外部坐标Xᵢ的D维梯度。由ℓᵢ↦i∂ᵢ，定义去掉这个i的一次算子V_L、V_R，并有

$$\boxed{\mathcal O=\mathcal L(i\partial)=-\mathcal V_L\cdot\mathcal V_R.}$$

负号来自两个一次顶点共同的i²，四个块都必须使用同一个符号。于是

$$\mathcal W_\rho=g^6\mathcal C_\rho
\left[-\frac{\Gamma(\nu)}{4\pi^{D/2}}\right]^5
\int_{\mathcal D_\rho}d^4s\;\mathcal O\Phi.$$

在所有X导数中，uᵢ保持固定。不能先把Xᵢ限制到一维路径，再把整个D维梯度都替换成一个s导数。

对真正的切向组合，链式法则则给出合法等式

$$\left.(u_i\cdot\partial_i f)\right|_{X=X(s)}
=\partial_{s_i}f(X(s)).$$

因此前一节的动量分解变成

$$\boxed{\begin{aligned}
\mathcal T_L&=u_2\partial_{s_1}-u_1\partial_{s_2},&
\mathcal R_L&=u_{12}(\partial_2-\partial_1)+2W_L(\partial_1+\partial_2),\\
\mathcal T_R&=u_4\partial_{s_3}-u_3\partial_{s_4},&
\mathcal R_R&=u_{34}(\partial_4-\partial_3)+2W_R(\partial_3+\partial_4).
\end{aligned}}$$

V_L=T_L+R_L、V_R=T_R+R_R。参数化之后，IBP仍由这两个恒等式产生，不必重新拟合分子。

## 5. 五α参数化：推导U、B及其幂次

### 5.1 一个单纯形约束

Feynman合并公式为

$$\prod_{r=1}^5\Delta_r^{-\nu}
=\frac{\Gamma(5\nu)}{\Gamma(\nu)^5}
\int_{\alpha_r\ge0}d^5\alpha\,\delta(1-\sum\alpha_r)
\frac{\prod_r\alpha_r^{\nu-1}}{(\sum_r\alpha_r\Delta_r)^{5\nu}}.$$

上述Γ(5ν)可以由五个Schwinger参数τᵢ=Tαᵢ直接看出：换元Jacobian为T⁴，五个参数幂连同测度给T^(5ν−1)，积分总尺度T便产生Γ(5ν)。

所以有五个非负α、四个独立α变量，而不是五个独立单位区间。可消去α₅=1−α₁−α₂−α₃−α₄，按α₁、α₂、α₃、α₄外到内积分，上界依次为1、1−α₁、1−α₁−α₂、1−α₁−α₂−α₃。

### 5.2 二元高斯完成平方

设

$$\Sigma_L=\alpha_1+\alpha_2+\alpha_5,\quad
\Sigma_R=\alpha_3+\alpha_4+\alpha_5,\quad
\Sigma=\begin{pmatrix}\Sigma_L&-\alpha_5\\-\alpha_5&\Sigma_R\end{pmatrix},
\quad U=\det\Sigma.$$

再记

$$Q_L=\alpha_1X_1+\alpha_2X_2,\quad
Q_R=\alpha_3X_3+\alpha_4X_4,\quad
K=\sum_{i=1}^4\alpha_iX_i^2.$$

若在坐标分母中暂时使用共同有限调节+iκ，则

$$\sum_r\alpha_r\Delta_r=
-Y^T\Sigma Y+2Y^TQ-K+i\kappa\sum_r\alpha_r$$

$$=-(Y-\Sigma^{-1}Q)^T\Sigma(Y-\Sigma^{-1}Q)
+\frac{B_\kappa}{U},$$

其中Y=(Y_L,Y_R)、Q=(Q_L,Q_R)，矩阵指标与Lorentz指标分别收缩，且

$$B_0=Q^T\operatorname{adj}(\Sigma)Q-UK,\qquad
B_\kappa=B_0+i\kappa U\sum_r\alpha_r.$$

以Φκ和Wρ,κ表示加入这个坐标辅助调节后的母核与图贡献；物理边界值仍按原始处方取κ↓0。κ在这里是坐标分母的辅助调节，具有长度平方量纲；不能将它当作动量传播子中有限质量型调节后仍适用同一个简单幂律FT公式。

把B₀展开成点间距离dᵢⱼ=−(Xᵢ−Xⱼ)²，得到

$$\boxed{\begin{aligned}
U={}&\Sigma_L\Sigma_R-\alpha_5^2,\\
B_0={}&\alpha_1\alpha_2\Sigma_R d_{12}
+\alpha_3\alpha_4\Sigma_L d_{34}\\
&+\alpha_5(\alpha_1\alpha_3d_{13}+\alpha_1\alpha_4d_{14}
+\alpha_2\alpha_3d_{23}+\alpha_2\alpha_4d_{24}).
\end{aligned}}$$

U在α中为二次，B₀为三次。在单纯形上Bκ=B₀+iκU；这才是共同坐标处方的继承形式。

### 5.3 积掉两个内部坐标

在允许Wick旋转和高斯积分的共同正规化解析区域，两个内部坐标积分给出i²=−1。移到Euclidean积分变量后，设M²=Bκ/U。利用一个Laplace参数，

$$\frac1{(Y_E^T\Sigma Y_E+M^2)^{5\nu}}
=\frac1{\Gamma(5\nu)}\int_0^\infty d\tau\,
\tau^{5\nu-1}e^{-\tau M^2}e^{-\tau Y_E^T\Sigma Y_E}.$$

两个D维Gaussian积分及剩余τ积分分别为

$$\int d^{2D}Y_E\,e^{-\tau Y_E^T\Sigma Y_E}
=\pi^D\tau^{-D}U^{-D/2},$$

$$\int_0^\infty d\tau\,\tau^{5\nu-D-1}e^{-\tau B_\kappa/U}
=\Gamma(5\nu-D)(B_\kappa/U)^{D-5\nu}.$$

这些公式先在相应实部条件下成立，再延拓至共同处方。因此U的指数为−D/2+(5ν−D)=5ν−3D/2；连同Feynman合并系数，

$$\Phi_\kappa=-\frac{\pi^D\Gamma(5\nu-D)}{\Gamma(\nu)^5}
\int d\mu_\alpha\,
\prod_r\alpha_r^{-\epsilon}
U^{5\nu-3D/2}B_\kappa^{D-5\nu},$$

$$d\mu_\alpha=d^5\alpha\,\delta(1-\sum\alpha_r),\qquad
5\nu-D=1-3\epsilon=\lambda,quad
5\nu-3D/2=-1-2\epsilon.$$

五个坐标传播子的Γ(ν)⁵与分母Γ(ν)⁵消去；其(−1)⁵与两个内部积分的i²相乘为+。因此

$$\boxed{\mathcal W_{\rho,\kappa}=
\frac{g^6\mathcal C_\rho\Gamma(\lambda)}{1024\pi^{6-3\epsilon}}
\int d\mu_\alpha\,w_\alpha\int_{\mathcal D_\rho}d^4s\,
\mathcal O B_\kappa^{-\lambda},\qquad
w_\alpha=\prod_r\alpha_r^{-\epsilon}U^{-1-2\epsilon}.}$$

O中的动量替换负号仍然保留，不能再加入一次i²。上述解析表示随后按维数正则作亚纯延拓，不据此宣称每个无调节积分在ε=0收敛。

## 6. 从共同标量核得到raw C、J

### 6.1 统一尺度提取

定义

$$\widehat B=B_0/A_0,quad \eta=\kappa/A_0,quad
\mathscr B_\eta=\widehat B+i\eta U,quad
\Psi_\eta=\mathscr B_\eta^{-\lambda}.$$

在正尺度起始区Bκ⁻λ=A₀⁻λΨη。A₀与uᵢ均不随独立Xi导数变化。再用g⁶=(4παs)³以及Γ(λ+1)=λΓ(λ)，可把上一节系数改写成

$$\boxed{C_{\rm raw}=\mathcal C_\rho
\left(\frac{\alpha_s}{\pi}\right)^3
\frac{\pi^{3\epsilon}\Gamma(2-3\epsilon)}{16}
A_{\rm ref}^{3\epsilon},}$$

$$\boxed{J_{\rm raw}=
\frac1{\lambda A_0}\int d\mu_\alpha\,w_\alpha
\int_{\mathcal D_\rho}d^4s\;\mathcal O\Psi,\qquad
\mathcal W_\rho=C_{\rm raw}J_{\rm raw}.}$$

因为−λ=3ε−1，剩下的尺度正是1/A₀。因子1/(λA₀)不是可随意丢弃的记号。这里物理J取共同处方的边界值；为便于核验，下文保留有限η，未逐式书写共同的η↓0极限。保留调节时，两边均指同一个辅助调节表达。

### 6.2 原始两分子形式作为连接检查

为把算子式连接到现有montopar多项式，定义

$$P_i=\partial_i B_0=-2\alpha_i\sum_{j\ne i}\alpha_jc_{ij}(X_i-X_j),$$

其中c₁₂=ΣR、c₃₄=ΣL，所有跨顶点cᵢⱼ=α₅。Hessian为

$$\partial_i^\mu\partial_j^\nu B_0=
\begin{cases}
2\alpha_i\alpha_jc_{ij}g^{\mu\nu},&i\ne j,\\
-2\alpha_i\sum_{k\ne i}\alpha_kc_{ik}g^{\mu\nu},&i=j.
\end{cases}$$

于是i≠j时∂ᵢ·∂ⱼB₀=2Dαᵢαⱼcᵢⱼ。定义

$$T_1=-\mathbf V_L(P_1,P_2)\cdot\mathbf V_R(P_3,P_4),\qquad T_2=\mathcal O B_0.$$

这里仅把V中的动量换成P，外部uᵢ不变。二阶链式法则直接给出

$$\mathcal O B_\kappa^{-\lambda}
=\lambda(\lambda+1)T_1B_\kappa^{-\lambda-2}
-\lambda T_2B_\kappa^{-\lambda-1}.$$

令T̂₁=T₁/vref³、T̂₂=T₂/vref²，则

$$\boxed{J_{\rm raw}=\int d\mu_\alpha\,w_\alpha\int_{\mathcal D_\rho} d^4s
\left[-\frac{\lambda+1}{8}\frac{\widehat T_1}{\mathscr B_\eta^{\lambda+2}}
-\frac14\frac{\widehat T_2}{\mathscr B_\eta^{\lambda+1}}\right].}$$

两个常数分别来自vref³/A₀³=−1/8和vref²/A₀²=1/4。T₂可以独立写成

$$T_2=-2\alpha_5\big[(\alpha_1\alpha_3+\alpha_2\alpha_4)K_+
+(\alpha_1\alpha_4+\alpha_2\alpha_3)K_-\big],$$

$$K_+=(D-2)u_{12}u_{34}+5u_{13}u_{24}-4u_{14}u_{23},$$

$$K_-=-(D-2)u_{12}u_{34}+4u_{13}u_{24}-5u_{14}u_{23}.$$

至此已完整得到raw C/J。下面只改变J的表示方式，C的全局规范继续保持。

## 7. 为什么两层IBP恰好产生四个算子块

将第4节的V=T+R代回O，得到

$$\boxed{\mathcal O=-\left(
\mathcal T_L\cdot\mathcal T_R+
\mathcal T_L\cdot\mathcal R_R+
\mathcal R_L\cdot\mathcal T_R+
\mathcal R_L\cdot\mathcal R_R\right).}$$

uᵢ是常量；左右算子作用在不同外腿坐标，常系数导数可交换。故可以在保留共同α标量核的同时逐块做IBP：

| 算子块 | 可以提取的路径导数 | 预期产生的对象 |
|---|---|---|
| T_L·T_R | 两次 | 双边界，以及接触面可能留下的法向导数 |
| T_L·R_R | 左顶点一次 | 左边界上的右顶点余项 |
| R_L·T_R | 右顶点一次 | 右边界上的左顶点余项 |
| R_L·R_R | 本规范尚未提取 | 双余项体积分 |

这不是把两个已经积分的三胶子结果相乘。两个顶点仍作用于同一个B与U，内部指标仍通过Lorentz内积相连。

## 8. 第一层IBP：边界和接触项从哪里来

### 8.1 独立区间

若s独立地在[0,1]积分，则

$$\int_0^1 ds\,c\,\partial_s f=c[f(1)-f(0)].$$

这给出两个cusp端点。沿路径的导数只在真正的独立区间上才可以这样单独使用。

### 8.2 两点有序三角域

设0≤a≤b≤1。对cₐ∂ₐf先积a，直接得到接触面a=b与下端a=0。对c_b∂bf，需要保留a积分上限b的导数：

$$\frac{d}{db}\int_0^b f(a,b)da
=f(b,b)+\int_0^b\partial_b f(a,b)da.$$

两部分合并，得到

$$\boxed{\begin{aligned}
\int_0^1db\int_0^bda\,(c_a\partial_a+c_b\partial_b)f
={}&-c_a\int_0^1db\,f(0,b)
+c_b\int_0^1da\,f(a,1)\\
&+(c_a-c_b)\int_0^1dt\,f(t,t).
\end{aligned}}$$

接触项来自依赖积分上限，不是人为添加的边界。若漏掉它，积分域已被改变。

### 8.3 任意有序链

对0≤s₁≤…≤sₘ≤1，配合继承的参数测度，三个面类型的通量协向量为

$$n_f=-e_1\quad(s_1=0),\qquad
n_f=+e_m\quad(s_m=1),\qquad
n_f=e_j-e_{j+1}\quad(s_j=s_{j+1}).$$

因此

$$\int_{\mathcal D_\rho}\sum_i c_i\partial_{s_i}f
=\sum_f\int_{D_f}d^3t\,
\left(\sum_i n_{fi}c_i\right)f(\phi_f(t)).$$

φf把面上的独立参数t映回原s；D_f仍保留未消失的排序约束。接触面使用这个参数测度时不另乘√2；单位法向与Euclidean面面积中的因子已经相抵。

### 8.4 应用于两个开放顶点

用固定Lorentz向量记录切向系数：

$$\boldsymbol a_1=u_2,\quad\boldsymbol a_2=-u_1,
\quad\boldsymbol a_3=\boldsymbol a_4=0,$$

$$\boldsymbol b_3=u_4,\quad\boldsymbol b_4=-u_3,
\quad\boldsymbol b_1=\boldsymbol b_2=0.$$

定义面通量向量

$$\omega_{L,f}=\sum_i n_{fi}\boldsymbol a_i,\qquad
\omega_{R,f}=\sum_j n_{fj}\boldsymbol b_j.$$

两个混合块的一层IBP就是

$$\int_{\mathcal D_\rho}\mathcal T_L\cdot\mathcal R_R\Psi
=\sum_f\int_{D_f}\left.(\omega_{L,f}\cdot\mathcal R_R\Psi)\right|_f,$$

$$\int_{\mathcal D_\rho}\mathcal R_L\cdot\mathcal T_R\Psi
=\sum_f\int_{D_f}\left.(\omega_{R,f}\cdot\mathcal R_L\Psi)\right|_f.$$

余算子必须先作用于父空间中的核，再做面限制。这两项各剩三个Wilson参数。

## 9. 第二层IBP：为什么接触面还需要法向修正

### 9.1 第一次IBP后的导数仍是父空间导数

对TT块先做左层，得到

$$\int_{\mathcal D_\rho}\mathcal T_L\cdot\mathcal T_R\Psi
=\sum_f\int_{D_f}\left.(b_f\cdot\partial_s\Psi)\right|_f,
\qquad (b_f)_j=\omega_{L,f}\cdot\boldsymbol b_j.$$

如果f是sᵣ=sₖ，面上的公共参数导数满足

$$\frac{d}{dt}\Psi(\ldots,s_r=t,\ldots,s_k=t,\ldots)
=\left.(\partial_{s_r}+\partial_{s_k})\Psi\right|_f.$$

因此不能将父导数∂sᵣ或∂sₖ中的任意一个直接替换成d/dt。

### 9.2 通用线性面投影

令s=φf(t)，E_f=∂s/∂t。每个余维1面消去一个原变量sᵣ，其余变量作为面内坐标保留。取βf等于bf在保留坐标上的分量，并定义

$$\chi_f=(b_f)_r-(E_f\beta_f)_r.$$

由于其余各行已按原坐标保留，向量恒等式为

$$b_f=E_f\beta_f+\chi_f e_r.$$

再用链式法则∂t(Ψ∘φf)=E_fᵀ(∂sΨ)∘φf，得到

$$\boxed{\left.(b_f\cdot\partial_s\Psi)\right|_f
=\beta_f\cdot\partial_t(\Psi\circ\phi_f)
+\chi_f\left.(\partial_{s_r}\Psi)\right|_f.}$$

χf在r002中记作νF，本版改记χ以免与传播子指数ν混淆。具体地，固定面sᵣ=0或1有χf=bfr；接触面sᵣ=sₖ有χf=bfr−bfk。

第一项可在D_f上再按第8节积分，第二项则保留为面上的法向导数。法向项不增加新的体积分，但不能丢弃。

### 9.3 第二层通量与复合面映射

对D_f的子面g，记其通量协向量为n_g^f，令

$$k_{fg}=n_g^f\cdot\beta_f.$$

用φ_g^f把两个独立参数r映到D_f的边界，用

$$\psi_{fg}=\phi_f\circ\phi_g^f$$

映回原四个s。于是

$$\boxed{\begin{aligned}
\int_{\mathcal D_\rho}\mathcal T_L\cdot\mathcal T_R\Psi
={}&\sum_{f,g}k_{fg}\int_{D_{fg}}d^2r\,
\Psi(\psi_{fg}(r))\\
&+\sum_f\chi_f\int_{D_f}d^3t\,
(\partial_{s_r}\Psi)(\phi_f(t)).
\end{aligned}}$$

不同路径可能到达同一个几何双边界。可以合并，但要先核对复合约束、继承域和参数对应，不能只凭面名合并。

### 9.4 一个不涉及发散的漏项反例

01的链为0<s41<s42<s43<1，另一个参数s31独立。取测试函数f=s41·s42，g=u₁·u₂，则TT=g∂s41(∂s43−∂s42)。在原有序域上直接积分得−g/6。

若在接触面s41=s42上只保留第二层切向项，所得为−g/3；法向项补回+g/6，才恢复−g/6。这个反例说明第二层投影是必要的链式法则，而非处理发散时的可选技巧。

## 10. 把IBP结果写成不含未执行微分算子的标量核

### 10.1 余项的显式分子

定义

$$G_L=\mathcal R_L\widehat B
=\frac1{A_0}\left[u_{12}(P_2-P_1)+2W_L(P_1+P_2)\right],$$

$$G_R=\mathcal R_R\widehat B
=\frac1{A_0}\left[u_{34}(P_4-P_3)+2W_R(P_3+P_4)\right],$$

$$H_{RR}=(\mathcal R_L\cdot\mathcal R_R)\widehat B.$$

Pᵢ已在第6节显式给出，所以GL、GR可直接生成，无需把母核当作未知函数。父参数导数也明确为

$$\partial_{s_i}\widehat B=\frac{u_i\cdot P_i}{A_0}.$$

对于H_RR，跨顶点Hessian仅用cᵢⱼ=α₅。利用

$$\langle\mathbf1,\mathbf1\rangle=D,\qquad
\langle\mathbf1,u\otimes v\rangle=u\cdot v,\qquad
\langle u\otimes v,w\otimes t\rangle=(u\cdot w)(v\cdot t),$$

逐项收缩得到

$$\boxed{\begin{aligned}
H_{RR}=\frac{2\alpha_5}{A_0}\Big[&D\,u_{12}u_{34}
(\alpha_1-\alpha_2)(\alpha_3-\alpha_4)\\
&+8(u_{13}u_{24}-u_{14}u_{23})
(\alpha_1+\alpha_2)(\alpha_3+\alpha_4)\Big].
\end{aligned}}$$

这与r001记号的T₂,RR=−A₀H_RR一致。显式D迹位于RR块；TT、TR、RT中没有两个单位张量彼此收缩产生的D迹。

### 10.2 对共同母核应用链式法则

由于ηU与Xi无关，有限η下也有

$$\mathcal R_L\Psi_\eta=-\lambda G_L\mathscr B_\eta^{-\lambda-1},\qquad
\mathcal R_R\Psi_\eta=-\lambda G_R\mathscr B_\eta^{-\lambda-1},$$

$$\partial_{s_r}\Psi_\eta=-\lambda(\partial_{s_r}\widehat B)
\mathscr B_\eta^{-\lambda-1},$$

$$\mathcal R_L\cdot\mathcal R_R\Psi_\eta
=\lambda(\lambda+1)(G_L\cdot G_R)\mathscr B_\eta^{-\lambda-2}
-\lambda H_{RR}\mathscr B_\eta^{-\lambda-1}.$$

后一式的Hessian项正是“两层”不能当成两份独立一阶结果相乘的原因之一。

定义每个余维1面的父空间分子

$$\mathcal N_f=
\omega_{L,f}\cdot G_R+\omega_{R,f}\cdot G_L
+\chi_f\partial_{s_r}\widehat B.$$

这里三个来源分别是TR、RT和TT接触面法向修正；代入面约束发生在求导之后。

### 10.3 最终通用IBP J

把第7至第10节代回J=−⟨(T+R)_L·(T+R)_RΨ⟩/(λA₀)，得到

$$\boxed{\begin{aligned}
J_{\rm IBP}=\int d\mu_\alpha\,w_\alpha\Bigg\{&
-\frac1{\lambda A_0}\sum_{f,g}k_{fg}
\int_{D_{fg}}d^2r\,
[\mathscr B_\eta\circ\psi_{fg}]^{-\lambda}\\
&+\frac1{A_0}\sum_f\int_{D_f}d^3t\,
\frac{\mathcal N_f\circ\phi_f}
{[\mathscr B_\eta\circ\phi_f]^{\lambda+1}}\\
&+\frac1{A_0}\int_{\mathcal D_\rho}d^4s\left[
-(\lambda+1)\frac{G_L\cdot G_R}{\mathscr B_\eta^{\lambda+2}}
+\frac{H_{RR}}{\mathscr B_\eta^{\lambda+1}}\right]\Bigg\}.
\end{aligned}}$$

这就是完整的双边界、单边界余项、双余项体积分表示。它与J_raw有相同的共同Feynman处方和全局C；在共同正规化下J_IBP=J_raw。

现在正文中的每个对象均已定义：B、U来自高斯参数化；GL、GR、H_RR来自开放顶点余项；kfg来自两次有序面通量；Nf另外保留面投影的χf修正。

### 10.4 可直接输出的C、J分配

也可以把上式逐项写成Wρ=ΣₜCₜJₜ。一种不改变规范的分配是：

| 项 | Cₜ/C_raw | Jₜ的核（均含α测度与相应s域） |
|---|---|---|
| 双边界fg | −kfg/(λA₀) | Bη⁻λ |
| 单边界f | 1 | (Nf/A₀) Bη⁻λ⁻¹ |
| RR梯度乘积 | −(λ+1) | (GL·GR/A₀) Bη⁻λ⁻² |
| RR Hessian | 1 | (H_RR/A₀) Bη⁻λ⁻¹ |

表中各项J都是无量纲的；kfg/A₀也是无量纲。Bη在这里指相应域上的无量纲B的共同调节值。两个RR行可保持为一个体积分组。

建议输出同时保存共同C_raw、局部比例和完整combined_C，而不要把Γ(λ+1)/λ的递推或几何尺度隐含在项名里。IBP不是新的Euler补偿，也不是新的颜色投影。

## 11. 两层退化、S/N区别以及项数

### 11.1 通用退化层级

| 顶点外腿条件 | 四块的简化 | 本批非零图 |
|---|---|---|
| 两个顶点各自同边 | RL=RR=0，只剩TT | 04、19 |
| 恰有一个顶点同边 | 双余项体积分消失，仍有另一侧面余项 | 01、03、07、09、10、12、13、18、22 |
| 两个顶点均不退化 | 不能预设任何余项消失 | 06、15、16、21、24、25 |

两个R均为零不自动取消有序域的接触面规则；本批04、19恰好左右两组排序链分离，所以其TT最终是纯双边界。仅一个R为零时，另一个顶点余项与跨顶点接触面修正仍要保留。

### 11.2 DT1-S：为什么只剩九个双边界

04的配对是(s31,s32;s41,s42)，u₁=u₂=v₃、u₃=u₄=v₄。于是

$$\mathcal T_L=v_3(\partial_{s31}-\partial_{s32}),\qquad
\mathcal T_R=v_4(\partial_{s41}-\partial_{s42}),\qquad
\mathcal R_L=\mathcal R_R=0.$$

每个有序三角域中，内参数导数的面系数是(−1,1,0)，外参数导数是(0,−1,1)，面顺序为下端、接触、上端。两者之差是

$$f=(-1,2,-1).$$

左右各三个面，双边界权重与f⊗f成比例，九个位置均非零。退化在顶点级已经成立，不依赖后来如何求α积分。

### 11.3 DT1-N：两个顶点都跨边

06的配对是(s31,s42;s32,s41)。为比较两个三角域，记a=s31、b=s32、c=s41、d=s42，a<b、c<d，g=v₃·v₄。直接展开切向块可得

$$\mathcal T_L\cdot\mathcal T_R
=-g(\partial_a\partial_c+\partial_b\partial_d).$$

合并同一几何双边界后，其矩阵与

$$\begin{pmatrix}1&-1&0\\-1&2&-1\\0&-1&1\end{pmatrix}$$

成比例，因此有7个非零双边界。两个混合块各有4个初始面贡献，共享的两个接触面合并后有6个非零单边界余项；BB在本图接触面上的χf恰好为零。RR仍非零。

所以本规范下N是“7个双边界＋6个单边界余项＋1组体余项”。这不是把它替换成7个已知S积分，也不证明体余项已经容易求值。

### 11.4 与3e-3g的6个边界项相比

单TGV的非退化三外腿情形可以提取三个路径导数，各有两个端点，得到6个边界项及一个余项。双TGV的每个开放顶点只有两条Wilson外腿，内部q没有Wilson参数。

若四个s彼此独立，当前规范形式上给4×4=16个双边界、4+4=8个单边界余项和一组RR体积分。存在同边排序时必须重新按域合并，不能将16+8+1作为所有图不加条件的项数。

原双TGV有4个s加4个独立α，共8维。双边界降为6维，单边界余项为7维，RR体积分仍为8维。计数依赖所选IBP规范，不是独立主积分数，也不是不可约解空间维数。

## 12. 正规化、边界值及不可省略的数学条件

先固定有限η>0，并将α限制在单纯形内部αᵢ≥δα>0。此时U>0，共同母核在紧致有序s域上光滑，换序及两层IBP可以按普通微积分解释。δα截断不依赖s，因此不会产生额外的路径边界。

随后，对全部IBP项采用同一调节、先组合，再沿与原raw相同的边界值/亚纯延拓过程移去调节。单独边界项可能比总和更奇异；若要分项积分和展开，必须另给共同收敛区、辅助解析调节或相消连接。

所有Xi导数先在父空间执行；所有面限制后再明确其内在导数。D维迹保留D=4−2ε，不以轮廓Gram秩代替。不同非整数幂不能独立重选复分支。

特别是04的几何关系T̂₁=−2B̂T̂₂，在有限η下不是T̂₁=−2BηT̂₂。直接使用本篇算子IBP时处方被完整保留；若另用04特殊合并核，须按解析区化简后的共同边界值解释。

## 13. 最终理论与后续实现的对应

本版固定的理论接口为：

1. 原始对象是第2节五传播子Feynman积分，含两个δ、有向顶点与有序域。
2. IBP来源是第3节自腿动量乘Wilson相位的全导数恒等式；T/R选择在参数化前确定。
3. 第4—6节把同一个对象变为共同标量母核与raw C/J，完整保留Gamma、尺度、相位和D维Hessian。
4. 第8—10节将路径导数变为带明确面映射的IBP C/J；两次IBP之间保留χf修正。
5. 原始C/J作为可回代基准；J_ibp每项保存来源块、域、面约束、父导数与拉回顺序、共同调节分组及局部C比例。

同一父积分的颜色不因边界位置相等而重算。默认只维护一个规范代表，轨道重数另存。向唯一I转换时再做固定颜色基底投影、一次三圈Euler补偿及旋转组装，不能把这些工作与IBP混同。

若以后加入额外pure传播子，路径导数还会作用于它们的s依赖因子，必须新增相应项。本版定稿范围是四外腿、两个三胶子顶点、一条内部连接线的4e-3g+3g模块。

## 14. 验证记录与当前状态

r002已有258项计算侧检查覆盖通用算子、D维迹、全部26图原始分子连接、面投影和漏项负例。本版保留这些对象和规范，复用其原字节证据；新增检查针对初始动量顶点、相位导数、完成平方、Gamma/尺度连接及本篇使用的项数。

检查与公式渲染结果见`workflow/checks/double-tgv-theory-r003/README.md` (Wilson workspace)。它们不替代Lean形式证明，也不自动证明每个拆分项独立收敛。

**文档定稿完成后，生产montopar仍保持原raw输出。按本轮最新要求，本次不启动代码改造、不改旧数据、不增加Laurent系数，也不提交Git。**
