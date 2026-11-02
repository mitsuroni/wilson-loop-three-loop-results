> **2026-09-19实现更新：montopar现已按本文修复，并重新生成26个单一规范C/J（schema 1.3）。** 见[当前逐图C/J](4e-3g+3g-CJ-formulas.md)。下文为理论r001的原始正文，其中“代码尚未修复”等状态语句保留为该推导完成时的历史记录；数学公式未改变。当前仍未求参数积分。

# 4e-3g+3g：从两个三胶子顶点到正确的 C×J

2026-09-19，修订r001。适用于四边形Wilson轮廓上、四条外腿通过一条内部胶子连接两个三胶子顶点的图。本文从动量式、张量收缩、Fourier变换和Feynman参数化逐步推导；四条腿可以落在重复的边上。

本文替代旧版同名通用推导。旧S/T/U命名、旧review正文及现有montopar输出不能作为本版公式的替代；当前Python实现尚待按本文修复。本页不计算参数积分的Laurent系数。

## 1. 输入必须同时保存腿、顶点配对和轮廓次序

取 $D=4-2\epsilon$，度规为 $(+,-,\ldots,-)$。四条轮廓边满足

$$v_e^2=0,\qquad \sum_{e=1}^4v_e=0,\qquad
V_1=0,\quad V_e=\sum_{r<e}v_r.$$

以 **i=1,2,3,4标记胶子腿**，以 $e_i\in\{1,2,3,4\}$ 标记该腿所在边，定义

$$u_i=v_{e_i},\qquad X_i=V_{e_i}+s_i u_i,\qquad 0\le s_i\le1,
\qquad u_{ij}=u_i\cdot u_j.$$

顶点配对规定为 $(1,2)$ 接左顶点、$(3,4)$ 接右顶点。一般输入可先通过**显式腿置换**整理成这个形式；同一置换必须作用于动量、位置、切向量、颜色和α标签。不能按边号排序，也不能按JSON中propagator数组的偶然位置分组。

在每条边上，还要保存其端点从小到大的腿序 $\rho_{e,1},\ldots,\rho_{e,m_e}$。原积分域是

$$\mathcal D_\rho=\prod_{e=1}^4
\{0\le s_{\rho_{e,1}}\le\cdots\le s_{\rho_{e,m_e}}\le1\}.$$

对应的外到内积分顺序为每条链的最大参数先、最小参数最后。不同边之间的次序可以沿用源文件；任何换序另有可积性/正则化义务。

运动学用

$$v_{\rm ref}=v_1\cdot v_4\ne0,\qquad z=\frac{v_1\cdot v_2}{v_{\rm ref}},
\qquad A_0=-2v_{\rm ref},\qquad A_{\rm ref}=A_0+i0,$$

$$\frac{v_e\cdot v_f}{v_{\rm ref}}=
\begin{pmatrix}
0&z&-1-z&1\\z&0&1&-1-z\\-1-z&1&0&z\\1&-1-z&z&0
\end{pmatrix}_{ef}$$

表示。先在 $A_0>0,z>0$ 的解析区域作无量纲分解，再沿同一Feynman处方延拓；不是给每个因子独立重选主分支。

本页将一张固定代表图的原始贡献记为 $\mathcal W_\rho=C_\rho^{\rm raw}J_\rho^{\rm raw}$。它尚未包含轨道重数、固定颜色基底投影及Euler补偿，也不等同于之后用于解表的唯一规范化I。

## 2. 颜色先于参数积分

每条腿的颜色指标为 $A_i$，即使 $e_i=e_j$，也不能令 $A_i=A_j$。按完整轮廓次序 $\pi$ 排列Wilson生成元，颜色系数为

$$\boxed{\mathcal C_\rho=
\frac1N\operatorname{tr}(T^{A_{\pi(1)}}T^{A_{\pi(2)}}T^{A_{\pi(3)}}T^{A_{\pi(4)}})
 f^{A_1A_2B}f^{A_3A_4B}.}$$

采用 $\operatorname{tr}(T^aT^b)=\delta^{ab}/2$、$[T^a,T^b]=if^{abc}T^c$。可用

$$\sum_a(T^a)_{ij}(T^a)_{kl}=\frac12
\left(\delta_{il}\delta_{jk}-\frac1N\delta_{ij}\delta_{kl}\right),
\quad C_A=N,\quad C_F=\frac{N^2-1}{2N}$$

进行精确收缩。对固定轮廓颜色顺序abcd，结果为

$$\begin{array}{c|ccc}
\text{顶点配对}&(ab)(cd)&(ac)(bd)&(ad)(bc)\\\hline
\mathcal C_\rho/(C_A^2C_F)&-1/4&0&+1/4
\end{array}$$

顶点内交换两条腿会使f和对应Lorentz顶点各变号，二者乘积不变；完整的腿方向信息不能丢失。

当前轨道04为 $(s_{31},s_{32};s_{41},s_{42})$，颜色为 $-C_A^2C_F/4$。05为 $(s_{31},s_{41};s_{32},s_{42})$，颜色为0；06为 $(s_{31},s_{42};s_{32},s_{41})$，颜色为 $+C_A^2C_F/4$。**05、06的边序同为3434，颜色却不同。** 颜色为零的图在本阶段筛除，不再计算一个可能奇异的标量J来与0相乘。

## 3. 正确的动量空间输入与两个δ

四条外腿动量为 $\ell_i$，内部动量为q。采用

$$\ell_1+\ell_2+q=0,\qquad \ell_3+\ell_4-q=0,$$

所以 $q=-(\ell_1+\ell_2)=\ell_3+\ell_4$。三胶子运动学顶点为

$$\Gamma_{\mu\nu\rho}(k_1,k_2,k_3)
=g_{\mu\nu}(k_1-k_2)_\rho+g_{\nu\rho}(k_2-k_3)_\mu
+g_{\rho\mu}(k_3-k_1)_\nu.$$

在左、右顶点分别代入 $(-\ell_1,-\ell_2,-q)$、$(-\ell_3,-\ell_4,q)$，并收缩四个Wilson切向量及内部度规。记所得标量为 $\mathcal L$。四个Wilson顶点的 $i^4=1$，两个三胶子顶点各含g，故整体耦合为 $g^6=(4\pi\alpha_s)^3$。

消去q以后，正确的紧凑动量式为

$$\boxed{\begin{aligned}
\mathcal W_\rho={}&g^6\mathcal C_\rho\int_{\mathcal D_\rho}d^4s
\int\prod_{i=1}^4\frac{d^D\ell_i}{(2\pi)^D}
(2\pi)^D\delta^D\!\left(\sum_i\ell_i\right)\\
&\times e^{-i\sum_i\ell_i\cdot X_i}
\prod_{i=1}^4\frac{-i}{\ell_i^2+i0}
\frac{-i}{(\ell_1+\ell_2)^2+i0}\,
\mathcal L(\ell_1,\ell_2,\ell_3,\ell_4).
\end{aligned}}$$

这里总体动量δ必须保留。L已经包含全部切向量收缩，标量传播子前不能再乘四个带自由指标的 $u_i^\mu$。等价的五动量写法包含两项 $(2\pi)^D\delta^D(\ell_1+\ell_2+q)$、$(2\pi)^D\delta^D(\ell_3+\ell_4-q)$，以及内部 $-i/(q^2+i0)$。

## 4. 两个顶点向量与完整二阶算子

直接缩并一个顶点得到

$$\boxed{V_L=u_{12}(\ell_2-\ell_1)
+u_1(2u_2\cdot\ell_1+u_2\cdot\ell_2)
-u_2(u_1\cdot\ell_1+2u_1\cdot\ell_2),}$$

$$\boxed{V_R=u_{34}(\ell_4-\ell_3)
+u_3(2u_4\cdot\ell_3+u_4\cdot\ell_4)
-u_4(u_3\cdot\ell_3+2u_3\cdot\ell_4).}$$

完整Lorentz分子就是 $\mathcal L=V_L\cdot V_R$。它有九个向量乘积；展开成标量单项式后是36项，未施加额外运动学等式之前不应与其他分组计数混淆。

为便于复核，定义

$$a_1=2u_2\cdot\ell_1+u_2\cdot\ell_2,\quad
 a_2=u_1\cdot\ell_1+2u_1\cdot\ell_2,$$

$$b_1=2u_4\cdot\ell_3+u_4\cdot\ell_4,\quad
 b_2=u_3\cdot\ell_3+2u_3\cdot\ell_4.$$

于是

$$\begin{aligned}
\mathcal L={}&u_{12}u_{34}(\ell_2-\ell_1)\cdot(\ell_4-\ell_3)\\
&+u_{12}\big[(u_4\cdot\ell_1-u_4\cdot\ell_2)b_2
 +(u_3\cdot\ell_2-u_3\cdot\ell_1)b_1\big]\\
&+u_{34}\big[(u_1\cdot\ell_4-u_1\cdot\ell_3)a_1
 +(u_2\cdot\ell_3-u_2\cdot\ell_4)a_2\big]\\
&+u_{13}a_1b_1-u_{14}a_1b_2-u_{23}a_2b_1+u_{24}a_2b_2.
\end{aligned}$$

因为 $\ell_i e^{-i\ell_i\cdot X_i}=i\partial_i e^{-i\ell_i\cdot X_i}$，二次分子统一给出

$$\boxed{\mathcal O=\mathcal L(i\partial_1,\ldots,i\partial_4)
=-\mathcal L(\partial_1,\ldots,\partial_4).}$$

此负号只计一次，已经进入O。L中的每个单项式都按同一规则变号，不能对不同组混用L和O。

## 5. Fourier变换到两个内部坐标

引入内部位置 $Y_L,Y_R$，用两个δ的Fourier表示拆开五条传播子。定义 $\nu=1-\epsilon$，采用本项目的FT约定

$$\int\frac{d^D\ell}{(2\pi)^D}\frac{-i}{\ell^2+i0}e^{-i\ell\cdot x}
=-\frac{\Gamma(\nu)}{4\pi^{D/2}}\frac1{(-x^2+i0)^\nu}.$$

五个坐标分母为

$$\begin{aligned}
\Delta_1&=-(X_1-Y_L)^2+i0,&\Delta_2&=-(X_2-Y_L)^2+i0,\\
\Delta_3&=-(X_3-Y_R)^2+i0,&\Delta_4&=-(X_4-Y_R)^2+i0,\\
\Delta_5&=-(Y_L-Y_R)^2+i0.
\end{aligned}$$

因此

$$\mathcal W_\rho=g^6\mathcal C_\rho
\left[-\frac{\Gamma(\nu)}{4\pi^{D/2}}\right]^5
\int_{\mathcal D_\rho}d^4s\;
\mathcal O\int d^DY_L\,d^DY_R\prod_{r=1}^5\Delta_r^{-\nu}.$$

O只对外部X求导，切向量 $u_i$ 固定；不能先把所有X限制为轮廓上的单变量后，再把D维偏导错误地当成 $d/ds_i$。

## 6. Feynman参数、高斯完成平方和处方

用五个非负α和**一个**单纯形约束，

$$\prod_{r=1}^5\Delta_r^{-\nu}
=\frac{\Gamma(5\nu)}{\Gamma(\nu)^5}
\int_{\alpha_r\ge0}d^5\alpha\,\delta\!\left(1-\sum_r\alpha_r\right)
\frac{\prod_r\alpha_r^{\nu-1}}{(\sum_r\alpha_r\Delta_r)^{5\nu}}.$$

设 $\Sigma_L=\alpha_1+\alpha_2+\alpha_5$、$\Sigma_R=\alpha_3+\alpha_4+\alpha_5$，则

$$\Sigma=\begin{pmatrix}\Sigma_L&-\alpha_5\\-\alpha_5&\Sigma_R\end{pmatrix},
\qquad\boxed{\mathcal U=\Sigma_L\Sigma_R-\alpha_5^2.}$$

定义 $Q_L=\alpha_1X_1+\alpha_2X_2$、$Q_R=\alpha_3X_3+\alpha_4X_4$、$K=\sum_{i=1}^4\alpha_iX_i^2$。完成平方给出几何多项式

$$B_0=(Q_L,Q_R)^T\operatorname{adj}(\Sigma)(Q_L,Q_R)-\mathcal U K,$$

其中二维矩阵指标和Lorentz指标分别收缩。写 $d_{ij}=-(X_i-X_j)^2$，可化为

$$\boxed{\begin{aligned}
B_0={}&\alpha_1\alpha_2\Sigma_R d_{12}
+\alpha_3\alpha_4\Sigma_L d_{34}\\
&+\alpha_5(\alpha_1\alpha_3d_{13}+\alpha_1\alpha_4d_{14}
+\alpha_2\alpha_3d_{23}+\alpha_2\alpha_4d_{24}).
\end{aligned}}$$

这里 $B_0$ 是α三次多项式，U为二次。若对原五个分母共同保留有限 $+i\eta$，完成平方后实际分母为

$$\boxed{B_\eta=B_0+i\eta\,\mathcal U\sum_r\alpha_r
=B_0+i\eta\,\mathcal U\quad\text{在单纯形上}.}$$

这是从原Feynman处方继承的形式。不能把给每个外部 $d_{ij}$ 随意补同一个 $+i0$ 当作有限η下完全相同的恒等式。符号 $B_F$ 以下表示 $\eta\downarrow0$ 的继承边界值；所有X导数只作用于 $B_0$。

在可以Wick旋转且积分收敛的初始解析区域，两次内部坐标积分贡献 $i^2=-1$：

$$\int d^DY_Ld^DY_R\prod_r\Delta_r^{-\nu}
=-\frac{\pi^D\Gamma(5\nu-D)}{\Gamma(\nu)^5}
\int d^5\alpha\,\delta\!\left(1-\sum\alpha\right)
\prod_r\alpha_r^{-\epsilon}\,
\mathcal U^{5\nu-3D/2}B_F^{D-5\nu}.$$

代入 $5\nu-D=1-3\epsilon$、$5\nu-3D/2=-1-2\epsilon$，令 $\lambda=1-3\epsilon$，得到

$$\boxed{\mathcal W_\rho=
\frac{g^6\mathcal C_\rho\Gamma(\lambda)}{1024\pi^{6-3\epsilon}}
\int_{\mathcal D_\rho}d^4s\int d^5\alpha\,
\delta\!\left(1-\sum\alpha\right)
\frac{\prod\alpha_r^{-\epsilon}}{\mathcal U^{1+2\epsilon}}
\mathcal O B_F^{-\lambda}.}$$

符号来源为：五个坐标传播子的 $(-1)^5$ 与两次Wick旋转的 $i^2$ 相乘为+；动量替换的 $i^2$ 已在O中，不能再次计入C。以上恒等式随后按维数正则作亚纯延拓；并未据此证明每个完整有序Wilson域在ε=0绝对收敛。

## 7. 二阶链式法则：两个分子都必须保留

定义梯度 $P_i^\mu=\partial_i^\mu B_0$。令

$$c_{ij}=\begin{cases}
\Sigma_R,&\{i,j\}=\{1,2\},\\
\Sigma_L,&\{i,j\}=\{3,4\},\\
\alpha_5,&i\in\{1,2\},\ j\in\{3,4\}\text{或相反}.
\end{cases}$$

则

$$P_i=-2\alpha_i\sum_{j\ne i}\alpha_jc_{ij}(X_i-X_j),$$

$$\partial_i^\mu\partial_j^\nu B_0=
\begin{cases}2\alpha_i\alpha_jc_{ij}g^{\mu\nu},&i\ne j,\\
-2\alpha_i\sum_{k\ne i}\alpha_kc_{ik}g^{\mu\nu},&i=j.
\end{cases}$$

因此 $\partial_i\cdot\partial_jB_0$ 在 $i\ne j$ 时为 $2D\alpha_i\alpha_jc_{ij}$。**D必须保留为 $4-2\epsilon$**；轮廓Gram矩阵的秩不是这个迹的维数。

把O写成常系数二阶算子，定义T1为两次导数分别打到两个B上所得的收缩，T2为两次都打到同一个B所得的收缩：

$$\boxed{\mathcal O B_F^{n}=n(n-1)B_F^{n-2}T_1+nB_F^{n-1}T_2.}$$

用上一节的顶点向量函数，可直接写成

$$\boxed{T_1=-V_L(P_1,P_2)\cdot V_R(P_3,P_4),\qquad T_2=\mathcal O B_0.}$$

这里 $V_L(P_1,P_2)$ 表示仅把顶点式的 $\ell_1,\ell_2$ 替换为梯度，切向量保持不变。这个T1自动包含原先遗漏的Type I，不能只保留切向量方向导数的乘积。

取 $n=-\lambda$，有

$$\boxed{\mathcal O B_F^{-\lambda}=
\lambda\left[(\lambda+1)\frac{T_1}{B_F^{\lambda+2}}
-\frac{T_2}{B_F^{\lambda+1}}\right],\qquad \lambda+1=2-3\epsilon.}$$

### 7.1 T2的紧凑显式公式

为避免再次展开36项，定义四个线性映射

$$\begin{aligned}
M_1&=-u_{12}\mathbf1+2u_1\otimes u_2-u_2\otimes u_1,&
M_2&=u_{12}\mathbf1+u_1\otimes u_2-2u_2\otimes u_1,\\
N_3&=-u_{34}\mathbf1+2u_3\otimes u_4-u_4\otimes u_3,&
N_4&=u_{34}\mathbf1+u_3\otimes u_4-2u_4\otimes u_3.
\end{aligned}$$

约定 $(u\otimes v)p=u(v\cdot p)$。于是 $V_L=M_1\ell_1+M_2\ell_2$、$V_R=N_3\ell_3+N_4\ell_4$。收缩两个映射与度规的结果只有两种：

$$K_+=(D-2)u_{12}u_{34}+5u_{13}u_{24}-4u_{14}u_{23},$$

$$K_-=-(D-2)u_{12}u_{34}+4u_{13}u_{24}-5u_{14}u_{23}.$$

这里的 $D-2$ 来自一项D维度规迹与两项交叉收缩之和，不是预先取了四维。因此

$$\boxed{T_2=-2\alpha_5\left[
(\alpha_1\alpha_3+\alpha_2\alpha_4)K_+
+(\alpha_1\alpha_4+\alpha_2\alpha_3)K_-\right].}$$

这个公式对重复边也成立，直接代入对应 $u_{ij}=0$ 即可。

## 8. 统一尺度提取与最终raw C、J

在 $v_e\mapsto r v_e,X_i\mapsto rX_i$ 下，$B_0\mapsto r^2B_0$、$P_i\mapsto rP_i$，故

$$T_1\mapsto r^6T_1,\qquad T_2\mapsto r^4T_2.$$

特别要区分：$u_p\cdot\partial_i(u_q\cdot\partial_jB_0)$ 带 $r^2$，而无切向量的 $\partial_i\cdot\partial_jB_0$ 带 $r^0$。算子前因子也有尺度，不能单看一个H就判定整个T2为尺度0。

定义

$$\widehat B=B_0/A_0,\qquad \widehat T_1=T_1/v_{\rm ref}^3,
\qquad \widehat T_2=T_2/v_{\rm ref}^2.$$

在初始正尺度区域，$v_{\rm ref}^3/A_0^3=-1/8$、$v_{\rm ref}^2/A_0^2=1/4$，于是最终统一分解为

$$\boxed{C_\rho^{\rm raw}=\mathcal C_\rho
\left(\frac{\alpha_s}{\pi}\right)^3
\frac{\pi^{3\epsilon}\Gamma(2-3\epsilon)}{16}
(-2v_{14}+i0)^{3\epsilon},}$$

$$\boxed{\begin{aligned}
J_\rho^{\rm raw}={}&\int_{\mathcal D_\rho}d^4s
\int_{\alpha_r\ge0}d^5\alpha\,\delta\!\left(1-\sum_{r=1}^5\alpha_r\right)
\frac{\prod_{r=1}^5\alpha_r^{-\epsilon}}{\mathcal U^{1+2\epsilon}}\\
&\times\left[-\frac{2-3\epsilon}{8}
\frac{\widehat T_1}{\widehat B_F^{3-3\epsilon}}
-\frac14\frac{\widehat T_2}{\widehat B_F^{2-3\epsilon}}\right].
\end{aligned}}$$

$\widehat B_F$继承第6节的处方；在正尺度区域为 $\widehat B+i0\,\mathcal U$。到复运动学时，须同时延拓尺度、z和核的边界值，不能把两因子独立取主分支。共同C中的 $\Gamma(2-3\epsilon)=\lambda\Gamma(\lambda)$ 来自链式法则；不能再遗漏J中的 $2-3\epsilon$ 或再次抽出同一个λ。

一个可相乘的C账本是：颜色 $\mathcal C_\rho$、耦合 $(\alpha_s/\pi)^3$、FT/Γ因子 $\pi^{3\epsilon}\Gamma(2-3\epsilon)/16$、尺度 $(-2v_{14}+i0)^{3\epsilon}$。此前相位的抵消存为来源说明；账本中不额外放 `I^6` 或 `v14^2`。

后续预处理再做固定颜色基底投影和一次 $e^{-3\gamma_E\epsilon}$ 补偿。本页的raw C不提前包含它们。轨道重数和crossing也应由图的独立映射记录决定，不能因为使用单一规范代表就删除其物理旋转贡献。

### 8.1 一个方便实现的无量纲向量配方

取形式基 $E_e$，内积 $E_e\cdot_*E_f=(v_e\cdot v_f)/v_{\rm ref}$，只使用内积，不求这个退化Gram矩阵的逆。令

$$\xi_i=\sum_{r<e_i}E_r+s_iE_{e_i},\qquad t_i=E_{e_i},$$

$$g_i=-2\alpha_i\sum_{j\ne i}\alpha_jc_{ij}(\xi_i-\xi_j).$$

由 $d_{ij}/A_0=(\xi_i-\xi_j)\cdot_*(\xi_i-\xi_j)/2$ 算出 $\widehat B$；在第4节的向量式中用 $u_i\to t_i$、$\ell_i\to g_i$、点积改为 $\cdot_*$，即可算出 $\widehat T_1$。T2的紧凑式则将 $u_{ij}$ 改为 $t_i\cdot_*t_j$，仍保留独立的 $D=4-2\epsilon$。

检查α齐次次数应得到：U为2，B及各g为3，T1为6，T2为3；两个J被积项的总α齐次次数都为−5。这是五参数投影测度的必要一致性检查。

## 9. 原有序域与可选超立方变换

默认输出保留第1节的原s变量、依赖上界和明确外到内顺序。若一条边的链为 $0\le s_1\le\cdots\le s_m\le1$，可另行使用

$$s_j=\prod_{r=j}^{m}x_r,\qquad 0\le x_r\le1,
\qquad \boxed{\left|\frac{\partial(s_1,\ldots,s_m)}{\partial(x_1,\ldots,x_m)}\right|
=\prod_{r=2}^{m}x_r^{r-1}.}$$

例如两点链为 $s_1=x_1x_2,s_2=x_2$，Jacobian为 $x_2$。不能一边声明升序s，一边使用降序映射 $s_j=x_1\cdots x_j$ 而不重标记。

α始终是五变量单纯形；转换Wilson参数为立方体并不会把α变成五个独立单位区间，也不会消去δ约束。

同边两点满足 $d_{ij}=-(s_i-s_j)^2v_{e_i}^2=0$。这只消去B的对应距离项，不自动证明整图为零；整图零性仍由颜色和完整Lorentz收缩判断。

## 10. 两个明确基准

### 10.1 04：3344，两同边对

取 $u_1=u_2=v_3$、$u_3=u_4=v_4$，原域 $s_{31}\le s_{32}$、$s_{41}\le s_{42}$。定义

$$A=\alpha_1(1-s_{31})+\alpha_2(1-s_{32}),\qquad
B=\alpha_3s_{41}+\alpha_4s_{42},\qquad
\delta_\alpha=(\alpha_1-\alpha_2)(\alpha_3-\alpha_4).$$

则

$$\widehat B_5=\alpha_5 zAB,\qquad
\widehat K_+=z^2,\qquad\widehat K_-=-z^2,$$

$$\widehat T_2=-2\alpha_5z^2\delta_\alpha,\qquad
\widehat T_1=4\alpha_5^2z^3AB\delta_\alpha
=-2\widehat B_5\widehat T_2.$$

所以正确的J核括号可进一步化成

$$\boxed{-\frac{(1-3\epsilon)\alpha_5z^2\delta_\alpha}
{2\widehat B_{5,F}^{2-3\epsilon}}.}$$

化简先在不穿越奇点的解析区域进行，再取继承的边界值。C中的共同Γ(2−3ε)保持不变；这里又出现的λ来自此图两项分子的相消，是特例结果，不能推广成所有图的额外退化因子。

与本轮修正文档之前的04导出相比，使用相同C和测度，正确核恰为旧核的 $(1-3\epsilon)/4$。这只是核连接检查，不表示已经求得04的I。

### 10.2 25：1234，必须有D维迹项

此时 $u_{12}/v_{\rm ref}=u_{34}/v_{\rm ref}=z$，$u_{13}/v_{\rm ref}=u_{24}/v_{\rm ref}=-(1+z)$，$u_{14}/v_{\rm ref}=u_{23}/v_{\rm ref}=1$。因此

$$\boxed{\widehat K_+=(7-2\epsilon)z^2+10z+1,
\qquad\widehat K_-=(2+2\epsilon)z^2+8z-1.}$$

$$\widehat T_2=-2\alpha_5\left[
(\alpha_1\alpha_3+\alpha_2\alpha_4)\big((7-2\epsilon)z^2+10z+1\big)
+(\alpha_1\alpha_4+\alpha_2\alpha_3)\big((2+2\epsilon)z^2+8z-1\big)\right].$$

T1按第8.1节向量配方给出，无需抄写数百项多项式。这两个显式的ε系数可直接检查实现是否漏掉D维迹。

## 11. 后续实现应满足的连接要求

每张图保存腿→边/端点/动量/颜色的完整映射、有向顶点配对、轮廓颜色次序和有序域；不能只存边占据数。由两个顶点向量生成O，或实现等价的T1/T2紧凑配方；若另存36项表，须逐项与顶点乘积比较。

最终J必须显式包含原s测度、五α非负单纯形、U、B、两类分子及全部ε系数。变换到超立方时另外记录映射和Jacobian。颜色零图保留零的原因与源身份。

本轮已将新紧凑T1/T2与审计中逐项二阶算子的独立结果对照全部26个轨道，并核对104条原相位到位置的连接。另检查了D=4/6分量矩阵迹、顶点内换腿与左右顶点互换、两种核的α次数和尺度、正序Jacobian、颜色以及04/25基准，共194项精确检查通过。

这些是计算侧公式核验，formal=not_run。没有修改montopar代码、重新发布正确C/J或开始求积分；也未完成一般物理复运动学的数值延拓验证。后续实现必须按本版公式重新检查，不能继承旧输出的状态。
