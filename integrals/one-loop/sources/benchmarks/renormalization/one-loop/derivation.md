# 一圈全链路：从图、C/J与I到重整化Wilson圈及cusp

日期：2026-09-24。计算侧一圈链路完成；原有两个I的定义和结果不变。[机器结果](result-r001.json) · [计算与来源检查](../../../workflow/checks/rg-one-loop/README.md) · [后续计划](../../../workflow/renormalization-cusp-plan.md)。

本次实际重放qgraf-pipeline的组合生成、分类和表达式导出API，以及montopar的纯交换C/J生产函数；不是重新调用外部QGRAF二进制。预处理也已执行，47条默认根的四个发布文件与执行前字节一致。本页依赖通常的维数正规化与Feynman传播子规则；没有新增Lean/MCP收据，formal=not_run。

## 1. 几何与对象

基本表示、无自交的闭合类光四边形满足
$$
v_i^2=0,\quad \sum_i v_i=0,\qquad
x_1=0,\ x_2=v_1,\ x_3=v_1+v_2,\ x_4=-v_4.
\tag{R1.1}
$$
沿边点为X_i(u)=x_i+u v_i，0≤u≤1。定义
$$
s=x_{13}^2=2v_{12},\quad t=x_{24}^2=2v_{14},\quad z=s/t=v_{34}/v_{14}.
\tag{R1.2}
$$

先在解析不变量区域s,t<0、z>0计算，不要求由实Minkowski四向量实现。初始取0<Re epsilon<1，再在固定z>0下作epsilon=0附近的亚纯延拓。一般复z和物理i0延拓不属于本次验证范围。

固定D=4-2epsilon、归一化迹Tr/Nc，g0为有质量维度的裸耦合：
$$
\lambda_0=\frac{g_0^2}{4\pi^2},\quad
\lambda=\frac{\alpha_s^{\overline{\rm MS}}(\mu)}{\pi},\quad
G(\epsilon)=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon).
\tag{R1.3}
$$
本轮只求lambda一次项，因此W-1与log W相同。普通QCD耦合反项首次贡献于lambda平方；本阶需要轮廓cusp减除，不需要Pi1。

## 2. qgraf-pipeline：六个边对与两个代表

Wilson有序指数展开并收缩两场后，利用二点核对称性，可写成
$$
W_B^{(1)}=-\frac{g_0^2C_F}{2}
\sum_{i,j=1}^4\int_0^1du\,dv\,
v_i^\mu v_j^\nu D_{\mu\nu}(X_i(u)-X_j(v)).
\tag{R1.4}
$$
直接使用有序路径积分时不再加这个1/2；不同边的贡献等于无序边对i<j各一次。同边项按类光维数正规化规则为零，其Lorentz收缩含v_i²。这是所采用的QFT规则，不是把分布奇点上的零乘无穷作为普通函数证明。

| 源ID | 代表 | 完整C4轨道 | m |
|---|---|---|---:|
| 2e-01 | (3,4) | (1,2)、(2,3)、(3,4)、(1,4) | 4 |
| 2e-02 | (2,4) | (1,3)、(2,4) | 2 |

实际来源：[冻结动量公式](../../../workflow/checks/rg-one-loop/sources/upstream/qgraf-pipeline/outputs/1loop/integrals.md)、[重放定义](../../../workflow/checks/rg-one-loop/generated/1loop/orbit_definitions.json)。独立枚举六边对并核对两轨道无遗漏、无交集。metadata中的n_diagrams=1不是C4重数。

## 3. montopar：从动量式得到C/J

代表图的原始贡献为
$$
\mathcal A_{ij}=\frac1{N_c}{\rm tr}
\int_0^1du\,dv\,(ig_0v_j^\nu T^b)(ig_0v_i^\mu T^a)
\int\frac{d^D\ell}{(2\pi)^D}
\frac{-ig_{\mu\nu}\delta^{ab}}{\ell^2+i0}
e^{-i\ell\cdot(X_i-X_j)}.
\tag{R1.5}
$$
源文件记Phi=-(Xi-Xj)并用exp(+i ell·Phi)，两图相位均已独立核对。颜色缩并delta_ab Tr(Ta Tb)/Nc=CF。

采用(+---)约定的坐标传播子
$$
D_{\mu\nu}^{ab}(X)=
-\delta^{ab}g_{\mu\nu}\frac{\Gamma(1-\epsilon)}{4\pi^{2-\epsilon}}
(-X^2+i0)^{-1+\epsilon}.
\tag{R1.6}
$$
系数由Schwinger热核与Feynman边界值给出，与[0709.2368式(8)](https://arxiv.org/pdf/0709.2368)一致。两个Wilson顶点的i²与此负号合并，得到
$$
\mathcal A_{ij}=
g_0^2C_F(v_i\cdot v_j)\frac{\Gamma(1-\epsilon)}{4\pi^{2-\epsilon}}
\int_0^1du\,dv\,[-(X_i-X_j)^2+i0]^{-1+\epsilon}.
\tag{R1.7}
$$

令b=v14并提出尺度-2b：
$$
C_{\rm src}=-\frac12C_F\lambda_0\Gamma(1-\epsilon)(-2\pi b+i0)^\epsilon,
\qquad J=\int d\boldsymbol s\,d\,f^{-1+\epsilon}.
\tag{R1.8}
$$

| 代表 | f=(Xi-Xj)²/(2b) | d=(vi·vj)/b | 原外到内顺序 |
|---|---|---|---|
| 相邻 | z(1-s31)s41 | z | s31→s41 |
| 对边 | (1-s21)(1-s41)+z s21 s41 | -(1+z) | s21→s41 |

各参数均取[0,1]。montopar实际重建的C、f、d、边界、顺序、C4信息与[prepared 2e](../../2ne/prepared/2e.json)相符。本阶cube变换只是变量改名，不产生额外Jacobian。预处理加入一次Euler补偿：
$$
C_{\rm prep}=e^{-\gamma_E\epsilon}C_{\rm src}.
\tag{R1.9}
$$
它是计算归一化，尚未做MS-bar重整化。

## 4. 从原J连接到现有唯一I

相邻图按原顺序：
$$
J_A=z^\epsilon\int_0^1ds_{31}(1-s_{31})^{\epsilon-1}
\int_0^1ds_{41}s_{41}^{\epsilon-1}
=\frac{z^\epsilon}{\epsilon^2}.
\tag{R1.10}
$$
对边u=s21、v=s41，先积v：
$$
J_O=-(1+z)\int_0^1du\,
\frac{(zu)^\epsilon-(1-u)^\epsilon}{\epsilon[(1+z)u-1]}
=z^{\epsilon/2}H_0(\epsilon,z).
\tag{R1.11}
$$
u=1/(1+z)是可去差商；须保留连续值，不拆成各自有内部极点的积分。

复用[2e-02已有推导](../../2ne/families/single-opposite/basis.md)：
$$
H_0=-\frac{(1+z)z^{-\epsilon/2}}{\Gamma(1-\epsilon)\Gamma(1+\epsilon)^2}
G_{3,3}^{3,3}\!\left(z\;\middle|\;
\begin{matrix}\epsilon,0,0\\0,\epsilon,\epsilon\end{matrix}\right).
\tag{R1.12}
$$
本次以原参数差商积分独立数值复核此全epsilon表达。

现有系数提取契约为
$$
\frac m2 C_{\rm prep}[J(z)+z^\epsilon J(1/z)]
=C_F\lambda_0[(-2\pi b)^\epsilon+(-2\pi bz)^\epsilon]I,
\tag{R1.13}
$$
从而
$$
\boxed{I_A=-\frac G{\epsilon^2},\qquad
I_O=-\frac{G H_0(\epsilon,z)}{2\cosh(\epsilon\log z/2)}.}
\tag{R1.14}
$$
此处lambda0是将原计算标签还原为裸耦合时的记号；(R1.13)仍带prepared归一化。下一节统一连接到真实裸图。I不改定义，也不再重复乘m。

来源：[2e-01/r002](../../2ne/results/2e-01/r002.json)、[2e-02/r002](../../2ne/results/2e-02/r002.json)。旧输入、代码及证据哈希已重查，未改旧结果或赋予新形式状态。

## 5. 固定MS-bar尺度字典

选定物理耦合关系
$$
\lambda_0(4\pi)^\epsilon e^{-\gamma_E\epsilon}
=\mu^{2\epsilon}Z_\alpha\lambda,\qquad Z_\alpha=1+O(\lambda).
\tag{R1.15}
$$
同一raw因子因此满足
$$
\lambda_0\Gamma(1-\epsilon)\pi^\epsilon(-t)^\epsilon
=\lambda G(\epsilon)
\left[\frac{e^{2\gamma_E}}4\mu^2(-t)\right]^\epsilon+O(\lambda^2).
\tag{R1.16}
$$
故尺度关系确定为
$$
\boxed{\mu_x^2=\frac{e^{2\gamma_E}}4\mu^2,\qquad
\mu_{\rm prep}^2=\frac{e^{2\gamma_E}}{4\pi}\mu^2,\qquad
\mu_x^2=\pi\mu_{\rm prep}^2.}
\tag{R1.17}
$$
mu是(R1.15)的MS-bar尺度；mu_x是论文式(9)重定义后的坐标尺度。先逆用(R1.9)恢复raw，再代入裸耦合，正好得到(R1.16)。这不修改I，也不再次补偿Euler。

定义
$$
L_s=\log[\mu_x^2(-s)],\quad L_t=\log[\mu_x^2(-t)],\quad
L=(L_s+L_t)/2,\quad \rho=L_s-L_t=\log z.
\tag{R1.18}
$$
恢复物理前因子并求和：
$$
\boxed{\left.\log W_B\right|_\lambda
=C_F\lambda(e^{\epsilon L_s}+e^{\epsilon L_t})(I_A+I_O)
=-C_F\lambda G e^{\epsilon L}
\left[\frac{2\cosh(\epsilon\rho/2)}{\epsilon^2}+H_0(\epsilon,z)\right].}
\tag{R1.19}
$$
W_B此时以lambda、mu表示，仍含轮廓极点。普通耦合反项影响O(lambda²)，本阶不加入。

## 6. 有限项与严格MS减除

已有母函数及Gamma展开：
$$
H_0(0,z)=M_{0,1}=-\frac{\log^2z+\pi^2}{2},
\qquad G=1+\frac{\pi^2}{12}\epsilon^2+O(\epsilon^3).
\tag{R1.20}
$$
相邻/对边贡献分别为
$$
\frac{W_A^{(1)}}{C_F\lambda}
=-\frac2{\epsilon^2}-\frac{L_s+L_t}{\epsilon}
-\frac{L_s^2+L_t^2}{2}-\frac{\pi^2}{6}+O(\epsilon),
\tag{R1.21}
$$
$$
\frac{W_O^{(1)}}{C_F\lambda}
=\frac{(L_s-L_t)^2+\pi^2}{2}+O(\epsilon).
\tag{R1.22}
$$
因此
$$
\boxed{\log W_B=C_F\lambda
\left[-\frac2{\epsilon^2}-\frac{L_s+L_t}{\epsilon}
-L_sL_t+\frac{\pi^2}{3}+O(\epsilon)\right]+O(\lambda^2).}
\tag{R1.23}
$$

固定W_B=Z4 W_R，只取Laurent负幂：
$$
\boxed{\log Z_4^{(1)}
=-C_F\lambda\left[\frac2{\epsilon^2}+\frac{L_s+L_t}{\epsilon}\right].}
\tag{R1.24}
$$
实际减除是乘Z4的逆；不可翻转此定义中的号。得到
$$
\boxed{\log W_R=C_F\frac{\alpha_s^{\overline{\rm MS}}}{\pi}
\left[-L_sL_t+\frac{\pi^2}{3}\right]+O(\alpha_s^2).}
\tag{R1.25}
$$
对边图虽无UV极点，对有限项仍必不可少。

## 7. 两条路线提取RG与cusp

固定几何时，mu导数作用于两个尺度对数：
$$
\mathcal D_\mu=\mu\partial_\mu+\beta_\lambda\partial_\lambda,\quad
\mu\partial_\mu L_s=\mu\partial_\mu L_t=2,\quad
\beta_\lambda=O(\lambda^2).
\tag{R1.26}
$$
从有限结果直接得到
$$
\boxed{\mathcal D_\mu\log W_R
=-2C_F\lambda(L_s+L_t)+O(\lambda^2).}
\tag{R1.27}
$$
独立地，对Z4使用D维耦合导数：
$$
\beta_\lambda^{(D)}=-2\epsilon\lambda+O(\lambda^2),\qquad
\mathcal D_\mu^{(D)}\log Z_4^{(1)}
=2C_F\lambda(L_s+L_t).
\tag{R1.28}
$$
这里显式mu导数与经典维数项产生的1/epsilon相消。对Z4仅用四维beta会出错。

两条路线给出相同结论：
$$
\boxed{K_F=C_F\frac{\alpha_s}{\pi}+O(\alpha_s^2),\qquad
\gamma_4^{(1)}=0,\qquad
\Gamma_0=4C_F\quad(a_s=\alpha_s/(4\pi)).}
\tag{R1.29}
$$
gamma4=0采用(R1.18)的对数定义；改变对数中的常数尺度会重分配非cusp常数，不能混用。cusp不含nf，与本阶自由交换的物理内容一致。

全epsilon式还精确满足
$$
\left(2\partial_L-2\epsilon\lambda\partial_\lambda\right)
\left.\log W_B\right|_\lambda=0.
\tag{R1.30}
$$
裸对象的D维尺度独立性，与重整化对象的非零反常演化并不冲突。

## 8. 与0709.2368的一圈匹配

将(R1.19)写成
$$
\log W_B=C_F\lambda\left[
-\frac{e^{\epsilon L_s}+e^{\epsilon L_t}}{\epsilon^2}
+\frac12\log^2z+\frac{\pi^2}{3}+O(\epsilon)\right]+O(\lambda^2),
\tag{R1.31}
$$
即与[0709.2368式(10)、(17)](https://arxiv.org/pdf/0709.2368)一致。本阶没有内部物种回路，QCD/N=4比较有效；不由此推断两圈自能相同。

该文保留尺度幂的整块减除给出的有限形状项是
$$
\mathcal F_4^{(1)}=C_F\lambda[\tfrac12\log^2z+\tfrac{\pi^2}{3}].
\tag{R1.32}
$$
与严格MS有限结果的关系为
$$
\boxed{\log W_R^{(1)}
=\mathcal F_4^{(1)}-\frac{C_F\lambda}{2}(L_s^2+L_t^2).}
\tag{R1.33}
$$
因此mathcal F4对整体尺度不变，而log W_R有(R1.27)的演化。两者不能直接互换。

论文planar变量a=g²N/(8pi²)对应C_F lambda→a，论文Gamma_cusp=2a，本页基本表示K_F=C_F lambda。对尚含cusp极点的对象，
$$
\frac{d}{d\log\mu_x^2}\left.\log W_B\right|_\lambda
=-C_F\lambda[2/\epsilon+L_s+L_t]+O(\epsilon,\lambda^2),
\tag{R1.34}
$$
与论文式(29)的一圈版本相符；减除Z4后才得到有限RG。

## 9. 供后续反项复用的正epsilon幂

旧r002不改。本次另保存完整组装系数至epsilon²，余项O(epsilon³)，复用已有M/Hlog，不新增未求值母积分。

在(R1.11)令beta=1-u，得到
$$
H_0(\epsilon,z)=\int_0^1d\beta\,
[\beta(1-\beta)]^{\epsilon/2}
\frac{2\sinh(\epsilon h/2)}{\epsilon[\beta-z/(1+z)]},
\quad h=\log z+\log\frac{1-\beta}{\beta}.
\tag{R1.35}
$$
内部点h=0，取连续值。按[原M定义及求值](../../2ne/families/single-opposite/basis.md)：
$$
H_0=M_{0,1}+\frac{\epsilon}{2}M_{1,1}
+\epsilon^2\left(\frac{M_{2,1}}8+\frac{M_{0,3}}{24}\right)+O(\epsilon^3).
\tag{R1.36}
$$
本次精确核对这与既有H0_q_plus/Hlog文件一致，并独立数值核对定义矩。

令
$$
E_j=(L_s^j+L_t^j)/j!,\quad
g_2=\pi^2/12,\quad g_3=\zeta(3)/3,\quad g_4=\pi^4/160.
\tag{R1.37}
$$
写log W_B/(C_F lambda)=sum c_k epsilon^k：

| k | c_k |
|---:|---|
| -2 | $-2$ |
| -1 | $-E_1$ |
| 0 | $-E_2-2g_2-M_{0,1}=-L_sL_t+\pi^2/3$ |
| 1 | $-E_3-g_2E_1-2g_3-LM_{0,1}-M_{1,1}/2$ |
| 2 | $-E_4-g_2E_2-g_3E_1-2g_4-(L^2/2+g_2)M_{0,1}-LM_{1,1}/2-M_{2,1}/8-M_{0,3}/24$ |

固定z>0时，对边D只在两个角点为零，局部与两个非负坐标之和可比。delta<1时D^(-1-delta)乘任意有限次log幂二维可积，允许邻域内逐阶求导。相邻部分由显式Gamma/epsilon²给出亚纯结构，因而这里不是在发散积分内直接令epsilon=0。

## 10. 验证与范围

116项来源/精确检查及29项数值检查通过。涵盖上游重放、六边对、相位/Gram、C/J到I、尺度字典、旧证据、Hlog连接、完整有限项、两种RG提取、DHKS减除转换及错误符号/重数/Euler负控制。

数值部分包括70位工作精度的矩、Meijer G和原参数降维积分，以及两组双精度SciPy原二维积分；epsilon²余项呈O(epsilon³)趋势。它们是辅助核对，不是区间认证。

计算侧已从两条实际图来源走到完整一圈W_R和cusp。QFT Fourier边界值、同边规则、一般复运动学及新增RG关系没有Lean整链；当前工具列表没有Wilson MCP，formal保持not_run。源输入、旧I、生产代码及其他积分任务未改写。
