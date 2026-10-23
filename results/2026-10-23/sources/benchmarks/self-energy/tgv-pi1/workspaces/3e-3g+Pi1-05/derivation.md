# 3e-3g+Pi1 05的完整解析解

## 当前物理积分的入口（2026-10-15，r003）

从共同作用量出发，U′=igAU确定真实Wilson词；FT/Dyson规则在提取C前生效。当前[原生C/J](../../sources/montopar-r002.json)给出下式，J就是本文实际积分的原核，域、分支和IBP全保留：

$$R_P=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad I_P=2R_P\frac{J(z)+z^{3\epsilon}J(1/z)}{1+z^{3\epsilon}}.$$

用下文已经求出的J代入此定义，得到[当前物理I](../../results/3e-3g+Pi1-05/r003.json)；组装系数为+1，Euler只补一次。以下带日期的求值过程与r001公式属于旧根历史，计算J的步骤可按同域积分线性复用，旧形式状态不继承。

2026-10-15表示修订：主I、原J和主Laurent系数保持r002的值；当前源哈希、配对引用及附属公式统一到新修订。[修订说明](../../revisions/r003.md)。以下2026-10-01/02的r001系数表保留历史符号，不能直接当作当前结果。

2026-10-02。05属于334核心，Π₁在边3后槽，r=2。本页从[实际raw包](../../sources/montopar-r001.json)连接到已经完整求出的02，给出全ε等式 $I_{05}=I_{02}$ 及完整结果。核验使用本图自己的标记权重、两个完整分子、原域和C，不由其他图的反射关系自动推断。

[结果r001](../../results/3e-3g+Pi1-05/r001.json) · [02原积分求值推导](../3e-3g+Pi1-02/derivation.md)。

## 1. 显式变量变换

将本图路径参数记为A=s₃₁、B=s₃₂、C=s₄₁，原域为0≤A≤B≤1、0≤C≤1，实际外到内次序dB dA dC。本图α参数另记β₁、β₂、β₃，均非负且和为1。

将参考02号的344路径参数记为s=s₃₁、t=s₄₁、u=s₄₂。定义从本图到参考图的映射

$$s=1-C,\qquad t=1-B,\qquad u=1-A,$$
$$\alpha_1=\beta_3,\qquad\alpha_2=\beta_2,\qquad\alpha_3=\beta_1.$$

于是t≤u恰好等于A≤B，区间端点都映到原区间；路径Jacobian为+1，消去第三个α之后的单纯形Jacobian绝对值为1。

参考图实际标记权重为
$$w_{02}=\alpha_1^{-\epsilon}\alpha_2^{-2\epsilon}\alpha_3^{-\epsilon},$$
经映射恰为本图的
$$w_{05}=\beta_1^{-\epsilon}\beta_2^{-2\epsilon}\beta_3^{-\epsilon}.$$

这同时核对了自能位置，不能只交换未标记的几何腿。z及参考尺度保持。

## 2. 原分母、完整分子与有序积分的连接

参考图的分母多项式拉回为

$$z\alpha_1(1-s)(\alpha_2t+\alpha_3u)
=z\beta_3 C[\beta_1(1-A)+\beta_2(1-B)],$$

即本图q₀。实际JSON中两个完整分子分别满足

$$K_{1,02}\circ\phi=K_{1,05},\qquad
K_{2,02}\circ\phi=K_{2,05}.$$

因此两个分母幂2−3ε、3−3ε独立相符；在相同有限iη下仍是完整带权核的逐点等式，没有删除调节项。全部全颜色C、Gamma、相位、尺度和raw归一化字段亦相同，C4重数同为4。

参考图的嵌套次序经变换成为C外层、A次层、B从A到1。本图采用B外层、A从0到B、C独立。先在有限共同η与α内部截断下换序；再利用参考图在z>0、A₀>0、0<Re ε<1/3的绝对可积控制，通过上述双射传给本图。因而域换序及共同边界值合法，之后按同一ε亚纯延拓解释。

也已从本图直接代入A=1−x、B=1−xp、C=u及β₃=a、β₂=(1−a)v、β₁=(1−a)(1−v)。路径Jacobian=x；提取两个端点尺度后，q和两个完整约化分子与参考图已经求解的核精确一致。此检查补足了原有序域到已求解析表达的连接。

于是
$$J_{05}(\epsilon,z)=J_{02}(\epsilon,z).$$

## 3. 完整I及颜色恢复

沿用单一归一化
$$R=\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}
{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},$$
$$I=2R\frac{J(\epsilon,z)+z^{3\epsilon}J(\epsilon,z^{-1})}{1+z^{3\epsilon}},$$
便得到全阶I相等，而不仅是若干Laurent系数相同。

记
$$P_2=36\epsilon^4+18\epsilon^3-28\epsilon^2-6\epsilon+4,\quad
Q_2=18\epsilon^4+13\epsilon^3-9\epsilon^2-2\epsilon-2.$$
则完整结果是
$$I_{05}=I_{02}=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2}
{2304\epsilon^5(1+\epsilon)(1+2\epsilon)(1-2\epsilon)(1-3\epsilon)(2-3\epsilon)}
\left[P_2\Gamma(1-3\epsilon)\Gamma(1+2\epsilon)+Q_2\Gamma(1-\epsilon)\right].$$

唯一I与z无关，最高五阶。这个定义包含一次Euler补偿和规定的C4组装；原raw和恢复的轨道贡献仍带共同外尺度。

两个颜色通道的恢复式为
$$\mathcal W_{05}^{prepared}=\lambda_b^3(S_{14}+S_{12})
\left[C_FC_A^2\frac{10-6\epsilon}{3-2\epsilon}-C_FC_AT_Fn_f\frac{8(1-\epsilon)}{3-2\epsilon}\right]I_{05}.$$

各通道的完整有理系数卷积保存在结果中；共同外尺度未展开，完整w₃/RG尚未组装。源轨道仍单独计权，不因I相等删除物理贡献。

## 4. 全部极点至有限项

$$I_{05}=\sum_{k=-5}^0c_k\epsilon^k+O(\epsilon).$$

| k | c_k |
|---:|---|
| -5 | $- \frac{1}{2304}$ |
| -4 | $\frac{1}{4608}$ |
| -3 | $\frac{65}{9216} - \frac{\pi^{2}}{1024}$ |
| -2 | $- \frac{31 \pi^{2}}{18432} - \frac{13 \zeta\left(3\right)}{2304} + \frac{55}{2048}$ |
| -1 | $- \frac{199 \pi^{4}}{221184} - \frac{7 \pi^{2}}{4096} - \frac{47 \zeta\left(3\right)}{4608} + \frac{3925}{36864}$ |
| 0 | $- \frac{259 \pi^{4}}{147456} - \frac{23 \pi^{2} \zeta\left(3\right)}{3072} - \frac{47 \zeta\left(5\right)}{1280} - \frac{127 \zeta\left(3\right)}{9216} - \frac{17 \pi^{2}}{73728} + \frac{25199}{73728}$ |

## 5. 证据及限制

本图实际完成33项连接/系数检查，其中有6个带共同调节的raw数值点。参考图的Gamma/Beta求值按同一原输入版本及明确等式复用，旧检查不重复计数。06的两个颜色通道十二个导数控制是本轮新增；05复用02已经核对的两颜色卷积。

[统一检查与发布记录](../../../../../workflow/checks/tgv-pi1-05-06/README.md)。没有可调用的Wilson数学MCP，formal=not_run；原输入、旧结果和生产/Kernel/Lean源保持。数值残差不作为严格误差包围，一般复运动学延拓和原高维独立数值检验仍未完成。
