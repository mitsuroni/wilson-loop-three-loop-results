# 4e-Pi1-07：I求解工作区

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1/workspaces/4e-Pi1-07/derivation.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。

以下保留迁移前的求解/准备过程记录；其中r001的C/I公式不作为现行数值入口。现行理论及解表已经更新，J与母函数的推导仍可复用。

所属[P5：两条对边线：同一通道、有序耦合](../../families/p5-opposite-same/README.md)。本图只使用[主推导](../../theory.md)定义的主稿式(31)的I。

## 固定输入

- 来源：[不可变JSON包](../../sources/4e-Pi1-r001.json)中的`4e-Pi1-07`。
- 源SHA-256：`3b2cab8babd9b1f4630d77d01dc64dffb9b6868c90ccab950a55c4c30230ee0a`。
- 总圈数：3；Wilson颜色迹：`CF*(CF-CA/2)`。
- 以专用prepared C进行a₁/Euler和C4系数提取；完整CA/NF颜色因子位于外部前因子，I不含颜色和尺度。
- 外到内：`s22 → s21 → s42 → s41`。

| 原变量 | 下限 | 上限 |
|---|---|---|
| s22 | 0 | 1 |
| s21 | 0 | s22 |
| s42 | 0 | 1 |
| s41 | 0 | s42 |

## 本图推导及同组连接

只使用[唯一I定义](../../theory.md)式(31)，$e=\varepsilon$、$L=\log z$、$G=e^{-\gamma_Ee}\Gamma(1-e)$。两图 $n=3,a_1=-1/2,\mathrm{sym}=2$。完整CA/NF算符在共同外部因子中，Euler补偿只计一次。

## 1. 原域与完整核

令 $x=s_{21},u=s_{22},y=s_{41},v=s_{42}$，原域为 $0<x<u<1,0<y<v<1$，外到内顺序 $u\to x\to v\to y$。记

$$
D(x,y)=(1-x)(1-y)+zxy,\qquad g_e(x,y)=-(1+z)D(x,y)^{-1+e},
$$
$$
K_e(x,y)=\frac{-2(1+z)D(x,y)+(1-2e)z}{2(1+e)}D(x,y)^{-2+2e}.
$$

07为 $K_e(x,y)g_e(u,v)$，08为 $g_e(x,y)K_e(u,v)$。两核共享排序域，不能各自积遍整个方形后相乘。反射 $(x,u,y,v)\mapsto(1-u,1-x,1-v,1-y)$ 保留原域且绝对Jacobian为1；距离/完整分子分别带z/z²因子，因此 $J_{08}(z)=z^{3e}J_{07}(1/z)$。

Re(e)>0时完整立方体上的SE核和裸核分别绝对可积，限制至排序域仍可积，故可按Fubini重组积分以使用完整核恒等式

$$
K_e=g_{2e}+\frac{\partial_x\partial_yD^{2e}}{4e(1+e)}.
$$

07的导数核积遍内矩形后给出 $D(u,v)^{2e}-(1-u)^{2e}-(1-v)^{2e}+1$；08则在固定内点后积外矩形，给出 $z^{2e}-(zx)^{2e}-(zy)^{2e}+D(x,y)^{2e}$。八个角项全部保留。

## 2. 已积完的不同指数双核级数

使用[原同通道级数](../../../../2ne/families/same-channel-two/method.md)的有限块，并保留两条线不同的指数。以下定义的是显式数值级数，没有未求值参数积分：

$$
\mathcal T(a,b;z)=\left(\frac{1+z}{\sqrt z}\right)^{a+b}
\sum_{k,l\ge0}\frac{(1-a)_k(1-b)_l}{k!l!}\,B_{kl}(z),
$$
$$
B_{kl}(z)=\sum_{i,j=0}^k\sum_{m,n=0}^l
\frac{p_k(i,j;z)p_l(m,n;z)}{(i+1)(i+m+2)(j+1)(j+n+2)},
$$
$$
p_k(i,j;z)=\sum_{h=\max(0,i+j-k)}^{\min(i,j)}
\frac{k!(-1)^h z^{k-i-j+h}}
{(k-i-j+h)!(i-h)!(j-h)!h!(1+z)^{k-h}}.
$$

这里p是 $[1-D/(1+z)]^k$ 的有限多项式系数，B的四个分母来自原排序域上的单项式积分。因此在共同收敛域，T等于 $z^{-(a+b)/2}(1+z)^2$ 乘以两条有序传播子的积分（指数分别a−1、b−1），并非另行假设的独立域乘积。

**求和顺序：先算每个有限B块，再作双重无穷和。** 实 $0<a,b<1$ 时，每块等于非负核的积分；近a=b=0的小复邻域用可积的 $D_1^{-1-\delta}D_2^{-1-\delta}$ 控制，故分块级数绝对且局部一致收敛。不能任意打散成有符号无穷单项式并换序。全ε求值为此收敛级数，不声称有限特殊函数闭式或快速数值收敛。

## 3. 完整I

H₀、H₂采用已有的[Meijer G定义](../../../../2ne/families/single-opposite/basis.md)，$H_2=\mathcal H_{00;20}$。具体地，$H_c(e,z)$为

$$
-\frac{(1+z)z^{-e/2}}{\Gamma(1-e)\Gamma(1+e)\Gamma(1+(c+1)e)}
G_{3,3}^{3,3}\!\left(z\;\middle|\;\begin{matrix}e,0,-ce\\0,e,e\end{matrix}\right),
\quad c=0,2.
$$

记 $H_2^\pm=[H_2(e,z)\pm H_2(e,1/z)]/2$，定义已知函数的组合

$$
A(e,z)=H_0(3e,z)+\cosh(eL)[H_0(e,z)-2H_2^+]-2\sinh(eL)H_2^-.
$$

八个角项之和为 $z^{3e/2}A/[2e(1+e)]$。因此

$$
J_{07}+J_{08}=z^{3e/2}\left[\mathcal T(2e,e;z)+\mathcal T(e,2e;z)+\frac{A(e,z)}{2e(1+e)}\right].
$$

本图sym=2，唯一I为 $G^3(J_{07}+J_{08})/[64e(1-2e)(1+z^{3e})]$，即

$$
\boxed{I_{07}=I_{08}=\frac{G^3}{128e(1-2e)\cosh(3eL/2)}
\left[\mathcal T(2e,e;z)+\mathcal T(e,2e;z)+\frac{A(e,z)}{2e(1+e)}\right]}.
$$

H₀与H₂的常数相等，故A(0,z)=0，显式A/e可去。没有额外二阶极点。

## 4. Laurent系数与旧结果复用

已有[双指标母族系数](../../../../2ne/families/same-channel-two/structural-basis-functions.md)为

$$
\mathcal T(e,e;z)=\rho_0+e\rho_1+O(e^2),\qquad
\rho_0=M_{0,1}^2/6,\quad \rho_1=\tau^{[2],\mathrm{id}}_{0,0;1}.
$$

ρ₁的完整Hlog表达早已求值，机器键为`rho_id_0_0_1`；本批绑定旧文件的哈希并保存每个最终系数的完整Hlog坐标，不重新定义新的未知函数。一般T的两参数在原点求导，直接给出

$$
\mathcal T(2e,e)+\mathcal T(e,2e)=2\rho_0+3e\rho_1+O(e^2).
$$

这是不同指数双核的一阶精确连接，不把它们全阶误认成相同指数的母积分。

用原τ/η的实际Hlog坐标约化边界系数，写 $I=\sum_{k=-6}^0c_ke^k+O(e)$：

| $k$ | $c_k$ |
|---|---|
| $-6,-5,-4,-3,-2$ | $0$ |
| $-1$ | $M_{0,1}^2/384-M_{1,1}/256$ |
| $0$ | $[3\rho_1-\tfrac12M_{2,1}-\tfrac13M_{0,3}+\tfrac14L^2M_{0,1}+\tfrac23M_{0,1}^2-\tfrac12M_{1,1}]/128$ |

最高简单极点，有限项无需新基函数。A须展开到e²，有序双核展开到e¹；G³从e²才开始修正，对本图有限项无贡献。固定正实z下，双核和单对边函数在调节参数零点附近解析，所列余项为O(e)，不声称联合极端z极限的一致估计。

## 5. 验证及范围

24项精确检查覆盖原核、反射、八角边界、sym=2组装、4个终止级数与原顺序多项式积分的对照，以及7个系数位置的Hlog约化。9组数值控制包括原07/08完整核与降维式、旧ρ₀/ρ₁与独立有序积分、有限项余项缩放；具体分辨率及差值见[检查记录](../../../../../workflow/checks/self-energy-pi1-p5/README.md)。终止级数测试只针对数学核，不在Gamma极点评价物理I。

非整数全ε数值检查用独立降维求积，没有将慢收敛混合级数高精度求和；浮点求积和分辨率差不是认证误差界。正实z、$0<\Re e<1/2$起步，再向e=0附近亚纯延拓；一般复z物理延拓未验证。无可调用Wilson MCP，未修改Kernel/Lean，formal=not_run。

[07不可变结果](../../results/4e-Pi1-07/r001.json) · [08不可变结果](../../results/4e-Pi1-08/r001.json)
