# 4e-Pi1-09：I求解工作区

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1/workspaces/4e-Pi1-09/derivation.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。

以下保留迁移前的求解/准备过程记录；其中r001的C/I公式不作为现行数值入口。现行理论及解表已经更新，J与母函数的推导仍可复用。

所属[P4：相邻SE＋对边裸线](../../families/p4-opposite-bare/README.md)。本图只使用[主推导](../../theory.md)定义的主稿式(31)的I。

## 固定输入

- 来源：[不可变JSON包](../../sources/4e-Pi1-r001.json)中的`4e-Pi1-09`。
- 源SHA-256：`3b2cab8babd9b1f4630d77d01dc64dffb9b6868c90ccab950a55c4c30230ee0a`。
- 总圈数：3；Wilson颜色迹：`CF*(CF-CA/2)`。
- 以专用prepared C进行a₁/Euler和C4系数提取；完整CA/NF颜色因子位于外部前因子，I不含颜色和尺度。
- 外到内：`s22 → s21 → s31 → s41`。

| 原变量 | 下限 | 上限 |
|---|---|---|
| s22 | 0 | 1 |
| s21 | 0 | s22 |
| s31 | 0 | 1 |
| s41 | 0 | 1 |

## 本图计算与同组连接

只采用[唯一主推导](../../theory.md)式(31)。令 $e=\varepsilon$、$L=\log z$、$G=e^{-\gamma_Ee}\Gamma(1-e)$。两成员均为 $n=3,a_1=-1/2,\mathrm{sym}=4$；完整CA/NF自能算符保持在共同前因子中，Euler因子不重复加入。

## 原域与端点降维

04的完整原核为

$$
\frac{z^{2e}}{1+e}(1-s_{31})^{2e-1}s_{42}^{2e-1}
[-(1+z)]D_z(s_{21},s_{41})^{-1+e},\qquad 0<s_{41}<s_{42}<1.
$$

09的完整原核为

$$
\frac{1}{1+e}(1-s_{21})^{2e-1}s_{31}^{2e-1}
[-(1+z)]D_z(s_{22},s_{41})^{-1+e},\qquad 0<s_{21}<s_{22}<1.
$$

其余参数为独立单位区间，$D_z(x,y)=(1-x)(1-y)+zxy$。原外到内顺序分别为 $s_{21}\to s_{31}\to s_{42}\to s_{41}$ 和 $s_{22}\to s_{21}\to s_{31}\to s_{41}$。

正实z、Re(e)>0时，完整核绝对值受相邻端点可积幂和裸对边角点核控制；放宽排序后的积域积分有限。因此可以用Fubini先积相邻SE端点，不改变原根的定义。两种有序权重分别来自

$$
\int_y^1t^{2e-1}dt=\frac{1-y^{2e}}{2e},\qquad
\int_0^x(1-t)^{2e-1}dt=\frac{1-(1-x)^{2e}}{2e}.
$$

另一个SE端点给出 $1/(2e)$。以既有[H₀/H₂](../../../../2ne/families/single-opposite/basis.md)的定义连接，得到

$$
J_{04}=\frac{z^{5e/2}}{4e^2(1+e)}[H_0(e,z)-H_2(e,z)],\qquad
J_{09}=\frac{z^{e/2}}{4e^2(1+e)}[H_0(e,z)-H_2(e,1/z)].
$$

这里 $H_2=\mathcal H_{00;20}$，其明确求值为

$$
H_2(e,z)=-\frac{(1+z)z^{-e/2}}{\Gamma(1-e)\Gamma(1+e)\Gamma(1+3e)}
G_{3,3}^{3,3}\!\left(z\;\middle|\;\begin{matrix}e,0,-2e\\0,e,e\end{matrix}\right).
$$

$H_0$将上式的 $-2e$ 改为0、$\Gamma(1+3e)$改为$\Gamma(1+e)$。它们是已有标准特殊函数表达；H₂含2e端点权重，但裸传播子的指数仍为−1+e。

## 完整I

记 $H_2^\pm=[H_2(e,z)\pm H_2(e,1/z)]/2$。本图前因子与轨道组装给出

$$
I=\frac{G^3}{32e(1-2e)}\frac{J(z)+z^{3e}J(1/z)}{1+z^{3e}}.
$$

故两成员全阶相等：

$$
\boxed{I_{04}=I_{09}=\frac{G^3}{128e^3(1-2e)(1+e)}
\frac{\cosh(eL)(H_0-H_2^+)-\sinh(eL)H_2^-}{\cosh(3eL/2)}}.
$$

完整原J关系为 $J_{09}(z)=z^{3e}J_{04}(1/z)$，I相等是在同一归一化组装之后成立。

## Laurent系数

复用[已求值系数表](../../../../2ne/families/single-opposite/basis-functions.json)，令

$$
D_k=\tau_{0,k}-\tau_{2,k}=\operatorname{coeff}_{e^k}(H_0-H_2^+),\qquad
E_k=\eta_{2,k}=L\operatorname{coeff}_{e^k}H_2^-.
$$

D、E只是已有函数的本页缩写，不增加基函数。常数项 $D_0=E_0=0$ 已在完整Hlog坐标中精确核对，因此实际最高极点从保守上界三阶降低为二阶。

写 $I=\sum_{k=-6}^0c_ke^k+O(e)$：

| $k$ | $c_k$ |
|---|---|
| $-6,-5,-4,-3$ | $0$ |
| $-2$ | $D_1/128$ |
| $-1$ | $(D_1+D_2-E_1)/128$ |
| $0$ | $[D_3+D_2-E_1-E_2+(3+\tfrac32\zeta_2-\tfrac58L^2)D_1]/128$ |

完整Hlog展开随每个不可变JSON结果保存。有限项使用H₀/H₂至e³及H₂奇部分至e²；没有漏掉cosh分母或G³的二阶贡献。对固定正实z，核的log矩在e=0附近可积，归一化分子解析，除以e³后截断余项为O(e)。不声称联合极端z极限的一致余项。

## 验证与范围

16项精确检查连接实际完整原核、端点原函数、两图尺度关系、C4组装及7个系数位置；另外断言源哈希、域、sym/a₁和Hlog常数抵消。12组数值控制包含两图降维原核与Meijer G对照、04剩余有序三角域128/256阶求积，以及独立H₀/H₂矩与有限项余项减半。浮点求积没有认证区间界。

起始域 $z>0,v_{14}<0,0<\Re e<1/2$，随后向e=0附近亚纯延拓；一般复z物理延拓未验证。当前无可调用Wilson MCP，未改Kernel/Lean，formal=not_run。

[04结果](../../results/4e-Pi1-04/r001.json) · [09结果](../../results/4e-Pi1-09/r001.json) · [检查记录](../../../../../workflow/checks/self-energy-pi1-p4/README.md)
