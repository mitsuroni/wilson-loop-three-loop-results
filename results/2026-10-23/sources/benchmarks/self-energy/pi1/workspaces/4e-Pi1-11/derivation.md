# 4e-Pi1-11：I求解工作区

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1/workspaces/4e-Pi1-11/derivation.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。

以下保留迁移前的求解/准备过程记录；其中r001的C/I公式不作为现行数值入口。现行理论及解表已经更新，J与母函数的推导仍可复用。

所属[P6：两条对边线：不同通道、域可分离](../../families/p6-opposite-distinct/README.md)。本图只使用[主推导](../../theory.md)定义的主稿式(31)的I。

## 固定输入

- 来源：[不可变JSON包](../../sources/4e-Pi1-r001.json)中的`4e-Pi1-11`。
- 源SHA-256：`3b2cab8babd9b1f4630d77d01dc64dffb9b6868c90ccab950a55c4c30230ee0a`。
- 总圈数：3；Wilson颜色迹：`CF*(CF-CA/2)`。
- 以专用prepared C进行a₁/Euler和C4系数提取；完整CA/NF颜色因子位于外部前因子，I不含颜色和尺度。
- 外到内：`s11 → s21 → s31 → s41`。

| 原变量 | 下限 | 上限 |
|---|---|---|
| s11 | 0 | 1 |
| s21 | 0 | 1 |
| s31 | 0 | 1 |
| s41 | 0 | 1 |

## 计算推导

沿用[唯一主推导](../../theory.md)式(31)，$e=\varepsilon$、$L=\log z$、$G(e)=e^{-\gamma_Ee}\Gamma(1-e)$。本图 $n=3$、$a_1=-1/2$、$\mathrm{sym}=2$；完整 $C_A/N_F$ 算符位于外部共同因子中，Euler补偿只使用prepared中的一次。

## 1. 原输入到乘积

[prepared](../../prepared/4e-Pi1.json)中的原域为四个独立单位区间，原外到内顺序为 $s_{11}\to s_{21}\to s_{31}\to s_{41}$。令

$$
F=z(1-s_{11})(1-s_{31})+s_{11}s_{31},\qquad
D=(1-s_{21})(1-s_{41})+zs_{21}s_{41}.
$$

完整原核恰为

$$
\frac{-2(1+z)F+(1-2e)z}{2(1+e)}F^{-2+2e}
\,[-(1+z)]D^{-1+e}.
$$

两因子分别仅含13、24端点。对 $z>0,\Re e>0$，SE核角点绝对值由距离的 $-2+2\Re e$ 次幂控制，二维径向测度给出可积的 $r^{-1+2\Re e}dr$；裸核为可积的 $r^{\Re e}dr$。全域绝对可积，允许在保留原根顺序的前提下用Fubini进行因子化。

$F=zD_{1/z}$，分子也等于标准BO分子的 $z^2$ 倍。因此本图的SE积分为 $z^{2e}J_{B_O}(1/z,e)$。利用已核验的BO自交叉关系及其[完整解](../../families/b-opposite/solutions.md)，

$$
\begin{aligned}
J_{11}(z,e)
&=z^{2e}J_{B_O}(1/z,e)\,z^{e/2}H_0(e,z)\\
&=z^{3e/2}H_0(e,z)
\left[H_0(2e,z)+\frac{\cosh(eL)}{2e(1+e)}\right].
\end{aligned}
$$

这里使用的是已求值的[标准Meijer G函数](../../../../2ne/families/single-opposite/basis.md)：

$$
H_0(\delta,z)=
-\frac{(1+z)z^{-\delta/2}}{\Gamma(1-\delta)\Gamma(1+\delta)^2}
G_{3,3}^{3,3}\!\left(z\;\middle|\;
\begin{matrix}\delta,0,0\\0,\delta,\delta\end{matrix}\right).
$$

$H_0(\delta,1/z)=H_0(\delta,z)$，所以 $J_{11}(z)=z^{3e}J_{11}(1/z)$。上述式子是J层的因子化，不能改成两个已经归一化I的直接乘积。

## 2. 本图前因子与完整 I

按本图的sym和a₁，

$$
I_{11}=\frac{G(e)^3}{64e(1-2e)}
\frac{J_{11}(z,e)+z^{3e}J_{11}(1/z,e)}{1+z^{3e}}.
$$

因此

$$
\boxed{I_{11}=
\frac{G(e)^3 H_0(e,z)}{64e(1-2e)\cosh(3eL/2)}
\left[H_0(2e,z)+\frac{\cosh(eL)}{2e(1+e)}\right]}.
$$

这给出同一个完整I，保留完整自能核。本图 $\mathrm{sym}=2$，不同于P1/P2的4，不能照搬它们的轨道系数。

## 3. 已有函数及 Laurent 系数

使用原有[已求值Hlog系数表](../../../../2ne/families/single-opposite/basis-functions.json)，记

$$
H_0(\delta,z)=T_0+\delta T_1+\delta^2T_2+O(\delta^3),\qquad
T_0=\tau_{0,0}=M_{0,1},\quad T_1=\tau_{0,1}=M_{1,1}/2,\quad T_2=\tau_{0,2}.
$$

这些均是已有函数，不是未求值积分的新命名。不可变结果文件同时保存下列每个系数的完整Hlog坐标。

写 $I_{11}=\sum_{k=-6}^0c_ke^k+O(e)$：

| $k$ | $c_k$ |
|---|---|
| $-6,-5,-4,-3$ | $0$ |
| $-2$ | $T_0/128$ |
| $-1$ | $(T_0+T_1+2T_0^2)/128$ |
| $0$ | $[T_2+T_1+(3+\tfrac32\zeta_2-\tfrac58L^2)T_0+4T_0^2+6T_0T_1]/128$ |

最高为二阶极点；$T_0=-(L^2+\pi^2)/2$，所以首项也可写为 $-(L^2+\pi^2)/(256e^2)$。有限项只需要既有 $\tau_{0,2}$、M01/M11及log/ζ。

展开预算：含显式 $e^{-2}$ 的项需将H₀保留至二阶；含 $e^{-1}$ 的乘积保留至一阶。补偿Gamma因子为 $G^3=1+\tfrac32\zeta_2e^2+O(e^3)$。对固定正实z，各解析因子在e=0附近解析，故所列余项为O(e)；不声称极端z与e联合极限的一致性。

## 4. 验证与边界

15项精确检查连接实际源分子、距离、多项式因子化、C4归一化和7个系数位置。数值检查在精确因子化之后，分别求积两个二维原核，使用128/256阶；并以裸核的一维除差表示独立计算H₀的前三个矩，检查完整I的O(e)与H₀的O(e³)余项缩放。浮点求积不构成认证区间证明。

共同起始域为 $z>0,v_{14}<0,0<\Re e<1/2$，结果亚纯延拓到e=0附近；一般复z物理延拓未验证。当前无可调用Wilson MCP，未改变Kernel或Lean，formal=not_run。

[不可变结果r001](../../results/4e-Pi1-11/r001.json) · [计算证据](../../../../../workflow/checks/self-energy-pi1-p6/README.md)
