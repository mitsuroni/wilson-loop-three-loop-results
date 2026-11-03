# 4e-Pi1-06：I求解工作区

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1/workspaces/4e-Pi1-06/derivation.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。

以下保留迁移前的求解/准备过程记录；其中r001的C/I公式不作为现行数值入口。现行理论及解表已经更新，J与母函数的推导仍可复用。

所属[P2：相邻SE＋相邻裸线：中边双端点](../../families/p2-adjacent-middle/README.md)。本图只使用[主推导](../../theory.md)定义的主稿式(31)的I。

## 固定输入

- 来源：[不可变JSON包](../../sources/4e-Pi1-r001.json)中的`4e-Pi1-06`。
- 源SHA-256：`3b2cab8babd9b1f4630d77d01dc64dffb9b6868c90ccab950a55c4c30230ee0a`。
- 总圈数：3；Wilson颜色迹：`CF*(CF-CA/2)`。
- 以专用prepared C进行a₁/Euler和C4系数提取；完整CA/NF颜色因子位于外部前因子，I不含颜色和尺度。
- 外到内：`s21 → s32 → s31 → s41`。

| 原变量 | 下限 | 上限 |
|---|---|---|
| s21 | 0 | 1 |
| s32 | 0 | 1 |
| s31 | 0 | s32 |
| s41 | 0 | 1 |


## 从本图输入到完整 I

准备输入：`benchmarks/self-energy/pi1/prepared/4e-Pi1.json`，SHA-256 `25cf28567a251b1c94a5f140ba6b38512e6334b98466352559deedf3bdc5fb56`；a₁=-1/2，sym=4。

本工作区对应 **4e-Pi1-06**。以下保留P2两成员的原核，逐一给出相同规范下的连接；本图结果为[不可变r001](../../results/4e-Pi1-06/r001.json)。

令 $e=\varepsilon$、$L=\log z$、$G(e)=e^{-\gamma_E e}\Gamma(1-e)$，其中 $z>0$。本页的 $I$ 严格采用[唯一主推导](../../theory.md)式(31)。四个成员均有 $n=3,\ a_1=-1/2,\ \mathrm{sym}=4$，因此

$$
I=\frac{G(e)^3}{32e(1-2e)}
\frac{J(z,e)+z^{3e}J(1/z,e)}{1+z^{3e}}.
$$

完整自能因子 $\mathcal B_1=[(10-6e)C_A-8(1-e)T_Fn_f]/(3-2e)$ 保留在共同外部前因子 $(\alpha_s/\pi)^3C_AC_F\mathcal B_1$ 中。这里只求同一标量 $I$；prepared 中的一次 Euler 补偿不再重复加入。

## 原域降维

05的完整原核为

$$
\frac{z^e}{1+e}(1-s_{21})^{2e-1}
s_{32}^{2e-1}(1-s_{31})^{e-1}s_{41}^{e-1}.
$$

06的完整原核为

$$
\frac{z^{2e}}{1+e}(1-s_{21})^{e-1}
s_{32}^{e-1}(1-s_{31})^{2e-1}s_{41}^{2e-1}.
$$

两者均有 $0<s_{31}<s_{32}<1$，$s_{21},s_{41}\in(0,1)$，原顺序为 $s_{21}\to s_{32}\to s_{31}\to s_{41}$（外到内）。先积 $s_{41}$，再积 $s_{31}$ 和 $s_{32}$，最后积 $s_{21}$，保持原迭代顺序。

共享边的有序矩为

$$
\begin{aligned}
\int_0^1db\,b^{q-1}\int_0^bda\,(1-a)^{p-1}
&=\frac1p\left[\frac1q-B(q,1+p)\right]\\
&=\frac{1}{pq}\left[1-\frac{\Gamma(1+p)\Gamma(1+q)}{\Gamma(1+p+q)}\right].
\end{aligned}
$$

这也直接说明 $p,q$ 互换不改变该矩。05取 $(p,q)=(e,2e)$，06取 $(p,q)=(2e,e)$；另两个外端点总贡献为 $1/(2e^2)$。令

$$
Q(e)=\frac{\Gamma(1+e)\Gamma(1+2e)}{\Gamma(1+3e)},
$$

则

$$
J_{05}=\frac{z^e[1-Q(e)]}{4e^4(1+e)},\qquad
J_{06}=\frac{z^{2e}[1-Q(e)]}{4e^4(1+e)}.
$$

两图各自直接求值；此结果亦核对了既有 $J_{06}(z)=z^{3e}J_{05}(1/z)$ 的连接。

## 完整结果

$$
\boxed{I_{05}=I_{06}=\frac{G(e)^3[1-Q(e)]}{128e^5(1-2e)(1+e)}
\frac{\cosh(eL/2)}{\cosh(3eL/2)}}.
$$

这里用到 $(z^e+z^{2e})/(1+z^{3e})=\cosh(eL/2)/\cosh(3eL/2)$。保留了两条线不同的尺度指数；完整 $I$ 满足 $z\leftrightarrow1/z$ 对称。

## Laurent 系数

写 $I=\sum_{k=-6}^0c_ke^k+O(e)$，并记 $\zeta_k=\zeta(k)$。05、06共用下表。

| $k$ | $c_k$ |
|---|---|
| $-6,-5,-4$ | $0$ |
| $-3$ | $\zeta_2/64$ |
| $-2$ | $(2\zeta_2-6\zeta_3)/128$ |
| $-1$ | $[2(3-L^2)\zeta_2-6\zeta_3+\tfrac{37}{2}\zeta_4]/128$ |
| $0$ | $[2(5-L^2)\zeta_2-6(3-L^2)\zeta_3+\tfrac{37}{2}\zeta_4+5\zeta_2\zeta_3-42\zeta_5]/128$ |

虽然闭式显式提出 $e^{-5}$，$1-Q(e)=2\zeta_2e^2-6\zeta_3e^3+11\zeta_4e^4+(12\zeta_2\zeta_3-42\zeta_5)e^5+O(e^6)$，所以实际最高极点为三阶。共享边矩在 $e\to0$ 时趋于 $\zeta_2$，没有增加端点极点。

有限项需要上述分子至 $e^5$；尺度比为 $1-L^2e^2+O(e^4)$。按固定正实 $z$ 展开，余项为 $O(e)$。当 $z\to0$ 或 $\infty$ 同时变化时，不声称该余项一致。

## 收敛、验证与交付边界

原参数积分在 $\Re e>0$ 绝对收敛；与前因子共同采用 $0<\Re e<1/2$。正实 $z$、$v_{14}<0$ 下所有距离为正，幂分解使用同一实对数；由此得到零点附近的亚纯表达。P1/P2均按原嵌套顺序求积，没有假定未经检查的换序桥。没有验证一般复 $z$ 的物理解析延拓。

本批共69项精确检查及46组数值控制：四个实际原核在3个正 $e$、3个正 $z$ 上独立四维 Gauss–Jacobi 求积（8/16阶）；4次实际 prepared C重构（含 $n_f=0,4$）；6组高精度余项减半检查。数值求积为浮点支持，不是认证区间证明。详见[计算与证据](../../../../../workflow/checks/self-energy-pi1-adjacent/README.md)。

当前工具集没有可调用的 Wilson MCP；未改动 Kernel/Lean，formal=not_run。计算结果不提升为形式整链完成。
