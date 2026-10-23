# 4e-Pi1-01：I求解工作区

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1/workspaces/4e-Pi1-01/derivation.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。

以下保留迁移前的求解/准备过程记录；其中r001的C/I公式不作为现行数值入口。现行理论及解表已经更新，J与母函数的推导仍可复用。

所属[P1：相邻SE＋相邻裸线：两边双端点](../../families/p1-adjacent-double/README.md)。本图只使用[主推导](../../theory.md)定义的主稿式(31)的I。

## 固定输入

- 来源：[不可变JSON包](../../sources/4e-Pi1-r001.json)中的`4e-Pi1-01`。
- 源SHA-256：`3b2cab8babd9b1f4630d77d01dc64dffb9b6868c90ccab950a55c4c30230ee0a`。
- 总圈数：3；Wilson颜色迹：`CF*(CF-CA/2)`。
- 以专用prepared C进行a₁/Euler和C4系数提取；完整CA/NF颜色因子位于外部前因子，I不含颜色和尺度。
- 外到内：`s32 → s31 → s42 → s41`。

| 原变量 | 下限 | 上限 |
|---|---|---|
| s32 | 0 | 1 |
| s31 | 0 | s32 |
| s42 | 0 | 1 |
| s41 | 0 | s42 |


## 从本图输入到完整 I

准备输入：`benchmarks/self-energy/pi1/prepared/4e-Pi1.json`，SHA-256 `25cf28567a251b1c94a5f140ba6b38512e6334b98466352559deedf3bdc5fb56`；a₁=-1/2，sym=4。

本工作区对应 **4e-Pi1-01**。以下保留P1两成员的原核，逐一给出相同规范下的连接；本图结果为[不可变r001](../../results/4e-Pi1-01/r001.json)。

令 $e=\varepsilon$、$L=\log z$、$G(e)=e^{-\gamma_E e}\Gamma(1-e)$，其中 $z>0$。本页的 $I$ 严格采用[唯一主推导](../../theory.md)式(31)。四个成员均有 $n=3,\ a_1=-1/2,\ \mathrm{sym}=4$，因此

$$
I=\frac{G(e)^3}{32e(1-2e)}
\frac{J(z,e)+z^{3e}J(1/z,e)}{1+z^{3e}}.
$$

完整自能因子 $\mathcal B_1=[(10-6e)C_A-8(1-e)T_Fn_f]/(3-2e)$ 保留在共同外部前因子 $(\alpha_s/\pi)^3C_AC_F\mathcal B_1$ 中。这里只求同一标量 $I$；prepared 中的一次 Euler 补偿不再重复加入。

## 完整结果

两条相邻线的原核在正域化为

$$
\frac{z^{3e}}{1+e}
(1-s_{31})^{2e-1}(1-s_{32})^{e-1}
s_{41}^{2e-1}s_{42}^{e-1}
$$

（01）；02交换两条线的 $e,2e$ 标签。两者的域均为 $0<s_{31}<s_{32}<1$、$0<s_{41}<s_{42}<1$，原顺序均为 $s_{32}\to s_{31}\to s_{42}\to s_{41}$（外到内）。

依原顺序从内到外积分。对 $\Re p,\Re q>0$，两个所需有序矩明确为

$$
\int_0^1 db\,(1-b)^{q-1}\int_0^b da\,(1-a)^{p-1}
=\frac1{q(p+q)},\qquad
\int_0^1 db\,b^{q-1}\int_0^b da\,a^{p-1}
=\frac1{p(p+q)}.
$$

01取 $(p,q)=(2e,e)$，02取 $(p,q)=(e,2e)$。两矩乘积均为 $1/(18e^4)$，从而

$$
J_{01}=J_{02}=\frac{z^{3e}}{18e^4(1+e)},\qquad
\boxed{I_{01}=I_{02}=\frac{G(e)^3}{576e^5(1-2e)(1+e)}}.
$$

这是完整 $I$ 的全阶等式；$z^{3e}$ 在 C4 系数组装后消去，因此 P1 全阶与 $z$ 无关。

## Laurent 系数

写 $I=\sum_{k=-6}^0 c_ke^k+O(e)$，并记 $\zeta_k=\zeta(k)$。01、02共用下表。

| $k$ | $c_k$ |
|---|---|
| $-6$ | $0$ |
| $-5$ | $1/576$ |
| $-4$ | $1/576$ |
| $-3$ | $(3+\tfrac32\zeta_2)/576$ |
| $-2$ | $(5+\tfrac32\zeta_2+\zeta_3)/576$ |
| $-1$ | $(11+\tfrac92\zeta_2+\zeta_3+\tfrac{57}{16}\zeta_4)/576$ |
| $0$ | $(21+\tfrac{15}{2}\zeta_2+3\zeta_3+\tfrac{57}{16}\zeta_4+\tfrac32\zeta_2\zeta_3+\tfrac35\zeta_5)/576$ |

最高五阶极点。展开采用 $\log G(e)^3=3\sum_{k\ge2}\zeta_k e^k/k$，保留解析分子至 $e^5$。其余项除以 $e^5$ 后为 $O(e)$。

## 收敛、验证与交付边界

原参数积分在 $\Re e>0$ 绝对收敛；与前因子共同采用 $0<\Re e<1/2$。正实 $z$、$v_{14}<0$ 下所有距离为正，幂分解使用同一实对数；由此得到零点附近的亚纯表达。P1/P2均按原嵌套顺序求积，没有假定未经检查的换序桥。没有验证一般复 $z$ 的物理解析延拓。

本批共69项精确检查及46组数值控制：四个实际原核在3个正 $e$、3个正 $z$ 上独立四维 Gauss–Jacobi 求积（8/16阶）；4次实际 prepared C重构（含 $n_f=0,4$）；6组高精度余项减半检查。数值求积为浮点支持，不是认证区间证明。详见[计算与证据](../../../../../workflow/checks/self-energy-pi1-adjacent/README.md)。

当前工具集没有可调用的 Wilson MCP；未改动 Kernel/Lean，formal=not_run。计算结果不提升为形式整链完成。
