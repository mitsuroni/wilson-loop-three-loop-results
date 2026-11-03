# 6e-21：[3] 族 I 的解析解

日期：2026-09-17。prepared 0.3；完整ε显式收敛级数与 ε⁻⁶…ε⁰ 已计算检查，formal=not_run。

## 原始输入与排序

- [输入](../../prepared/6e.json)：6e-21，包SHA-256：9699532d1ee3a4027c2eb93c961a5b467a9bceda229c0204c13775cfa2451406。
- 配对π=[1, 3, 2]，a1=1/4，sym=2，κ=-1/32。
- 外到内次序：s43 → s23 → s42 → s22 → s41 → s21。
- 正实z；完整I先取0<ε<1，再在ε=0邻域展开。C已含一次Euler补偿。

| 变量 | 下界 | 上界 |
|---|---|---|
| s43 | 0 | 1 |
| s23 | 0 | 1 |
| s42 | 0 | s43 |
| s22 | 0 | s23 |
| s41 | 0 | s42 |
| s21 | 0 | s22 |

原被积函数：

$$
\left(- z - 1\right)^{3} \left(s_{21} s_{41} z + \left(1 - s_{21}\right) \left(1 - s_{41}\right)\right)^{\varepsilon - 1} \left(s_{22} s_{43} z + \left(1 - s_{22}\right) \left(1 - s_{43}\right)\right)^{\varepsilon - 1} \left(s_{23} s_{42} z + \left(1 - s_{23}\right) \left(1 - s_{42}\right)\right)^{\varepsilon - 1}
$$

## 结构、完整 J 与 I

两条有序链各有三个端点，每个端点属于一个传播子，无邻边核。立方体上的绝对可积界允许按公共方法分组；不因此假定固定排列能拆成单核乘积。
母核R的尺度、符号和配对定义见[公共方法](../../families/same-channel-three/method.md)。

$$
J=R_{132} z^{\frac{3 \varepsilon}{2}}
$$

$$
I=- \frac{G^{3} \left(R_{132} + R_{213}\right)}{64 \cosh{\left(\frac{3 \ell \varepsilon}{2} \right)}}
$$

crossing把排列π变成321·π·321，因此123自镜像，132与213互换；6e-21与6e-22的完整I相等。
[结果r001](../../results/6e-21/r001.json)保存每个R的显式收敛三重级数，不止于定义积分。

## Laurent 系数

| ε幂 | 结构表达 |
|---:|---|
| -6 | $0$ |
| -5 | $0$ |
| -4 | $0$ |
| -3 | $0$ |
| -2 | $0$ |
| -1 | $0$ |
| 0 | $- \frac{M_{0,1}^{3}}{1152} + \frac{T_3}{32}$ |

余项为 $O(\varepsilon)$；全部负幂为零来自ε=0邻域全纯性。

有限项的紧凑M+E表达：

$$
c_0=\frac{E}{16} - \frac{M_{0,1}^{3}}{1152} - \frac{M_{0,1} M_{2,1}}{64} + \frac{\pi^{2} M_{0,3}}{288} - \frac{M_{0,5}}{384} + \frac{M_{1,1}^{2}}{64}
$$

[基函数说明](../../families/same-channel-three/basis-functions.md)定义全部M、E和T3。T3=R123(0,z)已实际求值，且可化成已有M与[2]的E；没有新增独立有限项函数。

## 实际检查与限制

[本轮数值/精确检查](../../../../workflow/checks/same-channel-three/runs/20260917T024056933831Z/checks.json)；[复算入口](../../../../workflow/checks/same-channel-three/README.md)。
原顺序整数J与终止三重级数、非整数原六维Sobol、独立降维有限项求积、crossing与O(ε)缩放均通过。
HyperInt在新会话并开启发散检查下计算R123与R132；完整Hlog、精确表示转换和不相容见证另存。
没有运行本批Lean/MCP形式整链，没有物理复z/i0或一般更高正幂交付；完整级数尚未压缩为有限特殊函数闭式。
