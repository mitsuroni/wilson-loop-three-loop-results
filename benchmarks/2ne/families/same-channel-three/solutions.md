# [3] 族：完整 I 与 Laurent 解析解

更新：2026-09-17。6e-20、6e-21、6e-22 均已完成正实 z 下的全 ε 显式收敛级数与 ε⁻⁶…ε⁰。formal=not_run；余项 O(ε)。

## 有限项

使用既有 M 与 [2] 的 $E=\eta^{[2],\mathrm{id}}_{0,1;1}$，定义已经求值的结构组合
$$
T_3(z)=R_{123}(0,z)
=2E-\frac12M_{0,1}M_{2,1}
+\frac12M_{1,1}^2+\frac{\pi^2}{9}M_{0,3}-\frac1{12}M_{0,5}.
$$

| ε 的幂次 | 6e-20 | 6e-21 | 6e-22 |
|---:|---|---|---|
| −6 | 0 | 0 | 0 |
| −5 | 0 | 0 | 0 |
| −4 | 0 | 0 | 0 |
| −3 | 0 | 0 | 0 |
| −2 | 0 | 0 | 0 |
| −1 | 0 | 0 | 0 |
| 0 | $-T_3/16$ | $T_3/32-M_{0,1}^3/1152$ | $T_3/32-M_{0,1}^3/1152$ |

T3 不是未求值函数，也不是新增的独立基底元素：它的完整 Hlog 已由本轮 HyperInt 计算，并精确化成上面的 M+E。所有符号见 [基函数说明](basis-functions.md)。

具体地，
$$
I_{6e-20}=
-\frac E8+\frac{M_{0,1}M_{2,1}}{32}
-\frac{M_{1,1}^2}{32}
-\frac{\pi^2M_{0,3}}{144}
+\frac{M_{0,5}}{192}+O(\varepsilon),
$$
$$
I_{6e-21}=I_{6e-22}=
\frac E{16}-\frac{M_{0,1}M_{2,1}}{64}
+\frac{M_{1,1}^2}{64}
+\frac{\pi^2M_{0,3}}{288}
-\frac{M_{0,5}}{384}
-\frac{M_{0,1}^3}{1152}+O(\varepsilon).
$$

三条均无极点：原六维域的每个角点核在 ε=0 邻域可积。完整归一化不会在零点引入极点，见 [方法 §2](method.md)。

## 完整 ε 表达

令 $G=e^{-\gamma_E\varepsilon}\Gamma(1-\varepsilon)$、$\ell=\ln z$。三标记有序母核为
$$
R_\pi=-(1+z)^3z^{-3\varepsilon/2}
\int_{\Omega_3}\prod_{i=1}^3D_z(t_i,s_{\pi(i)})^{\varepsilon-1}\,d^3t\,d^3s,
$$
其中 $\Omega_3$ 要求两条边的三个端点各自递增。完整求值采用 [显式收敛三重级数](method.md)，所有有限块均已积完。

$$
J_{20}=z^{3\varepsilon/2}R_{123},\quad
J_{21}=z^{3\varepsilon/2}R_{132},\quad
J_{22}=z^{3\varepsilon/2}R_{213},
$$
$$
\boxed{
I_{20}=-\frac{G^3R_{123}}{16\cosh(3\varepsilon\ell/2)},\qquad
I_{21}=I_{22}
=-\frac{G^3(R_{132}+R_{213})}{64\cosh(3\varepsilon\ell/2)}.
}
$$

这里20/21/22简写完整编号6e-20/21/22。21与22的完整I相等来自 crossing 交换两个配对，未把它们的原J当作同一个函数。

零阶的额外关系
$$R_{123}(0,z)+R_{132}^+(0,z)=M_{0,1}^3/36$$
使有限项只需一个组合T3；没有据此声称全 ε 母核也满足这一式。

## 表示结论与数据

所测固定权重6空间中，纯M、M加显式ln z、旧单对边扩展τ都不能覆盖；加入此前[2]已有的E后即可。无需再增加一个独立权重6方向。此结论限定于所声明的有理多项式空间；更高ε阶和任意函数闭包仍开放。

| 积分 | 配对 | 推导 | 当前首个结果 |
|---|---|---|---|
| 6e-20 | 123 | [derivation](../../workspaces/6e-20/derivation.md) | [r001](../../results/6e-20/r001.json) |
| 6e-21 | 132 | [derivation](../../workspaces/6e-21/derivation.md) | [r001](../../results/6e-21/r001.json) |
| 6e-22 | 213 | [derivation](../../workspaces/6e-22/derivation.md) | [r001](../../results/6e-22/r001.json) |

[当前清单](catalog.json)、[基函数完整坐标](basis-functions.json)、[实际检查](../../../../workflow/checks/same-channel-three/README.md)。

本批未修改Kernel/Lean，未进行形式整链或物理复z/i0延拓。全ε级数尚未压缩为有限特殊函数闭式；旧[2]基底表及迁移状态保持历史记录。
