# [1,2] / [2,1] 族：I 的解析解

更新：2026-09-17。6e-37、6e-38 的完整 ε 表达与 ε⁻⁶…ε⁰ 已计算并核验。取 $z>0$，完整 I 先取 $0<\varepsilon<1$；展开余项为 $O(\varepsilon)$。formal=not_run。

## Laurent 结果

| ε 的幂次 | 6e-37 | 6e-38 |
|---:|---|---|
| −6 | 0 | 0 |
| −5 | 0 | 0 |
| −4 | 0 | 0 |
| −3 | 0 | 0 |
| −2 | 0 | 0 |
| −1 | 0 | 0 |
| 0 | $-M_{0,1}^3/96$ | $-M_{0,1}^3/96$ |

这里 $M_{0,1}(z)=\tau_{0,0}(z)=-(\ln^2z+\pi^2)/2$，所以两图均有
$$
\boxed{I(z,\varepsilon)=\frac{(\ln^2z+\pi^2)^3}{768}+O(\varepsilon).}
$$

无极点来自原积分在 ε=0 邻域的全纯性：没有邻边核，每个传播子的角点奇性在其两个独立端点上可积。不是根据拓扑标签推断，也不是把未求系数记作零。

## 完整 ε 结果

令
$$
G=e^{-\gamma_E\varepsilon}\Gamma(1-\varepsilon),\quad
\ell=\ln z,\quad
Q_{\rm id}=Q^{[2]}_{{\rm id};0,0},\quad
Q_{\rm swap}=Q^{[2]}_{{\rm swap};0,0}.
$$

H0 是已有单对边母核；Q 是已有 [2] 母核的零插入项。其定义、Meijer G 与显式收敛级数见 [统一推导 §3–5](method.md)。

$$
J_{6e-37}=z^{3\varepsilon/2}H_0Q_{\rm id},\qquad
J_{6e-38}=z^{3\varepsilon/2}H_0Q_{\rm swap},
$$
$$
\boxed{
I_{6e-37}=-\frac{G^3H_0Q_{\rm id}}{16\cosh(3\varepsilon\ell/2)},\qquad
I_{6e-38}=-\frac{G^3H_0Q_{\rm swap}}{32\cosh(3\varepsilon\ell/2)}.
}
$$

这两条的有限项相同，完整函数一般不同。乘积域的排序分割给出
$$
Q_{\rm id}+Q_{\rm swap}=\frac{H_0^2}{2},
\qquad
I_{6e-37}+2I_{6e-38}
=-\frac{G^3H_0^3}{32\cosh(3\varepsilon\ell/2)}.
$$

因此完整解也可只用一个 H0 和一个 Qid 表示；没有新增 σ、E 或其他基函数。完整 Qid 采用已积完的收敛双级数，尚未得到有限特殊函数闭式；有限项已经是显式 log/π 多项式。

## 单图、数据与检查

| 积分 | 原通道 | 双核配对 | 推导 | 不可变结果 |
|---|---|---|---|---|
| 6e-37 | 13×24² | id | [derivation](../../workspaces/6e-37/derivation.md) | [r001](../../results/6e-37/r001.json) |
| 6e-38 | 13×24² | swap | [derivation](../../workspaces/6e-38/derivation.md) | [r001](../../results/6e-38/r001.json) |

[当前清单](catalog.json)保留修订指针；每个结果保存当前 C/J 身份、完整表达、全部系数、Hlog 坐标、条件与源码/复用数据哈希。[实际检查与复算](../../../../workflow/checks/mixed-opposite/README.md)包含原顺序精确积分、原六维域数值对照与余项缩放。

本批没有修改 Kernel/Lean；形式整链、物理复 z/i0 与更高正幂系数的交付仍未实施。
