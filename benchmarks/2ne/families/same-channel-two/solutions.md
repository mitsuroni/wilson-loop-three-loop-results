# [2] 族：I 的解析结果

7/7 个积分已得到 ε⁻⁶ 到 ε⁰ 的解析系数，余项 O(ε)。完整 ε 依赖使用[已逐项积完的收敛级数](method.md#4-完整-ε-的显式收敛级数)，不是有限特殊函数闭式。参数范围为 z>0；本批 formal: not_run。

所有 ε⁻⁶、ε⁻⁵、ε⁻⁴、ε⁻³、ε⁻² 系数为零。4e-04 的 ε⁻¹ 也为零。

令 A、B、C 为下列简单极点系数；ΣA/B/C 是[已求值的权重6 Hlog](basis-functions.md)，不含待求积分。

$$
A(z)=- \frac{M{\left(0,1,z \right)} \zeta\left(3\right)}{8} - \frac{M{\left(1,3,z \right)}}{96}
$$

$$
B(z)=\frac{M{\left(0,1,z \right)} M{\left(1,1,z \right)}}{64} + \frac{\pi^{2} M{\left(1,1,z \right)}}{96} + \frac{M{\left(1,3,z \right)}}{96} - \frac{N{\left(1,2,z \right)}}{64}
$$

$$
C(z)=- \frac{M{\left(0,1,z \right)} M{\left(1,1,z \right)}}{64} + \frac{M{\left(0,1,z \right)} \zeta\left(3\right)}{16} - \frac{\pi^{2} M{\left(1,1,z \right)}}{48} - \frac{M{\left(1,3,z \right)}}{64} + \frac{N{\left(1,2,z \right)}}{32}
$$

| 积分 | I 至有限项 |
|---|---|
| 4e-04 | $-M_{0,1}^2/48+O(\varepsilon)$ |
| 6e-10、6e-23 | $A(z)/\varepsilon+\Sigma_A(z)+O(\varepsilon)$ |
| 6e-11、6e-25 | $B(z)/\varepsilon+\Sigma_B(z)+O(\varepsilon)$ |
| 6e-12、6e-24 | $C(z)/\varepsilon+\Sigma_C(z)+O(\varepsilon)$ |

三对 I 的相等是完整 ε 恒等式。M(k,r,z) 表示原 M_{k,r}(z)，N(k,r,z) 表示原 N_{k,r}(z)，ζ(n) 是 Riemann ζ。N 的引入使用 ln z·M_{k,r}=M_{k,r+1}-N_{k,r}。

## 完整 ε 表达

令 G=e^{-γ_E ε}Γ(1-ε)，ℓ=ln z；F 的显式求和定义见 method.md 与 full-series.json。

$$I_{4e-04}=-\frac{G^2F_{\mathrm{id},0}(\varepsilon,z)}{8\cosh(\varepsilon\ell)}.$$

$$I_{\pi,X}=\frac{\kappa G^3}{2\varepsilon\cosh(3\varepsilon\ell/2)}\big[e^{\varepsilon\ell/2}F_{\pi,X}(\varepsilon,z)+e^{-\varepsilon\ell/2}F_{\pi,X}(\varepsilon,1/z)\big].$$

10/23：(π,X,κ)=(id,A,-1/8)；11/25：(id,B,-1/16)；12/24：(swap,A,-1/16)。

## 结果表示的边界

在已列出的固定权重、crossing 偶、有理 Hlog/周期坐标中，4e-04 有限项属于纯 M；10/23 极点属于纯 M，其余极点可用 M/N。三圈有限项都超出既有 M、M/N 和单对边扩展 τ 空间。加入 ΣA 后，ΣB、ΣC 可由旧 τ+ΣA 精确表达；这是本批权重6的一维扩充结论，不是全阶最小闭包定理。

## 逐图记录

- [4e-04](../../workspaces/4e-04/derivation.md)：[r002](../../results/4e-04/r002.json)。
- [6e-10](../../workspaces/6e-10/derivation.md)：[r002](../../results/6e-10/r002.json)。
- [6e-11](../../workspaces/6e-11/derivation.md)：[r002](../../results/6e-11/r002.json)。
- [6e-12](../../workspaces/6e-12/derivation.md)：[r002](../../results/6e-12/r002.json)。
- [6e-23](../../workspaces/6e-23/derivation.md)：[r002](../../results/6e-23/r002.json)。
- [6e-24](../../workspaces/6e-24/derivation.md)：[r002](../../results/6e-24/r002.json)。
- [6e-25](../../workspaces/6e-25/derivation.md)：[r002](../../results/6e-25/r002.json)。
