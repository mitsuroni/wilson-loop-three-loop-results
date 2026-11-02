# 4e-4g：单一规范代表的C/J输出

2026-09-18修订，schema 1.2。生成器为纯Python；每图仅保留原轮廓顺序的规范代表，不再输出A/B双份。完整推导见[4GV通用推导](4GV-general-derivation.md)，09的错误修正见[09几何](4e-4g-09-derivation.md)。

$$C_{\rm raw}=\frac{C_A^2C_F}{16}(\alpha_s/\pi)^3\Gamma(2-3\epsilon)
\pi^{3\epsilon}(-2v_{14}+i0)^{3\epsilon},$$

$$J_g=\int_{\mathcal D_g}d^4s\int_{\alpha_i\ge0}d^4\alpha\,\delta(1-\sum\alpha)
\ell_g(z)(\alpha_1\alpha_2\alpha_3\alpha_4)^{-\epsilon}S_g^{-2+3\epsilon}.$$

$S_g=\sum_{i<j}\alpha_i\alpha_j(X_i-X_j)^2/(2v_{14})$，位置、腿序和原始有序界由本图输入给出。以共同解析起始区固定分支。换成cube时须再乘对应Jacobian。

| 图 | 边 | lorentzZ | cube Jacobian | 非零贡献通道 | C4几何重数（未乘入C） |
|---|---|---|---|---|---|
| 01 | 3444 | 0 | x2²x3 | 无 | 4 |
| 02 | 3344 | −z²/4 | x1x2 | 3 | 4 |
| 03 | 3334 | 0 | x1²x3 | 无 | 4 |
| 04 | 2444 | 0 | x2²x3 | 无 | 4 |
| 05 | 2344 | z(1+z)/4 | x3 | 3 | 4 |
| 06 | 2334 | −z/4 | x2 | 1 | 4 |
| 07 | 2244 | −(1+z)²/4 | x1x2 | 3 | 2 |
| 08 | 2234 | (1+z)/4 | x1 | 3 | 4 |
| 09 | 1234 | −(1+4z+z²)/4 | 1 | 1、3 | 1 |

06的零距离配对是α2α3，09的首尾边项是α1α4 x1(1−x4)。C4几何重数只对当前完整Bose对称四胶子顶点的附着拓扑计数，不等于已完成Wick/重整化或指数化验收。

读取 `workspaces/3loop/4e-4g/montopar-initial.json` 可获得九图所有显式原s多项式、边界、腿序和cube式。`C.components` 按 `factor_order` 相乘等于combined；`J_raw_s`与`J_raw`均有四α单纯形。三零图也保留共同C和原始域，以Lorentz零原因标记J=0。

```powershell
python -m montopar.cli --source ../outputs/3loop/orbit_definitions.json --type 4e-4g --out-dir ../workspaces --package-initial --fgv-compat-out-dir montopar-outputs/3loop
```

该命令同时重建旧路径下的final/state/export三种文件，它们仅是同一新C/J的导出视图。α权重已在integrand内；δ测度不能再次当权重相乘，也不能将四α独立积成四维立方体。raw不含Euler补偿，FGV未执行IBP或Laurent求值。
