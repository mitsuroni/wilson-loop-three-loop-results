# 2e-Pi2-01：从实际动量源到全ε的I

2026-10-06审核。相邻边3/4，原次序s31→s41，两参数范围[0,1]，C4重数4。有效Π₂含两圈、完整三个颜色及g₀⁴；与两个Wilson端点组成三圈g₀⁶源。原raw负号经现行2.0导入规范化一次，当前I不再额外改号。

完整链见[本类独立推导](../../source-to-I-review.md)：从qgraf-pipeline实际输出重放Montopar，用作者全D数据和独立Fourier得到C；再由Φ=−X重建原J。令G=e^{−γ_Eε}Γ(1−ε)，U为现行Gamma比值，得到

$$J_{01}=\frac{2+\epsilon}{2(1+2\epsilon)}z^{3\epsilon}
\int_0^1ds\int_0^1dt\,[(1-s)t]^{-1+3\epsilon}
=\frac{(2+\epsilon)z^{3\epsilon}}{18\epsilon^2(1+2\epsilon)}.$$

初始Reε>0，两个端点各产生1/(3ε)，不是根据Laurent表反推。恢复C4与两个crossing半项后

$$I_{01}=\frac{G^3U(2+\epsilon)}{1152\epsilon^4(1-2\epsilon)^2(1+2\epsilon)}.$$

三个物理颜色系数为4ε²p_c I，p_c来自作者A_c+B_c/U；必须保留Ff。新展开与现行[r002](../../results/2e-Pi2-01/r002.json)的全部7个几何和21个颜色位置相同。原结果保持，formal=not_run。

[核验记录](../../../../../workflow/checks/pi2-chain-audit/README.md)。内部Π₂的原始两圈图未在本轮重新积分，范围沿用主审核报告。
