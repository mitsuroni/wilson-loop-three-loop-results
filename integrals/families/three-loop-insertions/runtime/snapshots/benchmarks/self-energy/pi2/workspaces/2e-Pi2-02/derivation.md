# 2e-Pi2-02：原对边核、边界项与2ne共同生成函数

2026-10-06审核。边2/4，原次序s21→s41，两参数范围[0,1]，C4重数2。符号及作者算符、C的来源见[本类独立推导](../../source-to-I-review.md)，原完整核为

$$J_{02}=\int_0^1ds\int_0^1dt\,
\frac{-(2+\epsilon)(1+z)F+(1-3\epsilon)z}{2(1+2\epsilon)}F^{-2+3\epsilon},
\qquad F=(1-s)(1-t)+zst.$$

令ν=3ε，完整核等于−(1+z)F^{ν−1}+∂s∂tF^ν/[2ν(1+2ε)]。四角值1、0、0、z给出完整边界质量1+z^ν。定义Hν=∫∫F^{ν−1}、H₀(ν,z)=−(1+z)z^{−ν/2}Hν，则

$$J_{02}=z^{\nu/2}H_0(\nu,z)+\frac{1+z^\nu}{2\nu(1+2\epsilon)},$$
$$I_{02}=\frac{G^3U}{128\epsilon^2(1-2\epsilon)^2}
\left[\frac{H_0(3\epsilon,z)}{\cosh(3\epsilon\log z/2)}
+\frac1{3\epsilon(1+2\epsilon)}\right].$$

初始0<Reν<1、z>0。H₀从原核精确化为原2ne的矩生成密度，并有

$$[\nu^n]H_0(\nu,z)=\sum_{j=0}^{\lfloor n/2\rfloor}
\frac{M_{n-2j,2j+1}(z)}{2^n(n-2j)!(2j+1)!}.$$

Gauss表示的相反关联支路与一重Mellin–Barnes直接推得现行Meijer G。全εI与[r002](../../results/2e-Pi2-02/r002.json)相同；原M生成的几何7位置及全颜色21位置重新匹配。9组非整数ε/z点的原二维核与闭式比较支持完整边界质量，不只核对ε=0的结果。

[核验记录](../../../../../workflow/checks/pi2-chain-audit/README.md)。内部Π₂算符仍采用作者裸全D输入，原多维内部图未在本轮再积分；数值控制非区间证明，formal=not_run。
