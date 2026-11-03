# −1/(1+ε) 的两个来源及旧 SE 推导的具体差异

2026-09-23。对照用户指定的[旧稿 I3_SE_因子分解推导.md](/C:/Users/Lenovo/Documents/learning/test/TGV/I3_SE_因子分解推导.md)，以及当前[理论r002](theory-r002.md)。旧文件、生产代码和C/J均未修改。

## 1. 先区分比较对象

记P为用户此前给出的论文度规/约化表达式，I_g为当前给定动量起点的度规项，I_full为同一个起点保留完整张量的结果。

$$
\boxed{\frac{I_{\rm full}}P
=\underbrace{\frac{I_g}P}_{-1}
\underbrace{\frac{I_{\rm full}}{I_g}}_{1/(1+\varepsilon)}
=-\frac1{1+\varepsilon}.}
$$

其中只有1/(1+ε)来自补齐−ℓμℓν导数项；常数−1是当前整图定义与所给论文表达式之间的符号关系。完整自能张量ℓ²g−ℓℓ本身横向，旧稿称作“横向项”的I_k专指其中−ℓℓ的贡献。

旧稿[第1节](/C:/Users/Lenovo/Documents/learning/test/TGV/I3_SE_因子分解推导.md:16)已经有最外层负号，[Π标量定义](/C:/Users/Lenovo/Documents/learning/test/TGV/I3_SE_因子分解推导.md:35)也已经含+i。因此新推导没有比旧动量起点多加一个负号。当前与论文的比较不能被误读为“新旧两版之差完全是横向修正”。

## 2. 第一处关键差异：Fourier 核的坐标符号

旧稿[第4b节](/C:/Users/Lenovo/Documents/learning/test/TGV/I3_SE_因子分解推导.md:139)写为

$$
\int\frac{d^D\ell}{(2\pi)^D}(-\ell^2-i0)^\alpha e^{-i\ell\cdot X}
\stackrel{\rm old}{=}
\frac{i}{(4\pi)^{D/2}}\frac{\Gamma(D/2+\alpha)}{\Gamma(-\alpha)}
\left(\frac{X^2}{4}-i0\right)^{-D/2-\alpha}.
$$

在同一个动量母式和Fourier相位下，应为

$$
\boxed{
\int\frac{d^D\ell}{(2\pi)^D}(-\ell^2-i0)^\alpha e^{-i\ell\cdot X}
=\frac{i}{(4\pi)^{D/2}}\frac{\Gamma(D/2+\alpha)}{\Gamma(-\alpha)}
\left(-\frac{X^2}{4}+i0\right)^{-D/2-\alpha}.}
$$

这不是给整张图乘同一个相位就能修复的问题：两个Fourier指数分别产生距离幂−1+2ε和2ε，两幂相差1。反转坐标底数会使两个核之间多出相对负号。旧稿在后续第4c、5c节沿用了同一个错误母式。

可以不依赖相邻积分就检验这一点。定义

$$
F_a=\int\frac{d^D\ell}{(2\pi)^D}e^{-i\ell\cdot X}(-\ell^2-i0)^{-a}.
$$

动量积分直接要求

$$
\boxed{\Box_XF_{2+\varepsilon}=F_{1+\varepsilon}.}
$$

将旧稿第4c和5c两式代入却得到负的F_(1+ε)。当前母式满足正号关系，故这是旧核与其自身动量定义之间的矛盾，不只是不同作者的命名差异。

## 3. 旧的正修正怎样变成当前的负修正

旧稿第5b节的ℓμℓν→−∂μ∂ν是正确的，第6c节的Gamma递推本身也正确。问题不是这两步，而是它们作用于第2节指出的错误Fourier核。

设两个Fourier系数的公共 `(4π)^(-D/2)` 已约去，则

$$
R=\frac{\Gamma(-2\varepsilon)/\Gamma(2+\varepsilon)}
{\Gamma(1-2\varepsilon)/\Gamma(1+\varepsilon)}
=-\frac1{2\varepsilon(1+\varepsilon)}.
$$

相邻几何满足(v₃·X)(v₄·X)=X²v34/2。对旧Y=X²/4−i0求导后，系数是+2ε²；对正确的Y=−X²/4+i0求导后，相应系数是−2ε²。因此

$$
\left(\frac{I_k}{I_g}\right)_{\rm old}=-R(2\varepsilon^2)
=+\frac{\varepsilon}{1+\varepsilon},
\qquad
\left(\frac{I_k}{I_g}\right)_{\rm current}=-R(-2\varepsilon^2)
=-\frac{\varepsilon}{1+\varepsilon}.
$$

这正对应旧稿[第6c–6d节](/C:/Users/Lenovo/Documents/learning/test/TGV/I3_SE_因子分解推导.md:351)和当前结果的差异：

| 同一动量起点内的比值 | 旧稿 | 当前完整核 |
|---|---|---|
| 导数项/度规项 | +ε/(1+ε) | −ε/(1+ε) |
| 完整项/度规项 | (1+2ε)/(1+ε) | 1/(1+ε) |

这些都是同一个整图符号下的相对关系。它们不包含第1节与论文P比较时的额外常数−1。

旧稿第6d节把ε→0时的相对比值趋零写成“D=4时I_k消失”，对发散单图并不成立：I_g具有ε⁻³极点，因而I_k仍有ε⁻²极点。删去导数项必须另有完整图组的相消依据。

## 4. 第二类差异：旧前因子还有独立的归一化问题

### 4.1 第7a节的2的幂次计算不成立

旧稿[396–407行](/C:/Users/Lenovo/Documents/learning/test/TGV/I3_SE_因子分解推导.md:396)合并FT₃与公共Fourier因子时，将下列分母化为128：

$$
2\times2^{4-2\varepsilon}\times4^{2-\varepsilon}
=2^{9-4\varepsilon}\ne128.
$$

仅看这个幂次因子，在ε=0已经是512而不是128。正确的该中间步骤为

$$
{\rm FT}_3\frac{i}{(4\pi)^{2-\varepsilon}}
=-\frac{\Gamma(\varepsilon)\Gamma^2(1-\varepsilon)}
{2^{9-4\varepsilon}\pi^{4-2\varepsilon}\Gamma(2-2\varepsilon)}.
$$

随后还要乘坐标核的1/4尺度贡献，不能把这个中间分母直接当成最终C。

### 4.2 第7c–7d节重复保留了v14

取r=v14、A=−2r>0。正确坐标底数为Y=Af/4。Lorentz中的r与距离的r幂合并为

$$
r\left(\frac A4\right)^{-1+2\varepsilon}
=-2^{1-4\varepsilon}A^{2\varepsilon}.
$$

与上面的负Fourier系数相乘后，数值因子是

$$
\left(-2^{-(9-4\varepsilon)}\right)
\left(-2^{1-4\varepsilon}\right)=\frac1{256}.
$$

所以在未把g⁴改成αs的分解中，共同Gamma前因子的数值分母为256；尺度只剩A^(2ε)。旧稿[第7d节](/C:/Users/Lenovo/Documents/learning/test/TGV/I3_SE_因子分解推导.md:457)保留128并在A^(2ε)之外再乘一次r，是另一个独立问题。它使恢复尺度后的齐次次数从2ε变为1+2ε。开头设r=1会掩盖这项错误，但不能证明额外r应该存在。

旧稿第7d节注释中“负号与g⁴的正号相消”也不是有效的符号运算。当前正的C₂来自原动量相位链与r/A=−1/2的组合，详见r002式(18)、(24)、(28)，不能靠g⁴为正来消除负号。

## 5. 用同一个J比较最终两版，避免误判因子放置

旧稿把张量修正放进C，当前主输出把它保留在完整J中。因此单独比较两个C不够，应比较C×J，或先固定相同的J_g。

令公共因子为

$$
\mathcal N=C_FC_{\rm SE}(4\pi\alpha_s)^2
\frac{\Gamma(1-2\varepsilon)\Gamma(\varepsilon)\Gamma^2(1-\varepsilon)}
{\pi^{4-2\varepsilon}\Gamma(1+\varepsilon)\Gamma(2-2\varepsilon)}
(-2r+i0)^{2\varepsilon}.
$$

在给旧末式补上与当前一致的尺度处方后，按[第10节最终C/J](/C:/Users/Lenovo/Documents/learning/test/TGV/I3_SE_因子分解推导.md:531)字面比较：

$$
I_{\rm old}=\mathcal N\frac{r}{128}\frac{1+2\varepsilon}{1+\varepsilon}J_g,
\qquad
I_{\rm current}=\mathcal N\frac1{256}\frac1{1+\varepsilon}J_g.
$$

故

$$
\boxed{\frac{I_{\rm old}}{I_{\rm current}}=2r(1+2\varepsilon),
\qquad r=1:\quad\frac{I_{\rm old}}{I_{\rm current}}=2(1+2\varepsilon).}
$$

即使只将旧`transverse`改为1/(1+ε)，按其最终公式仍差2r，在r=1时仍差2。这个比值是统一分支后的末式比较；它不为旧稿中途从X²−i0切换到−X²+i0的操作提供合法性。

因此应分别记住三个比值：当前full/current g=1/(1+ε)，当前full/给定论文P=−1/(1+ε)，旧末式/current full=2(1+2ε)（旧r=1）。它们的分母是不同的对象。

本次[22项精确检查](../../../workflow/checks/self-energy-pi1-legacy/README.md)复现上述诊断。没有修改旧推导、当前生产或r001/r002；这些检查不是形式证明，也不把此前论文绝对符号验收改记为通过。
