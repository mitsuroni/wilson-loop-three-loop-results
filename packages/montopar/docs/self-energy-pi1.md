# 单个 Π₁ 的完整 montopar 参数化

日期：2026-09-22。适用 2e-Pi1 两图（总二圈）及 4e-Pi1 十一图（总三圈）。普通 Feynman 规范，D=4−2ε，裸无质量一圈 1PI 自能，TF=1/2，nf 为无质量夸克味数。主输出是原 s 有序域的完整 C×J，保留度规项和导数项。raw 不投影颜色、不补 Euler、不含 C4 重数。

## 1. 起始动量式与唯一符号账本

保持源 qgraf-pipeline/integrals.md 最外层的负号、Wilson 顶点 igv、传播子 −ig/(ℓ²+i0)。Π₁ 使用 [Davydychev–Osland–Tarasov 1996](https://arxiv.org/pdf/hep-ph/9605348) 附录 C (C.1)、(C.3)，其规范参数取零：

$$
\Pi^{ab}_{\mu\nu}(\ell)=i\delta^{ab}(\ell^2g_{\mu\nu}-\ell_\mu\ell_\nu)
\frac{g^2}{2(4\pi)^{2-\varepsilon}}
\frac{\Gamma(\varepsilon)\Gamma^2(1-\varepsilon)}{\Gamma(2-2\varepsilon)}
C_{\rm SE}(\varepsilon)(-\ell^2-i0)^{-\varepsilon},
$$

$$
C_{\rm SE}(\varepsilon)=\frac{(10-6\varepsilon)C_A-8(1-\varepsilon)T_Fn_f}{3-2\varepsilon}.
$$

JSON 中的 Φ 与位置差 X 的关系是 Φ=−X，指数 exp(+iℓ·Φ)=exp(−iℓ·X)。允许整条线反向路由，两个 X 收缩同步反号；核不变。SE 两侧内部颜色指标用 δ 连接，随后做 Tr/N 的完整 Wilson 颜色迹。

令 q=−ℓ²−i0、Q=−X²+i0。Fourier 母式是

$$
\mathcal F_a(X)=\int\frac{d^D\ell}{(2\pi)^D}e^{-i\ell\cdot X}q^{-a}
=i\frac{\Gamma(D/2-a)}{4^a\pi^{D/2}\Gamma(a)}Q^{a-D/2}.
$$

两侧传播子与 Π₁ 相乘得到 iA[g q^(−1−ε)+ℓℓ q^(−2−ε)]，A 是 Π₁ 去掉 i q^(−ε) 的标量系数。ℓμℓν 转成 −∂μ∂ν，故 Fourier 后为 −A 乘以下实系数核：

$$
H_{\mu\nu}(X)=h_1\left[g_{\mu\nu}Q^{-1+2\varepsilon}
+\frac{\partial_\mu\partial_\nu Q^{2\varepsilon}}{8\varepsilon(1+\varepsilon)}\right],
\qquad h_1=\frac{\Gamma(1-2\varepsilon)}{4^{1+\varepsilon}\pi^{2-\varepsilon}\Gamma(1+\varepsilon)}.
$$

这是因为 h₂/h₁=−1/[8ε(1+ε)]。实际求导后：

$$
H_{\mu\nu}/h_1=
\frac{1+2\varepsilon}{2(1+\varepsilon)}g_{\mu\nu}Q^{-1+2\varepsilon}
-\frac{1-2\varepsilon}{1+\varepsilon}X_\mu X_\nu Q^{-2+2\varepsilon}.
$$

这一核在 D=4−2ε 下横向。光锥奇点上的含义继承同一个 Feynman 边界值与解析延拓，不能把离奇点求导式当作任意普通函数等式。

## 2. 无量纲几何与完整 J

记 v14=v₁·v₄，z=v12/v14，A_ref=−2v14。取起始区域 v14<0、z>0，再共同解析延拓。类光闭合四边形满足 v12=v34=z v14、v23=v14、v13=v24=−(1+z)v14。对每条交换线定义

$$
f=\frac{X^2}{2v_{14}},\qquad d=\frac{v_i\cdot v_j}{v_{14}},
\qquad u=\frac{v_i\cdot X}{v_{14}},\qquad w=\frac{v_j\cdot X}{v_{14}},
\qquad X=x_i(s_i)-x_j(s_j).
$$

完整自能线的无量纲核为

$$
\mathcal S(f,d,u,w)=
\frac{(1+2\varepsilon)df+(1-2\varepsilon)uw}{2(1+\varepsilon)}(f+i0)^{-2+2\varepsilon}.
$$

逐图保存的原度规项是 d(f+i0)^(−1+2ε)，导数项是

$$
\mathcal S_{\partial}=\frac{-df+(1-2\varepsilon)uw}{2(1+\varepsilon)}(f+i0)^{-2+2\varepsilon}.
$$

无额外裸线时 J=∫Ω S；有额外裸线 b 时 J=∫Ω S d_b(f_b+i0)^(−1+ε)。Ω 保存所有 Wilson 端点的原有序约束，不把三圈四端点拆成两个独立二重积分。实际代码把上述导数全部执行，并导出显式多项式分子；不留下待作用的微分算子。

## 3. C 的完整公式

总圈数 L=2 或 3。令 T_W 是通过 SE 内部 δ 连接后的归一化 Wilson 颜色迹。则

$$
\boxed{C_L(T_W)=
\frac{(-1)^{L-2}T_W C_{\rm SE}(\varepsilon)}{2^{L+2}\varepsilon(1-2\varepsilon)}
\left(\frac{\alpha_s}{\pi}\right)^L
\left[\pi^{\varepsilon}\Gamma(1-\varepsilon)\right]^L
(-2v_{14}+i0)^{L\varepsilon}.}
$$

两端点 T_W=CF。四端点交错配对 T_W=CF(CF−CA/2)，非交错配对为 CF²；当前十一张来源均为交错配对。保留 C_SE 的两个颜色分量。上述 g²=4παs 为裸耦合记号转换；没有偷偷乘入 μ、MSbar、反项或 Euler 补偿。

总二圈在抽取尺度前的相位账本为：源负号 × I²(Wilson) × (−I)²(传播子) × I(Π) × (−1)(ℓ²=−q) × I(Fourier)=−1。再乘 v14/A_ref=−1/2，C₂ 符号为正。总三圈另乘一对 Wilson 顶点和一条裸传播子，其尺度抽取后的因子为 −(αs/π)π^ε Γ(1−ε)A_ref^ε/2；故 C₃ 的共同符号为负。最终颜色因子本身可以为负。

Gamma 简化只用递推：Γ(ε)/Γ(1+ε)=1/ε，Γ(2−2ε)=(1−2ε)Γ(1−2ε)。分解账本各分量的乘积严格等于 combined；没有额外 1/2!、n_diagrams 或旋转重数。

## 4. 两个控制例与适用范围

相邻边 3/4：f=z(1−s31)s41、d=z、uw=z f。因此 S=z(f+i0)^(−1+2ε)/(1+ε)，完整核/度规核=1/(1+ε)，不是旧错误因子 (1+2ε)/(1+ε)。

对边 2/4：f=(1−s21)(1−s41)+z s21 s41、d=−(1+z)，u=−1+(1+z)s41，w=1−(1+z)s21。即使在 z=−1 时 d=0，导数分子也不为零；不能按裸线规则删除 SE。

因为 f_s=u、f_t=−w、f_st=−d，导数项的积分是

$$
\int_0^1ds\int_0^1dt\,\mathcal S_{\partial}
=\frac{[f(1,1)]^{2\varepsilon}-[f(1,0)]^{2\varepsilon}-[f(0,1)]^{2\varepsilon}+[f(0,0)]^{2\varepsilon}}
{4\varepsilon(1+\varepsilon)}.
$$

在 Re ε>0 的端点解释下，六个无序边对的矩形边界贡献相消；CA 与 TF nf 通道分别成立。这只作为二点扇区控制，不删除任何单图的导数项，不假设四端点有序三圈的导数和消失。

## 5. 输出与运行

schema 1.4、dispatch=self_energy。C 是上述完整共同前因子；J_raw 和 J_raw_s 均是同一个原 s 表示，没有 rotated/cube/IBP 副本。operator、representative、color_contraction、c4_symmetry、branch_prescription 保存可追溯连接。translated 表示已参数化，evaluation_status=not_evaluated、formal=not_run。

```powershell
D:/py/python.exe -m montopar.cli --source ../qgraf-pipeline/outputs/2loop/orbit_definitions.json --type 2e-Pi1 --package-initial
D:/py/python.exe -m montopar.cli --source ../qgraf-pipeline/outputs/3loop/orbit_definitions.json --type 4e-Pi1 --package-initial
```

Π₂、多次插入、带外接 TGV、同边 SE 及其他未实现类型明确拒绝。CLI 先按用户选定类型筛选再解析，因此其他未支持类型不妨碍单独生成 Π₁。未知内部顶点不能再回退 pure。打包时从来源重建并核对完整 C/J，拒绝损坏、降级或来源不匹配的结果。

本轮未计算十三个 Wilson 积分的解析解或 Laurent 系数，未建立一般复运动学的解析延拓证明，也未扩大 Wilson Kernel 的默认 47 条根输入。
