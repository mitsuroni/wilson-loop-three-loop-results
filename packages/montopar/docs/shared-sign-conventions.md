# 共同作用量、路径序和 Fourier 符号（2026-10-06）

当前原生输出采用 `wilson-qcd-action-path-fourier/1.0`。同一个作用量和相位规则先确定完整动量图，再提取 C/J。运行时不读取图号变号表，也不使用 RG 校准图的符号。

取 Minkowski (+---)、路径积分 exp(iS)、Dμ=∂μ−igAμ，F=dA−dA+gfAA，场 Fourier 为 exp(−ipX)。Wilson 线满足 U′=igAU，轮廓上较晚的生成元放在迹的左边；腿的 μ/颜色/动量标签保持原编号。坐标标量传播子负号、动量圈 Wick +i、坐标中心 Wick −i 必须分别保留。

完整含 i 的 Π/Γ 已经包括内部费米／鬼闭圈负号。外部 Dyson 是 D₀+D₀ΠD₀+⋯，系数 +1。Γ₁沿用经过全 D 母积分核对的 V=iΓ_DOT 字典，不再补内部闭圈号。

由此直接得到：

| 类型 | 原生相位与旧输出的关系 |
|---|---|
| pure、树 TGV、mixed、双 TGV | 完整 C×J 保持；树/mixed 旧同序迹与旧星相位的成对号差同时统一 |
| 单 Π₁、Π₁串联、Π₂ | 原生正 Dyson 链直接重现已规范化 r002，不再由导入器补号 |
| TGV 腿 Π₁、Γ₁ | 路径反序在提取 C 前产生负迹；C 为旧固定早到晚源的相反数，J 保持 |
| FGV | i⁴·(−i)·(−1)⁴·(−i)=−1；原生 C=−CA²CF Γ(2−3ε)(αs/π)³(πA₀)^{3ε}/16，六个非零图改变、三个零图保持 |

`physical_convention` 声明完整路径词及 Fourier/Wick/Dyson 字典。CLI 要求此声明，拒绝未声明的旧振幅；打包器拒绝混合旧输出。顶点数组的存储顺序不决定路径序，端点 s 标签决定轮廓顺序。低层拓扑 API 重建物理图时明确标记其解释。

原生 C 的 `components`、`factor_order` 可乘，`assembly_sign=1`。raw 不含 Euler、C4 或镜像重数；IBP 的 raw/边界/recoil 共用同一个 C。准备时仅补一次 exp(−n EulerGamma ε)，完整颜色/web 投影和最终唯一 I 提取仍按各族约定。

历史 r001 结果继续证明旧根。新源通过相同 J/域/分支和积分线性连接新 I；新原生 I 直接求和，不再乘腿 Π₁／Γ₁的旧族级负号。不能把旧导入器的符号转换用于新原生输出。

完整推导及实际 163 源检查在 Wilson 工作区 `benchmarks/theory/montopar/shared-sign-conventions.md`、`workflow/checks/shared-sign-production/README.md`。文献依据为 [Grozin](https://arxiv.org/pdf/hep-ph/0508242)、[Davydychev–Osland–Tarasov](https://arxiv.org/pdf/hep-ph/9605348v3)、[Ruijl 等](https://arxiv.org/pdf/1703.08532v2)、[Korchemsky–Radyushkin](https://arxiv.org/pdf/hep-ph/9210281)。全域形式证明与重新积分全部内部 Π₂ 图不由此次生产迁移宣称完成。
