# 当前归一化：一个 C，一个 I

这是 2e-02 当前数学约定，取代日常推导中并行使用 Iraw/Ihat 的方式。旧 derivation.md、physical-normalization.md 和已验收证据保留为历史来源，不能按新名称直接解释其旧公式。

一圈取 N = C_F(alpha_s/pi)、S = (-2 pi v14)^epsilon，定义

\[
G(\varepsilon)=e^{-\gamma_E\varepsilon}\Gamma(1-\varepsilon),\qquad
C=-\frac12 N G(\varepsilon) S,
\]
\[
I(z,\varepsilon)=-\frac{G(\varepsilon)}{1+z^\varepsilon}J(z,\varepsilon).
\]

这就是旧稿中的 Ihat；不再新增日常使用的 Iraw 定义。补偿因子只放进 C 一次。旧 C 的表达式与旧系数表达式仅作为迁移定理的右端，用于证明新约定与旧 Ihat 公式完全对应。

首期目标保持 I^(0)=-M(0,1)/2、I^(1)=-M(1,1)/4。历史 Claim 计划中的 Ihat 读作本页的 I；未证明的系数结论不因改名而获得验证状态。

数学参数范围为 z>0、实 |epsilon|<1/4。Lean 定义中的 N、S 为抽象实数前因子，Gamma 和 Euler 常数使用 Mathlib 定义。定义与迁移证明见 lean/Wilson/Benchmarks/TwoE02Normalization.lean。

当 sym=2、S12=z^epsilon S 且本例 crossing 成立时，新 C 和独立积分定义的 I 满足

\[
\frac22 C[J(z,\varepsilon)+z^\varepsilon J(1/z,\varepsilon)]
=N(S+z^\varepsilon S)I(z,\varepsilon).
\]

这不是通过重新定义提取系数来代替积分对应证明：I 已由原始 J 定义，symmetrized_identity 证明它与对称化等式一致；coefficient_eq_I 另外在 N、S 非零时连接抽象系数提取定义。

完整物理尺度可能涉及复数分支。本轮实数前因子模型不证明 S 的物理来源或复幂分支，P-06 保持开放。一般 n 圈采用什么补偿、颜色通道和 crossing 权重仍须单独规定，不能从这个一圈定理直接推广。
