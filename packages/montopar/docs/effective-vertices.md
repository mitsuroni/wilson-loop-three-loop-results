# Wilson 圈有效顶点：TGV、FGV、双 TGV

> 有效顶点 = 颜色结构 × Lorentz 缩并（$v_i$ 向量与 QCD 顶点的 contraction）
> 日期：2026-08-11

---

## 1. 三胶子顶点（TGV）有效顶点

### 1a. QCD 三胶子顶点

$$
\Gamma^{\text{3g}}_{\mu\nu\rho}(k_1,k_2,k_3) = g_{\mu\nu}(k_1-k_2)_\rho + g_{\nu\rho}(k_2-k_3)_\mu + g_{\rho\mu}(k_3-k_1)_\nu,
\qquad k_1+k_2+k_3=0.
$$

### 1b. 有效顶点（与 Wilson 线缩并后）

三条腿分别附着在边 $a,b,c$ 上，缩并 $v_a^\mu v_b^\nu v_c^\rho$：

$$
\boxed{
\begin{aligned}
V_{\text{TGV}}(a,b,c;\ell_1,\ell_2,\ell_3)
&= f^{ABC}\;v_a^\mu v_b^\nu v_c^\rho\;\Gamma^{\text{3g}}_{\mu\nu\rho}(\ell_1,\ell_2,\ell_3) \\[4pt]
&= f^{ABC}\Bigl[
   v_{ab}\,v_c\!\cdot\!(\ell_1-\ell_2) \\
   &\qquad + v_{bc}\,v_a\!\cdot\!(\ell_2-\ell_3)
   + v_{ca}\,v_b\!\cdot\!(\ell_3-\ell_1)
   \Bigr].
\end{aligned}}
$$

其中 $v_{ij} \equiv v_i\!\cdot\!v_j$。三个 $\ell$ 满足 $\ell_1+\ell_2+\ell_3=0$。

**退化**：当 $a=b$ 时 $v_{ab}=0$，第一项自动消失。余下两项可合并为 $v_{ac}\,v_a\!\cdot\!(\ell_3-2\ell_2-\ell_1)$ 等形式。

### 1c. 完整积分（$C \times J$）

$$
I_{\text{TGV}} = \underbrace{\frac{i}{2}C_A C_F \times (4\pi\alpha_s)^2 \times \frac{\Gamma(1-2\varepsilon)}{64\pi^{4-2\varepsilon}} \times 2(1-2\varepsilon) \times (-2v_{\text{ref}}+i0)^{2\varepsilon}}_{C}
\times J_{\text{TGV}}.
$$

$$
J_{\text{TGV}} = \int_0^1 d^3x \int d^3\alpha\;\delta(1-{\textstyle\sum}\alpha)\;
\frac{\text{Jacobian} \times \text{Pnum}}
{(\alpha_1\alpha_2\alpha_3)^{\varepsilon}\;S^{\,2-2\varepsilon}},
$$

$S = \sum_{i<j} \alpha_i\alpha_j\Delta'_{ij}$，$\text{Pnum} = \frac{1}{4}\sum_{i<j} \alpha_i\alpha_j W_{ij}$（来自 $\partial_X$ 算子）。

---

## 2. 四胶子顶点（FGV）有效顶点

### 2a. QCD 四胶子顶点

$$
\begin{aligned}
\Gamma^{\text{4g}}_{\mu\nu\rho\sigma}
&= f^{ABE}f^{CDE}(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}) \\
&+ f^{ACE}f^{BDE}(g_{\mu\nu}g_{\rho\sigma} - g_{\mu\sigma}g_{\nu\rho}) \\
&+ f^{ADE}f^{BCE}(g_{\mu\nu}g_{\rho\sigma} - g_{\mu\rho}g_{\nu\sigma}).
\end{aligned}
$$

### 2b. 有效顶点（与 Wilson 线缩并后）

四条腿分别附着在边 $a,b,c,d$ 上，缩并 $v_a^\mu v_b^\nu v_c^\rho v_d^\sigma$：

$$
\boxed{
\begin{aligned}
V_{\text{FGV}}(a,b,c,d)
&= t^{A_1}t^{A_2}t^{A_3}t^{A_4}\;
   \Bigl[f^{A_1A_2E}f^{A_3A_4E}\,(v_{ac}v_{bd} - v_{ad}v_{bc}) \\
   &\qquad + f^{A_1A_3E}f^{A_2A_4E}\,(v_{ab}v_{cd} - v_{ad}v_{bc}) \\
   &\qquad + f^{A_1A_4E}f^{A_2A_3E}\,(v_{ab}v_{cd} - v_{ac}v_{bd})\Bigr].
\end{aligned}}
$$

**FGV 是纯接触顶点**——无动量依赖，无 $\ell$ 因子，无 $\partial_X$ 算子。

### 2c. 完整积分（$C \times J$）

$$
I_{\text{FGV}} = \underbrace{C_A^2 C_F \times (4\pi\alpha_s)^3 \times \frac{\Gamma(2-3\varepsilon)}{256\pi^{6-3\varepsilon}} \times \frac{1}{4} \times (-2v_{\text{ref}}+i0)^{3\varepsilon}}_{C}
\times J_{\text{FGV}}.
$$

$$
J_{\text{FGV}} = \int_0^1 d^4x \int d^4\alpha\;\delta(1-{\textstyle\sum}\alpha)\;
\frac{\text{Jacobian} \times \text{lorentzZ}}
{(\alpha_1\alpha_2\alpha_3\alpha_4)^{\varepsilon}\;Q'^{\,2-3\varepsilon}},
$$

$Q' = \sum_{i<j} \alpha_i\alpha_j\Delta'_{ij}$，$\text{lorentzZ} = \frac{1}{4}(-L_1 + L_3)$（$L_1, L_3$ 为 Lorentz 通道，$L_2$ 颜色恒为零）。

---

## 3. 双 TGV 有效顶点

### 3a. 结构

两个三胶子顶点通过内部胶子传播子耦合。有效顶点包含内部 $\ell_5$ 积分和传播子 $\frac{-i}{\ell_5^2}$——这是顶点内部结构的一部分。

$$
\boxed{
\begin{aligned}
V_{2\text{TGV}}(a,b,c,d;\ell_1,\ell_2,\ell_3,\ell_4)
&= f^{A_a A_b B}\;f^{A_c A_d B} \\
&\quad \times \int \frac{d^D\ell_5}{(2\pi)^D}\;
   \frac{-i\,g^{\nu_5\nu_6}}{\ell_5^2+i0} \\
&\quad \times (2\pi)^D\delta^D(\ell_1+\ell_2+\ell_5)\;
   (2\pi)^D\delta^D(\ell_3+\ell_4-\ell_5) \\
&\quad \times v_a^{\mu_1}v_b^{\mu_2}\;
   \Gamma^{\text{3g}}_{\mu_1\mu_2\nu_5}(-\ell_1,-\ell_2,-\ell_5) \\
&\quad \times v_c^{\mu_3}v_d^{\mu_4}\;
   \Gamma^{\text{3g}}_{\mu_3\mu_4\nu_6}(-\ell_3,-\ell_4,+\ell_5).
\end{aligned}}
$$

$\delta$ 函数强制 $\ell_5 = -(\ell_1+\ell_2) = \ell_3+\ell_4$，积分退化为代入：

$$
\boxed{
\begin{aligned}
V_{2\text{TGV}}(a,b,c,d;\ell_1,\ell_2,\ell_3,\ell_4)
&= f^{A_a A_b B}\;f^{A_c A_d B}\;
   \frac{-i}{(\ell_1+\ell_2)^2+i0}\;g^{\nu_5\nu_6} \\
&\quad \times v_a^{\mu_1}v_b^{\mu_2}\;
   \Gamma^{\text{3g}}_{\mu_1\mu_2\nu_5}(-\ell_1,-\ell_2,\;\ell_1+\ell_2) \\
&\quad \times v_c^{\mu_3}v_d^{\mu_4}\;
   \Gamma^{\text{3g}}_{\mu_3\mu_4\nu_6}(-\ell_3,-\ell_4,\;\ell_3+\ell_4).
\end{aligned}}
$$

> **与 FGV 的本质区别**：FGV 是纯接触顶点，无内部传播子。双 TGV 有效顶点含 $\frac{-i}{\ell_5^2}$ 极点——这是两个 TGV 之间的胶子传播子，是顶点内部结构的一部分。

### 3b. 完整化简表达式（FeynCalc 直接计算）

以下为 `lorentz_construct_general_V2.wl` 的 FeynCalc 输出经化简后的最终形式。9 项合并为 4 组：

$$
\boxed{
\begin{aligned}
V_{2\text{TGV}} &=
\underbrace{(v_a\!\cdot\!v_b)(v_c\!\cdot\!v_d)\;(\ell_2-\ell_1)\!\cdot\!(\ell_4-\ell_3)}_{\text{Group 1: } V_{\text{I}}} \\[4pt]
&\quad + \underbrace{(v_a\!\cdot\!v_b)\Bigl[
   2(v_c\!\cdot\!\ell_1)(v_d\!\cdot\!\ell_2) - 2(v_c\!\cdot\!\ell_2)(v_d\!\cdot\!\ell_1)
   + v_c\!\cdot\!(\ell_2-\ell_1)\,(v_d\!\cdot\!\ell_3)
   + v_d\!\cdot\!(\ell_1-\ell_2)\,(v_c\!\cdot\!\ell_4)
   \Bigr]}_{\text{Group 2: } V_{\text{II}}^{(ab)}} \\[4pt]
&\quad + \underbrace{2\,(v_c\!\cdot\!v_d)\Bigl[
   (-v_a\!\cdot\!\ell_1 - 2v_a\!\cdot\!\ell_2)(v_b\!\cdot\!\ell_4 - v_b\!\cdot\!\ell_3)
   + (2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)(v_a\!\cdot\!\ell_4 - v_a\!\cdot\!\ell_3)
   \Bigr]}_{\text{Group 3: } V_{\text{II}}^{(cd)}} \\[4pt]
&\quad + \underbrace{(v_a\!\cdot\!v_c)\,(2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)\;v_d\!\cdot\!(\ell_3-\ell_1-\ell_2)}_{\text{Group 4a}} \\
&\quad + \underbrace{(v_a\!\cdot\!v_d)\,(2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)\;v_c\!\cdot\!(\ell_1+\ell_2-\ell_4)}_{\text{Group 4b}} \\
&\quad - \underbrace{(v_b\!\cdot\!v_c)\,(v_a\!\cdot\!\ell_1 + 2v_a\!\cdot\!\ell_2)\;v_d\!\cdot\!(\ell_3-\ell_1-\ell_2)}_{\text{Group 4c}} \\
&\quad - \underbrace{(v_b\!\cdot\!v_d)\,(v_a\!\cdot\!\ell_1 + 2v_a\!\cdot\!\ell_2)\;v_c\!\cdot\!(\ell_1+\ell_2-\ell_4)}_{\text{Group 4d}}.
\end{aligned}}
$$

> **FeynCalc 验证**：`lorentz_construct_general_V2.wl` 直接计算 `Contract[Γ^A · Γ^B · g^{ν₅ν₆}]`。**此公式对任意 $a,b,c,d$ 成立，无需区分退化/非退化。**

### 3c. 完整表达式

将化简后的 4 组与颜色因子和内部传播子组合，得 $V_{2\text{TGV}}$ 的完整表达式：

$$
\boxed{
\begin{aligned}
V_{2\text{TGV}}(a,b,c,d;\ell_1,\ell_2,\ell_3,\ell_4)
&= f^{A_a A_b B}\,f^{A_c A_d B}\;
   \frac{-i}{(\ell_1+\ell_2)^2+i0} \\
&\quad \times \Bigl[
   (v_a\!\cdot\!v_b)(v_c\!\cdot\!v_d)\;(\ell_2-\ell_1)\!\cdot\!(\ell_4-\ell_3) \\[4pt]
&\quad + (v_a\!\cdot\!v_b)\Bigl(
   2(v_c\!\cdot\!\ell_1)(v_d\!\cdot\!\ell_2) - 2(v_c\!\cdot\!\ell_2)(v_d\!\cdot\!\ell_1) \\
   &\qquad\; + v_c\!\cdot\!(\ell_2-\ell_1)\,(v_d\!\cdot\!\ell_3)
   + v_d\!\cdot\!(\ell_1-\ell_2)\,(v_c\!\cdot\!\ell_4)
   \Bigr) \\[4pt]
&\quad + 2\,(v_c\!\cdot\!v_d)\Bigl(
   (-v_a\!\cdot\!\ell_1 - 2v_a\!\cdot\!\ell_2)(v_b\!\cdot\!\ell_4 - v_b\!\cdot\!\ell_3) \\
   &\qquad\; + (2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)(v_a\!\cdot\!\ell_4 - v_a\!\cdot\!\ell_3)
   \Bigr) \\[4pt]
&\quad + (v_a\!\cdot\!v_c)\,(2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)\;v_d\!\cdot\!(\ell_3-\ell_1-\ell_2) \\
&\quad + (v_a\!\cdot\!v_d)\,(2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)\;v_c\!\cdot\!(\ell_1+\ell_2-\ell_4) \\
&\quad - (v_b\!\cdot\!v_c)\,(v_a\!\cdot\!\ell_1 + 2v_a\!\cdot\!\ell_2)\;v_d\!\cdot\!(\ell_3-\ell_1-\ell_2) \\
&\quad - (v_b\!\cdot\!v_d)\,(v_a\!\cdot\!\ell_1 + 2v_a\!\cdot\!\ell_2)\;v_c\!\cdot\!(\ell_1+\ell_2-\ell_4)
   \Bigr].
\end{aligned}}
$$

> **注**：$\ell_5 = -(\ell_1+\ell_2) = \ell_3+\ell_4$ 被动量守恒确定，不是独立变量。$v_{ij}=0$ 时对应项自动消失——退化通过 $v_{ij}$ 的代数取值自动处理，无需引入额外退化因子。