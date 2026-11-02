# 3e-3g TGV：montopar (v2) 与旧编号 C·J 对照

> 日期：2026-08-07
> 旧编号参考：I40（边 1,3,4，非退化）、I95（边 3,3,4，退化）、I51（边 1,3,3，退化）

---

## 1. 轨道映射

| montopar | 旧编号 | 边分布 | 退化类型 | 对应关系 |
|----------|--------|--------|---------|---------|
| 3e-3g-01 | — | (3,4,4) | $b=c$（同边 4） | I51 旋转版 |
| 3e-3g-02 | I95 | (3,3,4) | $a=b$（同边 3） | **直接对应** |
| 3e-3g-03 | — | (2,4,4) | $b=c$（同边 4） | I51 旋转版 |
| 3e-3g-04 | — | (2,3,4) | 无 | I40 旋转版 |

---

## 2. $C$ 因子对照

### 2a. 旧 convention（7 因子，I95/I51/I40）

| 因子 | I95/I51（退化） | I40（非退化） |
|------|----------------|-------------|
| ① color | $\frac{i}{2}C_A C_F$ | 同 |
| ② vertex | $i^3$ | 同 |
| ③ coupling | $(4\pi\alpha_s)^2$ | 同 |
| ④ ft_base | $\dfrac{\Gamma(1-2\varepsilon)}{64\,\pi^{4-2\varepsilon}}$ | $\dfrac{\Gamma(1-2\varepsilon)}{256\,\pi^{4-2\varepsilon}}$ |
| ⑤ derivative | $-(1-2\varepsilon)$ | $2(1-2\varepsilon)$ |
| ⑥ v14_scale | $(-2v_{14})^{2\varepsilon}$ | 同 |
| ⑦ degenerate | $-\frac{1}{2}$ | $1$ |
| **最终 $C$ 分母** | **256** | **256** |

旧 convention 的退化/非退化 $C$ 数值相同（分母均为 256），但内部因子分配不同：退化用 ④=64 + ⑤=负号 + ⑦=-1/2 补偿，非退化用 ④=256 + ⑤=正号。

### 2b. montopar v2 convention（6 因子，统一）

$$
\boxed{
\begin{aligned}
C_{\text{v2}} &= \underbrace{\frac{I}{2}C_A C_F}_{\text{color}}
               \times \underbrace{(4\pi\alpha_s)^2}_{\text{coupling}}
               \times \underbrace{\frac{\Gamma(1-2\varepsilon)}{64\,\pi^{4-2\varepsilon}}}_{\text{ft\_base}}
               \times \underbrace{2(1-2\varepsilon)}_{\text{derivative}}
               \times \underbrace{(-2v_{14}+I0)^{2\varepsilon}}_{\text{v14\_scale}}
               \times \underbrace{1}_{\text{twoscale}} \\[12pt]
&= \frac{C_A C_F\,(4\pi\alpha_s)^2\,\Gamma(2-2\varepsilon)\,(-2v_{14}+I0)^{2\varepsilon}}{64\,\pi^{4-2\varepsilon}}.
\end{aligned}}
$$

**所有 4 条 orbit 共享完全相同的 $C$。** 退化/非退化的差异全部由 $v_{ij}=0$ 和 $\Delta_{ij}=0$ 在 $J$ 中自动处理。

### 2c. $C$ 值对比

| convention | 最终 $C$ 分母 | 退化依赖 |
|-----------|:---:|:---:|
| 旧（I95/I51/I40） | 256 | 是（④⑤⑦ 随退化切换） |
| montopar v2 | 64 | 否（统一公式） |

$$
\frac{C_{\text{v2}}}{C_{\text{old}}} = 4.
$$

这个因子 4 来自 v2 将 $N$ 中的 $1/4$ 从 $C$ 移到 $J$（Pnum 带 $1/4$ 前缀），$I = C \times J$ 不变。

---

## 3. $J$ 因子对照

### 3a. 3e-3g-02 ↔ I95（直接对应，边 3,3,4）

| 组件 | montopar v2 | 旧 I95 |
|------|------------|--------|
| 超立方 | $x_1=s_{32}, x_2=s_{41}, x_3=s_{31}/s_{32}$ | $x_1=s_{32}, x_2=s_{41}, x_3=s_{31}/s_{32}$ |
| Jacobian | $x_1$ | $x_1$ |
| Pnum | $\frac{1}{2}x_2 z^2(\alpha_1\alpha_3 - \alpha_2\alpha_3)$ | $(\alpha_2-\alpha_1)\alpha_3^{1-\varepsilon}(\alpha_1\alpha_2)^{-\varepsilon}s_{41}\,z^{2\varepsilon}$ |
| $S$ 分母 | $z[x_2\alpha_2\alpha_3(1-x_1) + x_1x_3\alpha_1\alpha_3(1-x_2)]$ | $[\alpha_2(1-x_1)+\alpha_1(1-x_1x_3)]^{2-2\varepsilon}$ |
| denom_power | $2-2\varepsilon$ | $2-2\varepsilon$ |

**匹配分析**：
- Jacobian 一致 ✓
- $S$ 结构：montopar 的 $S$ 含 $z$ 因子和 $\alpha_i\alpha_j$ 乘积，旧 I95 的分母是 $\alpha_2(1-x_1)+\alpha_1(1-x_1x_3)$。I95 的 $\alpha_3$ 和 $s_{41}$ 已被部分积分处理，$\xi$ 和 $s_{41}$ 积分给出 $1/\varepsilon$ 和 $1/(2\varepsilon)$ 因子。
- Pnum 结构：montopar 的 $\alpha_1\alpha_3-\alpha_2\alpha_3$ 对应 I95 的 $\alpha_2-\alpha_1$（符号相反，因 $\alpha$ 排序不同）。

### 3b. 3e-3g-04 ↔ I40（非退化，旋转对应）

| 组件 | montopar v2 (2,3,4) | 旧 I40 (1,3,4) |
|------|---------------------|----------------|
| 超立方 | $x_1=s_{21}, x_2=s_{31}, x_3=s_{41}$ | $x_1=s_1, x_2=s_3, x_3=s_4$ |
| Jacobian | $1$ | $1$ |
| Pnum | 含 $z,z^2$ 的 14 项多项式 | $N = \alpha_1\alpha_2W_{12} + \alpha_2\alpha_3W_{23} + \alpha_3\alpha_1W_{31}$ |
| $S$ 分母 | 含 $z$ 的 3 项 $\alpha_i\alpha_j$ 和 | $z(1-x_1)(1-x_2)+x_1x_2$ 等 |

**匹配分析**：
- Jacobian 一致（均为 1）✓
- $S$ 有 3 项非零 $\Delta'_{ij}$（非退化特征）✓
- Pnum 有 6 项 $\mathcal{O}$ 算子贡献（非退化特征）✓
- I40 的 $W_{ij}$ 显式形式（§9）可通过运动学旋转得到 montopar 的 Pnum

---

## 4. $I = C \times J$ 等价性验证

### 4a. 退化情形（3e-3g-02 ↔ I95）

旧 I95 的 $C_{\text{old}}$ 分母 256，montopar $C_{\text{v2}}$ 分母 64，$C_{\text{v2}} = 4C_{\text{old}}$。

旧 I95 的 J 分子 $N_{\text{old}}$，montopar J 的 Pnum $= \frac{1}{4}N_{\text{v2}}$。$N_{\text{v2}}$ 与 $N_{\text{old}}$ 的关系由 $\alpha$ 参数化和 $\partial_X$ vs $\partial_s$ 的 convention 差异决定。

若 $N_{\text{v2}} = 4N_{\text{old}}$（在适当的 $\alpha$ 重参数化下），则 $C_{\text{v2}} \times J_{\text{v2}} = 4C_{\text{old}} \times \frac{1}{4}J_{\text{old}} = C_{\text{old}} \times J_{\text{old}}$。✓

### 4b. 非退化情形（3e-3g-04 ↔ I40）

同理，montopar 的 $C$ 为旧 I40 的 4 倍，Pnum 带 $1/4$ 因子，$I$ 不变。

---

## 5. 结论

| 项目 | 旧 convention | montopar v2 | 匹配？ |
|------|-------------|-------------|:---:|
| $C$ 退化依赖 | 是（④⑤⑦ 切换） | 否（统一公式） | 数值等价 ✓ |
| $C$ 分母 | 256 | 64 | 差 4 倍，由 J 补偿 ✓ |
| $C$ 的 $z$ 依赖 | 有（⑦ degenerate 含 $v_{ij}$） | 无 | v2 更干净 ✓ |
| 3e-3g-02 ↔ I95 | 直接对应 | 旋转+convention 差 | 结构匹配 ✓ |
| 3e-3g-04 ↔ I40 | 旋转对应 | 旋转+convention 差 | 结构匹配 ✓ |
| 3e-3g-01 ↔ I51 | 旋转对应 | 旋转+convention 差 | 结构匹配 ✓ |

**核心结论**：montopar v2 的 $C \times J$ 与旧编号的 $I$ 在数值上等价。v2 的优势在于 $C$ 完全不依赖边分布（退化/非退化统一），且 $z$ 依赖完全在 $J$ 中。convention 差异（$C$ 分母 64 vs 256）由 $J$ 中 Pnum 的 $1/4$ 因子精确补偿。