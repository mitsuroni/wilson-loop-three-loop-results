# 5e-3g 混合型积分 $I = C \times J$ 公式

> 生成方式：montopar 3loop mixed dispatch，颜色因子由 wolframscript 计算
> 结构：5e-3g = 3e-3g TGV 底座 + 1 条裸胶子对（44 orbits）
> 日期：2026-08-12

---

## 1. 统一 $C$ 因子（所有 44 条 orbit 共享）

$$
\boxed{
\begin{aligned}
C &= \underbrace{-\frac{1}{4}I\cdot C_F}_{\text{① color}}
   \times \underbrace{I^5}_{\text{② vertex}}
   \times \underbrace{(4\pi\alpha_s)^3}_{\text{③ coupling}}
   \times \underbrace{\frac{\Gamma(1-\varepsilon)\Gamma(1-2\varepsilon)}{4\pi^{2-\varepsilon}\cdot 64\pi^{4-2\varepsilon}}}_{\text{④ ft\_base}}
   \times \underbrace{2(1-2\varepsilon)}_{\text{⑤ derivative}}
   \times \underbrace{(-2v_{\text{ref}}+i0)^{3\varepsilon}}_{\text{⑥ v14\_scale}}
   \times \underbrace{(-2)^{-1}}_{\text{⑦ twoscale}} \\[16pt]
&= \frac{C_F\,(4\pi\alpha_s)^3\,\Gamma(1-\varepsilon)\Gamma(1-2\varepsilon)\,(-2v_{\text{ref}}+i0)^{3\varepsilon}}
        {512\,\pi^{6-3\varepsilon}}.
\end{aligned}}
$$

**$C$ 为实数** $\checkmark$：$I^5 \times (-\frac{1}{4}I) = \frac{1}{4}I^6 = -\frac{1}{4}$，与 ⑦ $(-2)^{-1} = -\frac{1}{2}$ 合并得 $+\frac{1}{8}$，分母 $64\times 8 = 512$。

**所有 44 条 orbit 共享完全相同的 $C$。** $z$ 依赖完全在 $J$ 中。

---

## 2. $J$ 因子通式

$$
J = \int_0^1 d^5x \int d^3\alpha\;\delta(1-{\textstyle\sum}\alpha)\;
\frac{\text{Jacobian} \times \text{lorentzZ} \times f^{-1+\varepsilon} \times \text{Pnum}}
{(\alpha_1\alpha_2\alpha_3)^{\varepsilon}\;
S^{\,2-2\varepsilon}}.
$$

**结构**：$J = J_{\text{bare}} \times J_{\text{tgv}}$，裸胶子贡献 $\text{lorentzZ}$ 和 $f$，TGV 贡献 $\text{Pnum}$ 和 $S$。

> **关键**：$\partial_X$ 只作用于 TGV 的 $S$（产生 $\text{Pnum}$），不作用于裸胶子的 $f$。裸胶子缩并结果 $v_a\!\cdot\!v_b$ 是纯标量，无 $\ell$ 因子，无 $\partial_X$ 作用。

---

## 3. 4 个 TGV 底座

$S$ 和 $\text{Pnum}$ 来自 4 个 3e-3g 底座。在 5e-3g 中，底座 $x$ 变量 $\hat{x}_1,\hat{x}_2,\hat{x}_3$ 从 5 个 $x$ 中选取对应子集。$\alpha$ 是积分哑变量，排列不影响积分值。

### 底座 A：$(3,4,4)$ — 退化 $b=c$，Pnum 仅含 $z^2$

$$
\boxed{
\begin{aligned}
\text{Pnum}_A &= \frac{1}{4}z^2\,(1-\hat{x}_1)\,\alpha_1(\alpha_2 - \alpha_3), \\[4pt]
S_A &= \hat{x}_2\,z\,\alpha_1(1-\hat{x}_1)(\alpha_3 + \alpha_2\hat{x}_3).
\end{aligned}}
$$

### 底座 B：$(3,3,4)$ — 退化 $a=b$，Pnum 仅含 $z^2$

$$
\boxed{
\begin{aligned}
\text{Pnum}_B &= \frac{1}{4}\hat{x}_2\,z^2\,\alpha_3(\alpha_2 - \alpha_1), \\[4pt]
S_B &= \hat{x}_2\,z\,\alpha_3\bigl[\alpha_2(1-\hat{x}_1) + \alpha_1(1-\hat{x}_1\hat{x}_3)\bigr].
\end{aligned}}
$$

### 底座 C：$(1,3,3)$ — 退化 $b=c$，Pnum 含 $z+z^2$

$$
\boxed{
\begin{aligned}
\text{Pnum}_C &= \frac{1}{4}(1+z)\bigl(z - \hat{x}_1(1+z)\bigr)\,\alpha_1(\alpha_2 - \alpha_3), \\[4pt]
S_C &= \alpha_1\alpha_2\bigl[z(1-\hat{x}_1)(1-\hat{x}_2\hat{x}_3) + \hat{x}_1\hat{x}_2\hat{x}_3\bigr]
     + \alpha_1\alpha_3\bigl[z(1-\hat{x}_1)(1-\hat{x}_2) + \hat{x}_1\hat{x}_2\bigr].
\end{aligned}}
$$

### 底座 D：$(1,3,4)$ — 非退化，三边互异，Pnum 含常数+$z$+$z^2$

$$
\boxed{
\begin{aligned}
\text{Pnum}_D &= -\frac{1}{4}\alpha_1\alpha_2 + \frac{3}{4}\alpha_1\alpha_3 z
              - \frac{1}{4}\alpha_1\alpha_2 z - \frac{1}{4}\alpha_2\alpha_3 z \\
            &\quad + \frac{1}{4}\hat{x}_1\alpha_1\alpha_2 - \frac{1}{4}\hat{x}_1\alpha_1\alpha_3
              + \frac{1}{4}\hat{x}_1\alpha_1\alpha_2 z - \frac{1}{4}\hat{x}_1\alpha_1\alpha_3 z \\
            &\quad - \frac{1}{4}\hat{x}_2\alpha_1\alpha_2 z + \frac{1}{4}\hat{x}_2\alpha_2\alpha_3 z
              - \frac{1}{4}\hat{x}_3\alpha_1\alpha_3 z + \frac{1}{4}\hat{x}_3\alpha_2\alpha_3 z \\
            &\quad - \frac{1}{4}\hat{x}_2\alpha_1\alpha_2 z^2 + \frac{1}{4}\hat{x}_2\alpha_2\alpha_3 z^2, \\[8pt]
S_D &= \alpha_1\alpha_2\bigl[z(1-\hat{x}_1)(1-\hat{x}_2) + \hat{x}_1\hat{x}_2\bigr]
     + \alpha_2\alpha_3\,(1-\hat{x}_2)\,\hat{x}_3\,z
     + \alpha_1\alpha_3\,\hat{x}_1\,(1-\hat{x}_3).
\end{aligned}}
$$

---

## 4. 裸胶子 $\text{lorentzZ}$ 与 $f$

裸胶子连接两条边 $a,b$，$\text{lorentzZ} = v_a\!\cdot\!v_b/v_{\text{ref}}$：

| $\text{lorentzZ}$ | 边对 |
|:---:|------|
| $z$ | (1,2), (3,4) |
| $1$ | (2,3), (4,1) |
| $-(1+z)$ | (1,3), (2,4) |

$f$ 泛型（$s_a, s_b$ 超立方变换后化为 $x$ 表达式）：

| 拓扑 | 泛型 |
|------|------|
| 相邻边 | $(1-s_a)\,s_b\,(v_a\!\cdot\!v_b/v_{\text{ref}})$ |
| 对边 {1,3} | $z(1-s_a)(1-s_b) + s_a s_b$ |
| 对边 {2,4} | $(1-s_a)(1-s_b) + z\,s_a s_b$ |

---

## 5. Orbit 逐条公式

每条 orbit 列出：底座类型、$\hat{x}$ 映射、$\text{lorentzZ}$、$f$、Jacobian。
$\text{Pnum}$ 和 $S$ 由底座公式代入 $\hat{x}$ 得到。

> **注**：$\alpha$ 是积分哑变量，排列不影响积分值。所有 orbit 的 $\text{Pnum}$ 和 $S$ 均由底座公式直接代入 $\hat{x}$ 得到。

---

### 底座 A 族：$(3,4,4)$，退化 $b=c$（01-08）

**5e-3g-01**
- $\hat{x}_1 = x_1 x_3,\;\hat{x}_2 = x_2 x_4,\;\hat{x}_3 = x_5$
- $\text{lorentzZ} = z,\quad f = x_2 z(1-x_1),\quad\text{Jacobian} = x_1 x_2^2 x_4$

**5e-3g-02**
- $\hat{x}_1 = x_1 x_3,\;\hat{x}_2 = x_2,\;\hat{x}_3 = x_4 x_5$
- $\text{lorentzZ} = z,\quad f = x_2 x_4 z(1-x_1),\quad\text{Jacobian} = x_1 x_2^2 x_4$

**5e-3g-03**
- $\hat{x}_1 = x_1,\;\hat{x}_2 = x_2,\;\hat{x}_3 = x_4 x_5$
- $\text{lorentzZ} = z,\quad f = x_2 x_4 z(1-x_1 x_3),\quad\text{Jacobian} = x_1 x_2 x_4$

**5e-3g-04**
- $\hat{x}_1 = x_1,\;\hat{x}_2 = x_2,\;\hat{x}_3 = x_4$
- $\text{lorentzZ} = z,\quad f = x_2 x_4 x_5 z(1-x_1 x_3),\quad\text{Jacobian} = x_1 x_2 x_4 x_5$

**5e-3g-05**
- $\hat{x}_1 = x_1 x_3,\;\hat{x}_2 = x_2 x_4,\;\hat{x}_3 = x_5$
- $\text{lorentzZ} = z,\quad f = x_2 z(1-x_1),\quad\text{Jacobian} = x_1 x_2^2 x_4$

**5e-3g-06**
- $\hat{x}_1 = x_1 x_2,\;\hat{x}_2 = x_1 x_3,\;\hat{x}_3 = x_4 x_5$
- $\text{lorentzZ} = z,\quad f = x_1 z(1-x_2 x_4),\quad\text{Jacobian} = x_1^2 x_3 x_4$

**5e-3g-07**
- $\hat{x}_1 = x_1,\;\hat{x}_2 = x_2,\;\hat{x}_3 = x_3 x_5$
- $\text{lorentzZ} = z,\quad f = x_2 x_4 z(1-x_1 x_3),\quad\text{Jacobian} = x_1 x_2 x_4$

**5e-3g-08**
- $\hat{x}_1 = x_1,\;\hat{x}_2 = x_2,\;\hat{x}_3 = x_3$
- $\text{lorentzZ} = z,\quad f = x_2 x_4 z(1-x_1 x_3 x_5),\quad\text{Jacobian} = x_1 x_2 x_4 x_5$

---

### 底座 C 族：$(1,3,3)$，退化 $b=c$（09-12）

**5e-3g-09**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_2 x_4,\;\hat{x}_3 = x_5$
- $\text{lorentzZ} = z,\quad f = x_2 z(1-x_1),\quad\text{Jacobian} = x_1 x_2^2 x_4$

**5e-3g-10**
- $\hat{x}_1 = x_2,\;\hat{x}_2 = x_3,\;\hat{x}_3 = x_4 x_5$
- $\text{lorentzZ} = z,\quad f = x_3 x_4 z(1-x_1),\quad\text{Jacobian} = x_1 x_3^2 x_4$

**5e-3g-11**
- $\hat{x}_1 = x_2,\;\hat{x}_2 = x_3,\;\hat{x}_3 = x_4 x_5$
- $\text{lorentzZ} = -(1+z),\quad f = 1 - x_1 - x_3 x_4 + x_1 x_3 x_4 + x_1 x_3 x_4 z,\quad\text{Jacobian} = x_1 x_3 x_4$

**5e-3g-12**
- $\hat{x}_1 = x_2,\;\hat{x}_2 = x_3,\;\hat{x}_3 = x_4$
- $\text{lorentzZ} = -(1+z),\quad f = 1 - x_1 - x_3 x_4 x_5 + x_1 x_3 x_4 x_5 + x_1 x_3 x_4 x_5 z,\quad\text{Jacobian} = x_1 x_3 x_4 x_5$

---

### 底座 D 族：$(1,3,4)$，非退化（13-28）

**5e-3g-13**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_1 x_4,\;\hat{x}_3 = x_2 x_5$
- $\text{lorentzZ} = z,\quad f = x_2 z(1-x_1),\quad\text{Jacobian} = x_1 x_2 x_4 x_5$

**5e-3g-14**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_4,\;\hat{x}_3 = x_5$
- $\text{lorentzZ} = z,\quad f = x_1 x_2 z(1-x_1),\quad\text{Jacobian} = x_1^2 x_2$

**5e-3g-15**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_1 x_4,\;\hat{x}_3 = x_5$
- $\text{lorentzZ} = z,\quad f = x_2 z(1-x_1),\quad\text{Jacobian} = x_1 x_2 x_4$

**5e-3g-16**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_4,\;\hat{x}_3 = x_5$
- $\text{lorentzZ} = z,\quad f = x_1 x_2 z(1-x_1 x_4),\quad\text{Jacobian} = x_1^2 x_2 x_4$

**5e-3g-17**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_1 x_4,\;\hat{x}_3 = x_2 x_5$
- $\text{lorentzZ} = 1,\quad f = x_2(1-x_1),\quad\text{Jacobian} = x_1 x_2 x_4 x_5$

**5e-3g-18**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_4,\;\hat{x}_3 = x_5$
- $\text{lorentzZ} = 1,\quad f = x_1 x_2(1-x_1),\quad\text{Jacobian} = x_1^2 x_2$

**5e-3g-19**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_1 x_4,\;\hat{x}_3 = x_5$
- $\text{lorentzZ} = 1,\quad f = x_2(1-x_1),\quad\text{Jacobian} = x_1 x_2 x_4$

**5e-3g-20**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_4,\;\hat{x}_3 = x_5$
- $\text{lorentzZ} = 1,\quad f = x_1 x_2(1-x_1 x_4),\quad\text{Jacobian} = x_1^2 x_2 x_4$

**5e-3g-21**
- $\hat{x}_1 = x_4,\;\hat{x}_2 = x_1 x_5,\;\hat{x}_3 = x_2$
- $\text{lorentzZ} = z,\quad f = x_2 z(1-x_1),\quad\text{Jacobian} = x_1 x_2 x_5$

**5e-3g-22**
- $\hat{x}_1 = x_4,\;\hat{x}_2 = x_1 x_5,\;\hat{x}_3 = x_2$
- $\text{lorentzZ} = 1,\quad f = x_2(1-x_1),\quad\text{Jacobian} = x_1 x_2 x_5$

**5e-3g-23**
- $\hat{x}_1 = x_4,\;\hat{x}_2 = x_5,\;\hat{x}_3 = x_2$
- $\text{lorentzZ} = z,\quad f = x_1 x_2 z(1-x_1),\quad\text{Jacobian} = x_1^2 x_2$

**5e-3g-24**
- $\hat{x}_1 = x_4,\;\hat{x}_2 = x_5,\;\hat{x}_3 = x_2$
- $\text{lorentzZ} = 1,\quad f = x_1 x_2(1-x_1),\quad\text{Jacobian} = x_1^2 x_2$

**5e-3g-25**
- $\hat{x}_1 = x_4,\;\hat{x}_2 = x_1 x_5,\;\hat{x}_3 = x_2 x_3$
- $\text{lorentzZ} = 1,\quad f = x_2(1-x_1),\quad\text{Jacobian} = x_1 x_2 x_3 x_5$

**5e-3g-26**
- $\hat{x}_1 = x_4,\;\hat{x}_2 = x_1 x_5,\;\hat{x}_3 = x_2 x_3$
- $\text{lorentzZ} = z,\quad f = x_2 z(1-x_1),\quad\text{Jacobian} = x_1 x_2 x_3 x_5$

**5e-3g-27**
- $\hat{x}_1 = x_4,\;\hat{x}_2 = x_5,\;\hat{x}_3 = x_2 x_3$
- $\text{lorentzZ} = 1,\quad f = x_1 x_2(1-x_1),\quad\text{Jacobian} = x_1^2 x_2 x_3$

**5e-3g-28**
- $\hat{x}_1 = x_4,\;\hat{x}_2 = x_5,\;\hat{x}_3 = x_2 x_3$
- $\text{lorentzZ} = z,\quad f = x_1 x_2 z(1-x_1),\quad\text{Jacobian} = x_1^2 x_2 x_3$

---

### 底座 D 族续：$(1,3,4)$ 非退化（41-44）

**5e-3g-41**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_4,\;\hat{x}_3 = x_2 x_5$
- $\text{lorentzZ} = -(1+z),\quad f = 1 - x_1 - x_2 + x_1 x_2 + x_1 x_2 z,\quad\text{Jacobian} = x_1 x_2 x_5$

**5e-3g-42**
- $\hat{x}_1 = x_2,\;\hat{x}_2 = x_3,\;\hat{x}_3 = x_4$
- $\text{lorentzZ} = -(1+z),\quad f = 1 - x_1 - x_4 x_5 + x_1 x_4 x_5 + x_1 x_4 x_5 z,\quad\text{Jacobian} = x_1 x_4 x_5$

**5e-3g-43**
- $\hat{x}_1 = x_2,\;\hat{x}_2 = x_3,\;\hat{x}_3 = x_4$
- $\text{lorentzZ} = 1,\quad f = x_1(1-x_4 x_5),\quad\text{Jacobian} = x_1 x_4 x_5$

**5e-3g-44**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_4,\;\hat{x}_3 = x_5$
- $\text{lorentzZ} = -(1+z),\quad f = 1 - x_1 - x_4 + x_1 x_4 + x_1 x_4 z,\quad\text{Jacobian} = x_1 x_4$

---

### 底座 B 族：$(3,3,4)$，退化 $a=b$（29-40）

**5e-3g-29**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_2 x_4,\;\hat{x}_3 = x_5$
- $\text{lorentzZ} = z,\quad f = x_2 z(1-x_1),\quad\text{Jacobian} = x_1 x_2^2 x_4$

**5e-3g-30**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_2,\;\hat{x}_3 = x_4 x_5$
- $\text{lorentzZ} = z,\quad f = x_2 x_4 z(1-x_1),\quad\text{Jacobian} = x_1 x_2^2 x_4$

**5e-3g-31**
- $\hat{x}_1 = x_1,\;\hat{x}_2 = x_2,\;\hat{x}_3 = x_4 x_5$
- $\text{lorentzZ} = z,\quad f = x_2 x_4 z(1-x_1 x_3),\quad\text{Jacobian} = x_1 x_2 x_4$

**5e-3g-32**
- $\hat{x}_1 = x_1,\;\hat{x}_2 = x_2,\;\hat{x}_3 = x_4$
- $\text{lorentzZ} = z,\quad f = x_2 x_4 x_5 z(1-x_1 x_3),\quad\text{Jacobian} = x_1 x_2 x_4 x_5$

**5e-3g-33**
- $\hat{x}_1 = x_1 x_3,\;\hat{x}_2 = x_2 x_4,\;\hat{x}_3 = x_5$
- $\text{lorentzZ} = 1,\quad f = x_2(1-x_1),\quad\text{Jacobian} = x_1 x_2^2 x_4$

**5e-3g-34**
- $\hat{x}_1 = x_1 x_3,\;\hat{x}_2 = x_2,\;\hat{x}_3 = x_4$
- $\text{lorentzZ} = 1,\quad f = x_2 x_4(1-x_1),\quad\text{Jacobian} = x_1 x_2^2 x_4$

**5e-3g-35**
- $\hat{x}_1 = x_1 x_3,\;\hat{x}_2 = x_2 x_4,\;\hat{x}_3 = x_5$
- $\text{lorentzZ} = -(1+z),\quad f = 1 - x_1 - x_2 x_4 + x_1 x_2 x_4 + x_1 x_2 x_4 z,\quad\text{Jacobian} = x_1 x_2^2 x_4$

**5e-3g-36**
- $\hat{x}_1 = x_1,\;\hat{x}_2 = x_3,\;\hat{x}_3 = x_4 x_5$
- $\text{lorentzZ} = -(1+z),\quad f = 1 - x_1 - x_2 + x_1 x_2 + x_1 x_2 z,\quad\text{Jacobian} = x_1 x_2 x_3$

**5e-3g-37**
- $\hat{x}_1 = x_2,\;\hat{x}_2 = x_3,\;\hat{x}_3 = x_4 x_5$
- $\text{lorentzZ} = -(1+z),\quad f = 1 - x_1 - x_2 x_4 + x_1 x_2 x_4 + x_1 x_2 x_4 z,\quad\text{Jacobian} = x_1 x_2 x_3 x_4$

**5e-3g-38**
- $\hat{x}_1 = x_2,\;\hat{x}_2 = x_3,\;\hat{x}_3 = x_4 x_5$
- $\text{lorentzZ} = 1,\quad f = x_1(1-x_2 x_4),\quad\text{Jacobian} = x_1^2 x_3 x_4$

**5e-3g-39**
- $\hat{x}_1 = x_2,\;\hat{x}_2 = x_3,\;\hat{x}_3 = x_4$
- $\text{lorentzZ} = 1,\quad f = x_1(1-x_2 x_4 x_5),\quad\text{Jacobian} = x_1^2 x_3 x_4 x_5$

**5e-3g-40**
- $\hat{x}_1 = x_3,\;\hat{x}_2 = x_4,\;\hat{x}_3 = x_2 x_5$
- $\text{lorentzZ} = z,\quad f = x_2 z(1-x_1),\quad\text{Jacobian} = x_1 x_2 x_5$

---

## 6. 汇总

| 底座族 | 退化 | orbits | Pnum 特征 | lorentzZ 取值 |
|--------|------|:---:|----------|:---:|
| A $(3,4,4)$ | $b=c$ | 01-08 | 仅 $z^2$ | $z$ |
| C $(1,3,3)$ | $b=c$ | 09-12 | $z+z^2$ | $z$, $-(1+z)$ |
| D $(1,3,4)$ | 无 | 13-28, 41-44 | 常数+$z$+$z^2$ | $z$, $1$, $-(1+z)$ |
| B $(3,3,4)$ | $a=b$ | 29-40 | 仅 $z^2$（含 $\hat{x}_2$ 因子） | $z$, $1$, $-(1+z)$ |

**关键规律**：

1. 所有 44 条 orbit 共享相同的 $C$——边分布和配对不影响 $C$
2. TGV 部分（$S$, $\text{Pnum}$）由 4 个底座完全确定，每族内仅 $\hat{x}$ 变量索引不同
3. 裸胶子贡献 $\text{lorentzZ}$ 和 $f$，不改变 TGV 的 $S$ 和 $\text{Pnum}$ 结构
4. $z$ 依赖全部在 $J$ 中