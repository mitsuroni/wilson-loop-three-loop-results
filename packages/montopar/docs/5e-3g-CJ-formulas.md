# 5e-3g 混合型积分 $I = C \times J$ 公式

> 生成方式：montopar 3loop mixed dispatch，TGV 正规化 + C4 双翻译
> 结构：5e-3g = 3e-3g TGV 底座 + 1 条裸胶子（44 orbits → 4 种 TGV 类型）
> 更新日期：2026-08-15

---

## 1. 统一 $C$ 因子（所有 44 条 orbit 共享，A/B 组相同）

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

**因子来源**：`compute_c_mixed(n_pure=1)`，颜色因子由 wolframscript 计算。$C$ 为实数。

---

## 2. $J$ 因子通式

$$
J = \int_0^1 d^5x \int d^3\alpha\;\delta(1-{\textstyle\sum}\alpha)\;
\frac{\text{Jacobian} \times \text{lorentzZ} \times f^{-1+\varepsilon} \times \text{Pnum}}
{(\alpha_1\alpha_2\alpha_3)^{\varepsilon}\;
S^{\,2-2\varepsilon}}.
$$

**结构**：$J = J_{\text{bare}} \times J_{\text{tgv}}$。裸胶子贡献 $\text{lorentzZ}$ 和 $f$，TGV 贡献 $\text{Pnum}$ 和 $S$。所有 orbit 的 $x$ 变量均为 $x_1,x_2,x_3,x_4,x_5$。

### A/B 组约定

每条 orbit 经正规化后归入 4 种 TGV 底座之一，输出 A/B 双翻译：

| 组 | 定义 | lorentzZ 变换 |
|:---:|------|:---:|
| A | 正规化 canonical 形式 | 原始值 |
| B | A 组 C4+1 轮换 | $z \leftrightarrow 1$（$-(1+z)$ 不变） |

### 裸胶子 $\text{lorentzZ}$ 与 $f$ 泛型

| $\text{lorentzZ}$ | 边对 |
|:---:|------|
| $z$ | (1,2), (3,4) |
| $1$ | (2,3), (4,1) |
| $-(1+z)$ | (1,3), (2,4) |

| 拓扑 | $f(s)$ 泛型 |
|------|------|
| 相邻边 | $(1-s_a)\,s_b\,(v_a\!\cdot\!v_b/v_{\text{ref}})$ |
| 对边 {1,3} | $z(1-s_a)(1-s_b) + s_a s_b$ |
| 对边 {2,4} | $(1-s_a)(1-s_b) + z\,s_a s_b$ |

---

## 3. 四个 TGV 底座

$S$ 和 $\text{Pnum}$ 来自 4 个 3e-3g 底座。用 $\hat{x}_1,\hat{x}_2,\hat{x}_3$ 表示 TGV 的 $x$ 变量，每条 orbit 通过 $\hat{x} \to x$ 映射代入（$\hat{x}_i$ 是 $x_1,\dots,x_5$ 的单项式）。

### 底座 A：$(3,4,4)$ — 退化 $b=c$，邻边（对应 3e-3g-01）

$$
\boxed{
\begin{aligned}
\text{Pnum}_A &= \frac{1}{4}z^2\,(1-\hat{x}_1)\,\alpha_1(\alpha_2 - \alpha_3), \\[4pt]
S_A &= \hat{x}_2\,z\,\alpha_1(1-\hat{x}_1)(\alpha_3 + \alpha_2\hat{x}_3).
\end{aligned}}
$$

### 底座 B：$(3,3,4)$ — 退化 $a=b$，邻边（对应 3e-3g-02）

$$
\boxed{
\begin{aligned}
\text{Pnum}_B &= \frac{1}{4}\hat{x}_2\,z^2\,\alpha_3(\alpha_2 - \alpha_1), \\[4pt]
S_B &= \hat{x}_2\,z\,\alpha_3\bigl[\alpha_2(1-\hat{x}_1) + \alpha_1(1-\hat{x}_1\hat{x}_3)\bigr].
\end{aligned}}
$$

### 底座 C：$(2,4,4)$ — 退化 $b=c$，对边（对应 3e-3g-03）

$$
\boxed{
\begin{aligned}
\text{Pnum}_C &= \frac{1}{4}(1+z)\bigl(1 - \hat{x}_1(1+z)\bigr)\,\alpha_1(\alpha_2 - \alpha_3), \\[4pt]
S_C &= \alpha_1\alpha_2\bigl[z(1-\hat{x}_1)(1-\hat{x}_2\hat{x}_3) + \hat{x}_1\hat{x}_2\hat{x}_3\bigr]
     + \alpha_1\alpha_3\bigl[z(1-\hat{x}_1)(1-\hat{x}_2) + \hat{x}_1\hat{x}_2\bigr].
\end{aligned}}
$$

### 底座 D：$(2,3,4)$ — 非退化，连续三边（对应 3e-3g-04）

$$
\boxed{
\begin{aligned}
\text{Pnum}_D &= -\frac{1}{4}\alpha_1\alpha_2 + \frac{1}{4}\alpha_1\alpha_3
              - \frac{1}{4}\alpha_1\alpha_2 z + \frac{3}{4}\alpha_1\alpha_3 z + \frac{1}{4}\alpha_2\alpha_3 z \\
            &\quad + \frac{1}{4}\hat{x}_1\alpha_1\alpha_2 - \frac{1}{4}\hat{x}_1\alpha_1\alpha_3
              + \frac{1}{4}\hat{x}_1\alpha_1\alpha_2 z - \frac{1}{4}\hat{x}_1\alpha_1\alpha_3 z \\
            &\quad + \frac{1}{4}\hat{x}_2\alpha_1\alpha_2 z - \frac{1}{4}\hat{x}_2\alpha_2\alpha_3 z
              + \frac{1}{4}\hat{x}_3\alpha_1\alpha_3 z - \frac{1}{4}\hat{x}_3\alpha_2\alpha_3 z \\
            &\quad + \frac{1}{4}\hat{x}_3\alpha_1\alpha_3 z^2 - \frac{1}{4}\hat{x}_3\alpha_2\alpha_3 z^2, \\[8pt]
S_D &= \alpha_1\alpha_3
      - \hat{x}_1\alpha_1\alpha_3
      + \hat{x}_2\alpha_1\alpha_2
      - \hat{x}_3\alpha_1\alpha_3 \\
    &\quad - \hat{x}_1\hat{x}_2\alpha_1\alpha_2
      + \hat{x}_1\hat{x}_3\alpha_1\alpha_3
      + \hat{x}_3\alpha_2\alpha_3\,z \\
    &\quad + \hat{x}_1\hat{x}_3\alpha_1\alpha_3\,z
      - \hat{x}_2\hat{x}_3\alpha_2\alpha_3\,z.
\end{aligned}}
$$

---

## 4. Orbit 逐条数据

以下数据全部来自 montopar 实际输出。$\text{Pnum}$ 和 $S$ 由底座公式代入 $\hat{x}$ 得到（完整多项式见 `montopar.json`）。

### 底座 A 族：$(3,4,4)$，退化 $b=c$（10 条 orbit）

| Orbit | 组 | lorentzZ | $f_1$ | Jacobian | $\hat{x}_1$ | $\hat{x}_2$ | $\hat{x}_3$ |
|-------|:---:|:---:|------|------|:---:|:---:|:---:|
| 01 | A | $z$ | $x_2z(1-x_1)$ | $x_1x_2^2x_4$ | $x_1x_3$ | $x_2x_4$ | $x_5$ |
| 01 | B | $1$ | $x_2(1-x_1)$ | $x_1x_2^2x_4$ | $x_1x_3$ | $x_2x_4$ | $x_5$ |
| 02 | A | $z$ | $x_2x_4z(1-x_1)$ | $x_1x_2^2x_4$ | $x_1x_3$ | $x_2$ | $x_4x_5$ |
| 02 | B | $1$ | $x_2x_4(1-x_1)$ | $x_1x_2^2x_4$ | $x_1x_3$ | $x_2$ | $x_4x_5$ |
| 03 | A | $z$ | $x_2x_4z(1-x_1x_3)$ | $x_1x_2^2x_4$ | $x_1$ | $x_2$ | $x_4x_5$ |
| 03 | B | $1$ | $x_2x_4(1-x_1x_3)$ | $x_1x_2^2x_4$ | $x_1$ | $x_2$ | $x_4x_5$ |
| 04 | A | $z$ | $x_2x_4x_5z(1-x_1x_3)$ | $x_1x_2^2x_4$ | $x_1$ | $x_2$ | $x_4$ |
| 04 | B | $1$ | $x_2x_4x_5(1-x_1x_3)$ | $x_1x_2^2x_4$ | $x_1$ | $x_2$ | $x_4$ |
| 11 | A | $-(1+z)$ | $1-x_1-x_3x_4+x_1x_3x_4+x_1x_3x_4z$ | $x_3^2x_4$ | $x_2$ | $x_3$ | $x_4$ |
| 11 | B | $-(1+z)$ | $z-x_1z+x_1x_3x_4-x_3x_4z+x_1x_3x_4z$ | $x_3^2x_4$ | $x_2$ | $x_3$ | $x_4$ |
| 12 | A | $-(1+z)$ | $1-x_1-x_3x_4x_5+x_1x_3x_4x_5+x_1x_3x_4x_5z$ | $x_3^2x_4$ | $x_2$ | $x_3$ | $x_4$ |
| 12 | B | $-(1+z)$ | $z-x_1z+x_1x_3x_4x_5-x_3x_4x_5z+x_1x_3x_4x_5z$ | $x_3^2x_4$ | $x_2$ | $x_3$ | $x_4$ |
| 17 | A | $1$ | $x_2(1-x_1)$ | $x_2x_3$ | $x_2x_4$ | $x_2$ | $x_3$ |
| 17 | B | $z$ | $x_2z(1-x_1)$ | $x_2x_3$ | $x_2x_4$ | $x_2$ | $x_3$ |
| 18 | A | $1$ | $x_1(1-x_3x_4)$ | $x_3^2x_4$ | $x_2$ | $x_3$ | $x_4$ |
| 18 | B | $z$ | $x_1z(1-x_3x_4)$ | $x_3^2x_4$ | $x_2$ | $x_3$ | $x_4$ |
| 19 | A | $1$ | $x_1(1-x_3x_4x_5)$ | $x_3^2x_4$ | $x_2$ | $x_3$ | $x_4$ |
| 19 | B | $z$ | $x_1z(1-x_3x_4x_5)$ | $x_3^2x_4$ | $x_2$ | $x_3$ | $x_4$ |
| 32 | A | $-(1+z)$ | $z-x_1z+x_1x_2-x_2z+x_1x_2z$ | $x_1x_3$ | $x_1x_4$ | $x_1$ | $x_3$ |
| 32 | B | $-(1+z)$ | $1-x_1-x_2+x_1x_2+x_1x_2z$ | $x_1x_3$ | $x_1x_4$ | $x_1$ | $x_3$ |

---

### 底座 B 族：$(3,3,4)$，退化 $a=b$（10 条 orbit）

| Orbit | 组 | lorentzZ | $f_1$ | Jacobian | $\hat{x}_1$ | $\hat{x}_2$ | $\hat{x}_3$ |
|-------|:---:|:---:|------|------|:---:|:---:|:---:|
| 05 | A | $z$ | $x_2z(1-x_1)$ | $x_1^2x_2x_3$ | $x_1x_4$ | $x_2x_4$ | $x_5$ |
| 05 | B | $1$ | $x_2(1-x_1)$ | $x_1^2x_2x_3$ | $x_1x_4$ | $x_2x_4$ | $x_5$ |
| 06 | A | $z$ | $x_1z(1-x_2x_4)$ | $x_1x_2^2x_4$ | $x_1x_3$ | $x_1x_3$ | $x_4x_5$ |
| 06 | B | $1$ | $x_1(1-x_2x_4)$ | $x_1x_2^2x_4$ | $x_1x_3$ | $x_1x_3$ | $x_4x_5$ |
| 07 | A | $z$ | $x_2x_4z(1-x_1x_3)$ | $x_1^2x_2x_3$ | $x_1$ | $x_2$ | $x_3x_5$ |
| 07 | B | $1$ | $x_2x_4(1-x_1x_3)$ | $x_1^2x_2x_3$ | $x_1$ | $x_2$ | $x_3x_5$ |
| 08 | A | $z$ | $x_2x_4z(1-x_1x_3x_5)$ | $x_1^2x_2x_3$ | $x_1$ | $x_2$ | $x_3$ |
| 08 | B | $1$ | $x_2x_4(1-x_1x_3x_5)$ | $x_1^2x_2x_3$ | $x_1$ | $x_2$ | $x_3$ |
| 16 | A | $-(1+z)$ | $1-x_1-x_3x_5+x_1x_3x_5+x_1x_3x_5z$ | $x_2x_3$ | $x_3$ | $x_3$ | $x_5$ |
| 16 | B | $-(1+z)$ | $z-x_1z+x_1x_3x_5-x_3x_5z+x_1x_3x_5z$ | $x_2x_3$ | $x_3$ | $x_3$ | $x_5$ |
| 20 | A | $1$ | $x_2(1-x_1)$ | $x_2^2x_4$ | $x_3$ | $x_3$ | $x_5$ |
| 20 | B | $z$ | $x_2z(1-x_1)$ | $x_2^2x_4$ | $x_3$ | $x_3$ | $x_5$ |
| 21 | A | $1$ | $x_2x_4(1-x_1)$ | $x_2^2x_4$ | $x_3$ | $x_3$ | $x_5$ |
| 21 | B | $z$ | $x_2x_4z(1-x_1)$ | $x_2^2x_4$ | $x_3$ | $x_3$ | $x_5$ |
| 31 | A | $1$ | $x_1(1-x_3x_5)$ | $x_2x_3$ | $x_3$ | $x_3$ | $x_5$ |
| 31 | B | $z$ | $x_1z(1-x_3x_5)$ | $x_2x_3$ | $x_3$ | $x_3$ | $x_5$ |
| 36 | A | $-(1+z)$ | $z-x_1z+x_1x_2-x_2z+x_1x_2z$ | $x_1^2x_4$ | $x_3$ | $x_3$ | $x_5$ |
| 36 | B | $-(1+z)$ | $1-x_1-x_2+x_1x_2+x_1x_2z$ | $x_1^2x_4$ | $x_3$ | $x_3$ | $x_5$ |
| 37 | A | $-(1+z)$ | $z-x_1z+x_1x_2x_4-x_2x_4z+x_1x_2x_4z$ | $x_2^2x_4$ | $x_3$ | $x_3$ | $x_5$ |
| 37 | B | $-(1+z)$ | $1-x_1-x_2x_4+x_1x_2x_4+x_1x_2x_4z$ | $x_2^2x_4$ | $x_3$ | $x_3$ | $x_5$ |

---

### 底座 C 族：$(2,4,4)$，退化 $b=c$，对边（11 条 orbit）

| Orbit | 组 | lorentzZ | $f_1$ | Jacobian | $\hat{x}_1$ | $\hat{x}_2$ | $\hat{x}_3$ |
|-------|:---:|:---:|------|------|:---:|:---:|:---:|
| 09 | A | $1$ | $x_2(1-x_1)$ | $x_2^2x_4$ | $x_3$ | $x_2x_4$ | $x_5$ |
| 09 | B | $z$ | $x_2z(1-x_1)$ | $x_2^2x_4$ | $x_3$ | $x_2x_4$ | $x_5$ |
| 10 | A | $1$ | $x_3x_4(1-x_1)$ | $x_3^2x_4$ | $x_2$ | $x_3$ | $x_4x_5$ |
| 10 | B | $z$ | $x_3x_4z(1-x_1)$ | $x_3^2x_4$ | $x_2$ | $x_3$ | $x_4x_5$ |
| 22 | A | $-(1+z)$ | $z-x_1z+x_1x_2-x_2z+x_1x_2z$ | $x_1x_2^2x_4$ | $x_3$ | $x_2$ | $x_4$ |
| 22 | B | $-(1+z)$ | $1-x_1-x_2+x_1x_2+x_1x_2z$ | $x_1x_2^2x_4$ | $x_3$ | $x_2$ | $x_4$ |
| 23 | A | $-(1+z)$ | $z-x_1z+x_1x_2x_4-x_2x_4z+x_1x_2x_4z$ | $x_1x_2^2x_4$ | $x_3$ | $x_2$ | $x_4$ |
| 23 | B | $-(1+z)$ | $1-x_1-x_2x_4+x_1x_2x_4+x_1x_2x_4z$ | $x_1x_2^2x_4$ | $x_3$ | $x_2$ | $x_4$ |
| 24 | A | $-(1+z)$ | $z-x_1x_3z-x_2x_4z+x_1x_2x_3x_4+x_1x_2x_3x_4z$ | $x_1x_2^2x_4$ | $x_1$ | $x_2$ | $x_4$ |
| 24 | B | $-(1+z)$ | $1-x_1x_3-x_2x_4+x_1x_2x_3x_4+x_1x_2x_3x_4z$ | $x_1x_2^2x_4$ | $x_1$ | $x_2$ | $x_4$ |
| 25 | A | $-(1+z)$ | $z-x_1x_3z-x_2x_4x_5z+x_1x_2x_3x_4x_5+x_1x_2x_3x_4x_5z$ | $x_1x_2^2x_4$ | $x_1$ | $x_2$ | $x_4$ |
| 25 | B | $-(1+z)$ | $1-x_1x_3-x_2x_4x_5+x_1x_2x_3x_4x_5+x_1x_2x_3x_4x_5z$ | $x_1x_2^2x_4$ | $x_1$ | $x_2$ | $x_4$ |
| 26 | A | $1$ | $x_2(1-x_1)$ | $x_2x_3$ | $x_4$ | $x_2$ | $x_3$ |
| 26 | B | $z$ | $x_2z(1-x_1)$ | $x_2x_3$ | $x_4$ | $x_2$ | $x_3$ |
| 30 | A | $z$ | $x_1z(1-x_2x_4)$ | $x_2x_3$ | $x_2$ | $x_3$ | $x_5$ |
| 30 | B | $1$ | $x_1(1-x_2x_4)$ | $x_2x_3$ | $x_2$ | $x_3$ | $x_5$ |
| 38 | A | $z$ | $x_1z(1-x_2x_4)$ | $x_2^2x_4$ | $x_3$ | $x_2$ | $x_4$ |
| 38 | B | $1$ | $x_1(1-x_2x_4)$ | $x_2^2x_4$ | $x_3$ | $x_2$ | $x_4$ |
| 39 | A | $z$ | $x_1z(1-x_2x_4x_5)$ | $x_2^2x_4$ | $x_3$ | $x_2$ | $x_4$ |
| 39 | B | $1$ | $x_1(1-x_2x_4x_5)$ | $x_2^2x_4$ | $x_3$ | $x_2$ | $x_4$ |
| 44 | A | $-(1+z)$ | $1-x_1-x_2+x_1x_2+x_1x_2z$ | $x_4$ | $x_3$ | $x_4$ | $x_5$ |
| 44 | B | $-(1+z)$ | $z-x_1z+x_1x_2-x_2z+x_1x_2z$ | $x_4$ | $x_3$ | $x_4$ | $x_5$ |

---

### 底座 D 族：$(2,3,4)$，非退化连续三边（13 条 orbit）

| Orbit | 组 | lorentzZ | $f_1$ | Jacobian | $\hat{x}_1$ | $\hat{x}_2$ | $\hat{x}_3$ |
|-------|:---:|:---:|------|------|:---:|:---:|:---:|
| 13 | A | $1$ | $x_2(1-x_1)$ | $x_1x_2$ | $x_3$ | $x_1x_4$ | $x_2x_5$ |
| 13 | B | $z$ | $x_2z(1-x_1)$ | $x_1x_2$ | $x_3$ | $x_1x_4$ | $x_2x_5$ |
| 14 | A | $1$ | $x_1(1-x_3x_5)$ | $x_1x_3$ | $x_3$ | $x_4$ | $x_5$ |
| 14 | B | $z$ | $x_1z(1-x_3x_5)$ | $x_1x_3$ | $x_3$ | $x_4$ | $x_5$ |
| 15 | A | $1$ | $x_3x_5(1-x_2x_4)$ | $x_2x_3$ | $x_3$ | $x_1x_4$ | $x_5$ |
| 15 | B | $z$ | $x_3x_5z(1-x_2x_4)$ | $x_2x_3$ | $x_3$ | $x_1x_4$ | $x_5$ |
| 27 | A | $-(1+z)$ | $z-x_1z+x_1x_2-x_2z+x_1x_2z$ | $x_1x_2$ | $x_4$ | $x_5$ | $x_2x_3$ |
| 27 | B | $-(1+z)$ | $1-x_1-x_2+x_1x_2+x_1x_2z$ | $x_1x_2$ | $x_4$ | $x_5$ | $x_2x_3$ |
| 28 | A | $-(1+z)$ | $z-x_1z+x_1x_3x_5-x_3x_5z+x_1x_3x_5z$ | $x_1x_3$ | $x_4$ | $x_5$ | $x_2x_3$ |
| 28 | B | $-(1+z)$ | $1-x_1-x_3x_5+x_1x_3x_5+x_1x_3x_5z$ | $x_1x_3$ | $x_4$ | $x_5$ | $x_2x_3$ |
| 29 | A | $-(1+z)$ | $z-x_1x_4z-x_3x_5z+x_1x_3x_4x_5+x_1x_3x_4x_5z$ | $x_1x_3$ | $x_4$ | $x_5$ | $x_2x_3$ |
| 29 | B | $-(1+z)$ | $1-x_1x_4-x_3x_5+x_1x_3x_4x_5+x_1x_3x_4x_5z$ | $x_1x_3$ | $x_4$ | $x_5$ | $x_2x_3$ |
| 33 | A | $z$ | $x_2z(1-x_1)$ | $x_1x_2$ | $x_4$ | $x_1x_5$ | $x_2$ |
| 33 | B | $1$ | $x_2(1-x_1)$ | $x_1x_2$ | $x_4$ | $x_1x_5$ | $x_2$ |
| 34 | A | $z$ | $x_1z(1-x_2x_5)$ | $x_1x_2$ | $x_4$ | $x_5$ | $x_2$ |
| 34 | B | $1$ | $x_1(1-x_2x_5)$ | $x_1x_2$ | $x_4$ | $x_5$ | $x_2$ |
| 35 | A | $z$ | $x_2x_5z(1-x_1x_4)$ | $x_1x_2$ | $x_4$ | $x_5$ | $x_2$ |
| 35 | B | $1$ | $x_2x_5(1-x_1x_4)$ | $x_1x_2$ | $x_4$ | $x_5$ | $x_2$ |
| 40 | A | $1$ | $x_2(1-x_1)$ | $x_2$ | $x_3$ | $x_4$ | $x_2x_5$ |
| 40 | B | $z$ | $x_2z(1-x_1)$ | $x_2$ | $x_3$ | $x_4$ | $x_2x_5$ |
| 41 | A | $-(1+z)$ | $1-x_1-x_2+x_1x_2+x_1x_2z$ | $x_2$ | $x_3$ | $x_4$ | $x_2x_5$ |
| 41 | B | $-(1+z)$ | $z-x_1z+x_1x_2-x_2z+x_1x_2z$ | $x_2$ | $x_3$ | $x_4$ | $x_2x_5$ |
| 42 | A | $-(1+z)$ | $1-x_1-x_4x_5+x_1x_4x_5+x_1x_4x_5z$ | $x_4$ | $x_2$ | $x_3$ | $x_4$ |
| 42 | B | $-(1+z)$ | $z-x_1z+x_1x_4x_5-x_4x_5z+x_1x_4x_5z$ | $x_4$ | $x_2$ | $x_3$ | $x_4$ |
| 43 | A | $z$ | $x_1z(1-x_4x_5)$ | $x_4$ | $x_2$ | $x_3$ | $x_4$ |
| 43 | B | $1$ | $x_1(1-x_4x_5)$ | $x_4$ | $x_2$ | $x_3$ | $x_4$ |

---

## 5. 汇总

| 底座族 | 3e-3g | orbits | 数量 | $\hat{x}$ 模式 |
|--------|:---:|------|:---:|------|
| A $(3,4,4)$ | 01 | 01-04,11,12,17-19,32 | 10 | $\hat{x}_1\in\{x_1,x_1x_3,x_2,x_2x_4,x_1x_4\}$ |
| B $(3,3,4)$ | 02 | 05-08,16,20,21,31,36,37 | 10 | $\hat{x}_2$ 含 $x_3$ 或 $x_2x_4$ 等 |
| C $(2,4,4)$ | 03 | 09,10,22-26,30,38,39,44 | 11 | $\hat{x}_1\in\{x_1,x_2,x_3,x_4\}$ |
| D $(2,3,4)$ | 04 | 13-15,27-29,33-35,40-43 | 13 | 三项 $\alpha$ 结构 |

**A/B 组 lorentzZ 变换**：$z \leftrightarrow 1$（$-(1+z)$ 不变）

**关键规律**：

1. 所有 44 条 orbit 共享相同的 $C$
2. TGV 部分（$S$, $\text{Pnum}$）由 4 个底座完全确定，每族内仅 $\hat{x}$ 变量索引不同
3. 裸胶子贡献 $\text{lorentzZ}$ 和 $f$，不改变 TGV 的 $S$ 和 $\text{Pnum}$ 结构
4. $f$ 由裸胶子边对决定：相邻边 → $x_a(1-x_b)z$ 形式，对边 → $1-x_a-x_b+x_ax_b+x_ax_bz$ 形式
5. 超立方变换前的 $s$ 变量形式见 `montopar.json` 的 `J_raw_s` / `J_raw_s_rotated`