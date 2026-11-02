# 4e-3g+3g S 道 $(3,3,4,4)$：$J$ 积分化简与 Laurent 展开

> 轨道：`4e-3g+3g-04`  
> 边分布：$(3,3,4,4)$，配对 $(12)(34)$  
> 日期：2026-08-27

---

## 1. 初始条件：$C \times J$ 完整表达式

### 1a. $C$ 因子

$$
\boxed{
C = -\frac{1}{4}C_A^2 C_F \times (4\pi\alpha_s)^3 \times \frac{\Gamma(1-3\varepsilon)}{1024\,\pi^{6-3\varepsilon}} \times (1-3\varepsilon) \times (-2v_{14}+i0)^{3\varepsilon}.
}
$$

### 1b. $J$ 积分：原始形式

$$
\boxed{
\begin{aligned}
J &= \int_0^1 dx_1 dx_2 dx_3 dx_4 \;
    \underbrace{x_1 x_2}_{\text{Jacobian}}
    \int d^5\alpha\;\delta(1-{\textstyle\sum}\alpha)\;
    \frac{(\alpha_1\alpha_2\alpha_3\alpha_4\alpha_5)^{-\varepsilon}}
         {\mathcal{U}^{1+2\varepsilon}} \\
  &\quad \times \left[
    \frac{\text{Pnum}^{(1)}}{B_5^{\,3-3\varepsilon}}
    + \frac{\text{Pnum}^{(2)}}{B_5^{\,2-3\varepsilon}}
    \right].
\end{aligned}}
$$

其中各项：

| 量 | 表达式 |
|----|--------|
| $\mathcal{U}$ | $(\alpha_1+\alpha_2+\alpha_5)(\alpha_3+\alpha_4+\alpha_5) - \alpha_5^2$ |
| $B_5$ | $\alpha_5 z x_2 (\alpha_3 x_4 + \alpha_4)\bigl(\alpha_1(1-x_1x_3) + \alpha_2(1-x_1)\bigr)$ |
| $\text{Pnum}^{(1)}$ | $-4\alpha_5 z^2 (\alpha_1-\alpha_2)(\alpha_3-\alpha_4) \times B_5$ |
| $\text{Pnum}^{(2)}$ | $2\alpha_5 z^2 (\alpha_1-\alpha_2)(\alpha_3-\alpha_4)$ |
| Jacobian | $x_1 x_2$ |

> **超立方变换**：$s_{3,2}=x_1$, $s_{3,1}=x_1 x_3$, $s_{4,2}=x_2$, $s_{4,1}=x_2 x_4$

---

## 2. 关键简化：$\text{Pnum}^{(1)} \propto B_5$

**核心发现**：$\text{Pnum}^{(1)}$ 与 $B_5$ 成正比，使得两项可合并为单一 $B_5^{2-3\varepsilon}$ 分母。

$$
\frac{\text{Pnum}^{(1)}}{B_5^{3-3\varepsilon}} + \frac{\text{Pnum}^{(2)}}{B_5^{2-3\varepsilon}}
= \frac{1}{B_5^{2-3\varepsilon}}\left(\frac{\text{Pnum}^{(1)}}{B_5} + \text{Pnum}^{(2)}\right)
= \frac{-2\alpha_5 z^2 (\alpha_1-\alpha_2)(\alpha_3-\alpha_4)}{B_5^{2-3\varepsilon}}.
$$

### 简化后的 $J$

$$
\boxed{
\begin{aligned}
J &= -2z^2 \int_0^1 dx_1 dx_2 dx_3 dx_4 \; x_1 x_2
    \int d^5\alpha\;\delta(1-{\textstyle\sum}\alpha)\;
    \frac{(\alpha_1\alpha_2\alpha_3\alpha_4\alpha_5)^{-\varepsilon}}
         {\mathcal{U}^{1+2\varepsilon}} \\
  &\quad \times \frac{\alpha_5 (\alpha_1-\alpha_2)(\alpha_3-\alpha_4)}
         {B_5^{\,2-3\varepsilon}}.
\end{aligned}}
$$

---

## 3. $\delta(1-\Sigma\alpha)$ 消除

### 3a. 变量代换

利用 $\delta(1-\Sigma\alpha)$ 消去 $\alpha_5$，积分域变为 4 维单纯形：

$$
\alpha_5 = 1 - \Sigma\alpha,\qquad \Sigma\alpha \equiv \alpha_1+\alpha_2+\alpha_3+\alpha_4.
$$

$$
\int d^5\alpha\;\delta(1-{\textstyle\sum}\alpha) \to
\int_0^1 d\alpha_1\int_0^{1-\alpha_1} d\alpha_2\int_0^{1-\alpha_1-\alpha_2} d\alpha_3\int_0^{1-\alpha_1-\alpha_2-\alpha_3} d\alpha_4.
$$

### 3b. $\mathcal{U}$ 的化简

$$
\begin{aligned}
\mathcal{U} &= (\alpha_1+\alpha_2+\alpha_5)(\alpha_3+\alpha_4+\alpha_5) - \alpha_5^2 \\
&= (\alpha_1+\alpha_2)(\alpha_3+\alpha_4) + (\alpha_1+\alpha_2+\alpha_3+\alpha_4)\alpha_5 \\
&= (\alpha_1+\alpha_2)(\alpha_3+\alpha_4) + \Sigma\alpha\,(1-\Sigma\alpha).
\end{aligned}
$$

### 3c. 代换后的 $J$

$$
\boxed{
\begin{aligned}
J &= -2z^2 \int_0^1 dx_1 dx_2 dx_3 dx_4 \; x_1 x_2 \\
  &\quad \times \int_{\alpha_i \ge 0,\;\Sigma\alpha \le 1} d\alpha_1 d\alpha_2 d\alpha_3 d\alpha_4 \\
  &\quad \times \frac{(\alpha_1\alpha_2\alpha_3\alpha_4)^{-\varepsilon}\;(1-\Sigma\alpha)^{-\varepsilon}}
         {\bigl[(\alpha_1+\alpha_2)(\alpha_3+\alpha_4) + \Sigma\alpha(1-\Sigma\alpha)\bigr]^{1+2\varepsilon}} \\
  &\quad \times \frac{(1-\Sigma\alpha)\,(\alpha_1-\alpha_2)(\alpha_3-\alpha_4)}
         {B_5^{\,2-3\varepsilon}},
\end{aligned}}
$$

其中 $B_5 = (1-\Sigma\alpha)\,z\,x_2\,(\alpha_3 x_4 + \alpha_4)\bigl(\alpha_1(1-x_1x_3) + \alpha_2(1-x_1)\bigr)$。

---

## 4. $x$ 积分

### 4a. $B_5$ 的 $\alpha_5$ 因子合并

$$
(1-\Sigma\alpha) \cdot B_5^{-(2-3\varepsilon)} = (1-\Sigma\alpha)^{-(1-3\varepsilon)} \cdot z^{-(2-3\varepsilon)} \cdot x_2^{-(2-3\varepsilon)} \cdot (\alpha_3 x_4 + \alpha_4)^{-(2-3\varepsilon)} \cdot f_2^{-(2-3\varepsilon)}.
$$

其中 $f_2 = \alpha_1(1-x_1x_3) + \alpha_2(1-x_1)$。

### 4b. $x$ 积分分离

$$
\begin{aligned}
J &= -2z^{3\varepsilon} \int_{\alpha_i \ge 0,\;\Sigma\alpha \le 1} d^4\alpha\;
    \frac{(\alpha_1\alpha_2\alpha_3\alpha_4)^{-\varepsilon}\;(1-\Sigma\alpha)^{-(1-2\varepsilon)}}
         {\bigl[(\alpha_1+\alpha_2)(\alpha_3+\alpha_4) + \Sigma\alpha(1-\Sigma\alpha)\bigr]^{1+2\varepsilon}} \\
  &\quad \times (\alpha_1-\alpha_2)(\alpha_3-\alpha_4) \times I_x,
\end{aligned}
$$

$$
I_x = \int_0^1 dx_1 dx_2 dx_3 dx_4 \; x_1 x_2 \cdot x_2^{-(2-3\varepsilon)} \cdot (\alpha_3 x_4 + \alpha_4)^{-(2-3\varepsilon)} \cdot f_2^{-(2-3\varepsilon)}.
$$

### 4c. $x_2$, $x_4$, $x_3$ 逐次积分

| 变量 | 积分 | 结果 |
|:---:|------|------|
| $x_2$ | $\int_0^1 x_2^{-1+3\varepsilon} dx_2$ | $\displaystyle\frac{1}{3\varepsilon}$ |
| $x_4$ | $\int_0^1 (\alpha_3 x_4 + \alpha_4)^{-(2-3\varepsilon)} dx_4$ | $\displaystyle\frac{\alpha_4^{-(1-3\varepsilon)} - (\alpha_3+\alpha_4)^{-(1-3\varepsilon)}}{(1-3\varepsilon)\alpha_3}$ |
| $x_3$ | $\int_0^1 (A - \alpha_1 x_1 x_3)^{-(2-3\varepsilon)} dx_3$ | $\displaystyle\frac{A^{-(1-3\varepsilon)} - (A-\alpha_1 x_1)^{-(1-3\varepsilon)}}{(1-3\varepsilon)\alpha_1 x_1}$ |

其中 $A = \alpha_1 + \alpha_2(1-x_1)$，$A - \alpha_1 x_1 = (\alpha_1+\alpha_2)(1-x_1)$。

### 4d. $x_1$ 积分

$x_3$ 积分结果中的 $x_1$ 与 Jacobian 中的 $x_1$ 相消：

$$
\begin{aligned}
I_{x_1} &= \frac{1}{(1-3\varepsilon)\alpha_1} \int_0^1 dx_1\;
          \Bigl[(\alpha_1 + \alpha_2(1-x_1))^{-(1-3\varepsilon)} - ((\alpha_1+\alpha_2)(1-x_1))^{-(1-3\varepsilon)}\Bigr] \\[8pt]
        &= \frac{1}{(1-3\varepsilon)\alpha_1} \left[
            \frac{(\alpha_1+\alpha_2)^{3\varepsilon} - \alpha_1^{3\varepsilon}}{\alpha_2 \cdot 3\varepsilon}
            - \frac{(\alpha_1+\alpha_2)^{-(1-3\varepsilon)}}{3\varepsilon}
           \right].
\end{aligned}
$$

### 4e. $x$ 积分总结果

$$
I_x = \frac{1}{3\varepsilon} \cdot \frac{1}{(1-3\varepsilon)\alpha_3}\bigl[\alpha_4^{-(1-3\varepsilon)} - (\alpha_3+\alpha_4)^{-(1-3\varepsilon)}\bigr] \cdot I_{x_1}.
$$

---

## 5. 单纯形 → 超立方变换：$\alpha_i \to t_i \in [0,1]$

### 5a. 变量代换

将 4 维单纯形 $\{\alpha_i \ge 0,\;\Sigma\alpha \le 1\}$ 映射到超立方 $[0,1]^4$：

$$
\boxed{
\begin{aligned}
\alpha_1 &= t_1, \\
\alpha_2 &= t_2(1-t_1), \\
\alpha_3 &= t_3(1-t_2)(1-t_1), \\
\alpha_4 &= t_4(1-t_3)(1-t_2)(1-t_1).
\end{aligned}}
$$

Jacobian：$J_\alpha = (1-t_1)^3(1-t_2)^2(1-t_3)$。

### 5b. 关键因式分解

| 原表达式 | 变换后 |
|---------|--------|
| $1-\Sigma\alpha$ | $(1-t_1)(1-t_2)(1-t_3)(1-t_4)$ |
| $\alpha_1+\alpha_2$ | $t_1 + t_2 - t_1t_2$ |
| $\alpha_3+\alpha_4$ | $(1-t_1)(1-t_2)(t_3 + t_4 - t_3t_4)$ |
| $\alpha_1\alpha_2\alpha_3\alpha_4$ | $t_1t_2t_3t_4(1-t_1)^3(1-t_2)^2(1-t_3)$ |
| $(\alpha_1-\alpha_2)(\alpha_3-\alpha_4)$ | $(1-t_1)(1-t_2)(t_1t_2 + t_1 - t_2)(t_3t_4 + t_3 - t_4)$ |

### 5c. $\mathcal{U}$ 的化简

$$
\mathcal{U} = (\alpha_1+\alpha_2)(\alpha_3+\alpha_4) + \Sigma\alpha(1-\Sigma\alpha)
= (1-t_1)(1-t_2)\Bigl[(t_1+t_2-t_1t_2)(t_3+t_4-t_3t_4) + (1-t_3)(1-t_4)(t_1+t_2-t_1t_2 + \cdots)\Bigr].
$$

> **注**：$\mathcal{U}$ 的完整 $t$ 变量形式较复杂，但在后续积分中可逐变量处理。

### 5d. 代换后的纯 $t$ 积分

$$
\boxed{
\begin{aligned}
J &= -\frac{2z^{3\varepsilon}}{(3\varepsilon)^2(1-3\varepsilon)^2} 
    \int_0^1 dt_1 dt_2 dt_3 dt_4 \;
    (1-t_1)^3(1-t_2)^2(1-t_3) \\
  &\quad \times \frac{\bigl[t_1t_2t_3t_4(1-t_1)^3(1-t_2)^2(1-t_3)\bigr]^{-\varepsilon}}
         {\mathcal{U}^{1+2\varepsilon}} \\
  &\quad \times \bigl[(1-t_1)(1-t_2)(1-t_3)(1-t_4)\bigr]^{-(1-2\varepsilon)} \\
  &\quad \times \frac{(1-t_1)(1-t_2)(t_1t_2+t_1-t_2)(t_3t_4+t_3-t_4)}
         {t_1 \cdot t_3(1-t_2)(1-t_1)} \\
  &\quad \times \Bigl[(t_4(1-t_3)(1-t_2)(1-t_1))^{-(1-3\varepsilon)} - \bigl((1-t_1)(1-t_2)(t_3+t_4-t_3t_4)\bigr)^{-(1-3\varepsilon)}\Bigr] \\
  &\quad \times \left[
      \frac{(t_1+t_2-t_1t_2)^{3\varepsilon} - t_1^{3\varepsilon}}{t_2(1-t_1)}
      - (t_1+t_2-t_1t_2)^{-(1-3\varepsilon)}
     \right].
\end{aligned}}
$$

### 5e. $(1-t_i)$ 因子的合并

$(1-t_1)$ 和 $(1-t_2)$ 在 Jacobian、$\alpha_i^{-\varepsilon}$、$(1-\Sigma\alpha)$、分子中反复出现。合并各因子中 $(1-t_1)$ 和 $(1-t_2)$ 的幂次：

| 来源 | $(1-t_1)$ 幂次 | $(1-t_2)$ 幂次 | $(1-t_3)$ 幂次 | $(1-t_4)$ 幂次 |
|------|:---:|:---:|:---:|:---:|
| Jacobian | $+3$ | $+2$ | $+1$ | $0$ |
| $(\alpha_1\alpha_2\alpha_3\alpha_4)^{-\varepsilon}$ | $-3\varepsilon$ | $-2\varepsilon$ | $-\varepsilon$ | $0$ |
| $(1-\Sigma\alpha)^{-(1-2\varepsilon)}$ | $-(1-2\varepsilon)$ | $-(1-2\varepsilon)$ | $-(1-2\varepsilon)$ | $-(1-2\varepsilon)$ |
| 分子 $(1-t_1)(1-t_2)$ | $+1$ | $+1$ | $0$ | $0$ |
| 分母 $1/\alpha_1\alpha_3$ 中的 $(1-t_1)$ | $-1$ | $-1$ | $0$ | $0$ |

**净 $(1-t_1)$ 幂次**：$3 - 3\varepsilon - (1-2\varepsilon) + 1 - 1 = 2 - \varepsilon$

**净 $(1-t_2)$ 幂次**：$2 - 2\varepsilon - (1-2\varepsilon) + 1 - 1 = 1$

**净 $(1-t_3)$ 幂次**：$1 - \varepsilon - (1-2\varepsilon) = \varepsilon$

**净 $(1-t_4)$ 幂次**：$0 - (1-2\varepsilon) = -(1-2\varepsilon)$

所有变量积分域均为 $[0,1]$。$(1-t_i)$ 的幂次决定了 $t_i \to 1$ 端点的发散行为。

---

## 6. $\varepsilon$ 展开策略

### 6a. 极点结构

| 来源 | 极点 | 说明 |
|------|:---:|------|
| $x$ 积分 prefactor | $\varepsilon^{-3}$ | $\frac{1}{(3\varepsilon)^2(1-3\varepsilon)^2}$ |
| $t_4 \to 1$ 端点 | $\varepsilon^{-1}$ | $(1-t_4)^{-(1-2\varepsilon)}$ |
| $t$ 积分发散 | $\ge \varepsilon^{-1}$ | $\mathcal{U}^{-(1+2\varepsilon)}$ 和 $\alpha_i^{-\varepsilon}$ |

**预期总极点**：$\mathcal{O}(\varepsilon^{-4})$ 或更高。

### 6b. 后续步骤

| 步骤 | 操作 |
|:---:|------|
| 1 | 展开 $(1-t_i)$ 的 $(1-t_i)^{a+b\varepsilon}$ 形式 → $(1-t_i)^a[1 + b\varepsilon\ln(1-t_i) + \cdots]$ |
| 2 | 逐变量积分（$t_4 \to t_3 \to t_2 \to t_1$） |
| 3 | 提取 $\varepsilon$ 极点 + 有限部分 |
| 4 | 与 $C$ 因子合并得 $I = C \times J$ |

---

## 7. 与 TGV 的对比

| 特征 | TGV S 道 $(3,4,4)$ | 双 TGV S 道 $(3,3,4,4)$ |
|------|:---:|:---:|
| $x$ 维度 | 3 | 4 |
| $\alpha$ 维度（消除 $\delta$ 后） | 2 | 4 |
| $\mathcal{U}$ | $=1$ | $= (\alpha_1+\alpha_2)(\alpha_3+\alpha_4) + \Sigma\alpha(1-\Sigma\alpha)$ |
| 分母幂次 | $2-2\varepsilon$ | $2-3\varepsilon$（合并后） |
| 分子 | $\frac{1}{4}\sum\alpha_i\alpha_j W_{ij}$ | $-2\alpha_5 z^2(\alpha_1-\alpha_2)(\alpha_3-\alpha_4)$ |
| 预期极点 | $\mathcal{O}(\varepsilon^{-2})$ | $\mathcal{O}(\varepsilon^{-4})$ |