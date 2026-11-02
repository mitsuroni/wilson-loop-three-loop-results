# 4e-3g+3g S 道：从通用 $V_{2\text{TGV}}$ 推导 $C \times J$

> 风格对标 `test/TGV/I40_因子分解推导.md`——将 $V_{2\text{TGV}}$ 当作"单位顶点"，给出完整 $C, J$
> 轨道：4e-3g+3g-04 — 边分布 $(3,3,4,4)$，S 道配对 $(12)(34)$
> $D=4-2\varepsilon$，$v_i^2=0$，$v_{\text{ref}}=v_{14}$，$A_{\text{ref}}\equiv -2v_{\text{ref}}+i0$
> $z = v_{34}/v_{\text{ref}}$，$v_{13}=v_{24}=-(1+z)v_{\text{ref}}$，$v_{14}=v_{23}=v_{\text{ref}}$

---

## 1. 动量空间起点：$V_{2\text{TGV}}$ 通用解

双 TGV 有效顶点的通用解（FeynCalc `lorentz_construct_general_V2.wl` 直接计算，$\ell_5$ 已由动量守恒 $\ell_5=-(\ell_1+\ell_2)=\ell_3+\ell_4$ 消去）：

$$
\boxed{
\begin{aligned}
V_{2\text{TGV}}(a,b,c,d;\ell_1,\ell_2,\ell_3,\ell_4)
&= f^{A_a A_b B}\,f^{A_c A_d B}\;
   \frac{-i}{(\ell_1+\ell_2)^2+i0} \\
&\quad \times \Bigl[
   (v_a\!\cdot\!v_b)(v_c\!\cdot\!v_d)\;(\ell_2-\ell_1)\!\cdot\!(\ell_4-\ell_3) \\
&\quad + (v_a\!\cdot\!v_b)\Bigl(
   2(v_c\!\cdot\!\ell_1)(v_d\!\cdot\!\ell_2) - 2(v_c\!\cdot\!\ell_2)(v_d\!\cdot\!\ell_1) \\
   \quad + v_c\!\cdot\!(\ell_2-\ell_1)(v_d\!\cdot\!\ell_3)
   + v_d\!\cdot\!(\ell_1-\ell_2)(v_c\!\cdot\!\ell_4)\Bigr) \\
&\quad + 2(v_c\!\cdot\!v_d)\Bigl(
   (-v_a\!\cdot\!\ell_1-2v_a\!\cdot\!\ell_2)(v_b\!\cdot\!\ell_4-v_b\!\cdot\!\ell_3) \\
   \quad + (2v_b\!\cdot\!\ell_1+v_b\!\cdot\!\ell_2)(v_a\!\cdot\!\ell_4-v_a\!\cdot\!\ell_3)\Bigr) \\
&\quad + (v_a\!\cdot\!v_c)(2v_b\!\cdot\!\ell_1+v_b\!\cdot\!\ell_2)\;v_d\!\cdot\!(\ell_3-\ell_1-\ell_2) \\
&\quad + (v_a\!\cdot\!v_d)(2v_b\!\cdot\!\ell_1+v_b\!\cdot\!\ell_2)\;v_c\!\cdot\!(\ell_1+\ell_2-\ell_4) \\
&\quad - (v_b\!\cdot\!v_c)(v_a\!\cdot\!\ell_1+2v_a\!\cdot\!\ell_2)\;v_d\!\cdot\!(\ell_3-\ell_1-\ell_2) \\
&\quad - (v_b\!\cdot\!v_d)(v_a\!\cdot\!\ell_1+2v_a\!\cdot\!\ell_2)\;v_c\!\cdot\!(\ell_1+\ell_2-\ell_4)
   \Bigr].
\end{aligned}}
$$

完整积分 = 外部传播子 $\times$ Wilson 线因子 $\times$ $V_{2\text{TGV}}$：

$$
\boxed{
\begin{aligned}
I_S &= \int_0^1 ds_{3,2}\int_0^{s_{3,2}} ds_{3,1}\int_0^1 ds_{4,2}\int_0^{s_{4,2}} ds_{4,1} \int \frac{d^D\ell_1\cdots d^D\ell_4}{(2\pi)^{4D}} \\
&\quad \times \prod_{i=1}^4 (ig\,v_{e_i}^{\mu_i} t^{A_i})\,
   \frac{-i\,g_{\mu_i\nu_i}\,\delta^{A_i B_i}}{\ell_i^2+i0}\,
   e^{-i\ell_i\cdot X_i} \\
&\quad \times V_{2\text{TGV}}(a,b,c,d;\ell_1,\ell_2,\ell_3,\ell_4),
\end{aligned}}
$$

**S 道配对**：$(a,b,c,d) = (3,3,4,4)$，$e_i \in \{3,3,4,4\}$。

N=4 方形 Wilson 圈的角点：$z_1=0$，$z_2=v_1$，$z_3=v_1+v_2$，$z_4=v_1+v_2+v_3$。边的参数化：

$$
X_3(s) = z_3 + s\,v_3 = v_1+v_2+s\,v_3,\qquad
X_4(s) = z_4 + s\,v_4 = v_1+v_2+v_3+s\,v_4.
$$

四条外腿的附着点位置：

$$
\boxed{
\begin{aligned}
X_1 &= X_3(s_{3,1}) = v_1+v_2+s_{3,1}\,v_3, \\
X_2 &= X_3(s_{3,2}) = v_1+v_2+s_{3,2}\,v_3, \\
X_3 &= X_4(s_{4,1}) = v_1+v_2+v_3+s_{4,1}\,v_4, \\
X_4 &= X_4(s_{4,2}) = v_1+v_2+v_3+s_{4,2}\,v_4.
\end{aligned}}
$$

Schwinger 参数的积分域：$s_{3,1}\in[0,s_{3,2}]$，$s_{3,2}\in[0,1]$，$s_{4,1}\in[0,s_{4,2}]$，$s_{4,2}\in[0,1]$。

---

## 2. 颜色因子

$V_{2\text{TGV}}$ 中的颜色 $f^{A_3 A_3 B}f^{A_4 A_4 B}$ 由 FeynCalc `SUNSimplify` 求值：

$$
\text{① color} = -\frac{1}{4}C_A^2 C_F.
$$

提取颜色后（$t^{A_i}$ 与 $\delta^{A_i B_i}$ 缩并后消去）：

$$
\boxed{
\begin{aligned}
I_S &= -\frac{1}{4}C_A^2 C_F
   \int_0^1 ds_{3,2}\!\int_0^{s_{3,2}}\! ds_{3,1}\int_0^1 ds_{4,2}\!\int_0^{s_{4,2}}\! ds_{4,1} \int \frac{d^D\ell_1\cdots d^D\ell_4}{(2\pi)^{4D}} \\
&\quad \times \prod_{i=1}^4 \frac{-i\,v_{e_i}^{\mu_i}\,g_{\mu_i\nu_i}}{\ell_i^2+i0}\,
   e^{-i\ell_i\cdot X_i} \\
&\quad \times \frac{-i}{(\ell_1+\ell_2)^2+i0}\;
   \mathcal{L}(3,3,4,4),
\end{aligned}}
$$

其中 $\mathcal{L}(a,b,c,d)$ 是 $V_{2\text{TGV}}$ 的 Lorentz 部分（§1 中的 $[\cdots]$ 括号）。

---

## 3. Lorentz 缩并：代入 S 道

### 3a. 传播子缩并

$g_{\mu_i\nu_i}$ 缩并 $v_{e_i}^{\mu_i}$，将 Wilson 线向量传递到内部：

$$
v_{e_i}^{\mu_i}\,g_{\mu_i\nu_i} \to v_{e_i}^{\nu_i}.
$$

### 3b. 代入 $(a,b,c,d) = (3,3,4,4)$

$$
v_a = v_b = v_3,\quad v_c = v_d = v_4,
$$

$$
v_a\!\cdot\!v_b = v_{33} = 0,\quad v_c\!\cdot\!v_d = v_{44} = 0,\quad
v_a\!\cdot\!v_c = v_a\!\cdot\!v_d = v_b\!\cdot\!v_c = v_b\!\cdot\!v_d = v_{34} = z.
$$

### 3c. 逐组消失

| 组 | 因子 | S 道 | 状态 |
----|------|------|-----
| G1 | $(v_a\!\cdot\!v_b)(v_c\!\cdot\!v_d)$ | $0 \times 0$ | 消失 |
| G2 | $(v_a\!\cdot\!v_b)$ | $0$ | 消失 |
| G3 | $(v_c\!\cdot\!v_d)$ | $0$ | 消失 |
| G4a-d | $(v_a\!\cdot\!v_c)$ 等 $= z$ | $z$ | 存活 |

### 3d. G4a-d 合并

$$
\begin{aligned}
\text{G4a}+\text{G4c} &= z\,(v_3\!\cdot\!\ell_1-v_3\!\cdot\!\ell_2)\;v_4\!\cdot\!(\ell_3-\ell_1-\ell_2), \\
\text{G4b}+\text{G4d} &= z\,(v_3\!\cdot\!\ell_1-v_3\!\cdot\!\ell_2)\;v_4\!\cdot\!(\ell_1+\ell_2-\ell_4).
\end{aligned}
$$

两项相加（$v_4\!\cdot\!(\ell_3-\ell_1-\ell_2) + v_4\!\cdot\!(\ell_1+\ell_2-\ell_4) = v_4\!\cdot\!(\ell_3-\ell_4)$）：

$$
\boxed{\mathcal{L}_S = z\,\big(v_3\!\cdot\!(\ell_1-\ell_2)\big)\,\big(v_4\!\cdot\!(\ell_3-\ell_4)\big).}
$$

> **FeynCalc 验证**：`lorentz_construct02_S.wl` $\checkmark$

代入后：

$$
\boxed{
\begin{aligned}
I_S &= -\frac{1}{4}C_A^2 C_F
   \int_0^1 ds_{3,2}\!\int_0^{s_{3,2}}\! ds_{3,1}\int_0^1 ds_{4,2}\!\int_0^{s_{4,2}}\! ds_{4,1} \int \frac{d^D\ell_1\cdots d^D\ell_4}{(2\pi)^{4D}} \\
&\quad \times \prod_{i=1}^4 \frac{-i\,e^{-i\ell_i\cdot X_i}}{\ell_i^2+i0}
   \times \frac{-i}{(\ell_1+\ell_2)^2+i0} \\
&\quad \times z\,\big(v_3\!\cdot\!(\ell_1-\ell_2)\big)\,\big(v_4\!\cdot\!(\ell_3-\ell_4)\big).
\end{aligned}}
$$

> **与 I40 对照**：I40 的 $V$ 有 3 项，每项 1 个 $v\!\cdot\!\ell$。S 道 $V_S$ 仅 1 项但有 **2 个** $v\!\cdot\!\ell$——$\partial_X$ 是二次型。

---

## 4. 顶点 $i$ 因子

- 4 个 Wilson 线顶点：$(ig)^4 = g^4$（$i^4=1$）
- 2 个 TGV 顶点：$g^2$
- 合计 $g^6 = (4\pi\alpha_s)^3$

5 条传播子分子 $(-i)^5$ 暂留积分内（归入 §6 的 FT prefactor）。提取 $g^6$：

$$
\boxed{
\begin{aligned}
I_S &= -\frac{1}{4}C_A^2 C_F\;g^6
   \int_0^1 ds_{3,2}\!\int_0^{s_{3,2}}\! ds_{3,1}\int_0^1 ds_{4,2}\!\int_0^{s_{4,2}}\! ds_{4,1} \int \frac{d^D\ell_1\cdots d^D\ell_4}{(2\pi)^{4D}} \\
&\quad \times \prod_{i=1}^4 \frac{-i\,e^{-i\ell_i\cdot X_i}}{\ell_i^2+i0}
   \times \frac{-i}{(\ell_1+\ell_2)^2+i0} \\
&\quad \times z\,\big(v_3\!\cdot\!(\ell_1-\ell_2)\big)\,\big(v_4\!\cdot\!(\ell_3-\ell_4)\big).
\end{aligned}}
$$

---

## 5. 标量点乘 → 偏导算子

### 5a. 核心恒等式

$V_S$ 中的 $v\!\cdot\!\ell$ 因子作用于相位 $e^{-i\ell\cdot X}$：

$$
(v\!\cdot\!\ell)\,e^{-i\ell\cdot X} = i\,v\!\cdot\!\partial_X\,\bigl[e^{-i\ell\cdot X}\bigr].
$$

因此 $v\!\cdot\!\ell$ 可提到积分外，变为作用在 FT 结果上的偏导。

### 5b. 逐项替换

$$
v_3\!\cdot\!(\ell_1-\ell_2)\,e^{-i\ell_1\cdot X_1}e^{-i\ell_2\cdot X_2}
= i\,v_3\!\cdot\!(\partial_{X_1}-\partial_{X_2})\,e^{-i\ell_1\cdot X_1}e^{-i\ell_2\cdot X_2},
$$

$$
v_4\!\cdot\!(\ell_3-\ell_4)\,e^{-i\ell_3\cdot X_3}e^{-i\ell_4\cdot X_4}
= i\,v_4\!\cdot\!(\partial_{X_3}-\partial_{X_4})\,e^{-i\ell_3\cdot X_3}e^{-i\ell_4\cdot X_4}.
$$

### 5c. 定义 $\mathcal{O}_S$ 算子

两个替换的乘积给出 $i^2 = -1$，与 $z$ 前的符号合并：

$$
\boxed{\mathcal{O}_S = -z\;v_3\!\cdot\!(\partial_{X_1}-\partial_{X_2})\;v_4\!\cdot\!(\partial_{X_3}-\partial_{X_4}).}
$$

### 5d. 替换后的 $I_S$

$\ell$ 从分子中消失，只留在传播子分母中。$\mathcal{O}_S$ 提到 $\ell$ 积分外，作用在"标量 FT"上：

$$
\boxed{
\begin{aligned}
I_S &= -\frac{1}{4}C_A^2 C_F\;g^6
   \int_0^1 ds_{3,2}\!\int_0^{s_{3,2}}\! ds_{3,1}\int_0^1 ds_{4,2}\!\int_0^{s_{4,2}}\! ds_{4,1} \\
&\quad \times \mathcal{O}_S\!\left[
   \int \frac{d^D\ell_1\cdots d^D\ell_4}{(2\pi)^{4D}}
   \prod_{i=1}^4 \frac{-i\,e^{-i\ell_i\cdot X_i}}{\ell_i^2+i0}
   \times \frac{-i}{(\ell_1+\ell_2)^2+i0}
\right].
\end{aligned}}
$$

---

## 6. Fourier 变换：引入 $z_A, z_B$

### 6a. 问题：$1/(\ell_1+\ell_2)^2$ 耦合 $\ell_1, \ell_2$

方括号中的 5 条传播子：$1/\ell_1^2, 1/\ell_2^2, 1/\ell_3^2, 1/\ell_4^2, 1/(\ell_1+\ell_2)^2$。前 4 条各含独立 $\ell_i$，但第 5 条 $1/(\ell_1+\ell_2)^2$ 耦合 $\ell_1$ 和 $\ell_2$，无法逐个 FT。

### 6b. 解耦：重新引入 $\ell_5$ 和 $\delta$ 函数

$1/(\ell_1+\ell_2)^2$ 来自内部传播子 $1/\ell_5^2$（$\ell_5 = \ell_1+\ell_2$）。重新引入 $\ell_5$ 和两个 $\delta$ 函数解耦：

$$
\frac{1}{(\ell_1+\ell_2)^2} = \int \frac{d^D\ell_5}{(2\pi)^D}\;\frac{(2\pi)^D\delta^D(\ell_1+\ell_2-\ell_5)}{\ell_5^2}.
$$

TGV B 的动量守恒给出第二个 $\delta$：

$$
(2\pi)^D\delta^D(\ell_3+\ell_4+\ell_5).
$$

> **注**：TGV B 的 $\delta(\ell_3+\ell_4+\ell_5)$ 在 $V_{2\text{TGV}}$ 通用解中已被隐式使用（$\ell_5 = \ell_3+\ell_4$ 用于 Lorentz 缩并），但 FT 时需要显式恢复。

### 6c. $\delta$ 函数的 Fourier 表示 → $z_A, z_B$

$$
\delta^D(\ell_1+\ell_2-\ell_5) = \int \frac{d^Dz_A}{(2\pi)^D}\;e^{i(\ell_1+\ell_2-\ell_5)\cdot z_A},
$$

$$
\delta^D(\ell_3+\ell_4+\ell_5) = \int \frac{d^Dz_B}{(2\pi)^D}\;e^{i(\ell_3+\ell_4+\ell_5)\cdot z_B}.
$$

### 6d. 指数重组

$z_A, z_B$ 将每条 $\ell$ 的相位重新组合：

$$
\begin{aligned}
\ell_1:&\quad e^{-i\ell_1\cdot X_1}\,e^{i\ell_1\cdot z_A} = e^{-i\ell_1\cdot(X_1-z_A)}, \\
\ell_2:&\quad e^{-i\ell_2\cdot X_2}\,e^{i\ell_2\cdot z_A} = e^{-i\ell_2\cdot(X_2-z_A)}, \\
\ell_5:&\quad e^{-i\ell_5\cdot z_A}\,e^{i\ell_5\cdot z_B} = e^{-i\ell_5\cdot(z_A-z_B)}, \\
\ell_3:&\quad e^{-i\ell_3\cdot X_3}\,e^{i\ell_3\cdot z_B} = e^{-i\ell_3\cdot(X_3-z_B)}, \\
\ell_4:&\quad e^{-i\ell_4\cdot X_4}\,e^{i\ell_4\cdot z_B} = e^{-i\ell_4\cdot(X_4-z_B)}.
\end{aligned}
$$

5 条 $\ell$ 完全解耦——每条含一个传播子 $1/\ell_i^2$ 和一个相位 $e^{-i\ell_i\cdot x_i}$。

### 6e. 逐个 FT

标准 Fourier 变换：

$$
\int \frac{d^D\ell}{(2\pi)^D}\,\frac{-i}{\ell^2+i0}\,e^{-i\ell\cdot x}
= -\frac{\Gamma(a)}{4\pi^{D/2}}\,\frac{1}{(-x^2+i0)^a},
\qquad a = \tfrac{D}{2}-1 = 1-\varepsilon.
$$

5 个积分各给出一条位置空间传播子：

$$
\boxed{
\begin{aligned}
\Delta_1 &= -(X_1-z_A)^2+i0, &\quad
\Delta_2 &= -(X_2-z_A)^2+i0, \\
\Delta_5 &= -(z_A-z_B)^2+i0, &\quad
\Delta_3 &= -(X_3-z_B)^2+i0, &\quad
\Delta_4 &= -(X_4-z_B)^2+i0.
\end{aligned}}
$$

5-点拓扑：

```
X_1 ──Δ₁── z_A ──Δ₅── z_B ──Δ₃── X_3
              │                    │
             Δ₂                   Δ₄
              │                    │
             X_2                  X_4
```

FT prefactor：$\left(-\dfrac{\Gamma(a)}{4\pi^{D/2}}\right)^5$。

### 6f. FT 后的 $I_S$

$$
\boxed{
\begin{aligned}
I_S &= -\frac{1}{4}C_A^2 C_F\;g^6
   \left(-\frac{\Gamma(a)}{4\pi^{D/2}}\right)^{\!5} \\
&\quad \times \int_0^1 ds_{3,2}\!\int_0^{s_{3,2}}\! ds_{3,1}\int_0^1 ds_{4,2}\!\int_0^{s_{4,2}}\! ds_{4,1}
   \int d^Dz_A\,d^Dz_B \\
&\quad \times \mathcal{O}_S\!\left[
   \frac{1}{\Delta_1^a\,\Delta_2^a\,\Delta_5^a\,\Delta_3^a\,\Delta_4^a}
\right].
\end{aligned}}
$$

> $\mathcal{O}_S$ 的偏导 $\partial_{X_i}$ 作用于对应传播子 $\Delta_i^{-a}$，将指数从 $a$ 提升到 $a+1 = 2-\varepsilon$。

### 6g. 目标：5-点 $z_A, z_B$ 积分的初始形式

类比 I40 的 3-点 $z$ 积分 $\displaystyle\int d^Dz\,\prod_{i=1}^3\frac{1}{(-(X_i-z)^2+i0)^{D/2-1}}$，我们的 5-点标量形式为：

$$
\boxed{
\int d^Dz_A\,d^Dz_B\;
\frac{1}{(-(X_1-z_A)^2+i0)^a\,(-(X_2-z_A)^2+i0)^a\,(-(z_A-z_B)^2+i0)^a\,(-(X_3-z_B)^2+i0)^a\,(-(X_4-z_B)^2+i0)^a}
}
$$

其中 $a = D/2-1 = 1-\varepsilon$。

> **操作顺序**（与 I40 一致）：先 Feynman 参数化（将 5 条传播子合并为单一分母 $B$），再 $\mathcal{O}_S$ 作用于 $B^n$。$\mathcal{O}_S$ 的 $\partial_{X_i}$ 通过链式法则 $\partial_{X_i} B^n = n B^{n-1} \partial_{X_i} B$ 产生 Pnum，不改变 $\Gamma$ 结构。

---

## 7. Symanzik 参数化：5-点 $z_A, z_B$ 积分恒等式

### 7a. 一步 Feynman 参数化

对 5 条传播子做一次 Feynman 参数化（引入 $\alpha_1,...,\alpha_5$，$\delta(1\!-\!\sum\alpha)$）：

$$
\frac{1}{\Delta_1^a\cdots\Delta_5^a} = \frac{\Gamma(5a)}{\Gamma(a)^5}\int d\alpha\;\frac{\delta(1\!-\!\textstyle\sum\alpha)\,\prod\alpha_i^{a-1}}{Q^{5a}}
$$

### 7b. 矩阵形式

$Q = \sum\alpha_i\Delta_i$ 是 $z_A, z_B$ 的二次型：

$$
Q = \vec{z}^T(\Sigma\otimes I_D)\vec{z} - 2\vec{J}\cdot\vec{z} + C
$$

$$
\Sigma = \begin{pmatrix}\alpha_1\!+\!\alpha_2\!+\!\alpha_5 & -\alpha_5 \\ -\alpha_5 & \alpha_3\!+\!\alpha_4\!+\!\alpha_5\end{pmatrix},\quad
\vec{J} = \begin{pmatrix}\alpha_1 X_1\!+\!\alpha_2 X_2 \\ \alpha_3 X_3\!+\!\alpha_4 X_4\end{pmatrix},\quad
C = \sum_{i=1}^4\alpha_i X_i^2
$$

### 7c. Symanzik 多项式

**第一 Symanzik**（$\mathcal{U} = \det\Sigma$）：

$$
\boxed{\mathcal{U} = (\alpha_1\!+\!\alpha_2\!+\!\alpha_5)(\alpha_3\!+\!\alpha_4\!+\!\alpha_5) - \alpha_5^2 = \Sigma_A\Sigma_B - \alpha_5^2}
$$

**Feynman 分母**（$B_5 = -\mathcal{F}$，Python 验证所有 $X_i^2$ 消去）：

$$
\boxed{
B_5 = \alpha_1\alpha_2\Sigma_B\,\Delta_{12} + \alpha_3\alpha_4\Sigma_A\,\Delta_{34}
    + \alpha_5\!\!\sum_{\substack{i\in\{1,2\}\\j\in\{3,4\}}}\!\!\alpha_i\alpha_j\,\Delta_{ij}
}
$$

> **结构**：同组距离 $\Delta_{12},\Delta_{34}$ 系数 $\alpha_i\alpha_j\times\Sigma_{\text{对侧}}$；跨组距离系数 $\alpha_i\alpha_j\alpha_5$（桥接因子）。

### 7d. 高斯积分

配平方后 $2D$ 维高斯积分（$i^2=-1$，$\det(\Sigma\otimes I_D)=\mathcal{U}^D$）：

$$
\int d^{2D}z\;Q^{-5a} = \frac{-\pi^D\,\Gamma(5a\!-\!D)}{\Gamma(5a)}\;\frac{1}{\mathcal{U}^{1+2\varepsilon}\,B_5^{1-3\varepsilon}}
$$

### 7e. 5-点恒等式

$$
\boxed{
\int d^Dz_A\,d^Dz_B\;\prod_{i=1}^5\frac{1}{\Delta_i^a}
= \frac{-\pi^D\,\Gamma(1\!-\!3\varepsilon)}{\Gamma(1\!-\!\varepsilon)^5}
\int d\alpha\;\frac{\delta(1\!-\!\textstyle\sum\alpha)\,\prod\alpha_i^{a-1}}
{\mathcal{U}^{1+2\varepsilon}\;B_5^{1-3\varepsilon}}
}
$$

> **与 3-点对比**：3-点 $\mathcal{U}=1$（1-圈平凡），$B=\sum\alpha_i\alpha_j\Delta_{ij}$（3 项）。5-点 $\mathcal{U}\neq 1$（2-圈非平凡），$B_5$ 有 6 项。详见 `symanzik-parametrization.md`。

### 7f. 代入 $I_S$ 后

$$
\boxed{
\begin{aligned}
I_S &= -\frac{1}{4}C_A^2 C_F\;g^6
   \left(-\frac{\Gamma(a)}{4\pi^{D/2}}\right)^{\!5}
   \times \frac{-\pi^D\,\Gamma(1\!-\!3\varepsilon)}{\Gamma(1\!-\!\varepsilon)^5}
   \int_0^1 ds_{3,2}\!\cdots\! ds_{4,1}
   \int d\alpha\;\frac{\delta(1\!-\!\textstyle\sum\alpha)\,\prod\alpha_i^{a-1}}{1} \\
&\quad \times \mathcal{O}_S\!\left[
   \frac{1}{\mathcal{U}^{1+2\varepsilon}\;B_5^{1-3\varepsilon}}
\right].
\end{aligned}}
$$

FT prefactor 与 Symanzik prefactor 合并：$\Gamma(a)^5/\Gamma(1\!-\!\varepsilon)^5 = 1$，$(-1)^5\times(-1) = 1$，$4^5=1024$：

$$
\boxed{
I_S = -\frac{1}{4}\,\frac{C_A^2 C_F\;g^6\;\Gamma(1\!-\!3\varepsilon)}{1024\,\pi^{6-3\varepsilon}}
   \int_0^1 ds_{3,2}\!\cdots\! ds_{4,1}
   \int d\alpha\;\delta(1\!-\!\textstyle\sum\alpha)\,\prod\alpha_i^{a-1}
   \;\mathcal{O}_S\!\left[
   \frac{1}{\mathcal{U}^{1+2\varepsilon}\;B_5^{1-3\varepsilon}}
\right].
}
$$

---

## 8. $\mathcal{O}_S$ 作用于 Symanzik 分母

### 8a. 链式法则（二次型）

$\mathcal{O}_S = -z\,v_3\!\cdot\!(\partial_{X_1}\!-\!\partial_{X_2})\,v_4\!\cdot\!(\partial_{X_3}\!-\!\partial_{X_4})$ 作用于 $B_5^{-(1-3\varepsilon)}$（$\mathcal{U}$ 不含 $X$，不变）。

记 $n = -(1\!-\!3\varepsilon)$，分两步：

**Step 1**：$v_4\!\cdot\!(\partial_{X_3}\!-\!\partial_{X_4})\,B_5^n = n\,B_5^{n-1}\,G$，其中 $G = v_4\!\cdot\!(\partial_{X_3}\!-\!\partial_{X_4})\,B_5$。

**Step 2**：$v_3\!\cdot\!(\partial_{X_1}\!-\!\partial_{X_2})\,[n\,B_5^{n-1}\,G] = n\bigl[(n\!-\!1)B_5^{n-2}FG + B_5^{n-1}H\bigr]$，其中 $F = v_3\!\cdot\!(\partial_{X_1}\!-\!\partial_{X_2})\,B_5$，$H = v_3\!\cdot\!(\partial_{X_1}\!-\!\partial_{X_2})\,G$。

### 8b. 运动学恒等式（$v_{\text{ref}}=1$）

$$
\begin{aligned}
v_3\!\cdot\!X_1 &= v_3\!\cdot\!X_2 = -z, &\quad v_4\!\cdot\!X_3 &= v_4\!\cdot\!X_4 = 0, \\
v_3\!\cdot\!X_3 &= z(s_{4,1}\!-\!1), &\quad v_4\!\cdot\!X_1 &= z(s_{3,1}\!-\!1), \\
v_3\!\cdot\!X_4 &= z(s_{4,2}\!-\!1), &\quad v_4\!\cdot\!X_2 &= z(s_{3,2}\!-\!1).
\end{aligned}
$$

### 8c. $F$ 和 $G$ 计算（Python 验证）

$$
\boxed{F = v_3\!\cdot\!(\partial_{X_1}\!-\!\partial_{X_2})\,B_5 = 2\alpha_5\,z\,(\alpha_1\!-\!\alpha_2)\,(\alpha_3 s_{4,1}+\alpha_4 s_{4,2})}
$$

$$
\boxed{G = v_4\!\cdot\!(\partial_{X_3}\!-\!\partial_{X_4})\,B_5 = 2\alpha_5\,z\,(\alpha_4\!-\!\alpha_3)\,(\alpha_1(1\!-\!s_{3,1})+\alpha_2(1\!-\!s_{3,2}))}
$$

$F \propto z$，$G \propto z$ — 两者均含 $z$ 因子。

### 8d. 混合导数 $H \neq 0$

$$
\boxed{H = v_3\!\cdot\!(\partial_{X_1}\!-\!\partial_{X_2})\,G = 2\alpha_5\,z\,(\alpha_3\!-\!\alpha_4)\,(\alpha_1\!-\!\alpha_2) \neq 0}
$$

**推导**：$G$ 含 $v_4\!\cdot\!X_1$（D 维标量，经运动学恒等式化简为 $z(s_{3,1}\!-\!1)$）。$v_3\!\cdot\!\partial_{X_1}$ 对 $v_4\!\cdot\!X_1$ 求导给出 $v_3\!\cdot\!v_4 = z \neq 0$——**即使 $v_3^2=0$**。

> **关键**：$v_3\!\cdot\!\partial_{X_1}$ 是 D 维方向导数，**不等于** $\partial/\partial s_{3,1}$（Wilson 线参数导数，给出 $v_3^2=0$）。因为 $G$ 依赖 $v_4\!\cdot\!X_1$（$v_4$ 方向投影），沿 $v_3$ 方向移动 $X_1$ 会改变 $v_4\!\cdot\!X_1$（变化量 $v_3\!\cdot\!v_4 = z$）。

### 8e. $\mathcal{O}_S$ 作用结果：两项

$H \neq 0$ 使链式法则产生两项：

$$
\boxed{
\mathcal{O}_S\!\left[B_5^n\right] = -z\bigl[n(n\!-\!1)\,B_5^{n-2}\,FG + n\,B_5^{n-1}\,H\bigr]
}
$$

代入 $n = -(1\!-\!3\varepsilon)$，$n(n\!-\!1) = (1\!-\!3\varepsilon)(2\!-\!3\varepsilon)$：

**Term 1**（链式法则第一项，$FG$）：

$$
-z\,(1\!-\!3\varepsilon)(2\!-\!3\varepsilon)\,B_5^{-(3-3\varepsilon)}\;\times\;4\alpha_5^2 z^2\,\text{Pnum}_0
$$

**Term 2**（混合导数项，$H$）：

$$
+z\,(1\!-\!3\varepsilon)\,B_5^{-(2-3\varepsilon)}\;\times\;2\alpha_5 z\,(\alpha_3\!-\!\alpha_4)(\alpha_1\!-\!\alpha_2)
$$

其中：

$$
\text{Pnum}_0 = (\alpha_1\!-\!\alpha_2)(\alpha_4\!-\!\alpha_3)\,(\alpha_3 s_{4,1}+\alpha_4 s_{4,2})\,(\alpha_1(1\!-\!s_{3,1})+\alpha_2(1\!-\!s_{3,2}))
$$

### 8f. $B_5$ 的 $A_{\text{ref}}$ 因子

S 道中 $B_5 = A_{\text{ref}}\,z\,\hat{B}_5$（其中 $A_{\text{ref}} = -2v_{\text{ref}}$，$\hat{B}_5$ 不含 $z$ 和 $A_{\text{ref}}$）：

$$
B_5^{-(3-3\varepsilon)} = A_{\text{ref}}^{-3}\,A_{\text{ref}}^{3\varepsilon}\,z^{-(3-3\varepsilon)}\,\hat{B}_5^{-(3-3\varepsilon)},\quad
A_{\text{ref}}^{-3} = \frac{1}{(-2)^3} = -\frac{1}{8}
$$

$$
B_5^{-(2-3\varepsilon)} = A_{\text{ref}}^{-2}\,A_{\text{ref}}^{3\varepsilon}\,z^{-(2-3\varepsilon)}\,\hat{B}_5^{-(2-3\varepsilon)},\quad
A_{\text{ref}}^{-2} = \frac{1}{(-2)^2} = \frac{1}{4}
$$

### 8g. 代入后的 $I_S$

合并两项，提取 $z^{3\varepsilon}$ 和 $A_{\text{ref}}^{3\varepsilon}$（归入 ⑥ v14\_scale），剩余数值因子 $A_{\text{ref}}^{-3}=-1/8$ 和 $A_{\text{ref}}^{-2}=1/4$：

$$
\boxed{
\begin{aligned}
I_S &= -\frac{1}{4}C_A^2 C_F\,g^6\,\frac{\Gamma(1\!-\!3\varepsilon)}{1024\,\pi^{6-3\varepsilon}}
   \times z^{3\varepsilon}\,(-2v_{\text{ref}}\!+\!i0)^{3\varepsilon} \\
&\quad \times \frac{1}{2}(1\!-\!3\varepsilon)\int ds\,d\alpha\;\delta(1\!-\!\textstyle\sum\alpha)\,\prod\alpha_i^{a-1}\;\frac{1}{\mathcal{U}^{1+2\varepsilon}} \\
&\quad \times \biggl[\,(2\!-\!3\varepsilon)\,\alpha_5^2\,\frac{\text{Pnum}_0}{\hat{B}_5^{3-3\varepsilon}}
   \;+\;\alpha_5\,\frac{(\alpha_3\!-\!\alpha_4)(\alpha_1\!-\!\alpha_2)}{\hat{B}_5^{2-3\varepsilon}}\,\biggr]
\end{aligned}}
$$

> **⑤ derivative** $= \frac{1}{2}(1\!-\!3\varepsilon)$（不再是一个简单因子——含两项结构）。$\varepsilon\to 0$：$1/2$。
>
> **z 因子化**：两项均给出 $z^{3\varepsilon}$（Term 1: $z^{1+2-3}=z^0$, $+3\varepsilon$；Term 2: $z^{1+1-2}=z^0$, $+3\varepsilon$）。

---

## 9. 超立方变换

$$
x_1=s_{3,2},\;x_2=s_{3,1}/s_{3,2},\;x_3=s_{4,2},\;x_4=s_{4,1}/s_{4,2},\quad
\text{Jacobian}=x_1 x_3.
$$

代入后（Python 验证）：

$$
\boxed{
\begin{aligned}
\text{Pnum}_0(x) &= (\alpha_1\!-\!\alpha_2)(\alpha_4\!-\!\alpha_3)\;x_3\,(\alpha_3 x_4\!+\!\alpha_4)\;\bigl[\alpha_1(1\!-\!x_1 x_2)+\alpha_2(1\!-\!x_1)\bigr], \\[4pt]
\hat{B}_5(x) &= \alpha_5\;x_3\,(\alpha_3 x_4\!+\!\alpha_4)\;\bigl[\alpha_1(1\!-\!x_1 x_2)+\alpha_2(1\!-\!x_1)\bigr].
\end{aligned}}
$$

两者共享因子 $x_3(\alpha_3 x_4\!+\!\alpha_4)[\alpha_1(1\!-\!x_1 x_2)+\alpha_2(1\!-\!x_1)]$。

---

## 10. $C$ 因子分解

$$
\boxed{
C = \underbrace{-\frac{1}{4}C_A^2 C_F}_{\text{① color}}
   \times \underbrace{(4\pi\alpha_s)^3}_{\text{③ coupling}}
   \times \underbrace{\frac{\Gamma(1\!-\!3\varepsilon)}{1024\,\pi^{6-3\varepsilon}}}_{\text{④ ft\_base}}
   \times \underbrace{\frac{1}{2}(1\!-\!3\varepsilon)}_{\text{⑤ derivative}}
   \times \underbrace{(-2v_{\text{ref}}\!+\!i0)^{3\varepsilon}}_{\text{⑥ v14\_scale}}.
}
$$

**$C$ 为实数** $\checkmark$：$i^4=1$（Wilson 线），$i^2=-1$（$2D$ 高斯），$i^2=-1$（$\ell\to i\partial_X$），$i^4\times i^2\times i^2=1$。

> **与 I40 对照**：
>
> | 因子 | I40（单 TGV） | S 道（双 TGV） | 差异来源 |
> |------|:---:|:---:|---------|
> | ① color | $\frac{i}{2}C_A C_F$ | $-\frac{1}{4}C_A^2 C_F$ | 双 TGV 多一个 $C_A$ |
> | ③ coupling | $(4\pi\alpha_s)^2$ | $(4\pi\alpha_s)^3$ | 3-loop vs 2-loop |
> | ④ ft_base | $\Gamma(1\!-\!2\varepsilon)/(256\pi^{4-2\varepsilon})$ | $\Gamma(1\!-\!3\varepsilon)/(1024\pi^{6-3\varepsilon})$ | 5-pt vs 3-pt |
> | ⑤ derivative | $2(1\!-\!2\varepsilon)$ | $\frac{1}{2}(1\!-\!3\varepsilon)$ | 二次型（两项） |
> | ⑥ v14_scale | $(-2v_{14})^{2\varepsilon}$ | $(-2v_{\text{ref}})^{3\varepsilon}$ | 3-loop vs 2-loop |

---

## 11. $J$ 被积函数与汇总

$$
\boxed{
\begin{aligned}
I_S &= C \times J, \\[8pt]
C &= -\frac{1}{4}C_A^2 C_F\,(4\pi\alpha_s)^3\,
   \frac{\Gamma(1\!-\!3\varepsilon)}{1024\,\pi^{6-3\varepsilon}}\,
   \frac{1}{2}(1\!-\!3\varepsilon)\,
   (-2v_{\text{ref}}\!+\!i0)^{3\varepsilon}, \\[12pt]
J &= z^{3\varepsilon}\int_0^1 d^4x\;x_1 x_3
   \int d\alpha\;\delta(1\!-\!\textstyle\sum\alpha)\,\prod\alpha_i^{a-1}\;\frac{1}{\mathcal{U}^{1+2\varepsilon}} \\
  &\quad \times \biggl[\,(2\!-\!3\varepsilon)\,\alpha_5^2\,\frac{\text{Pnum}_0(x)}{\hat{B}_5(x)^{3-3\varepsilon}}
   \;+\;\alpha_5\,\frac{(\alpha_3\!-\!\alpha_4)(\alpha_1\!-\!\alpha_2)}{\hat{B}_5(x)^{2-3\varepsilon}}\,\biggr].
\end{aligned}}
$$

其中 $\mathcal{U} = (\alpha_1\!+\!\alpha_2\!+\!\alpha_5)(\alpha_3\!+\!\alpha_4\!+\!\alpha_5)-\alpha_5^2$，$a=1\!-\!\varepsilon$。

$J$ 含**两项**（对应链式法则的 $FG$ 项和 $H$ 项），$z^{3\varepsilon}$ 全局提取，kernel 不含 $z$。
