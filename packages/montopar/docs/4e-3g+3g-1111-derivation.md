# 4e-3g+3g [1,1,1,1] 道：完整 $C \times J$ 推导

> 风格对标 `4e-3g+3g-S-derivation.md`——将 $V_{2\text{TGV}}$ 当作"单位顶点"，给出完整 $C, J$
> 轨道：4e-3g+3g-25 — 边分布 $(1,2,3,4)$，配对 $(12)(34)$
> 四边互异，双非退化——4e-3g+3g 最一般情形
> $D=4-2\varepsilon$，$v_i^2=0$，$v_{\text{ref}}=v_{14}$，$A_{\text{ref}}\equiv -2v_{\text{ref}}+i0$
> $z = v_{12}/v_{\text{ref}} = v_{34}/v_{\text{ref}}$，$v_{13}=v_{24}=-(1+z)v_{\text{ref}}$，$v_{14}=v_{23}=v_{\text{ref}}$
> 代数计算：`compute_pnum_1111.py`（Python + sympy）

---

## 1. 动量空间起点：$V_{2\text{TGV}}$ 通用解

双 TGV 有效顶点的通用解（FeynCalc `lorentz_construct02_general.wl`）：

$$
\boxed{
\begin{aligned}
V_{2\text{TGV}}(a,b,c,d;\ell_1,\ell_2,\ell_3,\ell_4)
&= f^{A_a A_b B}\,f^{A_c A_d B}\;
   \frac{-i}{(\ell_1+\ell_2)^2+i0} \\
&\quad \times \Bigl[
   \underbrace{(v_a\!\cdot\!v_b)(v_c\!\cdot\!v_d)\;(\ell_2-\ell_1)\!\cdot\!(\ell_4-\ell_3)}_{\text{G1}} \\
&\quad + \underbrace{(v_a\!\cdot\!v_b)\Bigl[
   2(v_c\!\cdot\!\ell_1)(v_d\!\cdot\!\ell_2) - 2(v_c\!\cdot\!\ell_2)(v_d\!\cdot\!\ell_1)
   + v_c\!\cdot\!(\ell_2-\ell_1)\,(v_d\!\cdot\!\ell_3)
   + v_d\!\cdot\!(\ell_1-\ell_2)\,(v_c\!\cdot\!\ell_4)\Bigr]}_{\text{G2}} \\
&\quad + \underbrace{\begin{aligned}
   &(v_a\!\cdot\!v_b)\;4\Bigl[(v_d\!\cdot\!\ell_1)(v_c\!\cdot\!\ell_2) - (v_c\!\cdot\!\ell_1)(v_d\!\cdot\!\ell_2)\Bigr] \\
   &+ (v_c\!\cdot\!v_d)\Bigl[(v_a\!\cdot\!\ell_1+2v_a\!\cdot\!\ell_2)(v_b\!\cdot\!\ell_3-v_b\!\cdot\!\ell_4) \\
   &\qquad\qquad\; - (2v_b\!\cdot\!\ell_1+v_b\!\cdot\!\ell_2)(v_a\!\cdot\!\ell_3-v_a\!\cdot\!\ell_4)\Bigr] \\
   &+ (v_a\!\cdot\!v_c)\;2(2v_b\!\cdot\!\ell_1+v_b\!\cdot\!\ell_2)(v_d\!\cdot\!\ell_1+v_d\!\cdot\!\ell_2) \\
   &- (v_a\!\cdot\!v_d)\;2(2v_b\!\cdot\!\ell_1+v_b\!\cdot\!\ell_2)(v_c\!\cdot\!\ell_1+v_c\!\cdot\!\ell_2) \\
   &- (v_b\!\cdot\!v_c)\;2(v_a\!\cdot\!\ell_1+2v_a\!\cdot\!\ell_2)(v_d\!\cdot\!\ell_1+v_d\!\cdot\!\ell_2) \\
   &+ (v_b\!\cdot\!v_d)\;2(v_a\!\cdot\!\ell_1+2v_a\!\cdot\!\ell_2)(v_c\!\cdot\!\ell_1+v_c\!\cdot\!\ell_2)
   \end{aligned}}_{\text{G3}} \\
&\quad + \underbrace{(v_a\!\cdot\!v_c)\,(2v_b\!\cdot\!\ell_1+v_b\!\cdot\!\ell_2)\;v_d\!\cdot\!(\ell_3-\ell_1-\ell_2)}_{\text{G4a}}
   + \underbrace{(v_a\!\cdot\!v_d)\,(2v_b\!\cdot\!\ell_1+v_b\!\cdot\!\ell_2)\;v_c\!\cdot\!(\ell_1+\ell_2-\ell_4)}_{\text{G4b}} \\
&\quad - \underbrace{(v_b\!\cdot\!v_c)\,(v_a\!\cdot\!\ell_1+2v_a\!\cdot\!\ell_2)\;v_d\!\cdot\!(\ell_3-\ell_1-\ell_2)}_{\text{G4c}}
   - \underbrace{(v_b\!\cdot\!v_d)\,(v_a\!\cdot\!\ell_1+2v_a\!\cdot\!\ell_2)\;v_c\!\cdot\!(\ell_1+\ell_2-\ell_4)}_{\text{G4d}}
   \Bigr].
\end{aligned}}
$$

完整积分 = 外部传播子 $\times$ Wilson 线因子 $\times$ $V_{2\text{TGV}}$：

$$
\boxed{
\begin{aligned}
I_{(1,2,3,4)} &= \int_0^1 ds_1 ds_2 ds_3 ds_4 \int \frac{d^D\ell_1\cdots d^D\ell_4}{(2\pi)^{4D}} \\
&\quad \times \prod_{i=1}^4 (ig\,v_{e_i}^{\mu_i} t^{A_i})\,
   \frac{-i\,g_{\mu_i\nu_i}\,\delta^{A_i B_i}}{\ell_i^2+i0}\,
   e^{-i\ell_i\cdot X_i} \\
&\quad \times V_{2\text{TGV}}(1,2,3,4;\ell_1,\ell_2,\ell_3,\ell_4),
\end{aligned}}
$$

其中 $e_i = i$（腿 $i$ 附着在边 $i$ 上），$X_i = V_i + s_i v_i$，$V_1=0, V_2=v_1, V_3=v_1+v_2, V_4=v_1+v_2+v_3$。

> **与 S 道对比**：S 道 $(3,3,4,4)$ 中 $v_{ab}=v_{33}=0$ 使 G1/G2/G3 全部消失，仅 G4a-d 存活（4 项合并为 1 项）。[1,1,1,1] 中 $v_{12}=z\neq0$ 且 $v_{34}=z\neq0$，**所有 4 组全部存活**——共 15 个独立 Lorentz 项（G1:1 + G2:4 + G3 的 $v_{ab}$ 子项:1 + G3 的 $v_{cd}$ 子项:2 + G3 的交叉项:4 + G4a-d:4 = 15，合并 G2 和 G3 的 $v_{ab}$ 子项后为 15 项）。

---

## 2. 颜色因子

配对 $(12)(34)$，颜色 $f^{A_1 A_2 B}f^{A_3 A_4 B}$ 缩并 Wilson 线生成元 $t^{A_1}t^{A_2}t^{A_3}t^{A_4}$：

$$
\boxed{\text{① color} = t^{A_1}t^{A_2}t^{A_3}t^{A_4}\;f^{A_1 A_2 B}f^{A_3 A_4 B}.}
$$

> FeynCalc `SUNSimplify` 待求值。由于四边互异（$A_1,A_2,A_3,A_4$ 为不同生成元），颜色结构与 S 道 $-\frac14 C_A^2 C_F$ 和 T 道 $+\frac14 C_A^2 C_F$ 不同。

提取颜色后：

$$
\boxed{
\begin{aligned}
I_{(1,2,3,4)} &= \text{① color}
   \int_0^1 ds_1 ds_2 ds_3 ds_4 \int \frac{d^D\ell_1\cdots d^D\ell_4}{(2\pi)^{4D}} \\
&\quad \times \prod_{i=1}^4 \frac{-i\,v_{e_i}^{\nu_i}}{\ell_i^2+i0}\,
   e^{-i\ell_i\cdot X_i}
   \times \frac{-i}{(\ell_1+\ell_2)^2+i0} \\
&\quad \times \mathcal{L}(1,2,3,4;\ell_1,\ell_2,\ell_3,\ell_4),
\end{aligned}}
$$

其中 $\mathcal{L}$ 是 $V_{2\text{TGV}}$ 的 Lorentz 部分（§1 中 15 项的 $[\cdots]$ 括号）。

---

## 3. Lorentz 缩并：$\ell$ 结构分类

### 3a. 传播子缩并

$g_{\mu_i\nu_i}$ 缩并 $v_{e_i}^{\mu_i}$，将 Wilson 线向量传递到内部。$v_{e_i}^{\nu_i}$ 即 $v_1^{\nu_1},v_2^{\nu_2},v_3^{\nu_3},v_4^{\nu_4}$。

### 3b. 代入 $v_{ij}$

$$
\boxed{
\begin{aligned}
v_a=v_1,\; v_b=v_2,\; v_c=v_3,\; v_d=v_4, \\[4pt]
v_{ab}=v_{12}=z,\quad v_{cd}=v_{34}=z,\quad
v_{ac}=v_{13}=-(1+z),\quad v_{ad}=v_{14}=1,\quad
v_{bc}=v_{23}=1,\quad v_{bd}=v_{24}=-(1+z).
\end{aligned}}
$$

**所有 $v_{ij}\neq 0$**——没有任何 Lorentz 组因退化消失。

### 3c. $\ell$ 结构三类

15 项按 $\ell$ 结构分为三类（与通用推导 §4c 一致）：

| 类型 | $\ell$ 结构 | 项数 | 前因子 | 位置空间算子 |
|------|-------------|:---:|--------|-------------|
| I | $\ell_i\!\cdot\!\ell_j$ | 1（G1） | $v_{ab}v_{cd}=z^2$ | $-\partial_{X_i}\!\cdot\!\partial_{X_j}$ |
| II | $(v_i\!\cdot\!\ell)(v_j\!\cdot\!\ell)$ | 6 | $v_{ab}=z$ 或 $v_{cd}=z$ | $-(v_i\!\cdot\!\partial_X)(v_j\!\cdot\!\partial_X)$ |
| III | $(v_i\!\cdot\!\ell)(v_j\!\cdot\!\ell)$ | 8 | $v_{ac},v_{ad},v_{bc},v_{bd}$ | $-(v_i\!\cdot\!\partial_X)(v_j\!\cdot\!\partial_X)$ |

> **与 S 道对比**：S 道仅 Type III 存活（1 项）。[1,1,1,1] 三类全活（15 项）。

---

## 4. 顶点 $i$ 因子

- 4 个 Wilson 线顶点：$(ig)^4 = g^4$（$i^4=1$）
- 2 个 TGV 顶点：$g^2$
- 合计 $g^6 = (4\pi\alpha_s)^3$

5 条传播子分子 $(-i)^5$ 暂留积分内（归入 §6 的 FT prefactor）。提取 $g^6$：

$$
\boxed{
\begin{aligned}
I_{(1,2,3,4)} &= \text{① color}\;g^6
   \int_0^1 ds_1 ds_2 ds_3 ds_4 \int \frac{d^D\ell_1\cdots d^D\ell_4}{(2\pi)^{4D}} \\
&\quad \times \prod_{i=1}^4 \frac{-i\,e^{-i\ell_i\cdot X_i}}{\ell_i^2+i0}
   \times \frac{-i}{(\ell_1+\ell_2)^2+i0} \\
&\quad \times \mathcal{L}(1,2,3,4;\ell_1,\ell_2,\ell_3,\ell_4).
\end{aligned}}
$$

---

## 5. 标量点乘 → 偏导算子

### 5a. 核心恒等式

$$
(v\!\cdot\!\ell)\,e^{-i\ell\cdot X} = i\,v\!\cdot\!\partial_X\bigl[e^{-i\ell\cdot X}\bigr].
$$

### 5b. 定义微分算子 $\mathcal{O}(1,2,3,4)$

将 $\mathcal{L}$ 的 15 项中 $\ell \to i\partial_X$，展开为 20 个子项（通用推导 §5c 的算子项表）。按通用推导 §2d 的编号：

**Type I**（$\partial_i\!\cdot\!\partial_j$，4 子项，$C_t = \pm z^2$）：

| $t$ | $C_t$ | $(i,j)$ |
|:---:|------|--------|
| I-1 | $-z^2$ | $(2,4)$ |
| I-2 | $+z^2$ | $(2,3)$ |
| I-3 | $+z^2$ | $(1,4)$ |
| I-4 | $-z^2$ | $(1,3)$ |

**Type II**（$v_{ab}=z$ 或 $v_{cd}=z$ 前因子，8 子项）：

| $t$ | $C_t$ | $v_p$ | $D^{(p)}$ | $v_q$ | $D^{(q)}$ |
|:---:|------|------|----------|------|----------|
| IIa-1 | $-2z$ | $v_3$ | $\partial_1$ | $v_4$ | $\partial_2$ |
| IIa-2 | $+2z$ | $v_3$ | $\partial_2$ | $v_4$ | $\partial_1$ |
| IIa-3 | $-z$ | $v_3$ | $\partial_2-\partial_1$ | $v_4$ | $\partial_3$ |
| IIa-4 | $-z$ | $v_4$ | $\partial_1-\partial_2$ | $v_3$ | $\partial_4$ |
| IIb-1 | $-4z$ | $v_4$ | $\partial_1$ | $v_3$ | $\partial_2$ |
| IIb-2 | $+4z$ | $v_3$ | $\partial_1$ | $v_4$ | $\partial_2$ |
| IIc-1 | $-z$ | $v_1$ | $\partial_1+2\partial_2$ | $v_2$ | $\partial_3-\partial_4$ |
| IIc-2 | $+z$ | $v_2$ | $2\partial_1+\partial_2$ | $v_1$ | $\partial_3-\partial_4$ |

**Type III**（$v_{ac},v_{ad},v_{bc},v_{bd}$ 前因子，8 子项）：

| $t$ | $C_t$ | $v_p$ | $D^{(p)}$ | $v_q$ | $D^{(q)}$ |
|:---:|------|------|----------|------|----------|
| IIIa-1 | $-2v_{13}$ | $v_2$ | $2\partial_1+\partial_2$ | $v_4$ | $\partial_1+\partial_2$ |
| IIIa-2 | $+2v_{14}$ | $v_2$ | $2\partial_1+\partial_2$ | $v_3$ | $\partial_1+\partial_2$ |
| IIIa-3 | $+2v_{23}$ | $v_1$ | $\partial_1+2\partial_2$ | $v_4$ | $\partial_1+\partial_2$ |
| IIIa-4 | $-2v_{24}$ | $v_1$ | $\partial_1+2\partial_2$ | $v_3$ | $\partial_1+\partial_2$ |
| IIIb-1 | $-v_{13}$ | $v_2$ | $2\partial_1+\partial_2$ | $v_4$ | $\partial_3-\partial_1-\partial_2$ |
| IIIb-2 | $-v_{14}$ | $v_2$ | $2\partial_1+\partial_2$ | $v_3$ | $\partial_1+\partial_2-\partial_4$ |
| IIIb-3 | $+v_{23}$ | $v_1$ | $\partial_1+2\partial_2$ | $v_4$ | $\partial_3-\partial_1-\partial_2$ |
| IIIb-4 | $+v_{24}$ | $v_1$ | $\partial_1+2\partial_2$ | $v_3$ | $\partial_1+\partial_2-\partial_4$ |

其中 $v_{13}=-(1+z)$，$v_{14}=1$，$v_{23}=1$，$v_{24}=-(1+z)$。

### 5c. 替换后的 $I$

$\mathcal{O}$ 提到 $\ell$ 积分外：

$$
\boxed{
\begin{aligned}
I_{(1,2,3,4)} &= \text{① color}\;g^6
   \int_0^1 ds_1 ds_2 ds_3 ds_4 \\
&\quad \times \mathcal{O}(1,2,3,4)\!\left[
   \int \frac{d^D\ell_1\cdots d^D\ell_4}{(2\pi)^{4D}}
   \prod_{i=1}^4 \frac{-i\,e^{-i\ell_i\cdot X_i}}{\ell_i^2+i0}
   \times \frac{-i}{(\ell_1+\ell_2)^2+i0}
\right].
\end{aligned}}
$$

---

## 6. Fourier 变换：引入 $z_A, z_B$

### 6a. 解耦 $1/(\ell_1+\ell_2)^2$

重新引入 $\ell_5$ 和 $\delta$ 函数（$\ell_5 = \ell_1+\ell_2$）：

$$
\frac{1}{(\ell_1+\ell_2)^2+i0} = \int \frac{d^D\ell_5}{(2\pi)^D}\;
\frac{(2\pi)^D\delta^D(\ell_1+\ell_2-\ell_5)}{\ell_5^2+i0}.
$$

TGV B 的动量守恒给出第二个 $\delta$：

$$
(2\pi)^D\delta^D(\ell_3+\ell_4+\ell_5).
$$

### 6b. $\delta$ 函数的 Fourier 表示

$$
\delta^D(\ell_1+\ell_2-\ell_5) = \int \frac{d^Dz_A}{(2\pi)^D}\;e^{i(\ell_1+\ell_2-\ell_5)\cdot z_A},
\qquad
\delta^D(\ell_3+\ell_4+\ell_5) = \int \frac{d^Dz_B}{(2\pi)^D}\;e^{i(\ell_3+\ell_4+\ell_5)\cdot z_B}.
$$

### 6c. 5-点拓扑

指数重组后，5 条 $\ell$ 完全解耦：

```
X_1 ──Δ₁── z_A ──Δ₅── z_B ──Δ₃── X_3
  │                    │
  Δ₂                   Δ₄
  │                    │
 X_2                  X_4
```

### 6d. 逐个 FT

标准标量 FT（$a = D/2-1 = 1-\varepsilon$）：

$$
\int\frac{d^D\ell}{(2\pi)^D}\,\frac{-i}{\ell^2+i0}\,e^{-i\ell\cdot x}
= -\frac{\Gamma(a)}{4\pi^{D/2}}\,\frac{1}{(-x^2+i0)^a}.
$$

FT prefactor：$\left(-\dfrac{\Gamma(a)}{4\pi^{D/2}}\right)^5$。

### 6e. FT 后的 $I$

$$
\boxed{
\begin{aligned}
I_{(1,2,3,4)} &= \text{① color}\;g^6
   \left(-\frac{\Gamma(a)}{4\pi^{D/2}}\right)^{\!5} \\
&\quad \times \int_0^1 ds_1 ds_2 ds_3 ds_4
   \int d^Dz_A\,d^Dz_B \\
&\quad \times \mathcal{O}(1,2,3,4)\!\left[
   \frac{1}{\Delta_1^a\,\Delta_2^a\,\Delta_5^a\,\Delta_3^a\,\Delta_4^a}
\right].
\end{aligned}}
$$

---

## 7. Symanzik 参数化：5-点 $z_A,z_B$ 积分恒等式

### 7a. 一步 Feynman 参数化

引入 $\alpha_1,\ldots,\alpha_5$，$\delta(1-\sum\alpha)$：

$$
\frac{1}{\Delta_1^a\Delta_2^a\Delta_5^a\Delta_3^a\Delta_4^a}
= \frac{\Gamma(5a)}{\Gamma(a)^5}\int d\alpha\;\frac{\delta(1-\sum\alpha)\,\prod\alpha_i^{a-1}}{Q^{5a}}.
$$

### 7b. 矩阵形式

$Q = \sum\alpha_i\Delta_i$ 是 $z_A,z_B$ 的二次型（同 S 道 §7b）：

$$
\Sigma = \begin{pmatrix}\alpha_1+\alpha_2+\alpha_5 & -\alpha_5 \\ -\alpha_5 & \alpha_3+\alpha_4+\alpha_5\end{pmatrix},\quad
\vec{J} = \begin{pmatrix}\alpha_1 X_1+\alpha_2 X_2 \\ \alpha_3 X_3+\alpha_4 X_4\end{pmatrix},\quad
C = \sum_{i=1}^4\alpha_i X_i^2.
$$

### 7c. $\mathcal{U}$ 和 $B_5$

$$
\boxed{\mathcal{U} = \det\Sigma = (\alpha_1+\alpha_2+\alpha_5)(\alpha_3+\alpha_4+\alpha_5) - \alpha_5^2.}
$$

$$
\boxed{
\begin{aligned}
B_5 &= \alpha_1\alpha_2(\alpha_3+\alpha_4+\alpha_5)\,\Delta_{12}
    + \alpha_3\alpha_4(\alpha_1+\alpha_2+\alpha_5)\,\Delta_{34} \\
    &\quad + \alpha_5\bigl(
        \alpha_1\alpha_3\Delta_{13}
      + \alpha_1\alpha_4\Delta_{14}
      + \alpha_2\alpha_3\Delta_{23}
      + \alpha_2\alpha_4\Delta_{24}
      \bigr).
\end{aligned}}
$$

### 7d. $\Delta_{ij}$ 显式（代码计算）

定义 $\Delta'_{ij} \equiv \Delta_{ij}/A_{\text{ref}}$。代入 $X_i$ 位置（§1），代码展开：

$$
\boxed{
\begin{aligned}
\Delta'_{12} &= s_2\,z\,(1-s_1), \\
\Delta'_{34} &= s_4\,z\,(1-s_3), \\[4pt]
\Delta'_{13} &= s_1 s_3 z + s_1 s_3 - s_1 z - s_3 z + z, \\
\Delta'_{14} &= s_1(1-s_4), \\
\Delta'_{23} &= s_3(1-s_2), \\
\Delta'_{24} &= s_2 s_4 z + s_2 s_4 - s_2 - s_4 + 1.
\end{aligned}}
$$

**全部 6 项非零**——这是 [1,1,1,1] 与 S 道（$\Delta_{12}=\Delta_{34}=0$）的本质区别。

> **几何结构**：$\Delta'_{12},\Delta'_{34}$ 是邻边型（$\propto z$，单因子 $s$）。$\Delta'_{13}$ 是对边型（$v_{13}=-(1+z)$，含 $z$ 和 $1$ 混合）。$\Delta'_{14},\Delta'_{23}$ 是邻边型（$v_{14}=v_{23}=1$，不含 $z$）。$\Delta'_{24}$ 是对边型（$v_{24}=-(1+z)$，含 $z$ 和 $1$）。

### 7e. 高斯积分

配平方后 $2D$ 维高斯积分：

$$
\int d^{2D}z\;Q^{-5a}
= \frac{-\pi^D\,\Gamma(1-3\varepsilon)}{\Gamma(5-5\varepsilon)}\;\frac{1}{\mathcal{U}^{1+2\varepsilon}\;B_5^{\,1-3\varepsilon}}.
$$

### 7f. 5-点恒等式

$$
\boxed{
\int d^Dz_A\,d^Dz_B\;\prod_{i=1}^5\frac{1}{\Delta_i^a}
= \frac{-\pi^D\,\Gamma(1-3\varepsilon)}{\Gamma(1-\varepsilon)^5}
\int d\alpha\;\frac{\delta(1-\sum\alpha)\,\prod\alpha_i^{-\varepsilon}}
{\mathcal{U}^{1+2\varepsilon}\;B_5^{\,1-3\varepsilon}}.
}
$$

### 7g. 代入 $I$

FT prefactor 与 Symanzik prefactor 合并（$\Gamma(a)^5/\Gamma(1-\varepsilon)^5=1$，$(-1)^5\times(-1)=1$，$4^5=1024$）：

$$
\boxed{
\begin{aligned}
I_{(1,2,3,4)} &= \text{① color}\;g^6\,
   \frac{\Gamma(1-3\varepsilon)}{1024\,\pi^{6-3\varepsilon}} \\
&\quad \times \int_0^1 ds_1 ds_2 ds_3 ds_4
   \int d\alpha\;\delta(1-\textstyle\sum\alpha)\,\prod\alpha_i^{-\varepsilon} \\
&\quad \times \mathcal{O}(1,2,3,4)\!\left[
   \frac{1}{\mathcal{U}^{1+2\varepsilon}\;B_5^{\,1-3\varepsilon}}
\right].
\end{aligned}}
$$

---

## 8. $\mathcal{O}$ 作用于 Symanzik 分母

### 8a. 二次型链式法则

$\mathcal{O}$ 的每一项为 $C_t \times (v_p\!\cdot\!D^{(p)})(v_q\!\cdot\!D^{(q)})$。作用在 $B_5^n$（$n=-(1-3\varepsilon)$）上：

$$
\boxed{
(v_p\!\cdot\!D^{(p)})(v_q\!\cdot\!D^{(q)})\,B_5^n
= n(n-1)B_5^{n-2}\,F^{(p)}G^{(q)} + nB_5^{n-1}\,H^{(pq)},
}
$$

其中 $n(n-1)=(1-3\varepsilon)(2-3\varepsilon)$，$n=-(1-3\varepsilon)$。

### 8b. Building block 公式

**一次导数**（通用推导 §8c）：

$$
F_i^{(p)} = -2\alpha_i\sum_{j\neq i}\alpha_j\,c_{ij}\;v_p\!\cdot\!(X_i-X_j),
\qquad
F^{(p)} = \sum_i d_i^{(p)} F_i^{(p)}.
$$

**混合导数**（通用推导 §8d）：

$$
H_{ij}^{(pq)} = \begin{cases}
2\alpha_i\alpha_j\,c_{ij}\,(v_p\!\cdot\!v_q), & i\neq j, \\
-2\alpha_i\,(v_p\!\cdot\!v_q)\sum_{k\neq i}\alpha_k\,c_{ik}, & i=j,
\end{cases}
\qquad
H^{(pq)} = \sum_{i,j} d_i^{(p)} d_j^{(q)} H_{ij}^{(pq)}.
$$

**系数矩阵 $c_{ij}$**（通用推导 §8b）：

$$
c_{ij} = \begin{cases}
\alpha_3+\alpha_4+\alpha_5, & i,j\in\{1,2\}, \\
\alpha_1+\alpha_2+\alpha_5, & i,j\in\{3,4\}, \\
\alpha_5, & i\in\{1,2\},\;j\in\{3,4\}.
\end{cases}
$$

**Type I 特例**（$\partial_i\!\cdot\!\partial_j$）：

$$
F^{(p)}G^{(q)} = (\partial_i B_5)\!\cdot\!(\partial_j B_5),\qquad
H^{(pq)} = \partial_i\!\cdot\!\partial_j B_5.
$$

其中 $H_{ij}^{(\text{I})}$ 将 $(v_p\!\cdot\!v_q)$ 替换为 $D=4-2\varepsilon$。

### 8c. $v_p\!\cdot\!(X_i-X_j)$ 计算（代码）

对每条腿 $i$ 和每个速度 $v_p$，代码计算 $v_p\!\cdot\!(X_i-X_j)$。所有 $4\times4\times4=64$ 个组合中，非零项共 48 个。完整表见 `compute_pnum_1111.py`。

### 8d. Pnum 组装

对 20 个子项逐项计算 $F^{(p_t)},G^{(q_t)},H^{(p_t q_t)}$，求和：

$$
\boxed{
\begin{aligned}
\text{Pnum}^{(1)} &= (2-3\varepsilon)\sum_{t=1}^{20} C_t\,F^{(p_t)}G^{(q_t)}, \\[8pt]
\text{Pnum}^{(2)} &= -\sum_{t=1}^{20} C_t\,H^{(p_t q_t)}.
\end{aligned}}
$$

### 8e. $\text{Pnum}^{(2)}$ 结果（闭合形式）

代码计算得：

$$
\boxed{
\begin{aligned}
\text{Pnum}^{(2)} &= \alpha_5\Bigl[
   \alpha_1\alpha_3\bigl(4\varepsilon z^2 - 2z^2 + 4z - 2\bigr)
  +\alpha_1\alpha_4\bigl(-4\varepsilon z^2 + 8z^2 + 8z + 2\bigr) \\
  &\quad +\alpha_2\alpha_3\bigl(-4\varepsilon z^2 + 8z^2 + 8z + 2\bigr)
  +\alpha_2\alpha_4\bigl(4\varepsilon z^2 - 2z^2 + 4z - 2\bigr)
   \Bigr].
\end{aligned}}
$$

**关键性质**：
- **不含 $s_i$**——$H_{ij}^{(pq)}\propto(v_p\!\cdot\!v_q)$，与 $X$ 无关
- **仅 $\alpha_5^1$**——$\alpha_5^0$ 项因行和为零（$\sum_j H_{ij}=0$）而抵消，$\alpha_5^2$ 项因 $c_{ij}$ 不含 $\alpha_5^2$ 而不存在

**$\varepsilon\to 0$ 极限**：

$$
\text{Pnum}^{(2)}\big|_{\varepsilon=0} = 2\alpha_5\Bigl[
   -\alpha_1\alpha_3(z-1)^2 + \alpha_1\alpha_4(2z+1)^2 + \alpha_2\alpha_3(2z+1)^2 - \alpha_2\alpha_4(z-1)^2
   \Bigr].
$$

> **与 S 道退化的对应**：S 道中 $v_1,v_2\to v_3$ 且 $v_3,v_4\to v_4$，此时 $\alpha_1\to\alpha_3^{(1)},\alpha_2\to\alpha_3^{(2)},\alpha_3\to\alpha_4^{(1)},\alpha_4\to\alpha_4^{(2)}$。四项合并为 $\alpha_5(\alpha_3-\alpha_4)(\alpha_1-\alpha_2)z$——与 S-derivation §8d 一致 $\checkmark$。

### 8f. $\text{Pnum}^{(1)}$ 结构

$\text{Pnum}^{(1)}$ 是 $s_i,\alpha_i,z$ 的多项式（$\sim$24k 字符），按 $\alpha_5$ 幂次分层：

| $\alpha_5$ 幂次 | 来源 | $\alpha$ 结构 |
|:---:|------|------|
| $\alpha_5^0$ | $c_{ij}$ 中 $\Sigma_A,\Sigma_B$ 的 $\alpha_3+\alpha_4$ 和 $\alpha_1+\alpha_2$ | $\alpha_i\alpha_j\alpha_k\alpha_\ell$（4 次） |
| $\alpha_5^1$ | $c_{ij}$ 中 $\alpha_5$ 与 $\Sigma_A,\Sigma_B$ 交叉 | $\alpha_5\alpha_i\alpha_j\alpha_k$（4 次） |
| $\alpha_5^2$ | 两个 $c_{ij}$ 都取 $\alpha_5$（跨组 $FG$） | $\alpha_5^2\alpha_i\alpha_j$（4 次） |

完整表达式见 `compute_pnum_1111.py` 输出。

> **与 S 道对比**：S 道 $\text{Pnum}^{(1)}$ 仅 $\alpha_5^2$ 项（$\sim$200 字符），因为 $\Delta_{12}=\Delta_{34}=0$ 使 $c_{12}$ 和 $c_{34}$ 的项消失，且仅 Type III 存活。这是双退化 vs 双非退化的极端差异。

### 8g. $A_{\text{ref}}$ 标度

$B_5\propto A_{\text{ref}}$，$F,G\propto A_{\text{ref}}$，$H\propto A_{\text{ref}}^0$。因此：

$$
\frac{\text{Pnum}^{(1)}}{B_5^{3-3\varepsilon}} \propto A_{\text{ref}}^{-1+3\varepsilon},\qquad
\frac{\text{Pnum}^{(2)}}{B_5^{2-3\varepsilon}} \propto A_{\text{ref}}^{-2+3\varepsilon}.
$$

$A_{\text{ref}}^{3\varepsilon}$ 归入 $C$ 的 ⑥。$A_{\text{ref}}^{-1}=-1/(2v_{\text{ref}})$ 和 $A_{\text{ref}}^{-2}=1/(4v_{\text{ref}}^2)$ 与 $(1-3\varepsilon)$ 组合产生轨道依赖的 ⑤ derivative 数值系数。

### 8h. 代入后的 $I$

$$
\boxed{
\begin{aligned}
I_{(1,2,3,4)} &= \text{① color}\;g^6\,
   \frac{\Gamma(1-3\varepsilon)}{1024\,\pi^{6-3\varepsilon}}
   \times (1-3\varepsilon) \\
&\quad \times \int_0^1 ds_1 ds_2 ds_3 ds_4
   \int d\alpha\;\delta(1-\textstyle\sum\alpha)\,
   \frac{\prod\alpha_i^{-\varepsilon}}{\mathcal{U}^{1+2\varepsilon}} \\
&\quad \times \left[
   \frac{(2-3\varepsilon)\,\text{Pnum}^{(1)}}{B_5^{\,3-3\varepsilon}}
   + \frac{\text{Pnum}^{(2)}}{B_5^{\,2-3\varepsilon}}
   \right].
\end{aligned}}
$$

---

## 9. 超立方变换

四边各一腿，无嵌套：$x_i = s_i$（$i=1,2,3,4$），Jacobian $=1$。此步平凡。

---

## 10. $C$ 因子分解

$$
\boxed{
\begin{aligned}
C &= \underbrace{t^{A_1}t^{A_2}t^{A_3}t^{A_4}\,f^{A_1A_2B}f^{A_3A_4B}}_{\text{① color（待 FeynCalc）}}
   \times \underbrace{(4\pi\alpha_s)^3}_{\text{③ coupling}} \\
   &\quad \times \underbrace{\frac{\Gamma(1-3\varepsilon)}{1024\,\pi^{6-3\varepsilon}}}_{\text{④ ft\_base}}
   \times \underbrace{(1-3\varepsilon)}_{\text{⑤ derivative}}
   \times \underbrace{(-2v_{\text{ref}}+i0)^{3\varepsilon}}_{\text{⑥ v14\_scale}}.
\end{aligned}}
$$

**六个因子**：② vertex（$i$ 因子）已吸收进 ③ coupling。④ ft_base 与 S/T 道相同。⑤ derivative 的 $(1-3\varepsilon)$ 是"未提取 $A_{\text{ref}}$ 数值因子"的形式——$A_{\text{ref}}^{-1},A_{\text{ref}}^{-2}$ 与 $(1-3\varepsilon)$ 组合后产生轨道依赖的数值系数。

> **与 S 道对照**：
>
> | 因子 | S 道 | [1,1,1,1] | 差异 |
> |------|:---:|:---:|------|
> | ① color | $-\frac14 C_A^2 C_F$ | 待 FeynCalc | 四边互异，生成元不同 |
> | ③ coupling | $(4\pi\alpha_s)^3$ | 同 | — |
> | ④ ft_base | $\Gamma(1-3\varepsilon)/(1024\pi^{6-3\varepsilon})$ | 同 | — |
> | ⑤ derivative | $(1-3\varepsilon)$ | 同 | 数值系数不同 |
> | ⑥ v14_scale | $(-2v_{\text{ref}})^{3\varepsilon}$ | 同 | — |

---

## 11. $J$ 被积函数与汇总

### 11a. 完整 $J$

$$
\boxed{
\begin{aligned}
J_{(1,2,3,4)} &= \int_0^1 d^4s
   \int d^5\alpha\;\delta(1-\textstyle\sum\alpha)\,
   \frac{(\alpha_1\alpha_2\alpha_3\alpha_4\alpha_5)^{-\varepsilon}}
        {\mathcal{U}^{1+2\varepsilon}} \\
   &\quad \times \left[
   \frac{(2-3\varepsilon)\,\text{Pnum}^{(1)}}{B_5^{\,3-3\varepsilon}}
   + \frac{\text{Pnum}^{(2)}}{B_5^{\,2-3\varepsilon}}
   \right].
\end{aligned}}
$$

其中：
- $\mathcal{U} = (\alpha_1+\alpha_2+\alpha_5)(\alpha_3+\alpha_4+\alpha_5)-\alpha_5^2$（§7c）
- $B_5$ 由 §7c 公式给出，$\Delta'_{ij}$ 由 §7d 显式给出（全部 6 项非零）
- $\text{Pnum}^{(2)}$ 由 §8e 显式给出（4 项 $\alpha_5$ 多项式，不含 $s$）
- $\text{Pnum}^{(1)}$ 由 §8f 给出结构，完整表达式见 `compute_pnum_1111.py`

**$z$ 无法全局提取**：$\Delta'_{ij}$ 中 $z$ 与 $s_i$ 混合（$\Delta'_{13}$ 含 $z$ 和 $1$ 的混合，$\Delta'_{14}$ 不含 $z$），$B_5$ 不是 $z$ 的单项式。与 S 道 $z^{3\varepsilon}$ 全局提取不同——这是双非退化情形的一般特征。

### 11b. 最终 $C \times J$

$$
\boxed{I_{(1,2,3,4)} = C \times J_{(1,2,3,4)}.}
$$

### 11c. 与框架的对比验证

| # | 验证项 | 框架预测 | 代码结果 | 状态 |
|---|--------|---------|---------|:---:|
| 1 | 存活算子项 | 20（全部 $v_{ij}\neq 0$） | 20 | ✅ |
| 2 | $\Delta'_{ij}$ 全部非零 | 6 项 | 6 项（§7d） | ✅ |
| 3 | $B_5$ 全部非零 | 6 项 | 6 项（§7c） | ✅ |
| 4 | $\text{Pnum}^{(2)}\propto\alpha_5$ | $\alpha_5^1$ 仅，不含 $s$ | $\alpha_5^1$ 仅，不含 $s$（§8e） | ✅ |
| 5 | $\text{Pnum}^{(2)}$ 行和消去 $\alpha_5^0$ | $\sum_j H_{ij}=0$ | $\alpha_5^0=0$ | ✅ |
| 6 | $\text{Pnum}^{(1)}$ 含 $s$ | $F,G\propto A_{\text{ref}}$ | $\alpha_5^0,\alpha_5^1,\alpha_5^2$ 均含 $s$ | ✅ |
| 7 | $z$ 无法全局提取 | $\Delta'_{ij}$ 含 $z$ 和 $1$ 混合 | 确认（§11a） | ✅ |
| 8 | S 道退化极限 | $\text{Pnum}^{(2)}$ 合并为 $\alpha_5(\alpha_3-\alpha_4)(\alpha_1-\alpha_2)z$ | 一致（§8e） | ✅ |

### 11d. 总结

[1,1,1,1] 配对 (12)(34) 是 4e-3g+3g 的最一般情形——20 算子项全活，$B_5$ 全 6 项，$\text{Pnum}^{(1)}$ 达 24k 字符。通用框架的算子项表 + building block 公式能系统地处理这种复杂度，产出结构正确的 $\text{Pnum}^{(1)},\text{Pnum}^{(2)}$。框架预测的 8 个结构特征全部得到代码验证，S 道退化极限一致。

**待完成**：① color 的 FeynCalc 求值。