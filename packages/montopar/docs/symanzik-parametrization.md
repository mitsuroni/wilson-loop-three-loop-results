# Symanzik 参数化方法：3-点验证 + 5-点推导

> 用一步 Feynman 参数化 + 矩阵高斯积分重新推导 3-点恒等式，验证方法正确性，再推广到 5-点
> 日期：2026-08-12

---

## 1. 矩阵形式回顾

$L$-圈 $N$-点积分，Feynman 参数化后分母 $Q = \sum\alpha_i\Delta_i$ 是 $z$ 的二次型：

$$
Q = \vec{z}^T(\Sigma\otimes I_D)\vec{z} - 2\vec{z}\cdot\vec{J} + C
$$

其中 $\vec{z} = (z_1,...,z_L)^T$，$\Sigma$ 是 $L\times L$ 矩阵，$\vec{J}$ 是 $L$ 维向量。

配平方：$Q = (\vec{z}-\vec{z}_0)^T(\Sigma\otimes I_D)(\vec{z}-\vec{z}_0) + \mathcal{F}/\mathcal{U}$

- **第一 Symanzik 多项式**：$\mathcal{U} = \det\Sigma$
- **第二 Symanzik 多项式**：$\mathcal{F} = C\cdot\mathcal{U} - \vec{J}^T\text{adj}(\Sigma)\vec{J}$

$LD$ 维高斯积分：

$$
\int d^{LD}z\;Q^{-Na} = \frac{(i\pi^{D/2})^L\,\Gamma(Na-LD/2)}{\Gamma(Na)}\;\frac{\mathcal{U}^{Na-LD/2-LD/2}}{\mathcal{F}^{Na-LD/2}}
$$

---

## 2. 3-点验证（$L=1, N=3$）

### 2a. 设定

3 条传播子 $\Delta_i = -(X_i-z)^2$，1 个 $z$ 积分。

Feynman 参数化：

$$
\frac{1}{\Delta_1^a\Delta_2^a\Delta_3^a} = \frac{\Gamma(3a)}{\Gamma(a)^3}\int d\alpha\;\frac{\delta(1\!-\!\textstyle\sum\alpha)\prod\alpha_i^{a-1}}{Q^{3a}}
$$

### 2b. 矩阵 $\Sigma$（$1\times 1$）

$Q = \alpha_1(X_1\!-\!z)^2 + \alpha_2(X_2\!-\!z)^2 + \alpha_3(X_3\!-\!z)^2$（注意符号：$\Delta_i = -(X_i-z)^2$，但 $Q = -\sum\alpha_i(X_i-z)^2$，这里用正平方展开后取 $Q = \sum\alpha_i\Delta_i$）

展开 $Q$：

$$
Q = -\underbrace{(\alpha_1+\alpha_2+\alpha_3)}_{\Sigma}\,z^2 + 2\underbrace{(\alpha_1 X_1+\alpha_2 X_2+\alpha_3 X_3)}_{\vec{J}}\cdot z - \underbrace{(\alpha_1 X_1^2+\alpha_2 X_2^2+\alpha_3 X_3^2)}_{C}
$$

$\Sigma$ 是 $1\times 1$ 矩阵：

$$
\Sigma = \alpha_1+\alpha_2+\alpha_3 = 1 \quad(\delta\text{ 约束})
$$

### 2c. Symanzik 多项式

$$
\mathcal{U} = \det\Sigma = 1 \quad\text{(平凡)}
$$

$$
\mathcal{F} = C\cdot\mathcal{U} - \vec{J}\cdot\text{adj}(\Sigma)\cdot\vec{J} = C - J^2
$$

计算 $C - J^2$（利用 $\sum\alpha_i = 1$）：

$$
C - J^2 = \sum\alpha_i X_i^2 - (\sum\alpha_i X_i)^2 = -\sum_{i<j}\alpha_i\alpha_j(X_i-X_j)^2 = \sum_{i<j}\alpha_i\alpha_j\Delta_{ij}
$$

因此：

$$
\boxed{\mathcal{U} = 1,\qquad \mathcal{F} = \alpha_1\alpha_2\Delta_{12}+\alpha_2\alpha_3\Delta_{23}+\alpha_3\alpha_1\Delta_{31}}
$$

### 2d. 高斯积分

$L=1$，$N=3$，$Na-LD/2 = 3a-D/2$，$Na-3LD/2+LD/2 = 3a-D$：

$$
\int d^Dz\;Q^{-3a} = \frac{i\pi^{D/2}\,\Gamma(3a-D/2)}{\Gamma(3a)}\;\frac{\mathcal{U}^{3a-D}}{\mathcal{F}^{3a-D/2}} = \frac{i\pi^{D/2}\,\Gamma(3a-D/2)}{\Gamma(3a)}\;\frac{1}{\mathcal{F}^{3a-D/2}}
$$

（$\mathcal{U}=1$ 消失）

### 2e. 合并

$$
\int d^Dz\;\prod_{i=1}^3\frac{1}{\Delta_i^a}
= \frac{i\pi^{D/2}\,\Gamma(3a-D/2)}{\Gamma(a)^3}\int d\alpha\;\frac{\delta(1\!-\!\textstyle\sum\alpha)\prod\alpha_i^{a-1}}{\mathcal{F}^{3a-D/2}}
$$

代入 $a=1\!-\!\varepsilon$，$3a\!-\!D/2 = 1\!-\!2\varepsilon$：

$$
\boxed{
\int d^Dz\;\prod_{i=1}^3\frac{1}{\Delta_i^a}
= \frac{i\pi^{D/2}\,\Gamma(1\!-\!2\varepsilon)}{\Gamma(1\!-\!\varepsilon)^3}
\int d\alpha\;\frac{\delta(1\!-\!\textstyle\sum\alpha)\prod\alpha_i^{a-1}}
{[\alpha_1\alpha_2\Delta_{12}+\alpha_2\alpha_3\Delta_{23}+\alpha_3\alpha_1\Delta_{31}]^{1-2\varepsilon}}
}
$$

**与已知 3-点恒等式完全一致** $\checkmark$

---

## 3. 5-点推导（$L=2, N=5$）

### 3a. 设定

5 条传播子，2 个 $z$ 积分（$z_A, z_B$），拓扑：

```
X_1 ── z_A ── z_B ── X_3
       │            │
      X_2          X_4
```

### 3b. 矩阵 $\Sigma$（$2\times 2$）

$Q = \alpha_1\Delta_1+\alpha_2\Delta_2+\alpha_5\Delta_5+\alpha_3\Delta_3+\alpha_4\Delta_4$

展开（$\Delta_i = -(X_i-z)^2$，$\Delta_5 = -(z_A-z_B)^2$）：

$$
\Sigma = \begin{pmatrix}\alpha_1+\alpha_2+\alpha_5 & -\alpha_5 \\ -\alpha_5 & \alpha_3+\alpha_4+\alpha_5\end{pmatrix}
$$

$$
\vec{J} = \begin{pmatrix}\alpha_1 X_1+\alpha_2 X_2 \\ \alpha_3 X_3+\alpha_4 X_4\end{pmatrix},\qquad
C = \alpha_1 X_1^2+\alpha_2 X_2^2+\alpha_3 X_3^2+\alpha_4 X_4^2
$$

### 3c. Symanzik 多项式

$$
\boxed{\mathcal{U} = (\alpha_1+\alpha_2+\alpha_5)(\alpha_3+\alpha_4+\alpha_5) - \alpha_5^2}
$$

$$
\mathcal{F} = C\cdot\mathcal{U} - \vec{J}^T\,\text{adj}(\Sigma)\,\vec{J}
$$

其中 $\text{adj}(\Sigma) = \begin{pmatrix}\alpha_3+\alpha_4+\alpha_5 & \alpha_5 \\ \alpha_5 & \alpha_1+\alpha_2+\alpha_5\end{pmatrix}$

**Python 展开**（`scripts/compute_F_polynomial.py` 验证：所有 $X_i^2$ 项消去，$\mathcal{F}$ 纯由 $\Delta_{ij}$ 表达）。

定义 Feynman 分母 $B_5 = -\mathcal{F}$（物理正定）：

$$
\boxed{
B_5 = \alpha_1\alpha_2\underbrace{(\alpha_3+\alpha_4+\alpha_5)}_{\Sigma_B}\,\Delta_{12}
    + \alpha_3\alpha_4\underbrace{(\alpha_1+\alpha_2+\alpha_5)}_{\Sigma_A}\,\Delta_{34}
    + \alpha_5\!\!\sum_{\substack{i\in\{1,2\}\\j\in\{3,4\}}}\!\!\alpha_i\alpha_j\,\Delta_{ij}
}
$$

**结构**：
- $\Delta_{12}$（$z_A$ 侧同组）：$\alpha_1\alpha_2 \times \Sigma_B$（组内积 × 对侧总和）
- $\Delta_{34}$（$z_B$ 侧同组）：$\alpha_3\alpha_4 \times \Sigma_A$（对称！）
- $\Delta_{ij}$（跨组 $i\in\{1,2\}, j\in\{3,4\}$）：$\alpha_i\alpha_j\alpha_5$（桥接 $\alpha_5$）

### 3d. 高斯积分

$Q = -(z-z_0)^T(\Sigma\otimes I_D)(z-z_0) + B_5/\mathcal{U}$，$2D$ 维高斯：

$$
\int d^{2D}z\;Q^{-5a} = \frac{-\pi^D\,\Gamma(5a\!-\!D)}{\Gamma(5a)}\;\frac{\mathcal{U}^{5a-3D/2}}{B_5^{5a-D}}
= \frac{-\pi^D\,\Gamma(1\!-\!3\varepsilon)}{\Gamma(5a)}\;\frac{1}{\mathcal{U}^{1+2\varepsilon}\;B_5^{1-3\varepsilon}}
$$

### 3e. 合并

$$
\boxed{
\int d^Dz_A\,d^Dz_B\;\prod_{i=1}^5\frac{1}{\Delta_i^a}
= \frac{-\pi^D\,\Gamma(1\!-\!3\varepsilon)}{\Gamma(1\!-\!\varepsilon)^5}
\int d\alpha\;\frac{\delta(1\!-\!\textstyle\sum\alpha)\,\prod\alpha_i^{a-1}}
{\mathcal{U}^{1+2\varepsilon}\;B_5^{1-3\varepsilon}}
}
$$

其中 $a=1\!-\!\varepsilon$，且：

$$
\mathcal{U} = (\alpha_1\!+\!\alpha_2\!+\!\alpha_5)(\alpha_3\!+\!\alpha_4\!+\!\alpha_5) - \alpha_5^2
$$

$$
B_5 = \alpha_1\alpha_2\Sigma_B\,\Delta_{12} + \alpha_3\alpha_4\Sigma_A\,\Delta_{34}
    + \alpha_5\!\sum_{\substack{i\in\{1,2\}\\j\in\{3,4\}}}\!\alpha_i\alpha_j\,\Delta_{ij}
$$

$$
\Sigma_A = \alpha_1\!+\!\alpha_2\!+\!\alpha_5,\quad
\Sigma_B = \alpha_3\!+\!\alpha_4\!+\!\alpha_5,\quad
\Delta_{ij} = -(X_i\!-\!X_j)^2\!+\!i0
$$

---

## 4. 3-点 vs 5-点对比

| | 3-点（$L=1$） | 5-点（$L=2$） |
|---|---|---|
| $\Sigma$ | $1\times 1$：$\Sigma=1$ | $2\times 2$ 矩阵 |
| $\mathcal{U}$ | $= 1$（平凡，不出现） | $= \det\Sigma$（非平凡） |
| $\mathcal{F}$（即 $B$） | $\sum_{i<j}\alpha_i\alpha_j\Delta_{ij}$（3 项） | $\alpha_1\alpha_2\Sigma_B\Delta_{12}+\alpha_3\alpha_4\Sigma_A\Delta_{34}+\alpha_5\sum\alpha_i\alpha_j\Delta_{ij}$（6 项） |
| 参数 | $\alpha_1,\alpha_2,\alpha_3$（1 组，1 个 $\delta$） | $\alpha_1,...,\alpha_5$（1 组，1 个 $\delta$） |
| prefactor | $\frac{i\pi^{D/2}\Gamma(1-2\varepsilon)}{\Gamma(1-\varepsilon)^3}$ | $\frac{-\pi^D\Gamma(1-3\varepsilon)}{\Gamma(1-\varepsilon)^5}$ |
| 分母 | $B^{1-2\varepsilon}$ | $\mathcal{U}^{1+2\varepsilon}\,B_5^{1-3\varepsilon}$ |

> **关键**：3-点的 $\mathcal{U}=1$（1-圈平凡），5-点 $\mathcal{U}\neq 1$（2-圈非平凡）。这是 $L=1$ 与 $L=2$ 的本质区别——$\mathcal{U}$ 的出现是 2-圈图的特征。