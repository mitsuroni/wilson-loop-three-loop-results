# 4e-3g+3g T 道推导：双 TGV 顶点耦合（交叉配对，双双非退化）

> 轨道：4e-3g+3g-05 — 边分布 $(3,3,4,4)$，T 道配对 $(13)(24)$，两个 TGV 均非退化
> 参考：`test/TGV/I_4wl2x3g-02_因子分解推导.md`，`test/TGV/TGV-TGV_通用因子分解推导.md`
> 日期：2026-08-11

---

## 1. 运动学

### 1a. N=4 四边形

同 S 道。$v_i^2=0$，$v_{\text{ref}} = v_{14}$：

$$
\boxed{
\begin{aligned}
v_{12} &= v_{34} = z\,v_{\text{ref}}, \\
v_{14} &= v_{23} = v_{\text{ref}}, \\
v_{13} &= v_{24} = -(1+z)\,v_{\text{ref}}.
\end{aligned}}
$$

以下设 $v_{\text{ref}}=1$。

### 1b. 边分布与 Schwinger 参数

同 S 道：边 3 有 2 个附着点（$s_{3,1} < s_{3,2}$），边 4 有 2 个附着点（$s_{4,1} < s_{4,2}$）。

超立方变换相同：$x_1 = s_{3,2},\; x_2 = s_{3,1}/s_{3,2},\; x_3 = s_{4,2},\; x_4 = s_{4,1}/s_{4,2}$，Jacobian $= x_1 x_3$。

### 1c. 拓扑（与 S 道的本质区别）

```
边3: X_3(s_{3,1}) ──→ [TGV A, z_A] ←── X_4(s_{4,2}) :边4
                         │
                        ℓ₅  (内部传播子)
                         │
边3: X_3(s_{3,2}) ──→ [TGV B, z_B] ←── X_4(s_{4,1}) :边4
```

**关键**：TGV A 的腿在边 3 和边 4（$v_{34} = z \neq 0$），TGV B 的腿也在边 3 和边 4。两个 TGV 均**非退化**——每个 TGV 的两条外腿来自不同边。

四个外部点：

$$
\begin{aligned}
X_1 &\equiv X_3(s_{3,1}) = v_1+v_2 + x_1 x_2\,v_3 \quad\text{（TGV A，边 3）}\\
X_2 &\equiv X_4(s_{4,2}) = v_1+v_2+v_3 + x_3\,v_4 \quad\text{（TGV A，边 4）}\\
X_3 &\equiv X_3(s_{3,2}) = v_1+v_2 + x_1\,v_3 \quad\text{（TGV B，边 3）}\\
X_4 &\equiv X_4(s_{4,1}) = v_1+v_2+v_3 + x_3 x_4\,v_4 \quad\text{（TGV B，边 4）}
\end{aligned}
$$

### 1d. S 道 vs T 道对比

| 特征 | S 道 | T 道 |
|------|------|------|
| 配对 | (12)(34) — 同侧 | (13)(24) — 交叉 |
| TGV A | 边 3, 边 3（退化） | 边 3, 边 4（非退化） |
| TGV B | 边 4, 边 4（退化） | 边 3, 边 4（非退化） |
| $v_{ab}$ | $v_{33}=0$ | $v_{34}=z$ |
| $v_{cd}$ | $v_{44}=0$ | $v_{34}=z$ |
| 颜色 | $-\frac{1}{4}C_A^2 C_F$ | $+\frac{1}{4}C_A^2 C_F$ |
| Lorentz 项数 | 1 | 14（9 类） |
| $\Delta_{12}$ | $0$（同边 3） | $\neq 0$（跨边 3-4） |
| $\Delta_{34}$ | $0$（同边 4） | $\neq 0$（跨边 3-4） |

---

## 2. 动量空间积分

### 2a. Feynman 规则

同 S 道 §2a。

### 2b. 完整积分

**TGV A**：{ℓ₁ (边3, $s_{3,1}$), ℓ₂ (边4, $s_{4,2}$)}
**TGV B**：{ℓ₃ (边3, $s_{3,2}$), ℓ₄ (边4, $s_{4,1}$)}

动量守恒：$\ell_5 = -(\ell_1+\ell_2) = \ell_3+\ell_4$

$$
\boxed{
\begin{aligned}
I_T &= \int_0^{1} ds_{3,2} \int_0^{s_{3,2}} ds_{3,1} \int_0^{1} ds_{4,2} \int_0^{s_{4,2}} ds_{4,1} \\
&\quad \times (ig\,v_{3}^{\mu_{1}} t^{a_{1}}) (ig\,v_{4}^{\mu_{2}} t^{a_{2}}) (ig\,v_{3}^{\mu_{3}} t^{a_{3}}) (ig\,v_{4}^{\mu_{4}} t^{a_{4}}) \\
&\quad \times \int \frac{d^D\ell_{1}}{(2\pi)^D}\frac{d^D\ell_{2}}{(2\pi)^D}\frac{d^D\ell_{3}}{(2\pi)^D}\frac{d^D\ell_{4}}{(2\pi)^D}\frac{d^D\ell_{5}}{(2\pi)^D} \\
&\quad \times \frac{-ig_{\mu_{1}\nu_{1}}\delta^{a_{1}b_{1}}}{\ell_{1}^2+i0}\,e^{-i\ell_{1}\cdot X_{1}}
   \times \frac{-ig_{\mu_{2}\nu_{2}}\delta^{a_{2}b_{2}}}{\ell_{2}^2+i0}\,e^{-i\ell_{2}\cdot X_{2}} \\
&\quad \times \frac{-ig_{\mu_{3}\nu_{3}}\delta^{a_{3}b_{3}}}{\ell_{3}^2+i0}\,e^{-i\ell_{3}\cdot X_{3}}
   \times \frac{-ig_{\mu_{4}\nu_{4}}\delta^{a_{4}b_{4}}}{\ell_{4}^2+i0}\,e^{-i\ell_{4}\cdot X_{4}} \\
&\quad \times \frac{-ig_{\nu_{5}\nu_{6}}\delta^{b_{5}b_{6}}}{\ell_{5}^2+i0} \\
&\quad \times g\,f^{b_{1}b_{2}b_{5}}\;\Gamma^{\text{kin}}_{\nu_{1}\nu_{2}\nu_{5}}(-\ell_{1},-\ell_{2},-\ell_{5})
   \times g\,f^{b_{3}b_{4}b_{6}}\;\Gamma^{\text{kin}}_{\nu_{3}\nu_{4}\nu_{6}}(-\ell_{3},-\ell_{4},+\ell_{5}) \\
&\quad \times (2\pi)^D\delta^D(\ell_{1}+\ell_{2}+\ell_{5})
   \times (2\pi)^D\delta^D(\ell_{3}+\ell_{4}-\ell_{5}).
\end{aligned}}
$$

---

## 3. 颜色因子

### 3a. 颜色结构

T 道交叉配对：$f^{a_1 a_2 b} \to f^{a_1 a_4 b}$（ℓ₁ 连边 3，ℓ₂ 连边 4），$f^{a_3 a_4 b} \to f^{a_2 a_3 b}$（ℓ₃ 连边 3，ℓ₄ 连边 4）：

$$
\boxed{\text{color} = t^{a_1}t^{a_2}t^{a_3}t^{a_4}\;f^{a_1 a_4 b}\,f^{a_2 a_3 b}.}
$$

### 3b. 颜色求值

FeynCalc `SUNSimplify` 给出：

$$
\boxed{\text{① color} = +\frac{1}{4}\,C_A^2\,C_F.}
$$

> **与 S 道对照**：S 道 $f^{a_1 a_2 b}f^{a_3 a_4 b}$ 给出 $-\frac{1}{4}C_A^2 C_F$。T 道差一个符号。U 道 $f^{a_1 a_3 b}f^{a_2 a_4 b}$ 给出 $0$（Jacobi 恒等式：$f^{a_1 a_3 b}f^{a_2 a_4 b} + f^{a_1 a_2 b}f^{a_3 a_4 b} + f^{a_1 a_4 b}f^{a_2 a_3 b} = 0$）。

---

## 4. 位置空间、Lorentz 缩并、$\partial_X$ 算子

### 4a. 位置空间形式

$\delta$ 函数 → $z_A, z_B$ 积分，与 S 道完全相同：

$$
\boxed{
\begin{aligned}
I_T &= \int d^4s \int d^Dz_A\,d^Dz_B \\
&\quad \times (ig)^4 g^2\;v_{3}^{\mu_{1}}v_{4}^{\mu_{2}}v_{3}^{\mu_{3}}v_{4}^{\mu_{4}}\;
     t^{a_{1}}t^{a_{2}}t^{a_{3}}t^{a_{4}}\;f^{b_{1}b_{2}b_{5}}f^{b_{3}b_{4}b_{6}} \\
&\quad \times \Delta_{\mu_{1}\nu_{1}}(X_1-z_A)\;\Delta_{\mu_{2}\nu_{2}}(X_2-z_A) \\
&\quad \times \Delta_{\mu_{3}\nu_{3}}(X_3-z_B)\;\Delta_{\mu_{4}\nu_{4}}(X_4-z_B) \\
&\quad \times \Delta_{\nu_{5}\nu_{6}}(z_A-z_B) \\
&\quad \times \Gamma^{\text{kin}}_{\nu_{1}\nu_{2}\nu_{5}}\;\Gamma^{\text{kin}}_{\nu_{3}\nu_{4}\nu_{6}}.
\end{aligned}}
$$

### 4b. Lorentz 张量缩并（9 项分类）

T 道 $v_{ab}=v_{34}=z \neq 0$，$v_{cd}=v_{34}=z \neq 0$。所有 9 项均非零。

缩并后的 Lorentz 标量：

$$
V_T = v_3^{\nu_1}v_4^{\nu_2}v_3^{\nu_3}v_4^{\nu_4}\;
g^{\nu_5\nu_6}\;
\Gamma^A_{\nu_1\nu_2\nu_5}\;\Gamma^B_{\nu_3\nu_4\nu_6}.
$$

展开为 3 类 9 项（$3\times 3$）：

#### 类型 I：$\ell_i\cdot\ell_j$（1 项，$A_1B_1$）

$$
\boxed{V_{\text{I}} = v_{34}^2\;(\ell_2-\ell_1)\!\cdot\!(\ell_4-\ell_3).}
$$

#### 类型 II：$v_{ij}\,v_k\!\cdot\!\ell$（4 项，$A_1B_{2,3}$ 和 $A_{2,3}B_1$）

$$
\boxed{
\begin{aligned}
V_{\text{II}} &= v_{34}\Big[v_3\!\cdot\!(-\ell_3-2\ell_4+\ell_1+\ell_2)\;v_4\!\cdot\!(\ell_2-\ell_1) \\
&\qquad + v_4\!\cdot\!(2\ell_3+\ell_4-\ell_1-\ell_2)\;v_3\!\cdot\!(\ell_2-\ell_1)\Big] \\
&\quad + v_{34}\Big[v_3\!\cdot\!(-\ell_1-2\ell_2)\;v_4\!\cdot\!(\ell_4-\ell_3) \\
&\qquad + v_4\!\cdot\!(2\ell_1+\ell_2)\;v_3\!\cdot\!(\ell_4-\ell_3)\Big].
\end{aligned}}
$$

#### 类型 III：$(v_i\!\cdot\!\ell)(v_j\!\cdot\!\ell)$（4 项，$A_{2,3}B_{2,3}$）

$$
\boxed{
\begin{aligned}
V_{\text{III}} &= v_3\!\cdot\!(-\ell_1-2\ell_2)\;v_3\!\cdot\!(-\ell_3-2\ell_4+\ell_1+\ell_2)\;v_{44} \\
&\quad + v_3\!\cdot\!(-\ell_1-2\ell_2)\;v_4\!\cdot\!(2\ell_3+\ell_4-\ell_1-\ell_2)\;v_{34} \\
&\quad + v_4\!\cdot\!(2\ell_1+\ell_2)\;v_3\!\cdot\!(-\ell_3-2\ell_4+\ell_1+\ell_2)\;v_{34} \\
&\quad + v_4\!\cdot\!(2\ell_1+\ell_2)\;v_4\!\cdot\!(2\ell_3+\ell_4-\ell_1-\ell_2)\;v_{33}.
\end{aligned}}
$$

利用 $v_{33}=v_{44}=0$，$V_{\text{III}}$ 中首末两项消失，仅剩 2 项 $\propto v_{34}$。

### 4c. FeynCalc 验证

`lorentz_construct02.wl` 给出完整展开（14 个独立项）：

$$
\boxed{
\begin{aligned}
V_T &= z^2\,(\ell_1\!\cdot\!\ell_3 - \ell_1\!\cdot\!\ell_4 - \ell_2\!\cdot\!\ell_3 + \ell_2\!\cdot\!\ell_4) \\
  &- z\Bigl[
    3(v_3\!\cdot\!\ell_1)(v_4\!\cdot\!\ell_1)
    + 2(v_4\!\cdot\!\ell_1)(v_3\!\cdot\!\ell_2)
    + 4(v_3\!\cdot\!\ell_1)(v_4\!\cdot\!\ell_2)
    + 3(v_3\!\cdot\!\ell_2)(v_4\!\cdot\!\ell_2) \\
  &\quad + 2(v_4\!\cdot\!\ell_1)(v_3\!\cdot\!\ell_3)
    + (v_4\!\cdot\!\ell_2)(v_3\!\cdot\!\ell_3)
    + (v_3\!\cdot\!\ell_1)(v_4\!\cdot\!\ell_3)
    - (v_3\!\cdot\!\ell_2)(v_4\!\cdot\!\ell_3) \\
  &\quad - (v_4\!\cdot\!\ell_1)(v_3\!\cdot\!\ell_4)
    + (v_4\!\cdot\!\ell_2)(v_3\!\cdot\!\ell_4)
    + (v_3\!\cdot\!\ell_1)(v_4\!\cdot\!\ell_4)
    + 2(v_3\!\cdot\!\ell_2)(v_4\!\cdot\!\ell_4)
  \Bigr].
\end{aligned}}
$$

> **Python 验证**：`scripts/verify_4e3g3g_T.py` 确认 $\Delta_{12}$ 和 $\Delta_{34}$ 均非零，且 $X_2-X_1$ 和 $X_4-X_3$ 同时含 $v_3$ 和 $v_4$ 分量，导致 $\partial_X$ 因式分解中 $z_A$ 和 $z_B$ 的因子混合——这是 T 道 Lorentz 分子有 14 项的根本原因。

### 4d. $\ell \to i\partial_X$：两类微分算子

T 道 Lorentz 分子含两类 $\ell$ 结构：

| 类型 | $\ell$ 结构 | 位置空间 | 算子类型 |
|------|------------|---------|---------|
| I | $\ell_i\!\cdot\!\ell_j$ | $-\partial_{X_i}\!\cdot\!\partial_{X_j}$ | 二阶交叉导数 |
| II/III | $v_k\!\cdot\!\ell_i$ | $i\,v_k\!\cdot\!\partial_{X_i}$ | 一阶导数 |

**类型 I**（$A_1B_1$，$\propto z^2$）：

$$
\boxed{\mathcal{O}_T^{\text{I}} = -z^2\,(\partial_{X_1}\!\cdot\!\partial_{X_3} - \partial_{X_1}\!\cdot\!\partial_{X_4} - \partial_{X_2}\!\cdot\!\partial_{X_3} + \partial_{X_2}\!\cdot\!\partial_{X_4})}.
$$

**类型 II/III**（$\propto z$）：14 个 $(v\!\cdot\!\partial_X)(v\!\cdot\!\partial_X)$ 项，不逐一列出。

### 4e. $\partial_X$ 作用于传播子

一阶导数（同 S 道）：

$$
v_j\!\cdot\!\partial_{X_i}\,\Delta_i^{-a} = 2a\,v_j\!\cdot\!(X_i-z)\,\Delta_i^{-a-1}.
$$

二阶交叉导数（T 道特有）：

$$
\partial_{X_i}\!\cdot\!\partial_{X_j}\,(\Delta_i^{-a}\Delta_j^{-a})
= 2a\,g^{\mu\nu}\,\partial_{X_i^\mu}\Delta_i^{-a}\;\partial_{X_j^\nu}\Delta_j^{-a}
= 4a^2\,(X_i-z)\!\cdot\!(X_j-z)\,\Delta_i^{-a-1}\Delta_j^{-a-1}.
$$

---

## 5. Fourier 变换

### 5a. 5-点位置空间结构

5 条传播子，指数 $a = D/2-1 = 1-\varepsilon$。$\partial_X$ 作用后外部传播子指数变为 $a+1 = 2-\varepsilon$（同 S 道）。

### 5b. 两步 Feynman 参数化

与 S 道结构相同，但 $\Delta_{12} \neq 0$ 和 $\Delta_{34} \neq 0$ 使 $\tilde{c}_0 \neq 0$ 和 $B_B$ 含 $\beta_3\beta_4\Delta_{34}$ 项。

**Step 1**（$z_A$ 积分）：引入 $\alpha_1,\alpha_2,\alpha_5$，$\Sigma_A = \alpha_1+\alpha_2+\alpha_5$。得 Feynman 分母 $B_A = c\,\Delta_{\bar{X},B} + \tilde{c}_0$，其中 $\tilde{c}_0 \propto \Delta_{12} \neq 0$。

**Step 2**（$z_B$ 积分）：引入 $\beta_3,\beta_4,\beta_5$，$\Sigma_B = \beta_3+\beta_4+\beta_5$。得有效分母 $D = \beta_5\tilde{c}_0 T + B_B$，其中：

$$
B_B = \beta_3\beta_4\Delta_{34} + \beta_5 c(\beta_3\Delta_{3\bar{X}} + \beta_4\Delta_{4\bar{X}}).
$$

### 5c. FT prefactor

两步 Feynman 参数化（含 $\partial_X$ 修正的传播子指数）给出：

$$
\boxed{\text{④ ft\_base} = \frac{\Gamma(3-2\varepsilon)\,\Gamma(5-3\varepsilon)}
                                {1024\,\pi^{6-3\varepsilon}\,\Gamma(2-\varepsilon)^2\,\Gamma(1-\varepsilon)}.
}
$$

> 与 S 道相同——$\partial_X$ 对传播子指数的修正不依赖于边分布。

---

## 6. T 道 $\Delta$ 结构

### 6a. 交叉 $\Delta$ 均非零

与 S 道不同，T 道四条腿交叉分布在边 3 和边 4，导致：

$$
\begin{aligned}
\Delta_{12} &= -(X_1-X_2)^2 = -(v_1+v_2+s_{3,1}v_3 - v_1-v_2-v_3-s_{4,2}v_4)^2 \\
            &= -((s_{3,1}-1)v_3 - s_{4,2}v_4)^2 \\
            &= -2s_{4,2}(1-s_{3,1})\,z \neq 0, \\[4pt]
\Delta_{34} &= -(X_3-X_4)^2 = -((s_{3,2}-1)v_3 - s_{4,1}v_4)^2 \\
            &= -2s_{4,1}(1-s_{3,2})\,z \neq 0.
\end{aligned}
$$

### 6b. 与 S 道的本质区别

| 量 | S 道 | T 道 |
|----|------|------|
| $\Delta_{12}$ | $0$ | $-2x_3(1-x_1x_2)z$ |
| $\Delta_{34}$ | $0$ | $-2(1-x_1)(1-x_3x_4)z$ |
| $\tilde{c}_0$ | $0$ | $\neq 0$ |
| $B_B$ 结构 | $\beta_5c(\beta_3\Delta_{3\bar{X}}+\beta_4\Delta_{4\bar{X}})$ | $+\beta_3\beta_4\Delta_{34}$ 项 |
| $z$ 全局提取 | 可（仅 $z^{3\varepsilon}$） | **不可**（$z$ 与 $x$ 混合） |

---

## 7. C×J 分解

### 7a. C 因子

与 S 道共享相同的④ ft_base、⑤ derivative、⑥ v14_scale，仅① color 符号相反：

$$
\boxed{
\begin{aligned}
C &= \underbrace{+\frac{1}{4}C_A^2 C_F}_{\text{① color}}
   \times \underbrace{g^6}_{\text{② vertex}}
   \times \underbrace{(4\pi\alpha_s)^3}_{\text{③ coupling}}
   \times \underbrace{\frac{\Gamma(3-2\varepsilon)\,\Gamma(5-3\varepsilon)}
                            {1024\,\pi^{6-3\varepsilon}\,\Gamma(2-\varepsilon)^2\,\Gamma(1-\varepsilon)}}_{\text{④ ft\_base}}
   \times \underbrace{\bigl(-16(1-\varepsilon)^2\bigr)}_{\text{⑤ derivative}}
   \times \underbrace{(-2v_{\text{ref}}+i0)^{3\varepsilon}}_{\text{⑥ v14\_scale}}.
\end{aligned}}
$$

| 因子 | S 道 | T 道 |
|------|------|------|
| ① color | $-\frac{1}{4}C_A^2 C_F$ | $+\frac{1}{4}C_A^2 C_F$ |
| ②-⑥ | 相同 | 相同 |

### 7b. J 因子

T 道的 J 因子比 S 道复杂得多：

1. **Pnum**：14 项，含一阶和二阶 $\partial_X$ 算子
2. **$S$ 分母**：$\tilde{c}_0 \neq 0$，$B_B$ 含 $\beta_3\beta_4\Delta_{34}$ 项
3. **$z$ 依赖**：$\Delta_{12}, \Delta_{34} \propto z$，但 $\tilde{c}_0$ 和 $B_B$ 中 $z$ 与 $x$ 变量混合，**$z$ 无法全局提取**

$$
\boxed{
J = \int_0^1 d^4x \int d\alpha\,\delta(1-\Sigma_A) \int d\beta\,\delta(1-\Sigma_B)\;
\frac{\text{Jacobian} \times \text{Pnum}}
{(\alpha\text{-measure})\;(\beta\text{-measure})\;
 \Sigma_A^{3-\varepsilon}\,T^{3-\varepsilon}\,
 (\beta_5\tilde{c}_0 T + B_B)^{5-3\varepsilon}}.
}
$$

### 7c. S 道 vs T 道总结

| 特征 | S 道 | T 道 |
|------|------|------|
| 颜色 | $-\frac{1}{4}C_A^2 C_F$ | $+\frac{1}{4}C_A^2 C_F$ |
| Lorentz 项数 | 1 | 14 |
| $\partial_X$ 类型 | 仅一阶 | 一阶 + 二阶交叉 |
| $\Delta_{12}$ | $0$ | $\neq 0$ |
| $\Delta_{34}$ | $0$ | $\neq 0$ |
| $\tilde{c}_0$ | $0$ | $\neq 0$ |
| $z$ 全局提取 | 可（$z^{3\varepsilon}$） | 不可 |
| $J$ 复杂度 | 低 | 高 |

---

**总结**：T 道是非退化双 TGV 的典型代表。S 道的极端简化（$\Delta_{12}=\Delta_{34}=0$，$z$ 全局提取）是相邻边双退化的特殊结果，不具有一般性。T 道代表了双 TGV 耦合的一般情况。