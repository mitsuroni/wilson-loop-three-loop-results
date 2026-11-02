# 从通用 $V_{2\text{TGV}}$ 推导 S 道和 T 道

> 使用 `effective-vertices.md` §3c 的通用表达式，代入具体边分布 $(a,b,c,d)$，直接得到各道的 Lorentz 分子
> 日期：2026-08-11

---

## 1. 通用表达式回顾

$$
V_{2\text{TGV}}(a,b,c,d;\ell_1,\ell_2,\ell_3,\ell_4)
= f^{A_a A_b B}f^{A_c A_d B}\;
\frac{-i}{(\ell_1+\ell_2)^2+i0} \times \mathcal{L}(a,b,c,d),
$$

其中 $\mathcal{L}$ 为 4 组 Lorentz 项：

$$
\boxed{
\begin{aligned}
\mathcal{L} &=
\underbrace{(v_a\!\cdot\!v_b)(v_c\!\cdot\!v_d)\;(\ell_2-\ell_1)\!\cdot\!(\ell_4-\ell_3)}_{\text{G1}} \\
&+ \underbrace{(v_a\!\cdot\!v_b)\Bigl[
   2(v_c\!\cdot\!\ell_1)(v_d\!\cdot\!\ell_2) - 2(v_c\!\cdot\!\ell_2)(v_d\!\cdot\!\ell_1)
   + v_c\!\cdot\!(\ell_2-\ell_1)\,(v_d\!\cdot\!\ell_3)
   + v_d\!\cdot\!(\ell_1-\ell_2)\,(v_c\!\cdot\!\ell_4)\Bigr]}_{\text{G2}} \\
&+ \underbrace{2\,(v_c\!\cdot\!v_d)\Bigl[
   (-v_a\!\cdot\!\ell_1 - 2v_a\!\cdot\!\ell_2)(v_b\!\cdot\!\ell_4 - v_b\!\cdot\!\ell_3)
   + (2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)(v_a\!\cdot\!\ell_4 - v_a\!\cdot\!\ell_3)\Bigr]}_{\text{G3}} \\
&+ \underbrace{(v_a\!\cdot\!v_c)\,(2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)\;v_d\!\cdot\!(\ell_3-\ell_1-\ell_2)}_{\text{G4a}}
   + \underbrace{(v_a\!\cdot\!v_d)\,(2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)\;v_c\!\cdot\!(\ell_1+\ell_2-\ell_4)}_{\text{G4b}} \\
&- \underbrace{(v_b\!\cdot\!v_c)\,(v_a\!\cdot\!\ell_1 + 2v_a\!\cdot\!\ell_2)\;v_d\!\cdot\!(\ell_3-\ell_1-\ell_2)}_{\text{G4c}}
   - \underbrace{(v_b\!\cdot\!v_d)\,(v_a\!\cdot\!\ell_1 + 2v_a\!\cdot\!\ell_2)\;v_c\!\cdot\!(\ell_1+\ell_2-\ell_4)}_{\text{G4d}}.
\end{aligned}}
$$

---

## 2. S 道：$(a,b,c,d) = (3,3,4,4)$

### 2a. 代入 $v_{ij}$

$$
\begin{aligned}
v_a = v_b = v_3,&\quad v_c = v_d = v_4, \\
v_a\!\cdot\!v_b = v_{33} = 0,&\quad v_c\!\cdot\!v_d = v_{44} = 0, \\
v_a\!\cdot\!v_c = v_a\!\cdot\!v_d = v_b\!\cdot\!v_c = v_b\!\cdot\!v_d = v_{34} = z.
\end{aligned}
$$

### 2b. 逐组简化

**G1**：$(v_a\!\cdot\!v_b)(v_c\!\cdot\!v_d)(\cdots) = 0 \times 0 \times (\cdots) = 0$

**G2**：$(v_a\!\cdot\!v_b)[\cdots] = 0 \times [\cdots] = 0$

**G3**：$2(v_c\!\cdot\!v_d)[\cdots] = 2 \times 0 \times [\cdots] = 0$

**G4a**：$z \times (2v_3\!\cdot\!\ell_1 + v_3\!\cdot\!\ell_2) \times v_4\!\cdot\!(\ell_3-\ell_1-\ell_2)$

**G4b**：$z \times (2v_3\!\cdot\!\ell_1 + v_3\!\cdot\!\ell_2) \times v_4\!\cdot\!(\ell_1+\ell_2-\ell_4)$

**G4c**：$-z \times (v_3\!\cdot\!\ell_1 + 2v_3\!\cdot\!\ell_2) \times v_4\!\cdot\!(\ell_3-\ell_1-\ell_2)$

**G4d**：$-z \times (v_3\!\cdot\!\ell_1 + 2v_3\!\cdot\!\ell_2) \times v_4\!\cdot\!(\ell_1+\ell_2-\ell_4)$

### 2c. 合并

G4a + G4c 的 $v_4$ 因子同为 $v_4\!\cdot\!(\ell_3-\ell_1-\ell_2)$：

$$
\begin{aligned}
\text{G4a}+\text{G4c} &= z \times \bigl[(2v_3\!\cdot\!\ell_1 + v_3\!\cdot\!\ell_2) - (v_3\!\cdot\!\ell_1 + 2v_3\!\cdot\!\ell_2)\bigr] \times v_4\!\cdot\!(\ell_3-\ell_1-\ell_2) \\
&= z \times (v_3\!\cdot\!\ell_1 - v_3\!\cdot\!\ell_2) \times v_4\!\cdot\!(\ell_3-\ell_1-\ell_2).
\end{aligned}
$$

G4b + G4d 的 $v_4$ 因子同为 $v_4\!\cdot\!(\ell_1+\ell_2-\ell_4)$：

$$
\begin{aligned}
\text{G4b}+\text{G4d} &= z \times \bigl[(2v_3\!\cdot\!\ell_1 + v_3\!\cdot\!\ell_2) - (v_3\!\cdot\!\ell_1 + 2v_3\!\cdot\!\ell_2)\bigr] \times v_4\!\cdot\!(\ell_1+\ell_2-\ell_4) \\
&= z \times (v_3\!\cdot\!\ell_1 - v_3\!\cdot\!\ell_2) \times v_4\!\cdot\!(\ell_1+\ell_2-\ell_4).
\end{aligned}
$$

两项相加：

$$
\begin{aligned}
\mathcal{L}_S &= z \times (v_3\!\cdot\!\ell_1 - v_3\!\cdot\!\ell_2) \times \bigl[v_4\!\cdot\!(\ell_3-\ell_1-\ell_2) + v_4\!\cdot\!(\ell_1+\ell_2-\ell_4)\bigr] \\
&= z \times (v_3\!\cdot\!\ell_1 - v_3\!\cdot\!\ell_2) \times v_4\!\cdot\!(\ell_3-\ell_4).
\end{aligned}
$$

### 2d. 结果

$$
\boxed{
\mathcal{L}_S = z\,(v_3\!\cdot\!(\ell_1-\ell_2))\,(v_4\!\cdot\!(\ell_3-\ell_4)).
}
$$

与 FeynCalc `lorentz_construct02_S.wl` 一致 $\checkmark$。

---

## 3. T 道：$(a,b,c,d) = (3,4,3,4)$

### 3a. 代入 $v_{ij}$

$$
\begin{aligned}
v_a = v_c = v_3,&\quad v_b = v_d = v_4, \\
v_a\!\cdot\!v_b = v_{34} = z,&\quad v_c\!\cdot\!v_d = v_{34} = z, \\
v_a\!\cdot\!v_c = v_{33} = 0,&\quad v_b\!\cdot\!v_d = v_{44} = 0, \\
v_a\!\cdot\!v_d = v_{34} = z,&\quad v_b\!\cdot\!v_c = v_{34} = z.
\end{aligned}
$$

### 3b. 逐组简化

**G1**：$z^2 \times (\ell_2-\ell_1)\!\cdot\!(\ell_4-\ell_3) = z^2(\ell_1\!\cdot\!\ell_3 - \ell_1\!\cdot\!\ell_4 - \ell_2\!\cdot\!\ell_3 + \ell_2\!\cdot\!\ell_4)$

**G2**：代入 $v_a\!\cdot\!v_b = z$，$v_c = v_3$，$v_d = v_4$：

$$
\begin{aligned}
\text{G2} &= z \times \Bigl[
   2(v_3\!\cdot\!\ell_1)(v_4\!\cdot\!\ell_2) - 2(v_3\!\cdot\!\ell_2)(v_4\!\cdot\!\ell_1) \\
   &\qquad + v_3\!\cdot\!(\ell_2-\ell_1)(v_4\!\cdot\!\ell_3) + v_4\!\cdot\!(\ell_1-\ell_2)(v_3\!\cdot\!\ell_4)
   \Bigr].
\end{aligned}
$$

**G3**：代入 $v_c\!\cdot\!v_d = z$，$v_a = v_3$，$v_b = v_4$：

$$
\begin{aligned}
\text{G3} &= 2z \times \Bigl[
   (-v_3\!\cdot\!\ell_1 - 2v_3\!\cdot\!\ell_2)(v_4\!\cdot\!\ell_4 - v_4\!\cdot\!\ell_3) \\
   &\qquad + (2v_4\!\cdot\!\ell_1 + v_4\!\cdot\!\ell_2)(v_3\!\cdot\!\ell_4 - v_3\!\cdot\!\ell_3)
   \Bigr].
\end{aligned}
$$

**G4a**：$(v_a\!\cdot\!v_c) = v_{33} = 0 \Rightarrow 0$

**G4b**：$z \times (2v_4\!\cdot\!\ell_1 + v_4\!\cdot\!\ell_2) \times v_3\!\cdot\!(\ell_1+\ell_2-\ell_4)$

**G4c**：$(v_b\!\cdot\!v_c) = z$：
$-z \times (v_3\!\cdot\!\ell_1 + 2v_3\!\cdot\!\ell_2) \times v_4\!\cdot\!(\ell_3-\ell_1-\ell_2)$

**G4d**：$(v_b\!\cdot\!v_d) = v_{44} = 0 \Rightarrow 0$

### 3c. 合并（与 FeynCalc 对照）

G2 + G3 + G4b + G4c 展开后共 14 个独立项，与 FeynCalc `lorentz_construct02.wl` 一致：

$$
\boxed{
\begin{aligned}
\mathcal{L}_T &= z^2\,(\ell_1\!\cdot\!\ell_3 - \ell_1\!\cdot\!\ell_4 - \ell_2\!\cdot\!\ell_3 + \ell_2\!\cdot\!\ell_4) \\
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

---

## 4. U 道：$(a,b,c,d) = (3,4,4,3)$

### 4a. 颜色

$f^{A_3 A_4 B}f^{A_4 A_3 B} = f^{A_3 A_4 B}(-f^{A_3 A_4 B}) = 0$（Jacobi 恒等式）。

$$
\boxed{V_{2\text{TGV}}(3,4,4,3) = 0.}
$$

---

## 5. 总结

| 道 | $(a,b,c,d)$ | $v_{ab}$ | $v_{cd}$ | 非零组 | 结果 |
|----|-------------|:---:|:---:|--------|------|
| S | (3,3,4,4) | 0 | 0 | G4a-d | $z(v_3\!\cdot\!(\ell_1-\ell_2))(v_4\!\cdot\!(\ell_3-\ell_4))$ |
| T | (3,4,3,4) | $z$ | $z$ | G1-4 | 14 项 |
| U | (3,4,4,3) | $z$ | $z$ | — | 颜色 $=0$ |

**关键**：通用表达式 $\mathcal{L}(a,b,c,d)$ 对任意边分布成立。只需代入 $v_{ij}$ 值，$v_{ij}=0$ 的项自动消失——退化自动处理，无需额外规则。