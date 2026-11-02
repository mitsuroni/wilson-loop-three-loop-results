# 3e-3g 四个 orbit 的 Contour IBP $\sum_i C_iJ_i$ 结果

> 本文展示 montopar 使用 `--tgv-ibp` 生成的四个 3e-3g orbit 结果。
> 所有系数、边界名称和 A/B 旋转关系均取自当前
> `workspaces/2loop/3e-3g/*/montopar.json`。
> 这里的 A、B 分别表示原始代表和 C4+1 旋转代表；本文不额外乘 orbit multiplicity。

---

## 1. 公共约定

记

$$
z\equiv z_1,\qquad
\lambda\equiv1-2\varepsilon.
$$

原始三维表示的公共前因子为

$$
\boxed{
C_{\rm raw}
=\frac{C_AC_F(4\pi\alpha_s)^2\Gamma(2-2\varepsilon)}
{64\pi^{4-2\varepsilon}}
(-2v_{14}+i0)^{2\varepsilon}.
}
$$

Contour IBP 边界核的公共前因子为

$$
\boxed{
C_\star
=\frac{C_AC_F(4\pi\alpha_s)^2\Gamma(1-2\varepsilon)}
{256\pi^{4-2\varepsilon}}
(-2v_{14}+i0)^{2\varepsilon}
=\frac{C_{\rm raw}}{4(1-2\varepsilon)}.
}
$$

对任一 orbit/group $G$，定义求导前星型核

$$
\mathcal K_G(\boldsymbol s)
=\int_{\alpha_i\geq0}d\alpha_1d\alpha_2d\alpha_3\,
\delta(1-\alpha_1-\alpha_2-\alpha_3)
(\alpha_1\alpha_2\alpha_3)^{-\varepsilon}
S_G(\boldsymbol s,\boldsymbol\alpha)^{-\lambda}.
$$

若 $B$ 是一个固定面或同边接触面，则

$$
\boxed{
J_G[B]
=\int_{D_G\cap B}d^2s\,\mathcal K_G(\boldsymbol s).
}
$$

因此下文的每个边界 $J$ 都是二维 Wilson 参数积分，并仍带三个
Feynman 参数及
$\delta(1-\alpha_1-\alpha_2-\alpha_3)$。这里 $d^2s$ 表示从原始迭代积分继承的
参数测度；对 contact 面不额外引入欧氏曲面元的 $\sqrt2$。

---

## 2. Orbit 01：$(3,4,4)$ 与 $(4,1,1)$

### 2.1 A 组 $(a,b,c)=(3,4,4)$

有序域为

$$
0\leq s_{41}\leq s_{42}\leq1,\qquad 0\leq s_{31}\leq1.
$$

IBP 导数向量为

$$
\boldsymbol c_{01,A}=(0,-z,z).
$$

最终结果为

$$
\boxed{
\begin{aligned}
I_{01,A}
={}&zC_\star J_{01,A}[s_{41}=0]
+zC_\star J_{01,A}[s_{42}=1]\\
&-2zC_\star J_{01,A}[s_{41}=s_{42}].
\end{aligned}
}
$$

即

$$
\boxed{
I_{01,A}
=C_\star\left(
zJ_{41,0}+zJ_{42,1}-2zJ_{41=42}
\right).
}
$$

| $J_i$ | 边界类型 | montopar term |
|---|---|---|
| $J_{41,0}$ | cusp | `boundary_s41_0` |
| $J_{42,1}$ | cusp | `boundary_s42_1` |
| $J_{41=42}$ | contact | `boundary_s41_eq_s42` |

### 2.2 B 组 $(a,b,c)=(4,1,1)$

有序域为

$$
0\leq s_{11}\leq s_{12}\leq1,\qquad 0\leq s_{41}\leq1,
$$

且

$$
\boldsymbol c_{01,B}=(0,-1,1).
$$

因此

$$
\boxed{
I_{01,B}
=C_\star\left(
J_{11,0}+J_{12,1}-2J_{11=12}
\right).
}
$$

| $J_i$ | 边界类型 | montopar term |
|---|---|---|
| $J_{11,0}$ | cusp | `boundary_s11_0` |
| $J_{12,1}$ | cusp | `boundary_s12_1` |
| $J_{11=12}$ | contact | `boundary_s11_eq_s12` |

---

## 3. Orbit 02：$(3,3,4)$ 与 $(4,4,1)$

### 3.1 A 组 $(a,b,c)=(3,3,4)$

有序域为

$$
0\leq s_{31}\leq s_{32}\leq1,\qquad 0\leq s_{41}\leq1,
$$

且

$$
\boldsymbol c_{02,A}=(-z,z,0).
$$

最终结果为

$$
\boxed{
I_{02,A}
=C_\star\left(
zJ_{31,0}+zJ_{32,1}-2zJ_{31=32}
\right).
}
$$

| $J_i$ | 边界类型 | montopar term |
|---|---|---|
| $J_{31,0}$ | cusp | `boundary_s31_0` |
| $J_{32,1}$ | cusp | `boundary_s32_1` |
| $J_{31=32}$ | contact | `boundary_s31_eq_s32` |

与 02 手算文档的记号对应时，需注意该文档把外部的 $z$ 吸收到三个
边界积分中。若记

$$
J_i^{\rm(manual)}=zJ_i^{\rm(here)},
$$

则手算结果是

$$
I_{02,A}
=C_{\rm IBP}
\left(2J_{\rm diag}^{\rm(manual)}
-J_0^{\rm(manual)}-J_1^{\rm(manual)}\right),
\qquad
C_{\rm IBP}=-C_\star.
$$

若像本文一样不把 $z$ 吸入 $J_i$，相应的有效公共系数才是
$-zC_\star$。两种写法给出完全相同的 $C_iJ_i$。

### 3.2 B 组 $(a,b,c)=(4,4,1)$

有序域为

$$
0\leq s_{41}\leq s_{42}\leq1,\qquad 0\leq s_{11}\leq1,
$$

且

$$
\boldsymbol c_{02,B}=(-1,1,0).
$$

因此

$$
\boxed{
I_{02,B}
=C_\star\left(
J_{41,0}+J_{42,1}-2J_{41=42}
\right).
}
$$

| $J_i$ | 边界类型 | montopar term |
|---|---|---|
| $J_{41,0}$ | cusp | `boundary_s41_0` |
| $J_{42,1}$ | cusp | `boundary_s42_1` |
| $J_{41=42}$ | contact | `boundary_s41_eq_s42` |

---

## 4. Orbit 03：$(2,4,4)$ 与 $(3,1,1)$

### 4.1 A 组 $(a,b,c)=(2,4,4)$

有序域为

$$
0\leq s_{41}\leq s_{42}\leq1,\qquad 0\leq s_{21}\leq1,
$$

且

$$
\boldsymbol c_{03,A}=(0,1+z,-1-z).
$$

最终结果为

$$
\boxed{
\begin{aligned}
I_{03,A}
={}&-(1+z)C_\star J_{03,A}[s_{41}=0]\\
&-(1+z)C_\star J_{03,A}[s_{42}=1]\\
&+2(1+z)C_\star J_{03,A}[s_{41}=s_{42}].
\end{aligned}
}
$$

即

$$
\boxed{
I_{03,A}
=(1+z)C_\star
\left(-J_{41,0}-J_{42,1}+2J_{41=42}\right).
}
$$

| $J_i$ | 边界类型 | montopar term |
|---|---|---|
| $J_{41,0}$ | cusp | `boundary_s41_0` |
| $J_{42,1}$ | cusp | `boundary_s42_1` |
| $J_{41=42}$ | contact | `boundary_s41_eq_s42` |

### 4.2 B 组 $(a,b,c)=(3,1,1)$

有序域为

$$
0\leq s_{11}\leq s_{12}\leq1,\qquad 0\leq s_{31}\leq1,
$$

且导数向量仍为

$$
\boldsymbol c_{03,B}=(0,1+z,-1-z).
$$

因此

$$
\boxed{
I_{03,B}
=(1+z)C_\star
\left(-J_{11,0}-J_{12,1}+2J_{11=12}\right).
}
$$

| $J_i$ | 边界类型 | montopar term |
|---|---|---|
| $J_{11,0}$ | cusp | `boundary_s11_0` |
| $J_{12,1}$ | cusp | `boundary_s12_1` |
| $J_{11=12}$ | contact | `boundary_s11_eq_s12` |

---

## 5. Orbit 04：$(2,3,4)$ 与 $(3,4,1)$

04 的三个 Wilson 参数相互独立，因此没有 contact 面。IBP 结果由六个
cusp 面和一个三维 recoil 积分组成。

### 5.1 A 组 $(a,b,c)=(2,3,4)$

积分域为

$$
(s_{21},s_{31},s_{41})\in[0,1]^3,
$$

导数向量为

$$
\boldsymbol c_{04,A}=(z,-1-z,1).
$$

六个有向面的贡献给出

$$
\boxed{
\begin{aligned}
I_{04,A}={}&
-zC_\star J_{21,0}
+zC_\star J_{21,1}\\
&+(1+z)C_\star J_{31,0}
-(1+z)C_\star J_{31,1}\\
&-C_\star J_{41,0}
+C_\star J_{41,1}
+C_{\rm raw}J_{A,\rm recoil}.
\end{aligned}
}
$$

换言之，

$$
\boxed{
\begin{aligned}
I_{04,A}
=C_\star\Bigl[
&z(J_{21,1}-J_{21,0})
-(1+z)(J_{31,1}-J_{31,0})\\
&+(J_{41,1}-J_{41,0})
\Bigr]
+C_{\rm raw}J_{A,\rm recoil}.
\end{aligned}
}
$$

| $J_i$ | montopar term |
|---|---|
| $J_{21,0}$ | `boundary_s21_0` |
| $J_{21,1}$ | `boundary_s21_1` |
| $J_{31,0}$ | `boundary_s31_0` |
| $J_{31,1}$ | `boundary_s31_1` |
| $J_{41,0}$ | `boundary_s41_0` |
| $J_{41,1}$ | `boundary_s41_1` |
| $J_{A,\rm recoil}$ | `recoil` |

取

$$
x_1=s_{21},\qquad x_2=s_{31},\qquad x_3=s_{41},
$$

recoil 分子为

$$
\boxed{
\begin{aligned}
P_{A,\rm recoil}
=\frac12\Bigl[
&(1+z)\alpha_1\alpha_2(x_1-1)\\
&+z\alpha_1\alpha_3\bigl(1+(1+z)x_3\bigr)\\
&+z\alpha_2\alpha_3(1-x_2)
\Bigr].
\end{aligned}
}
$$

相应的三维积分为

$$
\boxed{
J_{A,\rm recoil}
=\int_{[0,1]^3}d^3x
\int_{\alpha_i\geq0}d^3\alpha\,
\delta(1-\textstyle\sum_i\alpha_i)
\frac{(\alpha_1\alpha_2\alpha_3)^{-\varepsilon}
P_{A,\rm recoil}}
{S_A^{\,2-2\varepsilon}}.
}
$$

### 5.2 B 组 $(a,b,c)=(3,4,1)$

积分域为

$$
(s_{31},s_{41},s_{11})\in[0,1]^3,
$$

导数向量为

$$
\boldsymbol c_{04,B}=(1,-1-z,z).
$$

因此

$$
\boxed{
\begin{aligned}
I_{04,B}={}&
-C_\star J_{31,0}
+C_\star J_{31,1}\\
&+(1+z)C_\star J_{41,0}
-(1+z)C_\star J_{41,1}\\
&-zC_\star J_{11,0}
+zC_\star J_{11,1}
+C_{\rm raw}J_{B,\rm recoil}.
\end{aligned}
}
$$

等价地，

$$
\boxed{
\begin{aligned}
I_{04,B}
=C_\star\Bigl[
&(J_{31,1}-J_{31,0})
-(1+z)(J_{41,1}-J_{41,0})\\
&+z(J_{11,1}-J_{11,0})
\Bigr]
+C_{\rm raw}J_{B,\rm recoil}.
\end{aligned}
}
$$

| $J_i$ | montopar term |
|---|---|
| $J_{31,0}$ | `boundary_s31_0` |
| $J_{31,1}$ | `boundary_s31_1` |
| $J_{41,0}$ | `boundary_s41_0` |
| $J_{41,1}$ | `boundary_s41_1` |
| $J_{11,0}$ | `boundary_s11_0` |
| $J_{11,1}$ | `boundary_s11_1` |
| $J_{B,\rm recoil}$ | `recoil` |

取

$$
x_1=s_{31},\qquad x_2=s_{41},\qquad x_3=s_{11},
$$

recoil 分子为

$$
\boxed{
\begin{aligned}
P_{B,\rm recoil}
=\frac12\Bigl[
&z(1+z)\alpha_1\alpha_2(x_1-1)\\
&+\alpha_1\alpha_3\bigl(z+(1+z)x_3\bigr)\\
&+z\alpha_2\alpha_3(1-x_2)
\Bigr].
\end{aligned}
}
$$

相应积分为

$$
\boxed{
J_{B,\rm recoil}
=\int_{[0,1]^3}d^3x
\int_{\alpha_i\geq0}d^3\alpha\,
\delta(1-\textstyle\sum_i\alpha_i)
\frac{(\alpha_1\alpha_2\alpha_3)^{-\varepsilon}
P_{B,\rm recoil}}
{S_B^{\,2-2\varepsilon}}.
}
$$

---

## 6. A 组每个 $J_i$ 的显式形式

为了直接看清每个积分，先定义统一的二维边界积分泛函

$$
\boxed{
\mathcal B[F]
\equiv
\int_0^1dx_1\int_0^1dx_2
\int_{\alpha_i\geq0}d\alpha_1d\alpha_2d\alpha_3\,
\delta(1-\alpha_1-\alpha_2-\alpha_3)
(\alpha_1\alpha_2\alpha_3)^{-\varepsilon}
F(x_1,x_2,\boldsymbol\alpha;z)^{-(1-2\varepsilon)}.
}
$$

不同 $J_i$ 中的 $x_1,x_2$ 都是该项自己的局部积分变量；其与原始 Wilson
参数的对应关系在每项后明确给出。

### 6.1 Orbit 01 A 的三个积分

#### $J_{41,0}$：$s_{41}=0$

局部变量为

$$
(x_1,x_2)=(s_{31},s_{42}).
$$

因此

$$
\boxed{
J_{41,0}^{(01,A)}
=\mathcal B\!\left[
z\,x_2\alpha_1\alpha_3(1-x_1)
\right].
}
$$

#### $J_{42,1}$：$s_{42}=1$

局部变量为

$$
(x_1,x_2)=(s_{31},s_{41}),
$$

所以

$$
\boxed{
J_{42,1}^{(01,A)}
=\mathcal B\!\left[
z(1-x_1)\alpha_1
(\alpha_3+x_2\alpha_2)
\right].
}
$$

#### $J_{41=42}$：同边 contact 面

令

$$
s_{41}=s_{42}=x_2,\qquad x_1=s_{31}.
$$

则

$$
\boxed{
J_{41=42}^{(01,A)}
=\mathcal B\!\left[
z\,x_2(1-x_1)\alpha_1(\alpha_2+\alpha_3)
\right].
}
$$

因此 Orbit 01 A 的完整结果可以无歧义地写成

$$
\boxed{
I_{01,A}
=zC_\star\left(
J_{41,0}^{(01,A)}
+J_{42,1}^{(01,A)}
-2J_{41=42}^{(01,A)}
\right).
}
$$

### 6.2 Orbit 02 A 的三个积分

#### $J_{31,0}$：$s_{31}=0$

局部变量为

$$
(x_1,x_2)=(s_{32},s_{41}),
$$

故

$$
\boxed{
J_{31,0}^{(02,A)}
=\mathcal B\!\left[
z\,x_2\alpha_3
\bigl(\alpha_1+\alpha_2(1-x_1)\bigr)
\right].
}
$$

#### $J_{32,1}$：$s_{32}=1$

局部变量为

$$
(x_1,x_2)=(s_{31},s_{41}),
$$

故

$$
\boxed{
J_{32,1}^{(02,A)}
=\mathcal B\!\left[
z\,x_2\alpha_1\alpha_3(1-x_1)
\right].
}
$$

#### $J_{31=32}$：同边 contact 面

令

$$
s_{31}=s_{32}=x_1,\qquad x_2=s_{41}.
$$

则

$$
\boxed{
J_{31=32}^{(02,A)}
=\mathcal B\!\left[
z\,x_2(1-x_1)\alpha_3(\alpha_1+\alpha_2)
\right].
}
$$

从而

$$
\boxed{
I_{02,A}
=zC_\star\left(
J_{31,0}^{(02,A)}
+J_{32,1}^{(02,A)}
-2J_{31=32}^{(02,A)}
\right).
}
$$

### 6.3 Orbit 03 A 的三个积分

#### $J_{41,0}$：$s_{41}=0$

局部变量为

$$
(x_1,x_2)=(s_{21},s_{42}).
$$

对应的边界多项式为

$$
F_{41,0}^{(03,A)}
=\alpha_1\left[
\alpha_2(1-x_1)
+\alpha_3\bigl((1-x_1)(1-x_2)+zx_1x_2\bigr)
\right],
$$

所以

$$
\boxed{
J_{41,0}^{(03,A)}
=\mathcal B\!\left[F_{41,0}^{(03,A)}\right].
}
$$

#### $J_{42,1}$：$s_{42}=1$

局部变量为

$$
(x_1,x_2)=(s_{21},s_{41}).
$$

定义

$$
F_{42,1}^{(03,A)}
=\alpha_1\alpha_2
\bigl((1-x_1)(1-x_2)+zx_1x_2\bigr)
+z x_1\alpha_1\alpha_3.
$$

于是

$$
\boxed{
J_{42,1}^{(03,A)}
=\mathcal B\!\left[F_{42,1}^{(03,A)}\right].
}
$$

#### $J_{41=42}$：同边 contact 面

令

$$
s_{41}=s_{42}=x_2,\qquad x_1=s_{21}.
$$

此时

$$
F_{41=42}^{(03,A)}
=\alpha_1(\alpha_2+\alpha_3)
\bigl((1-x_1)(1-x_2)+zx_1x_2\bigr),
$$

因此

$$
\boxed{
J_{41=42}^{(03,A)}
=\mathcal B\!\left[F_{41=42}^{(03,A)}\right].
}
$$

Orbit 03 A 的完整组合为

$$
\boxed{
I_{03,A}
=(1+z)C_\star\left(
-J_{41,0}^{(03,A)}
-J_{42,1}^{(03,A)}
+2J_{41=42}^{(03,A)}
\right).
}
$$

### 6.4 Orbit 04 A 的七个积分

六个边界积分仍由 $\mathcal B[F]$ 定义。为了避免把不同边界上的局部变量混淆，
下面逐项给出映射与多项式。

#### $J_{21,0}$：$s_{21}=0$

$$
(x_1,x_2)=(s_{31},s_{41}),
$$

$$
F_{21,0}^{(04,A)}
=\alpha_1\alpha_3(1-x_2)
+x_1\alpha_1\alpha_2
+zx_2\alpha_2\alpha_3(1-x_1),
$$

$$
\boxed{
J_{21,0}^{(04,A)}
=\mathcal B\!\left[F_{21,0}^{(04,A)}\right].
}
$$

#### $J_{21,1}$：$s_{21}=1$

$$
(x_1,x_2)=(s_{31},s_{41}),
$$

$$
F_{21,1}^{(04,A)}
=zx_2\alpha_3
\bigl(\alpha_1+\alpha_2(1-x_1)\bigr),
$$

$$
\boxed{
J_{21,1}^{(04,A)}
=\mathcal B\!\left[F_{21,1}^{(04,A)}\right].
}
$$

#### $J_{31,0}$：$s_{31}=0$

$$
(x_1,x_2)=(s_{21},s_{41}),
$$

$$
F_{31,0}^{(04,A)}
=\alpha_1\alpha_3
\bigl((1-x_1)(1-x_2)+zx_1x_2\bigr)
+zx_2\alpha_2\alpha_3,
$$

$$
\boxed{
J_{31,0}^{(04,A)}
=\mathcal B\!\left[F_{31,0}^{(04,A)}\right].
}
$$

#### $J_{31,1}$：$s_{31}=1$

$$
(x_1,x_2)=(s_{21},s_{41}),
$$

$$
F_{31,1}^{(04,A)}
=\alpha_1\left[
\alpha_2(1-x_1)
+\alpha_3\bigl((1-x_1)(1-x_2)+zx_1x_2\bigr)
\right],
$$

$$
\boxed{
J_{31,1}^{(04,A)}
=\mathcal B\!\left[F_{31,1}^{(04,A)}\right].
}
$$

#### $J_{41,0}$：$s_{41}=0$

$$
(x_1,x_2)=(s_{21},s_{31}),
$$

$$
F_{41,0}^{(04,A)}
=\alpha_1(1-x_1)(\alpha_3+x_2\alpha_2),
$$

$$
\boxed{
J_{41,0}^{(04,A)}
=\mathcal B\!\left[F_{41,0}^{(04,A)}\right].
}
$$

#### $J_{41,1}$：$s_{41}=1$

$$
(x_1,x_2)=(s_{21},s_{31}),
$$

$$
F_{41,1}^{(04,A)}
=x_2\alpha_1\alpha_2(1-x_1)
+zx_1\alpha_1\alpha_3
+z\alpha_2\alpha_3(1-x_2),
$$

$$
\boxed{
J_{41,1}^{(04,A)}
=\mathcal B\!\left[F_{41,1}^{(04,A)}\right].
}
$$

#### $J_{A,\rm recoil}$

定义三维泛函

$$
\boxed{
\mathcal R[P,S]
\equiv
\int_0^1dx_1\int_0^1dx_2\int_0^1dx_3
\int_{\alpha_i\geq0}d^3\alpha\,
\delta(1-\textstyle\sum_i\alpha_i)
(\alpha_1\alpha_2\alpha_3)^{-\varepsilon}
P\,S^{-(2-2\varepsilon)}.
}
$$

其中

$$
\begin{aligned}
S_A={}&
\alpha_1\alpha_2x_2(1-x_1)\\
&+\alpha_1\alpha_3
\left[1-x_1-x_3+(1+z)x_1x_3\right]\\
&+z\alpha_2\alpha_3x_3(1-x_2),
\end{aligned}
$$

$$
\begin{aligned}
P_{A,\rm recoil}
=\frac12\Bigl[
&(1+z)\alpha_1\alpha_2(x_1-1)\\
&+z\alpha_1\alpha_3\bigl(1+(1+z)x_3\bigr)\\
&+z\alpha_2\alpha_3(1-x_2)
\Bigr].
\end{aligned}
$$

所以

$$
\boxed{
J_{A,\rm recoil}^{(04,A)}
=\mathcal R[P_{A,\rm recoil},S_A].
}
$$

将七项放在一起，

$$
\boxed{
\begin{aligned}
I_{04,A}={}&C_\star\Bigl[
-zJ_{21,0}^{(04,A)}
+zJ_{21,1}^{(04,A)}
+(1+z)J_{31,0}^{(04,A)}\\
&-(1+z)J_{31,1}^{(04,A)}
-J_{41,0}^{(04,A)}
+J_{41,1}^{(04,A)}
\Bigr]
+C_{\rm raw}J_{A,\rm recoil}^{(04,A)}.
\end{aligned}
}
$$

---

## 7. 四个 orbit 的紧凑汇总

$$
\boxed{
\begin{aligned}
I_{01,A}
&=zC_\star(J_{41,0}+J_{42,1}-2J_{41=42}),\\
I_{01,B}
&=C_\star(J_{11,0}+J_{12,1}-2J_{11=12}),\\[3pt]
I_{02,A}
&=zC_\star(J_{31,0}+J_{32,1}-2J_{31=32}),\\
I_{02,B}
&=C_\star(J_{41,0}+J_{42,1}-2J_{41=42}),\\[3pt]
I_{03,A}
&=(1+z)C_\star(-J_{41,0}-J_{42,1}+2J_{41=42}),\\
I_{03,B}
&=(1+z)C_\star(-J_{11,0}-J_{12,1}+2J_{11=12}),\\[3pt]
I_{04,A}
&=C_\star\!\left[
z(J_{21,1}-J_{21,0})
-(1+z)(J_{31,1}-J_{31,0})
+J_{41,1}-J_{41,0}
\right]
+C_{\rm raw}J_{A,\rm recoil},\\
I_{04,B}
&=C_\star\!\left[
J_{31,1}-J_{31,0}
-(1+z)(J_{41,1}-J_{41,0})
+z(J_{11,1}-J_{11,0})
\right]
+C_{\rm raw}J_{B,\rm recoil}.
\end{aligned}
}
$$

这里不同 orbit/group 中下标相同的 $J$ 不表示同一个函数；它始终使用本节所属
orbit/group 的 $S_G$。完整的局部 $s\to x$ 映射和展开后的边界多项式保存在对应
`montopar.json` 的 `J_ibp.terms` 或 `J_ibp_rotated.terms` 中。

---

## 8. 与机器输出的等价性

对退化 orbit 01、02、03，原始复杂分子被完全消除：

$$
\boxed{
C_{\rm raw}J_{\rm raw}
=\sum_{i=1}^{3}C_iJ_i.
}
$$

对非退化 orbit 04，

$$
\boxed{
C_{\rm raw}J_{\rm raw}
=\sum_{i=1}^{6}C_{\star,i}J_{{\rm face},i}
+C_{\rm raw}J_{\rm recoil}.
}
$$

montopar 对 A/B 两组均验证

$$
\boxed{
P_{\rm num}=P_{\rm boundary}+P_{\rm recoil},
}
$$

当前四个 orbit 的 `operator_identity`、`boundary_orientation` 和
`pnum_identity` 证书均为 `true`。
