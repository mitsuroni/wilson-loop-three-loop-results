# 3e-3g-02 的 Contour IBP 手动推导

> **目标**：以退化 TGV orbit `3e-3g-02` 为例，在引入复杂的
> $\mathcal O(Q^{3-D})$ 分子之前，沿 Wilson 线参数做分部积分，并把结果写成
> $\sum_i C_iJ_i$。
>
> **适用约定**：$D=4-2\varepsilon$，$A_{\rm ref}=-2v_{14}+i0$，
> $v_{34}=zv_{14}$。本文沿用 montopar 当前的 TGV 动量差、Fourier 变换和
> $i0$ convention。闭式参数积分先在收敛域解析计算，再按维数正规化作解析延拓。

---

## 1. Orbit、腿标号与有序积分域

`3e-3g-02` 的三条 TGV 腿为

$$
(a,b,c)=(3,3,4),
$$

并取

$$
s_1\equiv s_{31},\qquad
s_2\equiv s_{32},\qquad
u\equiv s_{41}.
$$

上游 orbit 给出的路径排序是

$$
0\leq s_1\leq s_2\leq 1,
\qquad 0\leq u\leq 1.
$$

三个附着点为

$$
X_1=V_3+s_1v_3,\qquad
X_2=V_3+s_2v_3,\qquad
X_3=V_4+uv_4.
$$

由于 $X_1$ 与 $X_2$ 位于同一条类光边，

$$
\Delta_{12}=-(X_1-X_2)^2+i0=0.
$$

另外两个无量纲距离为

$$
\Delta'_{23}=zu(1-s_2),\qquad
\Delta'_{31}=zu(1-s_1),
$$

所以星型积分的无量纲二次型是

$$
\boxed{
S(s_1,s_2,u)
=zu\alpha_3\bigl[\alpha_1(1-s_1)+\alpha_2(1-s_2)\bigr].
}
$$

---

## 2. 退化 TGV 分子是严格的全导数

采用当前 montopar 对应的三胶子顶点缩并：

$$
\begin{aligned}
V(a,b,c)
&=v_{ab}\,v_c\!\cdot\!(\ell_2-\ell_1)
 +v_{bc}\,v_a\!\cdot\!(\ell_3-\ell_2)
 +v_{ca}\,v_b\!\cdot\!(\ell_1-\ell_3).
\end{aligned}
$$

对 $(a,b,c)=(3,3,4)$，有 $v_{33}=0$，从而

$$
\begin{aligned}
V(3,3,4)
&=v_{34}v_3\!\cdot\!(\ell_3-\ell_2)
  +v_{43}v_3\!\cdot\!(\ell_1-\ell_3)\\
&=v_{34}v_3\!\cdot\!(\ell_1-\ell_2).
\end{aligned}
$$

相位满足

$$
v_3\!\cdot\!\ell_1\,e^{-i\ell_1\cdot X_1}
=i\partial_{s_1}e^{-i\ell_1\cdot X_1},
\qquad
v_3\!\cdot\!\ell_2\,e^{-i\ell_2\cdot X_2}
=i\partial_{s_2}e^{-i\ell_2\cdot X_2},
$$

因此

$$
\boxed{
V(3,3,4)=i v_{34}(\partial_{s_1}-\partial_{s_2}).
}
$$

这里没有 recoil remainder；这是退化构型相对于非退化构型的根本简化。

---

## 3. 求导前的星型核与公共前因子

定义

$$
\lambda\equiv D-3=1-2\varepsilon
$$

以及求导前的无量纲星型核

$$
\boxed{
\mathcal K(s_1,s_2,u)
=\int_{\alpha_i\geq0}d\alpha_1d\alpha_2d\alpha_3\,
\delta(1-\alpha_1-\alpha_2-\alpha_3)
(\alpha_1\alpha_2\alpha_3)^{-\varepsilon}
S(s_1,s_2,u)^{-\lambda}.
}
$$

将颜色、顶点、三个传播子的 Fourier 变换、$iv_{34}$ 和
$Q=A_{\rm ref}S$ 合并，有

$$
I_{3e\text{-}3g\text{-}02}
=C_{\rm IBP}\,z
\int_0^1du\int_0^1ds_2\int_0^{s_2}ds_1
(\partial_{s_1}-\partial_{s_2})\mathcal K(s_1,s_2,u),
$$

其中

$$
\boxed{
C_{\rm IBP}
=-\frac{C_AC_F(4\pi\alpha_s)^2\Gamma(1-2\varepsilon)}
{256\pi^{4-2\varepsilon}}
(-2v_{14}+i0)^{2\varepsilon}.
}
$$

负号来自

$$
\frac{v_{34}}{A_{\rm ref}}
=\frac{zv_{14}}{-2v_{14}+i0}
\longrightarrow -\frac{z}{2},
$$

其中箭头表示按当前的标度提取 convention 保留公共 $i0$ 相位。

### 3.1 与现有 $C_{\rm raw}J_{\rm raw}$ 的局部核对

现有 montopar 使用

$$
C_{\rm raw}
=\frac{C_AC_F(4\pi\alpha_s)^2\Gamma(2-2\varepsilon)}
{64\pi^{4-2\varepsilon}}
(-2v_{14}+i0)^{2\varepsilon}.
$$

因此

$$
\boxed{
C_{\rm IBP}=-\frac{C_{\rm raw}}{4(1-2\varepsilon)}.
}
$$

另一方面，

$$
(\partial_{s_1}-\partial_{s_2})S
=zu\alpha_3(\alpha_2-\alpha_1),
$$

故

$$
(\partial_{s_1}-\partial_{s_2})S^{-\lambda}
=-\lambda zu\alpha_3(\alpha_2-\alpha_1)S^{-(2-2\varepsilon)}.
$$

于是

$$
C_{\rm IBP}\,z\,[-\lambda zu\alpha_3(\alpha_2-\alpha_1)]
=C_{\rm raw}\left[\frac14z^2u\alpha_3(\alpha_2-\alpha_1)\right].
$$

右端方括号正是当前 `3e-3g-02` 的 $P_{\rm num}$。这在被积函数层面验证了
Contour IBP convention 与现有 $C_{\rm raw}J_{\rm raw}$ 完全一致。

---

## 4. 在有序三角域上做 Contour IBP

令 $F(s_1,s_2)=\mathcal K(s_1,s_2,u)$。先处理 $s_1$ 导数：

$$
\int_0^1ds_2\int_0^{s_2}ds_1\,\partial_{s_1}F
=\int_0^1dt\,[F(t,t)-F(0,t)].
$$

对 $s_2$ 导数使用 Leibniz 公式：

$$
\frac{d}{dt}\int_0^t ds\,F(s,t)
=F(t,t)+\int_0^t ds\,\partial_tF(s,t),
$$

从而

$$
-\int_0^1dt\int_0^t ds\,\partial_tF(s,t)
=-\int_0^1ds\,F(s,1)+\int_0^1dt\,F(t,t).
$$

两部分相加得到

$$
\boxed{
\int_0^1ds_2\int_0^{s_2}ds_1
(\partial_{s_1}-\partial_{s_2})F
=2\int_0^1dr\,F(r,r)
-\int_0^1dt\,F(0,t)
-\int_0^1ds\,F(s,1).
}
$$

这里出现三个边界：

| 边界 | 几何类型 | 系数 |
|---|---|---:|
| $s_1=s_2$ | 同边两腿重合的接触边界 | $+2$ |
| $s_1=0$ | 边 3 起点 cusp | $-1$ |
| $s_2=1$ | 边 3 终点 cusp | $-1$ |

特别地，$s_1=s_2$ 的对角边界出现两次，不能在代码中当作普通超立方体面而忽略。

---

## 5. $\sum_i C_iJ_i$ 的参数积分形式

定义三个不含外部 prefactor 的边界积分：

$$
\boxed{
\begin{aligned}
J_{\rm diag}
&=z\int_0^1dr\int_0^1du\,\mathcal K(r,r,u),\\
J_0
&=z\int_0^1dt\int_0^1du\,\mathcal K(0,t,u),\\
J_1
&=z\int_0^1ds\int_0^1du\,\mathcal K(s,1,u).
\end{aligned}
}
$$

则 IBP 后的最终线性组合为

$$
\boxed{
I_{3e\text{-}3g\text{-}02}
=C_{\rm IBP}\bigl(2J_{\rm diag}-J_0-J_1\bigr)
=\sum_{i\in\{\rm diag,0,1\}}C_iJ_i,
}
$$

其中

$$
\boxed{
C_{\rm diag}=2C_{\rm IBP},\qquad
C_0=-C_{\rm IBP},\qquad
C_1=-C_{\rm IBP}.
}
$$

这已经是 montopar 应输出的最小正确 IBP 结构：一个 orbit 对应多个积分项，且每项有独立系数和边界类型。

---

## 6. 三个边界积分的闭式计算

以下记

$$
B(x,y)=\frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)}.
$$

### 6.1 对角接触边界 $s_1=s_2=r$

此时

$$
S(r,r,u)=zu\alpha_3(\alpha_1+\alpha_2)(1-r).
$$

令 $\beta=\alpha_1+\alpha_2=1-\alpha_3$，则

$$
\int_0^\beta d\alpha_1\,
[\alpha_1(\beta-\alpha_1)]^{-\varepsilon}
=\beta^{1-2\varepsilon}B(1-\varepsilon,1-\varepsilon).
$$

$\beta$ 幂与 $S^{-\lambda}$ 中的 $\beta^{-1+2\varepsilon}$ 精确抵消，剩余

$$
\int_0^1d\alpha_3\,\alpha_3^{-1+\varepsilon}
=\frac1\varepsilon.
$$

再利用

$$
\int_0^1dr\,(1-r)^{-1+2\varepsilon}
=\int_0^1du\,u^{-1+2\varepsilon}
=\frac1{2\varepsilon},
$$

得到

$$
\boxed{
J_{\rm diag}
=\frac{z^{2\varepsilon}}{4\varepsilon^3}
B(1-\varepsilon,1-\varepsilon).
}
$$

### 6.2 起点 cusp 边界 $s_1=0$

此时

$$
S(0,t,u)=zu\alpha_3[\alpha_1+\alpha_2(1-t)].
$$

作变量替换

$$
\alpha_1=(1-\alpha_3)r,\qquad
\alpha_2=(1-\alpha_3)(1-r),
$$

则 $\alpha_3$ 和 $u$ 积分分别给出 $1/\varepsilon$ 与 $1/(2\varepsilon)$。
剩余的 $t$ 积分为

$$
\int_0^1dt\,[r+(1-r)(1-t)]^{-1+2\varepsilon}
=\frac{1-r^{2\varepsilon}}{2\varepsilon(1-r)}.
$$

因此

$$
\boxed{
J_0
=\frac{z^{2\varepsilon}}{4\varepsilon^3}
\left[
B(1-\varepsilon,-\varepsilon)
-B(1+\varepsilon,-\varepsilon)
\right].
}
$$

两个单独的 Beta 函数按解析延拓理解；它们的差对应原始收敛组合
$1-r^{2\varepsilon}$。

### 6.3 终点 cusp 边界 $s_2=1$

此时

$$
S(s,1,u)=zu\alpha_1\alpha_3(1-s).
$$

$s$ 和 $u$ 积分给出 $1/(4\varepsilon^2)$。参数单纯形积分是

$$
\begin{aligned}
&\int d^3\alpha\,\delta(1-\textstyle\sum\alpha)
\alpha_1^{-1+\varepsilon}
\alpha_2^{-\varepsilon}
\alpha_3^{-1+\varepsilon}\\
&\qquad
=\frac{\Gamma(\varepsilon)^2\Gamma(1-\varepsilon)}
{\Gamma(1+\varepsilon)}
=\frac{\Gamma(\varepsilon)\Gamma(1-\varepsilon)}{\varepsilon}.
\end{aligned}
$$

故

$$
\boxed{
J_1
=\frac{z^{2\varepsilon}}{4\varepsilon^3}
\Gamma(\varepsilon)\Gamma(1-\varepsilon).
}
$$

---

## 7. 完全闭式的最终结果

代入三个边界积分：

$$
\boxed{
\begin{aligned}
I_{3e\text{-}3g\text{-}02}
&=C_{\rm IBP}\frac{z^{2\varepsilon}}{4\varepsilon^3}
\Bigl[
2B(1-\varepsilon,1-\varepsilon)\\
&\qquad
-B(1-\varepsilon,-\varepsilon)
+B(1+\varepsilon,-\varepsilon)
-\Gamma(\varepsilon)\Gamma(1-\varepsilon)
\Bigr].
\end{aligned}
}
$$

利用

$$
B(1+\varepsilon,-\varepsilon)
=-\Gamma(\varepsilon)\Gamma(1-\varepsilon),
$$

也可以写成

$$
\boxed{
\begin{aligned}
I_{3e\text{-}3g\text{-}02}
&=C_{\rm IBP}\frac{z^{2\varepsilon}}{4\varepsilon^3}
\Bigl[
2B(1-\varepsilon,1-\varepsilon)
-B(1-\varepsilon,-\varepsilon)\\
&\qquad
-2\Gamma(\varepsilon)\Gamma(1-\varepsilon)
\Bigr].
\end{aligned}
}
$$

本文保留第一种形式作为主要结果，因为它逐项对应
$2J_{\rm diag}-J_0-J_1$，更适合作为代码测试 oracle。

---

## 8. 对 montopar 代码设计的直接启示

### 8.1 一个 IBP 结果必须允许多个 term

当前 `J_raw` 是单个 integrand；Contour IBP 后至少需要：

```json
{
  "J_ibp": {
    "method": "contour_ibp",
    "status": "complete",
    "terms": [
      {
        "term_id": "diag",
        "kind": "contact_boundary",
        "coefficient": "2",
        "constraints": {"s31": "s32"}
      },
      {
        "term_id": "s31_low",
        "kind": "cusp_boundary",
        "coefficient": "-1",
        "fixed_parameters": {"s31": "0"}
      },
      {
        "term_id": "s32_high",
        "kind": "cusp_boundary",
        "coefficient": "-1",
        "fixed_parameters": {"s32": "1"}
      }
    ]
  }
}
```

### 8.2 必须保留边界来源

只保存化简后的 integrand 不够。代码必须保存：

- 原始有序域；
- 边界方程及朝向；
- 边界 multiplicity；
- 是 cusp 面还是同边接触面；
- 从哪个全导数项产生。

否则很容易漏掉本文系数为 $2$ 的对角面。

### 8.3 `C_raw` 与 `C_IBP` 不能混用

求导后的旧表示含 $\Gamma(2-2\varepsilon)$ 和 derivative factor；IBP 边界核含
$\Gamma(1-2\varepsilon)$。建议每个表示保存自己的完整 prefactor，并额外保存

$$
C_{\rm IBP}/C_{\rm raw}=-\frac1{4(1-2\varepsilon)}
$$

作为验证信息。

### 8.4 建议保存两层 IBP 结果

建议同时输出：

1. `boundary_integrals`：本文第 5 节的未积分 $J_i$，供通用下游继续处理；
2. `closed_form`：本文第 6--7 节的 Beta/Gamma 结果，作为快速通道和回归测试。

### 8.5 最低限度的验证证书

对每个退化 orbit，生成器至少应验证：

$$
V-i v_{34}(\partial_{s_1}-\partial_{s_2})=0,
$$

$$
C_{\rm IBP}\,z\,(\partial_{s_1}-\partial_{s_2})\mathcal K
=C_{\rm raw}\,P_{\rm num}S^{-(2-2\varepsilon)},
$$

以及有序域的 Stokes/Leibniz 恒等式

$$
\int_{0\leq s_1\leq s_2\leq1}
(\partial_{s_1}-\partial_{s_2})F
=2\int F_{s_1=s_2}-\int F_{s_1=0}-\int F_{s_2=1}.
$$

这三个证书分别覆盖顶点代数、prefactor convention 和边界组合，是实现时不可省略的三层检查。

---

## 9. 结论

对于 `3e-3g-02`，Contour IBP 确实完全消除了原来的多项式分子
$P_{\rm num}$。但正确结果不是“两个 cusp 项”，而是

$$
\boxed{
I_{3e\text{-}3g\text{-}02}
=2C_{\rm IBP}J_{\rm diag}
-C_{\rm IBP}J_0
-C_{\rm IBP}J_1.
}
$$

对角接触边界是有序 Wilson 线积分的必要组成部分。这个样例说明，未来的 montopar
扩展应被设计为“带边界元数据的积分线性组合生成器”，而不是把单个 `J_raw.integrand`
简单替换成另一个单个 integrand。
