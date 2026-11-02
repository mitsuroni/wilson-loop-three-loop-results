# 3e-3g-04 非退化 TGV 的 Contour IBP 推导

> **目标**：对 `3e-3g-04` 的三条不同 Wilson 边构型建立严格的
> Contour IBP 分解，将原始 15 项 $P_{\rm num}$ 改写为六个二维 cusp
> 边界积分与一个具有紧凑分子的三维 recoil 积分。
>
> **约定**：$D=4-2\varepsilon$，$\lambda=D-3=1-2\varepsilon$，
> $A_{\rm ref}=-2v_{14}+i0$。本文沿用 montopar 当前的
> $e^{-i\ell\cdot X}$、TGV 动量差和 $i0$ convention。

---

## 1. 一般非退化 TGV 的精确分解

三胶子顶点与 Wilson 线切向量的缩并为

$$
\boxed{
\begin{aligned}
V(a,b,c)
&=v_{ab}\,v_c\!\cdot\!(\ell_2-\ell_1)
 +v_{bc}\,v_a\!\cdot\!(\ell_3-\ell_2)\\
&\quad
 +v_{ca}\,v_b\!\cdot\!(\ell_1-\ell_3).
\end{aligned}
}
$$

利用 $\ell_1+\ell_2+\ell_3=0$，可以写成

$$
\boxed{V=V_{\rm contour}+V_{\rm recoil},}
$$

其中

$$
\boxed{
V_{\rm contour}
=-\left(
v_{bc}v_a\!\cdot\!\ell_1
+v_{ca}v_b\!\cdot\!\ell_2
+v_{ab}v_c\!\cdot\!\ell_3
\right),
}
$$

$$
\boxed{
V_{\rm recoil}
=-2\left(
v_{ab}v_c\!\cdot\!\ell_1
+v_{bc}v_a\!\cdot\!\ell_2
+v_{ca}v_b\!\cdot\!\ell_3
\right).
}
$$

这一分解可以直接验证。定义

$$
K=v_{ab}v_c+v_{bc}v_a+v_{ca}v_b,
$$

则

$$
\boxed{
V-V_{\rm contour}-V_{\rm recoil}
=K\!\cdot\!(\ell_1+\ell_2+\ell_3)=0.
}
$$

因此该恒等式不依赖特殊四边形运动学；动量守恒已经足够。

由

$$
v_a\!\cdot\!\ell_1\longrightarrow i\partial_{s_1},\qquad
v_b\!\cdot\!\ell_2\longrightarrow i\partial_{s_2},\qquad
v_c\!\cdot\!\ell_3\longrightarrow i\partial_{s_3},
$$

得到

$$
\boxed{
\mathcal O_{\rm contour}
=-i\left(
v_{bc}\partial_{s_1}
+v_{ca}\partial_{s_2}
+v_{ab}\partial_{s_3}
\right).
}
$$

按照 $\ell_i\to i\partial_{X_i}$，recoil 算子为

$$
\boxed{
\mathcal O_{\rm recoil}
=-2i\left(
v_{ab}v_c\!\cdot\!\partial_{X_1}
+v_{bc}v_a\!\cdot\!\partial_{X_2}
+v_{ca}v_b\!\cdot\!\partial_{X_3}
\right).
}
$$

这里的负号由 $e^{-i\ell\cdot X}$ convention 固定；不能改成 $+2i$。

---

## 2. 公共星型核与两种 prefactor

定义求导前的星型核

$$
\boxed{
\mathcal K(\boldsymbol{x})
=\int_{\alpha_i\geq0}d\alpha_1d\alpha_2d\alpha_3\,
\delta(1-\alpha_1-\alpha_2-\alpha_3)
(\alpha_1\alpha_2\alpha_3)^{-\varepsilon}
S(\boldsymbol{x})^{-\lambda}.
}
$$

当前 montopar 的求导后表示使用

$$
\boxed{
C_{\rm raw}
=\frac{C_AC_F(4\pi\alpha_s)^2\Gamma(2-2\varepsilon)}
{64\pi^{4-2\varepsilon}}
(-2v_{14}+i0)^{2\varepsilon}.
}
$$

对求导前的边界星型核，定义统一的

$$
\boxed{
C_\star
=\frac{C_AC_F(4\pi\alpha_s)^2\Gamma(1-2\varepsilon)}
{256\pi^{4-2\varepsilon}}
(-2v_{14}+i0)^{2\varepsilon}
=\frac{C_{\rm raw}}{4(1-2\varepsilon)}.
}
$$

非退化 TGV 的六个边界项使用 $C_\star$，recoil 项继续使用
$C_{\rm raw}$。这样无需把不同的 Gamma 函数 convention 塞进同一个 $J$。

---

## 3. A 组：$(a,b,c)=(2,3,4)$

取

$$
x_1=s_{21},\qquad x_2=s_{31},\qquad x_3=s_{41},
$$

三个变量均独立地位于 $[0,1]$。

运动学比值为

$$
\frac{v_{bc}}{v_{14}}=\frac{v_{34}}{v_{14}}=z,
\qquad
\frac{v_{ca}}{v_{14}}=\frac{v_{42}}{v_{14}}=-(1+z),
\qquad
\frac{v_{ab}}{v_{14}}=\frac{v_{23}}{v_{14}}=1.
$$

因此无量纲 Contour IBP 算子是

$$
\boxed{
\mathcal D_A
=z\partial_{x_1}-(1+z)\partial_{x_2}+\partial_{x_3}.
}
$$

### 3.1 A 组的 $S$ 多项式

$$
\boxed{
\begin{aligned}
S_A={}&
\alpha_1\alpha_2x_2(1-x_1)\\
&+\alpha_1\alpha_3
\left[1-x_1-x_3+(1+z)x_1x_3\right]\\
&+z\alpha_2\alpha_3x_3(1-x_2).
\end{aligned}
}
$$

### 3.2 六个边界多项式

$$
\boxed{
\begin{aligned}
S_{A,1-}\equiv S_A|_{x_1=0}
&=\alpha_1\alpha_2x_2
 +\alpha_1\alpha_3(1-x_3)
 +z\alpha_2\alpha_3x_3(1-x_2),\\
S_{A,1+}\equiv S_A|_{x_1=1}
&=z\alpha_3x_3[\alpha_1+\alpha_2(1-x_2)],\\[4pt]
S_{A,2-}\equiv S_A|_{x_2=0}
&=\alpha_1\alpha_3
 [1-x_1-x_3+(1+z)x_1x_3]
 +z\alpha_2\alpha_3x_3,\\
S_{A,2+}\equiv S_A|_{x_2=1}
&=\alpha_1\alpha_2(1-x_1)
 +\alpha_1\alpha_3
 [1-x_1-x_3+(1+z)x_1x_3],\\[4pt]
S_{A,3-}\equiv S_A|_{x_3=0}
&=\alpha_1(1-x_1)(\alpha_2x_2+\alpha_3),\\
S_{A,3+}\equiv S_A|_{x_3=1}
&=\alpha_1\alpha_2x_2(1-x_1)
 +z\alpha_1\alpha_3x_1
 +z\alpha_2\alpha_3(1-x_2).
\end{aligned}
}
$$

定义

$$
J_{A,i\pm}
=\int_0^1d^2x
\int d^3\alpha\,\delta(1-\textstyle\sum\alpha)
(\alpha_1\alpha_2\alpha_3)^{-\varepsilon}
S_{A,i\pm}^{-\lambda},
$$

其中 $d^2x$ 表示除被固定的 $x_i$ 外其余两个 Wilson 参数。

边界部分为

$$
\boxed{
J_{A,\rm bdry}
=z(J_{A,1+}-J_{A,1-})
-(1+z)(J_{A,2+}-J_{A,2-})
+(J_{A,3+}-J_{A,3-}).
}
$$

### 3.3 A 组的 recoil 分子

定义

$$
P_{A,\rm bdry}=-\frac14\mathcal D_AS_A.
$$

从当前完整分子中减去该部分：

$$
P_{A,\rm recoil}=P_{A,\rm num}-P_{A,\rm bdry}.
$$

符号化简得到

$$
\boxed{
\begin{aligned}
P_{A,\rm recoil}=\frac12\Bigl[&
(1+z)\alpha_1\alpha_2(x_1-1)\\
&+z\alpha_1\alpha_3\bigl(1+(1+z)x_3\bigr)\\
&+z\alpha_2\alpha_3(1-x_2)
\Bigr].
\end{aligned}
}
$$

展开后它只有 9 项，而当前 $P_{A,\rm num}$ 有 15 项；紧凑形式只有三组。

相应积分为

$$
\boxed{
J_{A,\rm recoil}
=\int_0^1d^3x\int d^3\alpha\,
\delta(1-\textstyle\sum\alpha)
\frac{P_{A,\rm recoil}}
{(\alpha_1\alpha_2\alpha_3)^\varepsilon
 S_A^{2-2\varepsilon}}.
}
$$

### 3.4 A 组最终的 $\sum_iC_iJ_i$

$$
\boxed{
I_{04,A}
=C_\star J_{A,\rm bdry}
+C_{\rm raw}J_{A,\rm recoil}.
}
$$

完全展开为七项：

$$
\boxed{
\begin{aligned}
I_{04,A}={}&
(zC_\star)J_{A,1+}
-(zC_\star)J_{A,1-}\\
&-[(1+z)C_\star]J_{A,2+}
+[(1+z)C_\star]J_{A,2-}\\
&+C_\star J_{A,3+}
-C_\star J_{A,3-}
+C_{\rm raw}J_{A,\rm recoil}.
\end{aligned}
}
$$

---

## 4. B 组：$(a,b,c)=(3,4,1)$

B 组是 A 组作 C4+1 轮换：

$$
x_1=s_{31},\qquad x_2=s_{41},\qquad x_3=s_{11}.
$$

运动学比值为

$$
\frac{v_{41}}{v_{14}}=1,
\qquad
\frac{v_{13}}{v_{14}}=-(1+z),
\qquad
\frac{v_{34}}{v_{14}}=z,
$$

故

$$
\boxed{
\mathcal D_B
=\partial_{x_1}-(1+z)\partial_{x_2}+z\partial_{x_3}.
}
$$

B 组二次型可以紧凑写成

$$
\boxed{
\begin{aligned}
S_B={}&
z\alpha_1\alpha_2x_2(1-x_1)\\
&+\alpha_1\alpha_3
[z(1-x_1)(1-x_3)+x_1x_3]\\
&+\alpha_2\alpha_3x_3(1-x_2).
\end{aligned}
}
$$

六个边界由 $S_{B,i\pm}=S_B|_{x_i=0,1}$ 生成，边界组合为

$$
\boxed{
J_{B,\rm bdry}
=(J_{B,1+}-J_{B,1-})
-(1+z)(J_{B,2+}-J_{B,2-})
+z(J_{B,3+}-J_{B,3-}).
}
$$

recoil 分子为

$$
\boxed{
\begin{aligned}
P_{B,\rm recoil}=\frac12\Bigl[&
z(1+z)\alpha_1\alpha_2(x_1-1)\\
&+\alpha_1\alpha_3[z+(1+z)x_3]\\
&+z\alpha_2\alpha_3(1-x_2)
\Bigr].
\end{aligned}
}
$$

最终结果是

$$
\boxed{
I_{04,B}
=C_\star J_{B,\rm bdry}
+C_{\rm raw}J_{B,\rm recoil}.
}
$$

B 组应由同一个通用算法在旋转后的腿和运动学映射上自动生成，不应在代码中保存第二套手写公式。

---

## 5. 被积函数级等价证书

对 $G=A,B$，定义

$$
P_{G,\rm bdry}=-\frac14\mathcal D_GS_G.
$$

由于

$$
\partial_{x_i}S_G^{-\lambda}
=-\lambda(\partial_{x_i}S_G)S_G^{-(\lambda+1)},
$$

以及

$$
C_\star=\frac{C_{\rm raw}}{4\lambda},
$$

边界导数部分在未积分时满足

$$
C_\star\mathcal D_GS_G^{-\lambda}
=C_{\rm raw}P_{G,\rm bdry}S_G^{-(2-2\varepsilon)}.
$$

再加上

$$
\boxed{
P_{G,\rm num}
=P_{G,\rm bdry}+P_{G,\rm recoil},
}
$$

便得到

$$
\boxed{
C_{\rm raw}J_{G,\rm raw}
=C_\star J_{G,\rm bdry}
+C_{\rm raw}J_{G,\rm recoil}.
}
$$

本文已对 A/B 两组使用当前 montopar 生成的 $S_G$ 和 $P_{G,\rm num}$ 做符号检查，差均严格化为零。

---

## 6. 为什么不继续声称“全部化成边界”

非退化情形与同边退化情形不同。三个 recoil 方向分别是

$$
v_c\!\cdot\!\partial_{X_1},\qquad
v_a\!\cdot\!\partial_{X_2},\qquad
v_b\!\cdot\!\partial_{X_3},
$$

它们不是相应附着边的切向导数。平移不变性可以在 $X_i$ 之间移动一个公共方向导数，
但一般不能把三个方向同时变成
$v_a\cdot\partial_{X_1}$、$v_b\cdot\partial_{X_2}$、$v_c\cdot\partial_{X_3}$。

若继续进行参数空间 syzygy IBP，可能进一步缩短 $P_{\rm recoil}$，但也可能引入
$1/\alpha_i$、奇异边界或更多分支。因此第一版实现应以本文的
“六个边界 + 一个 recoil”作为完成条件。

---

## 7. 对代码实现的约束

1. 非退化结果固定产生 7 个 term：6 个 `cusp_boundary` 和 1 个 `recoil`。
2. 边界项引用 $C_\star$，分母幂为 $1-2\varepsilon$，只保留两个 Wilson 参数。
3. recoil 项引用 $C_{\rm raw}$，分母幂为 $2-2\varepsilon$，保留三个 Wilson 参数。
4. 必须保存边界的 `fixed_parameter`、方向和系数，不能只保存展开后的 integrand。
5. A/B 组使用同一个生成器；B 组只能来自 `rotate_orbit_full()` 后重新计算。
6. 每次生成必须验证
   $P_{\rm num}-P_{\rm bdry}-P_{\rm recoil}=0$。
7. `J_raw` 必须保留，作为兼容输出和等价性基准。

---

## 8. 最终结论

`3e-3g-04` 的安全 Contour IBP 终点不是纯边界表达式，而是

$$
\boxed{
I_{04,G}
=\sum_{k=1}^{6}C_{G,k}^{\rm cusp}J_{G,k}^{\rm cusp}
+C_{\rm raw}J_{G,\rm recoil},
\qquad G=A,B.
}
$$

该表示严格保留原积分，同时把复杂部分从完整 15 项分子缩减为一个三组紧凑的 recoil 分子；
其余六项都降为二维 Wilson 参数积分。这是适合第一版 montopar 改造的非退化输出。
