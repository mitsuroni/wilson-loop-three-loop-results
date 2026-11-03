---
title: "纯交换 Wilson 积分｜空集与单对边族的统一理论"
category: "通用理论"
status: "理论汇编"
tags:
  - wilson-loop
  - pure-exchange
  - generating-function
aliases:
  - "空集与单对边 [1] 族：全圈结构、统一母积分与低圈解析系数"
catalog_updated: 2026-09-15
---

# 纯交换 Wilson 积分｜空集与单对边族的统一理论

> **内容定位**：整理本轮关于空集、任意圈阶单对边 `[1]` 族、统一母积分 $\mathcal T_{u;v}$、它与原 $M,N$ 记号的关系，以及 `2e-02`、`4e-02`、`6e-04` 的新记号重算。本文是前述讨论与已有推导文件的汇编，不是对所有圈阶、所有图重新进行的一次独立计算。
>
> **适用范围**：一般闭合光锥四边形上的纯传播子交换；无内部三胶子、四胶子顶点，无自能插入。本文的 $n$ 是这一范围内的圈阶，也等于直接交换传播子的条数。
>
> **表达目标**：从通用 $J_G$ 出发保留端点排序，统一生成积分；在固定 Laurent 阶数下以常数系数的基函数表示结果，尽量不把 $z$ 或 $\ln z$ 留作基函数外的系数。
>
> **结论边界**：统一母族的覆盖、固定阶数的可计算性、给定表达中的基函数数目、最小独立基底维数是不同问题。下文不把母核覆盖数量解释成已经证明的主积分维数。

## 目录

- [1. 统一约定与物理归一化](#sec-1)
- [2. 两类问题共用的有序字积分](#sec-2)
- [3. 空集族：从通用初始 J 到全圈分解](#sec-3)
- [4. 空集的常数空间、极点与 Laurent 系数](#sec-4)
- [5. 空集的独立基函数及其与 M 的关系](#sec-5)
- [6. 单对边 [1] 族的通用初始 J](#sec-6)
- [7. 单对边族的对称化与奇偶通道](#sec-7)
- [8. 单对边族的解析空间上界](#sec-8)
- [9. 从四个端点区间压缩到两个字](#sec-9)
- [10. 只定义一次的统一母族与有限圈阶覆盖](#sec-10)
- [11. M、N 与统一母积分之间的转换](#sec-11)
- [12. 重复 A 子族与统一系数函数 tau](#sec-12)
- [13. 三张图从原始排序到精确 I](#sec-13)
- [14. 三张图直到有限项的统一解析结果](#sec-14)
- [15. 覆盖空间、关系空间与最小基底](#sec-15)
- [16. 确定性计算流程与核验边界](#sec-16)
- [17. 结论对照与尚未完成的任务](#sec-17)
- [附录 A. 当前 18 个 [1] 图的母核归约](#appendix-a)
- [附录 B. 系数及符号转换速查](#appendix-b)
- [附录 C. 来源与整理说明](#appendix-c)

---

<a id="sec-1"></a>
## 1. 统一约定与物理归一化

### 1.1 一般闭合光锥四边形

四个切向量满足

$$
v_i^2=0,\qquad v_1+v_2+v_3+v_4=0.
$$

一般运动学已经给出

$$
\boxed{
v_{12}=v_{34},\qquad
v_{14}=v_{23},\qquad
v_{13}=v_{24},\qquad
v_{12}+v_{13}+v_{14}=0.
}
$$

以 $v_{14}$ 为公共尺度，令

$$
z=\frac{v_{12}}{v_{14}},\qquad \ell=\ln z,
$$

则

$$
\boxed{
(v_{12},v_{13},v_{14},v_{23},v_{24},v_{34})
=v_{14}\bigl(z,-(1+z),1,1,-(1+z),z\bigr).
}
$$

这里不使用 $v_3=-v_1,\ v_4=-v_2$。此前已区分：额外采用对边反向条件会使非退化情况下 $z=-1$，不能把它作为一般 $z$ 推导的前提。上述尺度关系与源文件的运动学表一致。[S1, §1]

先在 $z>0$、合适的 $\operatorname{Re}\varepsilon>0$ 收敛区域建立积分等式，然后向 $\varepsilon=0$ 附近解析延拓。复数运动学最后沿原有 $i0$ 处方延拓，不在变换途中分别重新选择对数支路。

### 1.2 原始 C × J 与最终 I 的区别

源文件将单图贡献写成 $C_GJ_G$。本文的 $I_G$ 专指按下式提取出的、指定色通道的无量纲对称化系数：

$$
\frac{\mathrm{sym}_G}{2}\,C_G
\left[J_G(z)+z^{n\varepsilon}J_G(1/z)\right]
=
C_A^{n-1}C_F\left(\frac{\alpha_s}{\pi}\right)^n
\left[(-2\pi v_{14})^{n\varepsilon}+(-2\pi v_{12})^{n\varepsilon}\right]
I_{G,\mathrm{raw}}.
$$

该等式按约定在目标 $C_A^{n-1}C_F$ 通道内理解；不能左侧保留全部颜色，右侧却只保留这一通道。

设 $a_G$ 是完整颜色表达投影到 $C_A^{n-1}C_F$ 后的系数，定义

$$
\boxed{\kappa_G=\frac{\mathrm{sym}_G a_G}{2(-2)^n}.}
$$

源文件的完整前因子已经包含 $(-2)^{-n}$。把它计入 $\kappa_G$ 后，后续不得重复乘入。[S1, §6b；S2–S4 的 `C.combined`]

于是

$$
I_{G,\mathrm{raw}}
=\kappa_G\Gamma(1-\varepsilon)^n
\frac{J_G(z)+z^{n\varepsilon}J_G(1/z)}{1+z^{n\varepsilon}}.
$$

沿用本轮的去 Euler 常数约定：

$$
\boxed{
Z(\varepsilon)=e^{-\gamma_E\varepsilon}\Gamma(1-\varepsilon),\qquad
I_G=e^{-n\gamma_E\varepsilon}I_{G,\mathrm{raw}}.
}
$$

因此全文默认

$$
\boxed{
I_G=\kappa_GZ(\varepsilon)^n
\frac{J_G(z)+z^{n\varepsilon}J_G(1/z)}{1+z^{n\varepsilon}}.
}
\tag{1.1}
$$

这里

$$
Z(\varepsilon)^n
=\exp\left[n\sum_{k\ge2}\frac{\zeta_k}{k}\varepsilon^k\right].
$$

这不等于在最终表达式中任意丢掉 $\gamma_E$；两种方案由明确的指数因子联系。[S5, 范围与约定]

### 1.3 色投影不是大 N 截断

使用用户指定的

$$
C_A^2-2C_AC_F=1,\qquad \frac1{C_A}=C_A-2C_F,
$$

将颜色因子转成对应圈阶的齐次色基，再取 $C_A^{n-1}C_F$ 系数。三圈使用

$$
\{C_A^2C_F,\ C_AC_F^2,\ C_F^3\}.
$$

后文所有 $a_G$、$\kappa_G$ 都遵循这一投影；不先固定 $N$，也不将这项操作替换为大 $N$ 近似。

---

<a id="sec-2"></a>
## 2. 两类问题共用的有序字积分

### 2.1 邻边传播子的端点分解

沿循环方向，邻边尺度为

$$
(q_1,q_2,q_3,q_4)=(z,1,z,1).
$$

一条位于 $i\to i+1$ 通道的邻边传播子，连同 Lorentz 缩并，贡献

$$
q_i\,[q_i(1-u)v]^{-1+\varepsilon}
=q_i^\varepsilon(1-u)^{\varepsilon-1}v^{\varepsilon-1}.
$$

定义

$$
\boxed{
A_\varepsilon(x)=x^{\varepsilon-1},\qquad
B_\varepsilon(x)=(1-x)^{\varepsilon-1}.
}
$$

每条邻边在前一条边留下一个 $B$ 端点，在后一条边留下一个 $A$ 端点。循环通道 $41$ 按这一循环方向处理，而不是仅按边编号大小处理。[S1, §5c]

### 2.2 有序字、区间与正规化

若 $w=w_1\cdots w_m$，$w_j\in\{A,B\}$，定义

$$
\boxed{
F_w(a,b;\varepsilon)
=\int_{a<x_1<\cdots<x_m<b}
\prod_{j=1}^m f_{w_j}(x_j)\,d^mx,
\qquad F_\varnothing=1.
}
\tag{2.1}
$$

字序对应从小到大的端点参数顺序。递归为

$$
\boxed{
F_{wc}(a,b;\varepsilon)
=\int_a^b f_c(x)F_w(a,x;\varepsilon)\,dx.
}
\tag{2.2}
$$

定义

$$
\boxed{\widehat F_w(a,b;\varepsilon)=\varepsilon^{|w|}F_w(a,b;\varepsilon).}
\tag{2.3}
$$

例如

$$
\widehat F_{A^m}(0,s;\varepsilon)=\frac{s^{m\varepsilon}}{m!},\qquad
\widehat F_{B^m}(0,1;\varepsilon)=\frac1{m!}.
$$

$\varnothing$ 表示空字；图的“空集族”则表示没有对边传播子。两种空的含义不同。

### 2.3 三条结构恒等式

同区间乘积满足 shuffle：

$$
\widehat F_u(a,b)\widehat F_v(a,b)
=\sum_{w\in u\shuffle v}\widehat F_w(a,b).
\tag{2.4}
$$

求和保留重数，例如 $A\shuffle A$ 有两个相同的 $AA$。

路径合成满足

$$
\widehat F_w(a,c)
=\sum_{w=uv}\widehat F_u(a,b)\widehat F_v(b,c).
\tag{2.5}
$$

设 $\rho(w)$ 是“先反转字序，再交换 $A,B$”，则反射满足

$$
\widehat F_w(0,1-t)=\widehat F_{\rho(w)}(t,1).
\tag{2.6}
$$

这些是前述统一母族构造使用的固定规则。它们不依赖具体图号。[D2；S6]

---

<a id="sec-3"></a>
## 3. 空集族：从通用初始 J 到全圈分解

### 3.1 抽象初始积分

令 $n_i$ 为循环邻边 $i\leftrightarrow i+1$ 的传播子数。空集族满足

$$
\sum_{i=1}^4n_i=n.
$$

边 $i$ 的端点数为 $m_i=n_{i-1}+n_i$，且

$$
\sum_i m_i=2n,
\qquad
\Delta_G=\prod_{i=1}^4
\{0<s_{i,1}<\cdots<s_{i,m_i}<1\}.
$$

任意空集图的初始积分为

$$
\boxed{
J_G(z,\varepsilon)
=\int_{\Delta_G}d^{2n}\mathbf s\;
\prod_{p=1}^n q_{i_p}
\left[q_{i_p}(1-u_p)v_p\right]^{-1+\varepsilon}.
}
\tag{3.1}
$$

仅给圈数 $n$ 不足以确定积分；端点顺序仍由 $G$ 指定。

### 3.2 提出全部运动学尺度

令

$$
r_G=n_1+n_3,
$$

即 $12$、$34$ 两个带 $z$ 尺度通道的传播子数。逐端点分解后

$$
J_G=z^{r_G\varepsilon}
\int_{\Delta_G}d^{2n}\mathbf s\;
\prod_{p=1}^n(1-u_p)^{\varepsilon-1}v_p^{\varepsilon-1}.
$$

四条边上的字记为 $w_1,w_2,w_3,w_4$。不同边之间不再有参数耦合，因此

$$
\boxed{
J_G(z,\varepsilon)=z^{r_G\varepsilon}\mathcal B_G(\varepsilon),\qquad
\mathcal B_G(\varepsilon)=\prod_{i=1}^4F_{w_i}(0,1;\varepsilon).
}
\tag{3.2}
$$

由此直接得到

$$
z\partial_zJ_G=r_G\varepsilon J_G.
$$

**空集的排序复杂性全部位于不依赖 $z$ 的常数积分中；运动学只剩一个尺度幂。** 这是后续空集函数空间结论的起点，不依赖任何低圈例子。[D0]

### 3.3 对称化后的精确 I

将式 (3.2) 代入式 (1.1)：

$$
\boxed{
I_G(z,\varepsilon)
=\kappa_GZ(\varepsilon)^n\mathcal B_G(\varepsilon)
\frac{\cosh[(r_G-n/2)\varepsilon\ell]}
{\cosh(n\varepsilon\ell/2)}.
}
\tag{3.3}
$$

因此 $I_G(z)=I_G(1/z)$，每个 Laurent 系数的运动学部分只可能是 $\ell^2$ 的有限多项式。

在 $r_G=0$ 或 $r_G=n$ 时，双曲函数比值为 $1$，单图的最终 $I_G$ 不依赖 $\ell$。一般 $r_G$ 下仍可依赖偶次对数。

这说的是逐阶系数。有限 $\varepsilon$ 的完整双曲函数比值通常不是有限的对数多项式。

---

<a id="sec-4"></a>
## 4. 空集的常数空间、极点与 Laurent 系数

### 4.1 任意字的端点极点分离

将非空字写成

$$
w=A^a\,v\,B^b,
$$

其中 $a$ 是开头连续 $A$ 的个数，$b$ 是末尾连续 $B$ 的个数。中间字 $v$ 要么为空，要么以 $B$ 开头、以 $A$ 结尾。

若 $|v|=d>0$，直接完成前后端点段的积分，得到

$$
\boxed{
F_w(0,1;\varepsilon)
=\frac1{a!b!\varepsilon^{a+b}}
\int_{0<t_1<\cdots<t_d<1}
 t_1^{a\varepsilon}(1-t_d)^{b\varepsilon}
\prod_{j=1}^df_{v_j}(t_j)\,d^dt.
}
\tag{4.1}
$$

剩下的有序积分在 $\varepsilon=0$ 收敛。令

$$
S(t)=a\ln t_1+b\ln(1-t_d)
+\sum_{v_j=A}\ln t_j+\sum_{v_j=B}\ln(1-t_j),
$$

则可以在提取端点极点后展开：

$$
F_w
=\frac1{a!b!\varepsilon^{a+b}}
\sum_{k\ge0}\frac{\varepsilon^k}{k!}
\int_{\Delta_d}\prod_j\omega_{v_j}(t_j)\,S(t)^k,
\quad
\omega_A=\frac{dt}{t},\quad\omega_B=\frac{dt}{1-t}.
\tag{4.2}
$$

利用

$$
\ln t=-\int_t^1\frac{dx}{x},\qquad
\ln(1-t)=-\int_0^t\frac{dx}{1-x},
$$

并对新增参数按大小分区，式 (4.2) 变成有限个收敛的 $0,1$ 字母迭代积分。

### 4.2 多重 zeta 与普通 zeta 必须区分

定义多重 zeta 常数代数

$$
\mathcal Z_{\mathrm{MZV}}
=\mathbb Q[\zeta(s_1,\ldots,s_d)],
\quad
\zeta(s_1,\ldots,s_d)
=\sum_{k_1>\cdots>k_d\ge1}
\frac1{k_1^{s_1}\cdots k_d^{s_d}},
\quad s_1\ge2.
$$

前述端点迭代积分结构给出

$$
\boxed{
F_w(0,1;\varepsilon),\ \mathcal B_G(\varepsilon)
\in\mathcal Z_{\mathrm{MZV}}((\varepsilon)).
}
\tag{4.3}
$$

这里保留此前讨论的范围：一般定理给出的是多重 zeta 常数，不自动证明任意圈、任意阶都可进一步化为普通 $\zeta(k)$ 的乘积。[D0]

若 $v=\varnothing$，即 $w=A^aB^b$，则有更特殊的 Gamma 比值通式：

$$
\boxed{
F_{A^aB^b}(0,1;\varepsilon)
=\frac1{a!b!\varepsilon^{a+b}}
\frac{\Gamma(1+a\varepsilon)\Gamma(1+b\varepsilon)}
{\Gamma(1+(a+b)\varepsilon)}.
}
\tag{4.4}
$$

此类字只需普通 zeta 展开。$a=0$ 或 $b=0$ 时得到单字母阶乘公式；空字按 $F_\varnothing=1$ 理解。

### 4.3 极点阶数

按本轮的端点分离论证，若目标色投影非零，单图实际极点阶数为

$$
\boxed{P_G=\sum_{i=1}^4(a_i+b_i)\le2n.}
\tag{4.5}
$$

剩余收敛积分在正实初始域有非零系数，式 (3.3) 的双曲比值在 $\varepsilon=0$ 为 $1$。因此这个结论针对单图成立；图间求和可以进一步抵消极点。

### 4.4 全部 Laurent 系数的有限生成式

定义

$$
H_G(\varepsilon)
=\varepsilon^{2n}Z(\varepsilon)^n\mathcal B_G(\varepsilon)
=\sum_{j\ge0}h_{G,j}\varepsilon^j,
\qquad h_{G,j}\in\mathcal Z_{\mathrm{MZV}},
$$

以及

$$
\mathcal C_{n,r}(x)
=\frac{\cosh[(r-n/2)x]}{\cosh(nx/2)}
=\sum_{m\ge0}c_{n,r,m}x^{2m}.
$$

令 $t=r(n-r)$，最前几项为

$$
c_{n,r,0}=1,\qquad c_{n,r,1}=-\frac t2,\qquad
c_{n,r,2}=\frac{t(n^2+t)}{24}.
$$

写

$$
I_G=\sum_{p=-2n}^{\infty}I_G^{(p)}\varepsilon^p,
$$

则

$$
\boxed{
I_G^{(p)}(z)
=\kappa_G\sum_{m=0}^{\lfloor(2n+p)/2\rfloor}
c_{n,r_G,m}\,h_{G,2n+p-2m}\,\ell^{2m}.
}
\tag{4.6}
$$

所以

$$
\boxed{I_G^{(p)}\in\mathcal Z_{\mathrm{MZV}}[\ell^2].}
\tag{4.7}
$$

用实际极点阶数可以加强次数限制：

$$
\deg_\ell I_G^{(p)}
\le2\left\lfloor\frac{P_G+p}{2}\right\rfloor.
$$

每项可按自然总权重 $2n+p$ 组织；这一计数不以不同权重的数值常数代数独立为前提。

---

<a id="sec-5"></a>
## 5. 空集的独立基函数及其与 M 的关系

### 5.1 空集可以有自己的单字母积分族

定义

$$
H_0(z)=1,\qquad H_{j+1}(z)=\int_1^z\frac{dt}{t}H_j(t).
$$

则 $H_j=\ell^j/j!$。取 crossing 偶部分

$$
\boxed{E_m(z)=H_{2m}(z)=\frac{\ell^{2m}}{(2m)!}.}
\tag{5.1}
$$

它满足

$$
E_m(1/z)=E_m(z),\qquad
E_mE_k=\binom{2m+2k}{2m}E_{m+k},
$$

$$
(z\partial_z)^2E_m=E_{m-1}\quad(m\ge1).
$$

因此空集结果可以独立写成

$$
I_G^{(p)}(z)=\sum_m d_{G,p,m}\,E_m(z),\qquad
d_{G,p,m}\in\mathcal Z_{\mathrm{MZV}}.
$$

### 5.2 与原 M 的表达等价，不是拓扑等价

采用原 $M$ 定义中的最低矩，前述讨论已建立

$$
\boxed{M_{0,1}(z)=-\frac12\ell^2-3\zeta_2.}
\tag{5.2}
$$

因此

$$
\ell^{2m}=(-2M_{0,1}-6\zeta_2)^m,
$$

$$
\boxed{
\mathcal Z_{\mathrm{MZV}}[\ell^2]
=\mathcal Z_{\mathrm{MZV}}[M_{0,1}],\qquad
E_m=\frac{(-2M_{0,1}-6\zeta_2)^m}{(2m)!}.
}
\tag{5.3}
$$

这允许 $M_{0,1}$ 的有限乘积，不需要 $M$ 分母，也不需要新增 $N$。将空集对数换写成 $M_{0,1}$，不表示图中出现了对边传播子；只是两个函数表达之间的恒等转换。[D0]

**空集的全圈结论**：运动学部分完全由偶次对数解决；一般常数空间是 MZV；可以使用自己的 $E_m$ 族，也可以统一输出为 $M_{0,1}$ 的多项式。

---

<a id="sec-6"></a>
## 6. 单对边 [1] 族的通用初始 J

### 6.1 唯一对边核与邻边尺度

经循环重标号，将唯一对边连接固定在 $24$，保留两个端点 $t,s$：

$$
\boxed{D_z(t,s)=(1-t)(1-s)+zts.}
$$

其余 $n-1$ 条传播子全部连接邻边，满足

$$
\sum_i n_i=n-1,\qquad r_G=n_1+n_3.
$$

通用初始积分为

$$
\boxed{
J_G=-(1+z)
\int_{\Delta_G}d^{2n}\mathbf s\;
D_z(t,s)^{-1+\varepsilon}
\prod_{p=1}^{n-1}q_{i_p}
[q_{i_p}(1-u_p)v_p]^{-1+\varepsilon}.
}
\tag{6.1}
$$

提出邻边尺度后

$$
\boxed{
J_G=-(1+z)z^{r_G\varepsilon}
\int_{\Delta_G}d^{2n}\mathbf s\;
D_z(t,s)^{-1+\varepsilon}
\prod_{\text{邻边端点 }e}f_{w_e}(s_e).
}
\tag{6.2}
$$

与空集不同，$D_z$ 保留了两个边参数之间的耦合，所以剩余积分不能全部变成常数。[D1；S1, §5c]

### 6.2 带标记 O 的有序字分解

用 $O$ 标记唯一对边的端点，四条边写成

$$
w_1,\qquad u_2Ov_2,\qquad w_3,\qquad u_4Ov_4.
$$

所有 $A/B$ 字母总数为 $2n-2$。定义

$$
\mathcal B_G=F_{w_1}(0,1)F_{w_3}(0,1),
$$

$$
P_G(t)=F_{u_2}(0,t)F_{v_2}(t,1),\qquad
Q_G(s)=F_{u_4}(0,s)F_{v_4}(s,1).
$$

省略的参数均为 $\varepsilon$。定义带权核

$$
\boxed{
\mathscr K_z[w]
=-(1+z)\int_0^1dt\,ds\;
D_z(t,s)^{-1+\varepsilon}w(t,s;\varepsilon).
}
\tag{6.3}
$$

整个 `[1]` 族有精确分解

$$
\boxed{
J_G=z^{r_G\varepsilon}\mathcal B_G(\varepsilon)
\mathscr K_z[P_G(t;\varepsilon)Q_G(s;\varepsilon)].
}
\tag{6.4}
$$

正规化所有字后

$$
\boxed{
J_G=\frac{z^{r_G\varepsilon}}{\varepsilon^{2n-2}}
\widehat{\mathcal B}_G(\varepsilon)
\mathscr K_z[\widehat P_G\widehat Q_G].
}
\tag{6.5}
$$

**裸 `2e-02` 是这个结构中两端权重均为 $1$ 的情况，不是所有 `[1]` 图的完整生成式。**

### 6.3 极点上界与实际端点计数

对固定 $z>0$，$D_z$ 仅可能在 $(t,s)=(0,1),(1,0)$ 的角点为零。例如第一角点附近

$$
D_z\sim zx+y,\qquad x=t,\quad y=1-s.
$$

局部二维积分 $dx\,dy/(zx+y)$ 乘有限对数幂仍可积。正规化后的邻边系数只有对数型端点奇性，不新增该角点的 $\varepsilon=0$ 极点。因此

$$
\boxed{\operatorname{pole\ order}I_G\le2n-2.}
\tag{6.6}
$$

本轮的进一步端点分析给出：在四条**原始完整边字**中，将 $O$ 视为中断符，分别数开头连续 $A$ 和末尾连续 $B$ 的个数 $a_i,b_i$。目标色系数非零时，单图实际极点阶数为

$$
\boxed{P_G=\sum_i(a_i+b_i).}
\tag{6.7}
$$

这个计数要在原始端点排序上使用，不能先将图拆成若干母核后分别计数、再忽略母核组合中的抵消。它也不排除不同图的线性组合进一步消极点。[D1]

---

<a id="sec-7"></a>
## 7. 单对边族的对称化与奇偶通道

### 7.1 核反射与尺度指数

直接有

$$
D_{1/z}(t,s)=z^{-1}D_z(1-t,1-s),
$$

$$
\boxed{
\mathscr K_{1/z}[w(t,s)]
=z^{-\varepsilon}\mathscr K_z[w(1-t,1-s)].
}
\tag{7.1}
$$

代入归一化公式，得到

$$
\boxed{
\begin{aligned}
I_G={}&\frac{\kappa_GZ^n\mathcal B_G}{1+z^{n\varepsilon}}\,
\mathscr K_z\bigl[
 z^{r_G\varepsilon}P_G(t)Q_G(s)\\
&\hspace{39mm}
+z^{(n-1-r_G)\varepsilon}P_G(1-t)Q_G(1-s)
\bigr].
\end{aligned}
}
\tag{7.2}
$$

第二项的尺度指数是 $n-1-r_G$，不是空集的 $n-r_G$。差别来自唯一对边核的 $z^{-\varepsilon}$。

### 7.2 中心化核与宇称投影

记

$$
\widetilde{\mathscr K}_z[w]=z^{-\varepsilon/2}\mathscr K_z[w],
\qquad
\Phi_z(t,s)=\ln D_z(t,s)-\frac\ell2.
$$

那么

$$
\widetilde{\mathscr K}_z[w]
=-(1+z)\int_0^1\frac{dt\,ds}{D_z(t,s)}
 e^{\varepsilon\Phi_z(t,s)}w(t,s;\varepsilon).
$$

令 $\mathsf Rw(t,s)=w(1-t,1-s)$，对
$\widehat w_G=\widehat P_G\widehat Q_G$ 分解

$$
\widehat w_{G,\pm}=\frac12(\widehat w_G\pm\mathsf R\widehat w_G),
$$

$$
\mathcal V_G^\pm
=\widehat{\mathcal B}_G\,
\widetilde{\mathscr K}_z[\widehat w_{G,\pm}],
\qquad
\mathcal V_G^\pm(1/z)=\pm\mathcal V_G^\pm(z).
$$

定义

$$
\delta_G=r_G-\frac{n-1}{2}.
$$

最终通用式为

$$
\boxed{
I_G=\frac{\kappa_GZ(\varepsilon)^n}{\varepsilon^{2n-2}}
\frac{
\cosh(\delta_G\varepsilon\ell)\mathcal V_G^+
+\sinh(\delta_G\varepsilon\ell)\mathcal V_G^-
}{\cosh(n\varepsilon\ell/2)}.
}
\tag{7.3}
$$

这与空集的全圈双曲公式相对应，但多出一个非平凡奇通道：

$$
\ell^{2m}\times\text{偶函数},\qquad
\ell^{2m+1}\times\text{奇函数}.
$$

两者整体都为 crossing 偶。故不能由 $I(z)=I(1/z)$ 推出其任意表达中都只有偶次显式 $\ell$。[D1]

### 7.3 固定 Laurent 阶数的生成公式

写

$$
Z(\varepsilon)^n=\sum_{a\ge0}g_{n,a}\varepsilon^a,\qquad
\mathcal V_G^\pm=\sum_{j\ge0}V_{G,j}^\pm\varepsilon^j,
$$

$$
\frac{\cosh(\delta_Gx)}{\cosh(nx/2)}
=\sum_{m\ge0}C_{G,m}x^{2m},\qquad
\frac{\sinh(\delta_Gx)}{\cosh(nx/2)}
=\sum_{m\ge0}S_{G,m}x^{2m+1}.
$$

令 $h=2n-2+p$，则

$$
\boxed{
\begin{aligned}
I_G^{(p)}={}&\kappa_G
\sum_{a+j+2m=h}g_{n,a}C_{G,m}\ell^{2m}V_{G,j}^+\\
&+\kappa_G
\sum_{a+j+2m+1=h}g_{n,a}S_{G,m}\ell^{2m+1}V_{G,j}^-.
\end{aligned}
}
\tag{7.4}
$$

每个目标阶数只有有限个求和项。若 $h<0$，系数为零。

---

<a id="sec-8"></a>
## 8. 单对边族的解析空间上界

### 8.1 邻边权重只使用 0、1 字母

经端点正规化后，$\widehat P_G,\widehat Q_G$ 的系数为
$dt/t$、$dt/(1-t)$ 生成的迭代积分。经过乘积展开，待求对象是有限个

$$
-(1+z)\int_0^1dt\,ds\,
\frac{H_{\mathbf u}(t)H_{\mathbf v}(s)\ln^qD_z(t,s)}{D_z(t,s)}.
\tag{8.1}
$$

此处 $H_{\mathbf u},H_{\mathbf v}$ 的字母是 $0,1$；图特有常数属于 MZV 空间。

### 8.2 两次积分的线性字母结构

令

$$
\rho(s)=\frac{1-s}{1-(1+z)s},\qquad a(z)=\frac1{1+z}.
$$

有精确的对数型二形式

$$
\boxed{
-\frac{(1+z)\,dt\,ds}{D_z(t,s)}
=-\frac{dt}{t-\rho(s)}\frac{ds}{s-a(z)}.
}
\tag{8.2}
$$

第一次积分的字母为 $0,1,\rho(s)$。其端点值与线性系数为

$$
D_z(0,s)=1-s,\qquad
D_z(1,s)=zs,\qquad
\partial_tD_z=(1+z)s-1.
$$

所以第二次积分只需 $0,1,1/(1+z)$。外部变量的字母碰撞位置为

$$
z=0,\qquad z=-1,\qquad z=\infty.
$$

前述讨论据此使用线性可约超对数求值及端点正规化，得到

$$
\boxed{
[\varepsilon^p]I_G
\in
\left(\mathcal Z_{\mathrm{MZV}}\mathcal H_{\{0,-1\}}\right)^+
_{\mathrm{wt}\le2n+p}.
}
\tag{8.3}
$$

这里 $\mathcal H_{\{0,-1\}}$ 的微分形式为 $dz/z$ 和 $dz/(1+z)$。上标 $+$ 表示 crossing 偶部分。依据是两个步骤保持线性字母和对数型积分，而不是单凭外部奇点数目猜测函数类型。[D1；S8, §3]

中间的表观奇项必须保持组合相消，不能将它们分别解释为未经声明的主值积分。

### 8.3 这不是最终要采用的任意 Hlog ansatz

式 (8.3) 只提供一个保证足够的上界，没有证明其中每个字都由某个合法 `[1]` 图产生。自然权重上界意味着 $n$ 圈到有限项至多需要权重 $2n$，但不意味着要为这个权重以下的所有 Hlog 分配独立系数。

后续统一母族的目的正是进一步保留拓扑和插入次数限制，减少相对于整个 Hlog 空间的冗余。

---

<a id="sec-9"></a>
## 9. 从四个端点区间压缩到两个字

### 9.1 来自四边形拓扑的字长限制

将四条边分为

$$
\{2,4\}\quad\text{和}\quad\{1,3\}.
$$

固定对边为 $24$ 后，每一条邻边传播子恰好在每一组中各有一个端点。因此留在两个对边端点所在边上、能形成带变量权重的邻边端点只有

$$
\boxed{n-1\text{ 个}.}
\tag{9.1}
$$

其他 $n-1$ 个端点在边 $1,3$ 上，积分后贡献常数。这比仅知道邻边端点总数 $2n-2$ 更强。[D2]

### 9.2 消去右侧区间

对于含对边端点的边字 $pOq$，权重为

$$
\widehat F_p(0,t)\widehat F_q(t,1).
$$

由路径合成与逆路径，有

$$
\boxed{
\widehat F_q(t,1)
=\sum_{q=ab}(-1)^{|a|}
\widehat F_{\operatorname{rev}(a)}(0,t)
\widehat F_b(0,1).
}
\tag{9.2}
$$

这里遍历字 $q$ 的所有切割位置，保留空前缀和空后缀。

再应用 shuffle：

$$
\boxed{
\widehat F_p(0,t)\widehat F_q(t,1)
=\sum_{q=ab}(-1)^{|a|}\widehat F_b(0,1)
\sum_{w\in p\shuffle\operatorname{rev}(a)}\widehat F_w(0,t).
}
\tag{9.3}
$$

所有依赖变量的积分都变成 $[0,t]$，而 $\widehat F_b(0,1)$ 进入常数系数。新字长度不会增加。

两端一起处理，得到

$$
\boxed{
\widehat{\mathcal B}_G\widehat P_G(t)\widehat Q_G(s)
=\sum_{|u|+|v|\le n-1}
c^G_{u,v}(\varepsilon)\,
\widehat F_u(0,t)\widehat F_v(0,s),
}
\tag{9.4}
$$

其中 $c^G_{u,v}$ 不依赖 $z$，由完整区间有序积分和精确有理系数组成。

---

<a id="sec-10"></a>
## 10. 只定义一次的统一母族与有限圈阶覆盖

### 10.1 统一定义

固定

$$
\boxed{
\mathcal T_{u;v}(\varepsilon,z)
=-(1+z)z^{-\varepsilon/2}
\int_0^1dt\,ds\,
D_z(t,s)^{-1+\varepsilon}
\widehat F_u(0,t;\varepsilon)
\widehat F_v(0,s;\varepsilon).
}
\tag{10.1}
$$

这是同一个积分族；核、端点递归、正规化均不随图号改变。所有纯交换 `[1]` 图满足

$$
\boxed{
J_G(z,\varepsilon)
=\frac{z^{(r_G+1/2)\varepsilon}}{\varepsilon^{2n-2}}
\sum_{|u|+|v|\le n-1}
c^G_{u,v}(\varepsilon)\mathcal T_{u;v}(\varepsilon,z).
}
\tag{10.2}
$$

这正是当前统一归约数据文件采用的定义和归一化。[S6]

### 10.2 固定的非交换生成系统

引入不交换符号 $X_A,X_B$，记

$$
\mathbf U_\varepsilon(t)
=\sum_w\widehat F_w(0,t;\varepsilon)X_w.
$$

递归等价于

$$
\boxed{
\frac{d\mathbf U_\varepsilon}{dt}
=\mathbf U_\varepsilon(t)\,
\varepsilon\left[t^{\varepsilon-1}X_A+(1-t)^{\varepsilon-1}X_B\right],
\qquad \mathbf U_\varepsilon(0)=1.
}
\tag{10.3}
$$

不交换性保存 $AB$ 与 $BA$ 的差别。高圈增加的是允许的字长，不是每次新增一套独立积分定义。

### 10.3 候选母核的数量上界

核对两端交换对称：

$$
\mathcal T_{u;v}=\mathcal T_{v;u}.
$$

总字长为 $m$ 的有序字对数是 $(m+1)2^m$。交换的不动点只在偶数 $m$ 出现，数目为 $2^{m/2}$。所以去掉交换重复后，字长不超过 $n-1$ 的候选核数为

$$
\boxed{
B_n=\sum_{m=0}^{n-1}
\frac{(m+1)2^m+\mathbf1_{m\,\mathrm{even}}2^{m/2}}2
=(n-1)2^{n-1}+2^{\lfloor(n-1)/2\rfloor}.
}
\tag{10.4}
$$

| 圈阶 $n$ | 总字长上限 | 候选母核覆盖数 $B_n$ |
|---:|---:|---:|
| 1 | 0 | 1 |
| 2 | 1 | 3 |
| 3 | 2 | 10 |
| 4 | 3 | 26 |
| 5 | 4 | 68 |
| 6 | 5 | 164 |

**这些数字不是最小主积分数，也不是固定 Laurent 深度下的独立系数函数数目。** 它们尚未除去更进一步的积分恒等式；也不声称每个候选核都必须由某个独立物理图单独出现。[D2；S6, `counts` 与 `rank_warning`]

三圈的十个覆盖核是

$$
\begin{gathered}
\mathcal T_{\varnothing;\varnothing},\quad
\mathcal T_{\varnothing;A},\quad
\mathcal T_{\varnothing;B},\\
\mathcal T_{\varnothing;AA},\quad
\mathcal T_{\varnothing;AB},\quad
\mathcal T_{\varnothing;BA},\quad
\mathcal T_{\varnothing;BB},\\
\mathcal T_{A;A},\quad\mathcal T_{A;B},\quad\mathcal T_{B;B}.
\end{gathered}
$$

### 10.4 Crossing 不增加字长

由反射恒等式 (2.6) 和右区间消去式 (9.2)，有

$$
\boxed{
\mathcal T_{u;v}(\varepsilon,1/z)
=\sum_{|a|+|b|\le|u|+|v|}
R_{u,v}^{a,b}(\varepsilon)\mathcal T_{a;b}(\varepsilon,z),
}
\tag{10.5}
$$

系数不依赖 $z$。例如

$$
\mathcal T_{\varnothing;A}(1/z)
=\mathcal T_{\varnothing;\varnothing}(z)-\mathcal T_{\varnothing;B}(z),
$$

$$
\mathcal T_{\varnothing;B}(1/z)
=\mathcal T_{\varnothing;\varnothing}(z)-\mathcal T_{\varnothing;A}(z).
$$

因此固定圈阶的母核列表在 crossing 下封闭，不必为镜像图定义新的积分种类。

### 10.5 最终 I 的统一母核式

代入物理归一化，令 $\delta_G=r_G-(n-1)/2$，得到

$$
\boxed{
\begin{aligned}
I_G={}&\frac{\kappa_GZ(\varepsilon)^n}{\varepsilon^{2n-2}}
\sum_{u,v}c^G_{u,v}(\varepsilon)\\
&\times
\frac{
 e^{\delta_G\varepsilon\ell}\mathcal T_{u;v}(\varepsilon,z)
+e^{-\delta_G\varepsilon\ell}\mathcal T_{u;v}(\varepsilon,1/z)
}{2\cosh(n\varepsilon\ell/2)}.
\end{aligned}
}
\tag{10.6}
$$

完整图的系数必须按这一组合展开。图内极点、镜像抵消、常数序列和母核序列共同决定目标 Laurent 系数。

### 10.6 新体系真正统一了什么？

$\mathcal T$ 是生成层，不是对原来所有函数的简单重命名。它同时保留唯一核、固定两种插入、字长预算和全 $\varepsilon$ 生成关系。若最终只写成未展开 $\mathcal T$，只能称完成母核归约；要交付固定阶数的解析结果，还必须定义并计算其系数函数。

---

<a id="sec-11"></a>
## 11. M、N 与统一母积分之间的转换

### 11.1 生成层和系数层不是同一个计数

记

$$
b_z=\frac{z}{1+z},\qquad
W(\beta)=\ln[\beta(1-\beta)],\qquad
R(\beta)=\ln\frac{1-\beta}{\beta},\qquad
L(\beta,z)=\ell+R(\beta).
$$

原来的定义为

$$
\boxed{
M_{k,r}(z)=\int_0^1\frac{W^kL^r}{\beta-b_z}\,d\beta,
\qquad
N_{k,r}(z)=\int_0^1\frac{W^kRL^r}{\beta-b_z}\,d\beta,
\quad k\ge0,\ r\ge1.
}
\tag{11.1}
$$

由 $L-R=\ell$，精确有

$$
\boxed{\ell M_{k,r}=M_{k,r+1}-N_{k,r}.}
\tag{11.2}
$$

其宇称为

$$
M_{k,r}(1/z)=(-1)^{r+1}M_{k,r}(z),\qquad
N_{k,r}(1/z)=(-1)^rN_{k,r}(z).
$$

$\mathcal T_{u;v}$ 的下标保存邻边插入顺序；$M,N$ 的下标保存展开后的对数幂次。前者仍含 $\varepsilon$，后者不含。因此它们之间不是逐项一一改名。[D3；S7]

### 11.2 单端空字母核的一维桥梁

对 $\mathcal T_{\varnothing;v}$ 完成 $t$ 的线性积分，令 $\beta=1-s$，得到

$$
\boxed{
\mathcal T_{\varnothing;v}(\varepsilon,z)
=\int_0^1\frac{d\beta}{\beta-b_z}
 e^{\varepsilon W/2}
 \frac{\sinh(\varepsilon L/2)}{\varepsilon/2}
 \widehat F_v(0,1-\beta;\varepsilon).
}
\tag{11.3}
$$

其中前两项正是裸 `2e-02` 的生成核，最后一项保存邻边排序权重。

### 11.3 空字母核产生原 M 的特定组合

取 $v=\varnothing$：

$$
\boxed{
\mathcal T_{\varnothing;\varnothing}
=\sum_{k,j\ge0}
\frac{\varepsilon^{k+2j}}{2^{k+2j}k!(2j+1)!}
M_{k,2j+1}.
}
\tag{11.4}
$$

所以

$$
\begin{aligned}
\mathcal T_{\varnothing;\varnothing}
={}&M_{0,1}+\frac\varepsilon2M_{1,1}
+\varepsilon^2\left(\frac18M_{2,1}+\frac1{24}M_{0,3}\right)\\
&+\varepsilon^3\left(\frac1{48}M_{3,1}+\frac1{48}M_{1,3}\right)
+O(\varepsilon^4).
\end{aligned}
$$

若 $t_j=[\varepsilon^j]\mathcal T_{\varnothing;\varnothing}$，利用 $M_{0,3}=-M_{0,1}^2$，有

$$
M_{0,1}=t_0,\qquad M_{1,1}=2t_1,\qquad
M_{2,1}=8t_2+\frac13t_0^2.
\tag{11.5}
$$

一个全 $\varepsilon$ 母核可以产生多个固定阶函数；某一展开阶也可能只给出几个矩的组合，而不能分别区分每个矩。

原双源生成式

$$
\mathscr V(p,q;z)=\int_0^1\frac{d\beta}{\beta-b_z}
 e^{pW}\frac{\sinh(qL)}q
$$

与新母核的关系是

$$
\boxed{
\mathcal T_{\varnothing;\varnothing}(\varepsilon,z)
=\mathscr V(\varepsilon/2,\varepsilon/2;z).
}
\tag{11.6}
$$

两个源参数分别读取矩，物理母核则把它们限制在同一条 $\varepsilon$ 参数线上。[D3]

### 11.4 N 是带权母核展开后的矩组合

最简单非空字满足 $\widehat F_A(0,s)=s^\varepsilon$。其 crossing 偶部分为

$$
\boxed{
\mathcal T_{\varnothing;A}^+
=\int_0^1\frac{d\beta}{\beta-b_z}
 e^{\varepsilon W}\cosh(\varepsilon R/2)
 \frac{\sinh(\varepsilon L/2)}{\varepsilon/2}.
}
\tag{11.7}
$$

取三阶系数得

$$
[\varepsilon^3]\mathcal T_{\varnothing;A}^+
=\frac16M_{3,1}+\frac16M_{1,3}
-\frac14\ell M_{1,2}+\frac18\ell^2M_{1,1}.
$$

用原转换字典消去对数：

$$
\boxed{
\begin{aligned}
[\varepsilon^3]\mathcal T_{\varnothing;A}^+
={}&\frac16M_{3,1}-\frac1{12}M_{1,3}+\frac14N_{1,2}\\
&-\frac14M_{0,1}M_{1,1}-\frac34\zeta_2M_{1,1}.
\end{aligned}
}
\tag{11.8}
$$

因此 $N$ 有统一来源，但这不等于它本来就是原 $M$ 多项式中的冗余项。前一阶段的闭包分析以固定的原 $M$ 定义为前提，区分了有限多项式、允许对数和允许分母的不同空间；新母族没有抹去这些区别。[S8；D3]

### 11.5 一般字可能产生比普通对数幂更丰富的权重

例如

$$
\widehat F_{BA}(0,s;\varepsilon)
=\varepsilon^2\operatorname{Li}_2(s)+O(\varepsilon^3),
$$

故

$$
[\varepsilon^2]\mathcal T_{\varnothing;BA}
=\int_0^1\frac{L(\beta,z)\operatorname{Li}_2(1-\beta)}{\beta-b_z}\,d\beta.
$$

这个式子本身不证明出现了新独立函数；低权重还可能存在归约。它说明字下标保存的是产生权重的积分次序，而不只是普通对数幂次。

因此不能由几个低阶母核成功转成 $M,N$，宣布任意字、任意阶都只需要固定的 $M,N$ 列表。

### 11.6 避免显式 z 的含义

目标不是使函数不再依赖 $z$，而是使结果具有

$$
I_G^{(p)}(z)=\sum_a c_{G,p,a}(\zeta)\,\mathcal F_a(z)
$$

的形式：运动学依赖放入明确的基函数，外部系数不含 $z$ 或 $\ell$。

对于已写成 $M,\ell,\zeta$ 的项，可以用 $N$ 或统一母族的系数表达。特殊商

$$
\ell=\frac{3M_{0,2}}{2M_{0,1}}
$$

提供另一种格式，但它不证明任意 `[1]` 结果都属于原 $M$ 的有理函数域；边界处还需维护分子分母共同趋零时的可去极限。新体系优先保留积分定义与常数系数，而非开放任意 $M$ 分母。[D3；S8]

---

<a id="sec-12"></a>
## 12. 重复 A 子族与统一系数函数 tau

### 12.1 本次三图使用的子族

令

$$
\boxed{T_m(\varepsilon,z)=\mathcal T_{\varnothing;A^m}(\varepsilon,z),\qquad m\ge0.}
\tag{12.1}
$$

由 $\widehat F_{A^m}=s^{m\varepsilon}/m!$：

$$
\boxed{
m!T_m=-(1+z)z^{-\varepsilon/2}
\int_0^1dt\,ds\,D_z(t,s)^{-1+\varepsilon}s^{m\varepsilon}.
}
\tag{12.2}
$$

$m$ 是重复 $A$ 插入字长，不是圈数。只有 $m=0,1,2$ 就足以表示下一节的三张图。[S5, §1]

### 12.2 固定系数规范

定义

$$
\boxed{
\tau_{m,q}(z)
=m!\,[\varepsilon^q]
\frac{T_m(\varepsilon,z)+T_m(\varepsilon,1/z)}2.
}
\tag{12.3}
$$

$q$ 指普通 Taylor 系数，不额外乘 $q!$。$m!$ 正好抵消重复 $A$ 的有序阶乘。

这些 $\tau$ 是统一母族的系数函数，不是按图号重新定义的余项。另一方面，$\tau_{m,q}$ 目前只对应重复 $A$ 子族的偶部分，不能将它擅自当作全部双字母母族系数的完整分类。

### 12.3 明确的一维积分定义

完成一个参数积分后

$$
\boxed{
m!T_m^+
=\int_0^1\frac{d\beta}{\beta-b_z}
 e^{(m+1)\varepsilon W/2}
 \cosh(m\varepsilon R/2)
 \frac{\sinh(\varepsilon L/2)}{\varepsilon/2}.
}
\tag{12.4}
$$

因而

$$
\boxed{
\tau_{m,q}(z)=\int_0^1\frac{P_{m,q}(W,R,L)}{\beta-b_z}\,d\beta,
}
\tag{12.5}
$$

$$
\boxed{
P_{m,q}
=\sum_{\substack{a,b,c\ge0\\a+2b+2c=q}}
\frac{(m+1)^a m^{2b}}
{2^q a!(2b)!(2c+1)!}
W^aR^{2b}L^{2c+1}.
}
\tag{12.6}
$$

在 $m=0,b=0$ 的项中，$m^{2b}$ 按生成式的常数项取 $1$。每个分子都有 $L$ 因子，内部表观奇点可去。这给出了不含未展开 $\varepsilon$ 的系数函数定义。[S5, §4]

所有 $\tau_{m,q}$ 均为 crossing 偶，自然权重为 $q+2$。

三圈需要的多项式为

$$
\boxed{
P_{m,3}
=\frac{(m+1)^3}{48}W^3L
+\frac{(m+1)m^2}{16}WR^2L
+\frac{m+1}{48}WL^3,
}
\tag{12.7}
$$

$$
\boxed{
\begin{aligned}
P_{m,4}={}&\frac{(m+1)^4}{384}W^4L
+\frac{(m+1)^2m^2}{64}W^2R^2L
+\frac{(m+1)^2}{192}W^2L^3\\
&+\frac{m^4}{384}R^4L
+\frac{m^2}{192}R^2L^3
+\frac1{1920}L^5.
\end{aligned}
}
\tag{12.8}
$$

### 12.4 低矩消对数规则

最低系数满足

$$
\boxed{\tau_{0,0}=M_{0,1}=-\frac12\ell^2-3\zeta_2,\qquad
\ell^2=-2\tau_{0,0}-6\zeta_2.}
\tag{12.9}
$$

其余 $k=0$ 矩的关系来自同一个生成式

$$
\int_0^1\frac{e^{\lambda L}-1}{\beta-b_z}\,d\beta
=\ell+\pi\cot(\pi\lambda)
-\frac{\pi e^{\lambda\ell}}{\sin(\pi\lambda)}.
\tag{12.10}
$$

这里使用 $\lambda$ 作辅助源参数，避免与 $\tau_{m,q}$ 的整数下标 $q$ 混淆。式 (12.10) 在源参数零点的可去展开按既定分支理解。[S5, §4]

---

<a id="sec-13"></a>
## 13. 三张图从原始排序到精确 I

本节与下一节整理最近一次统一符号重算，解析系数保持 [S5] 的规范和数值，不将历史核验记录说成本文件新增的计算。

### 13.1 端点排序与颜色表

| 图 | 圈阶 $n$ | 四条边字 | $\mathrm{sym}_G$ | $a_G$ | $\kappa_G$ |
|---|---:|---|---:|---:|---:|
| `2e-02` | 1 | $(\varnothing,O,\varnothing,O)$ | 2 | $1$ | $-1/2$ |
| `4e-02` | 2 | $(\varnothing,O,B,OA)$ | 4 | $-1/2$ | $-1/4$ |
| `6e-04` | 3 | $(\varnothing,O,BB,OAA)$ | 4 | $1/2$ | $-1/8$ |

端点顺序来自原 JSON 的 `phases` 与 `bounds`。三圈颜色分解为

$$
\frac{(1+C_A^2)C_F}{4C_A^2}
=\frac12C_A^2C_F-\frac32C_AC_F^2+C_F^3.
$$

来源为 [S2, `2e-02`]、[S3, `4e-02`]、[S4, `6e-04`]；不是根据最终系数反推排序。

### 13.2 2e-02

没有邻边积分：

$$
J_{2e-02}=\mathscr K_z[1],
$$

故

$$
\boxed{J_{2e-02}=z^{\varepsilon/2}T_0.}
\tag{13.1}
$$

### 13.3 4e-02

边 3 的 $B$ 端点提供 $1/\varepsilon$；边 4 的 $A$ 端点位于对边端点 $s$ 之后：

$$
\int_s^1u^{\varepsilon-1}\,du=\frac{1-s^\varepsilon}{\varepsilon}.
$$

因此

$$
J_{4e-02}=\frac{z^\varepsilon}{\varepsilon^2}
\mathscr K_z[1-s^\varepsilon],
$$

$$
\boxed{J_{4e-02}=\frac{z^{3\varepsilon/2}}{\varepsilon^2}(T_0-T_1).}
\tag{13.2}
$$

### 13.4 6e-04

边 3 的两个 $B$ 有序端点给出

$$
\int_{0<y_1<y_2<1}
(1-y_1)^{\varepsilon-1}(1-y_2)^{\varepsilon-1}\,dy_1dy_2
=\frac1{2\varepsilon^2}.
$$

边 4 的两个 $A$ 端点都在对边端点 $s$ 之后，给出

$$
\int_{s<u_1<u_2<1}u_1^{\varepsilon-1}u_2^{\varepsilon-1}\,du_1du_2
=\frac{(1-s^\varepsilon)^2}{2\varepsilon^2}.
$$

所以

$$
J_{6e-04}=\frac{z^{2\varepsilon}}{4\varepsilon^4}
\mathscr K_z[(1-s^\varepsilon)^2],
$$

$$
\boxed{
J_{6e-04}=\frac{z^{5\varepsilon/2}}{4\varepsilon^4}
(T_0-2T_1+2T_2).
}
\tag{13.3}
$$

这里外部 $1/4$ 来自两组有序积分各自的 $1/2!$；括号内 $2T_2$ 来自 $\widehat F_{AA}=s^{2\varepsilon}/2$。这两个阶乘来源不同。

### 13.5 三图完整对称化

令

$$
x=\varepsilon\ell,\qquad
\overline T_m(\varepsilon,z)=T_m(\varepsilon,1/z).
$$

空字核满足 $\overline T_0=T_0$。代入式 (1.1)：

$$
\boxed{I_{2e-02}=-\frac{Z}{2\cosh(x/2)}T_0,}
\tag{13.4}
$$

$$
\boxed{
I_{4e-02}=-\frac{Z^2}{8\varepsilon^2\cosh x}
\left[e^{x/2}(T_0-T_1)+e^{-x/2}(T_0-\overline T_1)\right],
}
\tag{13.5}
$$

$$
\boxed{
\begin{aligned}
I_{6e-04}=-\frac{Z^3}{64\varepsilon^4\cosh(3x/2)}
\bigl[&e^x(T_0-2T_1+2T_2)\\
&+e^{-x}(T_0-2\overline T_1+2\overline T_2)\bigr].
\end{aligned}
}
\tag{13.6}
$$

不能把 `4e-02` 或 `6e-04` 的未对称化 $J$ 直接假设为 self-mirror。

二圈还可精确约为

$$
\boxed{
I_{4e-02}=\frac{Z(\varepsilon)^2}{4\varepsilon^2\cosh x}
\left[T_0(2\varepsilon,z)-\cosh(x/2)T_0(\varepsilon,z)\right].
}
\tag{13.7}
$$

$2\varepsilon$ 只作用于母核的参数；不同时替换 $Z$、颜色或外部归一化。[S5, §3]

---

<a id="sec-14"></a>
## 14. 三张图直到有限项的统一解析结果

以下全部系数来自最新统一重算文件 [S5, §5]。所有 $\tau$ 的自变量均为同一个 $z$，省略不写。最终公式不含显式 $z,\ell$ 系数，也不使用 $M$ 分母；并不表示基函数自身不再依赖运动学。

### 14.1 2e-02

$$
\boxed{I_{2e-02}=-\frac12\tau_{0,0}+O(\varepsilon).}
\tag{14.1}
$$

全部负幂系数为零。

### 14.2 4e-02

由式 (13.7)

$$
T_0(2\varepsilon)-\cosh(x/2)T_0(\varepsilon)
=\varepsilon\tau_{0,1}
+\varepsilon^2\left(3\tau_{0,2}-\frac{\ell^2}{8}\tau_{0,0}\right)
+O(\varepsilon^3).
$$

消去偶次对数后

$$
\boxed{
I_{4e-02}
=\frac{\tau_{0,1}}{4\varepsilon}
+\frac34\tau_{0,2}
+\frac1{16}\tau_{0,0}^2
+\frac3{16}\zeta_2\tau_{0,0}
+O(\varepsilon).
}
\tag{14.2}
$$

$\varepsilon^{-4},\varepsilon^{-3},\varepsilon^{-2}$ 的系数全为零。到有限项只使用空字核的前三个系数 $\tau_{0,0},\tau_{0,1},\tau_{0,2}$。

### 14.3 6e-04

权重 $(1-s^\varepsilon)^2=O(\varepsilon^2)$，因此原先的显式 $\varepsilon^{-4}$ 只留下二阶极点。

写成

$$
\boxed{
I_{6e-04}=\frac{A_3}{\varepsilon^2}+\frac{B_3}{\varepsilon}+C_3+O(\varepsilon).
}
\tag{14.3}
$$

二阶极点：

$$
\boxed{A_3=-\frac{4\tau_{0,2}+\zeta_2\tau_{0,0}}{64}.}
\tag{14.4}
$$

一阶极点：

$$
\boxed{
B_3=-\frac{17\tau_{0,3}-4\tau_{1,3}+\tau_{2,3}
+(\tau_{0,0}+3\zeta_2)\tau_{0,1}}{32}.
}
\tag{14.5}
$$

有限项：

$$
\boxed{
\begin{aligned}
C_3={}&-\frac{81\tau_{0,4}-7\tau_{1,4}+\tau_{2,4}}{32}\\
&-\frac{(15\tau_{0,0}+51\zeta_2)\tau_{0,2}}{64}\\
&-\frac{
163\tau_{0,0}^3
+2982\zeta_2\tau_{0,0}^2
+8559\zeta_2^2\tau_{0,0}
}{46080}.
\end{aligned}
}
\tag{14.6}
$$

$\varepsilon^{-6},\varepsilon^{-5},\varepsilon^{-4},\varepsilon^{-3}$ 的系数全为零。$A_3,B_3,C_3$ 的自然权重分别是 $4,5,6$。

这些是同一母族系数的有限多项式表示。尚未断言所列 $\tau$ 是最小独立基底，特别是同一阶不同 $m$ 的系数之间仍可能存在关系。

### 14.4 与旧 M、N 表达核对

最低转换为

$$
\boxed{
\tau_{0,0}=M_{0,1},\qquad
\tau_{0,1}=\frac12M_{1,1},\qquad
\tau_{0,2}=\frac18M_{2,1}-\frac1{24}M_{0,1}^2.
}
\tag{14.7}
$$

二圈式 (14.2) 因而恢复为

$$
I_{4e-02}
=\frac{M_{1,1}}{8\varepsilon}
+\frac{3M_{2,1}}{32}
+\frac{M_{0,1}^2}{32}
+\frac{3\zeta_2M_{0,1}}{16}
+O(\varepsilon).
$$

对三圈简单极点，需要

$$
\tau_{0,3}=\frac{M_{3,1}+M_{1,3}}{48},
$$

$$
\tau_{1,3}=\frac{M_{3,1}+M_{1,3}}6
-\frac{\ell M_{1,2}}4+\frac{\ell^2M_{1,1}}8,
$$

$$
\tau_{2,3}=\frac9{16}M_{3,1}+\frac{13}{16}M_{1,3}
-\frac32\ell M_{1,2}+\frac34\ell^2M_{1,1}.
$$

代入式 (14.5) 得

$$
\boxed{
B_3=-\frac{M_{3,1}}{128}-\frac{M_{1,3}}{64}
+\frac{\ell M_{1,2}}{64}
=-\frac{M_{3,1}}{128}-\frac{N_{1,2}}{64}.
}
\tag{14.8}
$$

因此新记号没有改变旧函数，也没有证明原 $M$ 多项式闭包。旧 $N$ 所记录的组合已经被同一个母核系数系统覆盖，不再需要逐图决定新名字。[S5, §6]

### 14.5 原始 Gamma 方案的恢复

对任意圈阶

$$
I_{G,\mathrm{raw}}=e^{n\gamma_E\varepsilon}I_G,
$$

故系数满足

$$
I_{G,\mathrm{raw}}^{(p)}
=\sum_{j=0}^{p+P_G}\frac{(n\gamma_E)^j}{j!}I_G^{(p-j)}.
$$

例如三圈最高为二阶极点时

$$
I_{\mathrm{raw}}
=\frac{A_3}{\varepsilon^2}
+\frac{B_3+3\gamma_EA_3}{\varepsilon}
+C_3+3\gamma_EB_3+\frac92\gamma_E^2A_3
+O(\varepsilon).
$$

这说明有限项会随方案而改变；不能只在一份结果中写含 $\gamma_E$ 的系数，另一份结果中不声明地删去它。

---

<a id="sec-15"></a>
## 15. 覆盖空间、关系空间与最小基底

### 15.1 必须区分四个计数对象

本轮讨论中，“积分有几个”可能指不同对象：

| 对象 | 记录的内容 | 本轮已有结论 |
|---|---|---|
| 原始 orbit | 端点配对、路径排序、颜色和旋转重数 | 当前输入有 18 个 `[1]` orbit，其中 15 个为三圈 |
| 全 $\varepsilon$ 母核 | 固定的 $\mathcal T_{u;v}(\varepsilon,z)$ | 在第 $n$ 圈，取 $\lvert u\rvert+\lvert v\rvert\le n-1$ 即足够覆盖 |
| 固定阶数的系数函数 | $[\varepsilon^q]\mathcal T_{u;v}$ 或其对称化组合 | 可以用统一递归和积分求值产生；本次三图已有明确表达 |
| 最小独立生成元 | 扣除恒等式、低权重乘积等后真正新增的对象 | 未得到全圈维数公式，不能从母核数量直接读取 |

因此，$1,3,10,26,\ldots$ 的母核覆盖数，不等于原始图数，也不等于最后的 $M$ 型基本生成元数。[D2；S6]

同样，$T_0$ 是一个全 $\varepsilon$ 对象，但它的不同展开系数可以分别提供 $M_{0,1}$、$M_{1,1}$、$M_{2,1}$ 等不同权重的函数。说“二圈只需三个基本矩”，与“二圈经过 crossing 可写成同一个母核的不同频率”，并不矛盾。

### 15.2 应当研究实际图的像空间，而不是完整 Hlog 的任意 ansatz

固定圈阶 $n$ 和展开上限 $p_{\max}$，定义实际 `[1]` 图产生的系数空间

$$
\boxed{
\mathcal V^{[1]}_{n,\le p_{\max}}
=
\operatorname{span}_{\text{不依赖 }z\text{ 的常数}}
\left\{
I_G^{(p)}(z):
G\text{ 是合法的 }[1]\text{ 图},\
-(2n-2)\le p\le p_{\max}
\right\}.
}
\tag{15.1}
$$

这里允许的常数范围须随任务固定，例如取多重 $\zeta$ 常数及其乘积。式 (15.1) 强调：只有真实端点插入规则产生的组合才进入空间，不把字母表 $\{0,-1\}$ 的所有 Hlog 系数都设成独立未知数。

对应的计算链是

```text
合法邻边配对与端点排序
        ↓
受 |u|+|v| ≤ n−1 限制的统一母核组合
        ↓
crossing、颜色投影和目标 Laurent 截断
        ↓
实际系数函数的有限列表
        ↓
完整恒等式归约，求其关系空间
```

Hlog/GPL 表达可以用于执行和核验这个流程；它是运算表示，不必直接充当最终的最小基底。

### 15.3 最小线性基底与最小代数生成元仍需分开

若目标是最小线性基底，可以把实际系数函数转换到同一个完整的函数表示中，求它们的精确线性关系。必须保留端点常数和低长度项，不能只比较最高 symbol。

若目标是用户更关心的“新增几个基本函数”，就应扣除已有低权重函数的乘积。设 $\mathcal V_w^{[1]}$ 为当前任务在权重 $w$ 产生的部分，$\mathcal D_w$ 为已知低权重对象、常数乘法和乘积生成的权重 $w$ 部分，则新增生成元对应

$$
\boxed{
\mathcal P_w^{[1]}
=
\mathcal V_w^{[1]}
\big/
\left(\mathcal V_w^{[1]}\cap\mathcal D_w\right).
}
\tag{15.2}
$$

式 (15.2) 是本轮提出的进一步归约目标，而不是已经求出的维数结果。[D2]

还有一个系数域区别：通常 IBP 归约可以允许 $z,\varepsilon$ 的有理函数作系数；用户希望的最终表达要求基函数外没有显式运动学系数。这两种要求下的“主积分数”不应直接混用。

### 15.4 二圈经验应怎样用于高圈？

用户给出的二圈经验是：直到有限项，各图可以用少数 $M$ 生成元及 $\zeta$ 常数表示。这是有价值的低权重约化信息，但不能仅由此推出任意三圈 `[1]` 图在原 $M$ 多项式空间中闭合。

在本轮的自然权重计数中，二圈有限项到权重四，三圈有限项到权重六。新母族保留了更多插入字的系数，所以能够统一容纳提高权重后出现的组合；它没有据此宣称那些组合都是原 $M$ 的多项式。

新框架与旧分析的兼容关系是：

$$
\boxed{
\text{固定的原 }M\text{ 空间可能不闭合}
\quad\text{与}\quad
\text{整个 }[1]\text{ 族可由一个统一母族覆盖}
\quad\text{同时成立。}
}
$$

---

<a id="sec-16"></a>
## 16. 确定性计算流程与核验边界

### 16.1 计算程序的公共前端

对空集和 `[1]`，前端都应从原始 JSON 的结构字段读取信息，而不是根据图编号猜测。

```text
读取原始 JSON
    ↓
核对纯交换适用范围：无内部顶点、无自能插入
    ↓
从 phases/propagators 确定邻边或对边
    ↓
从 bounds 恢复每条边的原始端点顺序
    ↓
生成 A/B 字；唯一对边端点标记为 O
    ↓
核对有序域到 hypercube 的 Jacobian 与 expanded 被积函数
```

`channel_partition=[]` 或 `[1]` 本身不是跳过其他结构检查的理由。该字段只统计直接对边连接；对未来包含其他内部结构的数据，必须重新检查算法的适用条件。

### 16.2 空集后端

空集后端执行

$$
(w_1,w_2,w_3,w_4;r_G)
\longrightarrow
\mathcal B_G(\varepsilon)
\longrightarrow
I_G(z,\varepsilon).
$$

具体操作是端点发散段剥离、常数有序积分展开、Gamma 级数卷积、双曲函数比值展开。最终可以选择 $E_m$ 或 $M_{0,1}^m$ 表示运动学。

普通 $\zeta$ 与多重 $\zeta$ 的进一步约化是常数后端的问题；不影响空集的 $z$ 依赖已经完全确定这一结论。

### 16.3 单对边后端

单对边后端执行

$$
\begin{aligned}
&\text{带 }O\text{ 的有序字}\\
&\quad\longrightarrow
\text{前后区间分解和路径合成}\\
&\quad\longrightarrow
\sum c^G_{u,v}(\varepsilon)\mathcal T_{u;v}(\varepsilon,z)\\
&\quad\longrightarrow
\text{crossing 与物理归一化}\\
&\quad\longrightarrow
\text{目标阶数的母核系数}\\
&\quad\longrightarrow
\text{完整函数求值与经验证的基底转换}.
\end{aligned}
$$

程序可以缓存同一个字积分和母核系数，避免不同 orbit 重复计算。没有覆盖的归约关系应当保留为明确的未归约项，而不是由 agent 临时猜一个系数补进输出。

### 16.4 展开阶数的管理

若某个精确表示写成

$$
I(\varepsilon)=\varepsilon^{-P}A(\varepsilon),
\qquad A\text{ 在零点正则},
$$

要得到 $I^{(p)}$，须保留 $A$ 到 $\varepsilon^{P+p}$。对于通用实现，可以保守使用空集的 $P=2n$、单对边的 $P=2n-2$；在证明额外抵消之后再降低所需阶数。

同一积分的不同表示可能具有不同的显式极点分母。例如三圈 `6e-04` 的母核组合式含 $\varepsilon^{-4}$，完成一个参数积分后的一维差分式可以出现 $\varepsilon^{-5}$。后者的指数差分同时带来一个零点，不能只复制分母阶数而忽略分子的起始阶。

因此“求到 $\varepsilon^0$”不意味着各个正则因子只展开到零阶；各因子必须按完整乘积所需的阶数截断。

### 16.5 纯代码、agent 与 Lean 的分工

本轮采用的分工是：

| 环节 | 合适的处理方式 |
|---|---|
| 已确定的有序字规则、颜色代数、级数卷积 | 固定代码执行 |
| 新积分权重的研究、新归约恒等式的发现 | 可以由人或 agent 参与，但结果须单独核验 |
| 符号转录、原始被积函数重建、反向代入 | 自动检查 |
| 多点高精度积分 | 数值一致性检查，不作为恒等式证明的替代 |
| 对通用定理及证书解释器的形式化证明 | 可选的 Lean 工作，不是当前固定求值流程的前置条件 |

纯代码不等于自动正确；它只是避免每次运行重新选择或编造规则。可靠性来自适用条件、可审计的归约、独立检查和对未支持情形的明确拒绝。

### 16.6 本文引用的既有检查记录

最新三图重算文件 [S5, §7] 记录了：原始 JSON 与降维积分在三个整数 $\varepsilon$ 参数点的 9 次精确多项式比较；6 个非零 Laurent 系数的完整符号反向代入；颜色及归一化核对；在 $z=0.4,1,2$ 上进行的 18 次、40 位工作精度数值比较，最大差约为 $9.69\times10^{-37}$；以及有限 $\varepsilon$ 余项的阶数检查。

统一母族数据 [S6] 记录了当前 18 个 `[1]` 图的有序字归约。其积分前系数矩阵的秩不能当作积分后的主函数数目。

**以上是所引用源文件已经记录的检查，不是本次 Markdown 整理重新执行的一套数值计算。** 本文没有新增“所有圈阶逐图验证完成”的声明，没有运行未提供源码的 `finalize.wl`，也没有进行 Lean 内核检查。

---

<a id="sec-17"></a>
## 17. 结论对照与尚未完成的任务

### 17.1 空集与 [1] 的并列结论

| 层面 | 空集族 | 单对边 `[1]` 族 |
|---|---|---|
| 初始结构 | $n$ 条邻边 | 一条对边加 $n-1$ 条邻边 |
| 尺度提出后的余项 | 四条边的常数有序积分 | 常数有序积分乘双端点带权对边核 |
| 通用 $J$ | $z^{r\varepsilon}\prod_iF_{w_i}$ | $z^{r\varepsilon}\mathcal B_G\mathscr K_z[P_GQ_G]$ |
| 最高可能极点 | $\varepsilon^{-2n}$ | $\varepsilon^{-(2n-2)}$ |
| 对称化的运动学结构 | 偶次对数多项式 | 偶通道与“奇对数乘奇函数”通道 |
| 常数保证范围 | 多重 $\zeta$ 代数 | 多重 $\zeta$ 代数 |
| 自然统一对象 | $F_w$ 与 $E_m$ | $F_w$ 与 $\mathcal T_{u;v}$ |
| 与原 $M$ 的关系 | 全部偶次对数可换成 $M_{0,1}$ 的幂 | 原 $M,N$ 是某些母核展开后的矩表达 |
| 固定圈阶限制 | $2n$ 个端点形成四个常数字 | 非平凡双端点字长 $\lvert u\rvert+\lvert v\rvert\le n-1$ |
| 最小基底状态 | 运动学只有一个偶对数生成元；常数约化另计 | 尚未给出全圈最小维数公式 |

### 17.2 已经保留的结构结论

空集有

$$
\boxed{
I_G^{(p)}\in
\mathcal Z_{\mathrm{MZV}}[\ln^2z]
=
\mathcal Z_{\mathrm{MZV}}[M_{0,1}].
}
$$

单对边图有统一母核覆盖

$$
\boxed{
J_G=
\frac{z^{(r_G+1/2)\varepsilon}}{\varepsilon^{2n-2}}
\sum_{|u|+|v|\le n-1}
 c^G_{u,v}(\varepsilon)\mathcal T_{u;v}(\varepsilon,z).
}
$$

Hlog 字母表 $\{0,-1\}$ 给出单对边的可靠求值上界；实际目标空间应取受合法插入规则限制的母族像空间，而不是无约束的全部 Hlog。

三张图 `2e-02`、`4e-02`、`6e-04` 的指定色通道结果到有限项，已经在同一套 $\tau_{m,q}$ 系数记号下表达。它们不含基函数外的显式 $z,\ln z$，但每个 $\tau$ 自身仍然是 $z$ 的函数。

### 17.3 不应从现有结果推出的结论

当前讨论和源文件没有建立以下更强结论：

1. 任意圈阶的所有多重 $\zeta$ 常数都能化成普通 $\zeta_k$ 的乘积。
2. 任意 `[1]` 图的全部系数都属于原 $M_{k,r}$ 的有限多项式空间。
3. 所有合法双字母核都能由重复 $A$ 子族的 $\tau_{m,q}$ 覆盖。
4. 母核候选数量 $B_n$ 等于最小主积分数量或最小代数生成元数量。
5. 两个三圈例子的 $M$–$N$ 表达足以证明全部三圈乃至全圈的最小闭包。
6. 单图的函数空间与全部图求和、取对数或重整化后的空间必然相同。

这些不是停止计算的理由，而是后续关系归约要精确定义和验证的问题。新母族使它们变成固定规则内的像空间和关系空间研究，不再需要逐图临时增加函数名称。

### 17.4 最终组织原则

$$
\boxed{
\text{拓扑与端点排序}
\longrightarrow
\text{有序字积分}
\longrightarrow
\text{统一母核}
\longrightarrow
\text{对称化后的 Laurent 系数}
\longrightarrow
\text{经验证的简约基底}.
}
$$

空集对应这一体系的常数与单字母运动学部分；`[1]` 对应一个固定双线性对边核上的有限次有序插入。原 $M,N$ 与新的 $\mathcal T,\tau$ 不是互相排斥的理论，而是生成层和系数表达层之间的不同坐标。


---

<a id="appendix-a"></a>
## 附录 A. 当前 18 个 [1] 图的母核归约

本表转录自最新统一母族数据 `single_opposite_universal_core/reductions.json`。所有行采用同一归一化：

$$
\boxed{
J_G(z,\varepsilon)=
\frac{z^{(r_G+1/2)\varepsilon}}{\varepsilon^{2n-2}}
\sum_{u,v}c^G_{u,v}(\varepsilon)\mathcal T_{u;v}(\varepsilon,z).
}
\tag{A.1}
$$

最后一列只列出式 (A.1) 的母核线性组合，不再重复公共尺度与极点前因子。表中

$$
Q=Q(\varepsilon)=\frac{\Gamma(1+\varepsilon)^2}{\Gamma(1+2\varepsilon)}
$$

是常数 Gamma 比值，不是第 6 节的端点权重函数 $Q_G(s;\varepsilon)$。所有核的参数都是 $(\varepsilon,z)$，$\mathcal T_{u;v}=\mathcal T_{v;u}$ 的交换已用于统一排序。[S6]

| 图 | 圈阶 $n$ | $r_G$ | 四条边的有序字 | 母核线性组合 |
|---|---:|---:|---|---|
| `2e-02` | 1 | 0 | $\bigl(\varnothing,O,\varnothing,O\bigr)$ | $\mathcal T_{\varnothing;\varnothing}$ |
| `4e-02` | 2 | 1 | $\bigl(\varnothing,O,B,OA\bigr)$ | $\mathcal T_{\varnothing;\varnothing} - \mathcal T_{\varnothing;A}$ |
| `4e-05` | 2 | 0 | $\bigl(\varnothing,BO,A,O\bigr)$ | $\mathcal T_{\varnothing;B}$ |
| `6e-04` | 3 | 2 | $\bigl(\varnothing,O,BB,OAA\bigr)$ | $\frac{1}{4}\mathcal T_{\varnothing;\varnothing} - \frac{1}{2}\mathcal T_{\varnothing;A} + \frac{1}{2}\mathcal T_{\varnothing;AA}$ |
| `6e-05` | 3 | 2 | $\bigl(\varnothing,O,BB,OAA\bigr)$ | $\frac{1}{4}\mathcal T_{\varnothing;\varnothing} - \frac{1}{2}\mathcal T_{\varnothing;A} + \frac{1}{2}\mathcal T_{\varnothing;AA}$ |
| `6e-06` | 3 | 2 | $\bigl(\varnothing,O,BB,AOA\bigr)$ | $\frac{1}{2}\mathcal T_{\varnothing;A} - \mathcal T_{\varnothing;AA}$ |
| `6e-13` | 3 | 1 | $\bigl(\varnothing,BO,AB,OA\bigr)$ | $Q\mathcal T_{\varnothing;B} - Q\mathcal T_{A;B}$ |
| `6e-14` | 3 | 1 | $\bigl(\varnothing,BO,BA,OA\bigr)$ | $-\left(Q - 1\right)\mathcal T_{\varnothing;B} + \left(Q - 1\right)\mathcal T_{A;B}$ |
| `6e-15` | 3 | 1 | $\bigl(\varnothing,BO,BA,AO\bigr)$ | $-\left(Q - 1\right)\mathcal T_{A;B}$ |
| `6e-16` | 3 | 1 | $\bigl(\varnothing,OB,BA,OA\bigr)$ | $-\left(Q - 1\right)\mathcal T_{\varnothing;\varnothing} + \left(Q - 1\right)\mathcal T_{\varnothing;A} + \left(Q - 1\right)\mathcal T_{\varnothing;B} - \left(Q - 1\right)\mathcal T_{A;B}$ |
| `6e-26` | 3 | 0 | $\bigl(\varnothing,BBO,AA,O\bigr)$ | $\frac{1}{2}\mathcal T_{\varnothing;BB}$ |
| `6e-27` | 3 | 0 | $\bigl(\varnothing,BOB,AA,O\bigr)$ | $\frac{1}{2}\mathcal T_{\varnothing;B} - \mathcal T_{\varnothing;BB}$ |
| `6e-28` | 3 | 0 | $\bigl(\varnothing,BBO,AA,O\bigr)$ | $\frac{1}{2}\mathcal T_{\varnothing;BB}$ |
| `6e-29` | 3 | 1 | $\bigl(A,O,B,BOA\bigr)$ | $-\mathcal T_{\varnothing;AB} + \mathcal T_{\varnothing;B} - \mathcal T_{\varnothing;BA}$ |
| `6e-30` | 3 | 1 | $\bigl(A,O,B,BAO\bigr)$ | $\mathcal T_{\varnothing;BA}$ |
| `6e-31` | 3 | 1 | $\bigl(A,O,B,OBA\bigr)$ | $-\left(Q - 1\right)\mathcal T_{\varnothing;\varnothing} + \mathcal T_{\varnothing;AB} - \mathcal T_{\varnothing;B}$ |
| `6e-36` | 3 | 2 | $\bigl(B,OA,B,OA\bigr)$ | $\mathcal T_{\varnothing;\varnothing} - 2\mathcal T_{\varnothing;A} + \mathcal T_{A;A}$ |
| `6e-39` | 3 | 0 | $\bigl(A,BO,A,BO\bigr)$ | $\mathcal T_{B;B}$ |

该表保留了原始 orbit；相同的 $J$ 结构不能在颜色或重数未核对前擅自合并。它也不是“18 个最终 $I$ 已经全部解析求出”的结果表：本轮最新统一系数重算详细完成的是正文中的三张图。

数据文件还记录了三圈积分前权重系数矩阵的秩为 10，并明确附有 `rank_warning`。这一数据只能描述积分前的线性组合，不能直接作为积分后主函数数目或 Laurent 系数空间维数。[S6]

---

<a id="appendix-b"></a>
## 附录 B. 系数及符号转换速查

### B.1 各个下标的意义

| 记号 | 下标含义 | 是否含未展开的 $\varepsilon$ |
|---|---|---|
| $F_w(a,b;\varepsilon)$ | $w$ 为端点排序字 | 是 |
| $E_m(z)$ | $2m$ 为偶次对数次数 | 否 |
| $M_{k,r}(z)$ | $W^kL^r$ 的对数幂次 | 否 |
| $N_{k,r}(z)$ | $W^kRL^r$ 的对数幂次 | 否 |
| $\mathcal T_{u;v}(\varepsilon,z)$ | 两个对边端点的邻边插入字 | 是 |
| $T_m(\varepsilon,z)$ | $\mathcal T_{\varnothing;A^m}$ 的简记 | 是 |
| $\tau_{m,q}(z)$ | 重复 $A$ 长度 $m$、Taylor 系数阶数 $q$ | 否 |
| $I_G^{(p)}(z)$ | 最终物理归一化的 Laurent 次数 $p$ | 否 |

$m$、$q$、$k$、$r$ 均不能未经说明解释为圈数；圈数统一用 $n$。

### B.2 端点与 Gamma 字典

$$
\widehat F_{A^m}(0,s)=\frac{s^{m\varepsilon}}{m!},
\qquad
F_{A^m}(0,1)=F_{B^m}(0,1)=\frac1{m!\varepsilon^m}.
$$

$$
Z=e^{-\gamma_E\varepsilon}\Gamma(1-\varepsilon),
\qquad
Z^n=\exp\left(n\sum_{k\ge2}\frac{\zeta_k}{k}\varepsilon^k\right).
$$

若

$$
Z^n=\sum_{j\ge0}g_{n,j}\varepsilon^j,
$$

则本轮使用的指数级数递推可写为

$$
g_{n,0}=1,
\qquad
g_{n,j}=\frac{n}{j}\sum_{k=2}^{j}\zeta_k g_{n,j-k}
\quad(j\ge1).
$$

和为空时结果为零。这是对已用 Gamma 指数生成式求导并逐阶比较的写法，不是额外的积分归约假设。

### B.3 低矩与空集字典

$$
\begin{aligned}
M_{0,1}&=-\frac12\ell^2-3\zeta_2,\\
\ell^2&=-2M_{0,1}-6\zeta_2,\\
M_{0,2}&=\frac23\ell M_{0,1},\\
M_{0,3}&=-M_{0,1}^2,\\
\ell M_{0,2}&=-\frac43M_{0,1}^2-4\zeta_2M_{0,1},\\
E_m&=\frac{(-2M_{0,1}-6\zeta_2)^m}{(2m)!}.
\end{aligned}
$$

允许分母时可写 $\ell=3M_{0,2}/(2M_{0,1})$，但它不是原多项式空间内的等式表达规则。在相关零点须保持解析延拓的可去极限，不能把它当作引入任意 $M$ 分母的许可。

### B.4 M/N 与母核系数字典

$$
\ell M_{k,r}=M_{k,r+1}-N_{k,r}.
$$

$$
\tau_{0,0}=M_{0,1},
\qquad
\tau_{0,1}=\frac12M_{1,1},
\qquad
\tau_{0,2}=\frac18M_{2,1}-\frac1{24}M_{0,1}^2.
$$

$$
M_{0,1}=\tau_{0,0},
\qquad
M_{1,1}=2\tau_{0,1},
\qquad
M_{2,1}=8\tau_{0,2}+\frac13\tau_{0,0}^2.
$$

$$
\mathcal T_{\varnothing;\varnothing}(\varepsilon,z)
=\mathscr V\left(\frac\varepsilon2,\frac\varepsilon2;z\right).
$$

这最后一条式子把原来的双源矩生成函数限制在物理母核的同一 $\varepsilon$ 参数线上；不能反过来把一个参数线的单个系数当成所有双源矩的独立值。

### B.5 最容易混淆的三个阶乘与前因子

| 因子 | 来源 | 本文处理 |
|---|---|---|
| $1/m!$ | 重复 $A^m$ 的路径有序区域 | 保留在 $T_m$ 中；定义 $\tau_{m,q}$ 时乘 $m!$ |
| $1/q!$ | 指数的 Taylor 展开 | 已含在生成式及 $P_{m,q}$ 的相应系数中；不对 $[\varepsilon^q]$ 再乘一个 $q!$ |
| $(-2)^{-n}$ | 源文件的公共尺度提取 | 只计入一次 $\kappa_G$ |

`6e-04` 的邻边积分还产生两个 $1/2!$ 的乘积 $1/4$。它既不是 orbit 重数，也不是上表中定义 $T_2$ 的阶乘可以任意抵消的对象。

---

<a id="appendix-c"></a>
## 附录 C. 来源与整理说明

### C.1 原始输入

| 标记 | 文件 | 本文使用范围 |
|---|---|---|
| **S1** | `N4-n-gluon-CJ-derivation_v2.md` | §1 一般四边形与端点排序；§5c 邻边/对边核；§6b–6c 公共 $C\times J$、$(-2)^{-n}$ 和 hypercube 结构 |
| **S2** | `montopar-initial.json` | `2e-02` 的 `topology.phases`、`bounds`、`c4_symmetry`、`C` 与 `J` |
| **S3** | `montopar-initial(1).json` | `4e-02` 及本轮二圈 `[1]` 归约所需输入 |
| **S4** | `montopar-initial(2).json` | `6e-04` 原始结构和全部三圈 `[1]` 图的原始数据 |

S1 的原题名含有 “N=4 四边形” 字样；本文按其明确列出的四边形几何与纯传播子交换条件引用，不由题名额外添加内部顶点、物质环或其他未给出的动力学内容。

原始参数化来源定位：S1 第 153–174 行为邻边与对边核，第 245–279 行为公共前因子和 $J$；S2 的 `2e-02` 位于第 127–228 行；S3 的 `4e-02` 位于第 172–317 行；S4 的 `6e-04` 位于第 601–791 行。行号对应本轮读取的文本版本，实际检索宜优先使用 orbit ID 和字段名。

### C.2 本轮已有的推导与程序数据

| 标记 | 文件或数据 | 本文使用范围与状态 |
|---|---|---|
| **S5** | `three_diagrams_T/derivation.md`；`three_diagrams_T/coefficients.wl` | 最新三图统一重算；$T_m$、$\tau_{m,q}$ 定义；精确 $J,I$；全部目标 Laurent 系数和既有核验记录 |
| **S6** | `single_opposite_universal_core/reductions.json`；`single_opposite_universal_core/check_word_reduction.py` | 统一双字母核定义；18 个 `[1]` 图的积分前精确归约；覆盖数及秩的适用范围警告 |
| **S7** | `6e29_MN_solution/derivation.md` | 原 $M,N$ 的固定定义、$\ell M=M-N$ 转换以及系数表达层的具体例证 |
| **S8** | `single_opposite_general_space/derivation.md`；`M_closure_check/theorem.md` | 一般单对边函数空间上界、与原 $M$ 多项式空间的区别及闭包障碍讨论；此处主要作为前述结论的来源记录 |

S5 中的历史数值记录按原记录保留其工作精度与误差解释；本文没有把数值差改称为严格误差界。S6 明确限定为有序字归约，不将其描述为全部 18 个 $I$ 的完整求值。

### C.3 此次新讨论的对应来源

下列内容来自本次对话中的连续推导，尚未单独保存为此前的源文件；本文将它们系统化收入正文。

| 标记 | 对应讨论 | 收入位置 |
|---|---|---|
| **D0** | “任意 $n$ 圈的空集族积分”的通用 $J$、对称化 $I$、MZV 常数结构、端点极点计数和空集自身 $E_m$ 族 | 第 2–5 节 |
| **D1** | “任意 $n$ 圈下的 `[1]` 族”的标记字分解、极点结构、奇偶通道及线性可约空间上界 | 第 6–8 节 |
| **D2** | “对于 `[1]` 通用、不反复增加函数名、避免完整 Hlog 冗余”的统一母族、字长限制、覆盖计数和最小空间问题 | 第 9–10、15 节 |
| **D3** | “原 $M,N$ 和新母积分的关系”的降维桥梁、双源/单源区别、带权核产生 $N$ 的具体展开 | 第 11–12 节 |

第 13–14 节以 S5 的最新结果为准，并与 D3 的系数解释衔接。本文不把早先一版较大的 Hlog 包含空间替换成已经求得的最小空间，也不将后来的统一命名理解为旧 $M$ 多项式闭包得到恢复。

### C.4 整理时保持的边界

本文主要完成章节组织、统一符号、去除重复和建立交叉索引，保留前述推导中的适用范围、有限阶结果和未完成任务。本文不是一次新的外部文献研究，也不是对所有一般性数学定理重新进行形式化验证。

前述推导中使用的外部数学工具包括 Euler Beta–Gamma 恒等式、Gamma Taylor 展开、迭代积分的路径合成与 shuffle、端点正规化和线性可约超对数求值。其参考文献已在 S5、S8 等原推导中记录；本次没有将未重新检索的文献包装成新的核验结果。

本次 Markdown 整理与源文件中既有数值/符号核验应分开理解。最小基底、全圈普通 $\zeta$ 常数约化以及全部合法字的最小矩表示，均保持为正文明确列出的后续问题。

---

**文档主线**：空集以常数有序字积分与单字母运动学 $E_m$ 组织；单对边以受拓扑字长约束的 $\mathcal T_{u;v}$ 组织；$M,N$ 是展开后的矩表达，$\tau$ 是同一母族的系数记号。统一母族解决“从哪里生成”，最小基底归约解决“最终还剩多少独立信息”。
