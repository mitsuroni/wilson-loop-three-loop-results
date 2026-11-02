> **2026-09-19替代说明：本页以下内容为历史审查记录，不能再作为当前实现依据。** 正确推导已统一写入[4e-3g+3g-general-derivation.md](../4e-3g+3g-general-derivation.md)，包含完整D维算子、链式法则、尺度和有序域。当前Python实现及旧C/J尚未修复；历史勾选状态不授予新实现正确性。

---

# 4e-3g+3g：$I(a,b,c,d) = C \times J$ 完整推导

> 独立重推导，基于 [effective-vertices.md](../effective-vertices.md) §3 的 $V_{2\text{TGV}}$ 起点。
> 验证脚本：`verify/v*.m`
> 更新日期：2026-08-27

---

## 推导流程

```
effective-vertices.md §3 (V_{2TGV})
  → §2: 动量空间 I(a,b,c,d)
  → §3: 颜色因子提取
  → §4-§5: Lorentz 缩并 + ℓ→i∂_X 替换 → 20 项算子表
  → §6: Fourier 变换 z_A,z_B 解耦
  → §7: Symanzik 参数化 → B_5, U
  → §8: 二次型链式法则 → Pnum^(1), Pnum^(2)
  → §9: C × J 因子分解
```

---

## §2c: $V_{2\text{TGV}}$ Lorentz 分解

> 验证：`verify/v01-lorentz-decomposition.m`
> 状态：✅ G1 通过，⚠️ G2/G3/G4a-d 与 effective-vertices.md 不一致，以 FeynCalc 直接结果为准

### Step 1：初始表达式

两个三胶子顶点通过内部胶子传播子 $g^{\nu_5\nu_6}$ 耦合。$\delta$ 函数强制 $\ell_5 = -(\ell_1+\ell_2) = \ell_3+\ell_4$，积分退化为代入：

$$
\boxed{
\begin{aligned}
V_{2\text{TGV}}(a,b,c,d;\ell_1,\ell_2,\ell_3,\ell_4)
&= f^{A_a A_b B}\,f^{A_c A_d B}\;
   \frac{-i}{(\ell_1+\ell_2)^2+i0}\;g^{\nu_5\nu_6} \\
&\quad \times v_a^{\mu_1}v_b^{\mu_2}\;
   \Gamma^{\text{3g}}_{\mu_1\mu_2\nu_5}(-\ell_1,-\ell_2,\;\ell_1+\ell_2) \\
&\quad \times v_c^{\mu_3}v_d^{\mu_4}\;
   \Gamma^{\text{3g}}_{\mu_3\mu_4\nu_6}(-\ell_3,-\ell_4,\;\ell_3+\ell_4).
\end{aligned}}
$$

三胶子运动学顶点：

$$
\Gamma^{\text{3g}}_{\mu\nu\rho}(k_1,k_2,k_3) = g_{\mu\nu}(k_1-k_2)_\rho + g_{\nu\rho}(k_2-k_3)_\mu + g_{\rho\mu}(k_3-k_1)_\nu.
$$

### Step 2：缩并 TGV A

缩并 $v_a^{\mu_1} v_b^{\mu_2}$ 与 $\Gamma^{\text{3g}}_{\mu_1\mu_2\nu_5}$：

$$
\boxed{
\begin{aligned}
\text{TGV A}_{\nu_5} &= v_{ab}\,(\ell_2-\ell_1)_{\nu_5}
   + v_{b,\nu_5}\,(-v_a\!\cdot\!\ell_1 - 2v_a\!\cdot\!\ell_2) \\
   &\quad + v_{a,\nu_5}\,(2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2).
\end{aligned}}
$$

| 分量 | 向量方向 | 标量系数 | 来源 |
|:---:|------|------|------|
| $A_1$ | $(\ell_2-\ell_1)_{\nu_5}$ | $v_{ab}$ | $g_{\mu_1\mu_2}$ 项 |
| $A_2$ | $v_{b,\nu_5}$ | $-v_a\!\cdot\!\ell_1 - 2v_a\!\cdot\!\ell_2$ | $g_{\nu_5\mu_1}$ 项 |
| $A_3$ | $v_{a,\nu_5}$ | $2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2$ | $g_{\mu_2\nu_5}$ 项 |

### Step 3：缩并 TGV B

缩并 $v_c^{\mu_3} v_d^{\mu_4}$ 与 $\Gamma^{\text{3g}}_{\mu_3\mu_4\nu_6}$：

$$
\boxed{
\begin{aligned}
\text{TGV B}_{\nu_6} &= v_{cd}\,(\ell_4-\ell_3)_{\nu_6}
   + v_{d,\nu_6}\,(-v_c\!\cdot\!\ell_3 - 2v_c\!\cdot\!\ell_4) \\
   &\quad + v_{c,\nu_6}\,(2v_d\!\cdot\!\ell_3 + v_d\!\cdot\!\ell_4).
\end{aligned}}
$$

| 分量 | 向量方向 | 标量系数 | 来源 |
|:---:|------|------|------|
| $B_1$ | $(\ell_4-\ell_3)_{\nu_6}$ | $v_{cd}$ | $g_{\mu_3\mu_4}$ 项 |
| $B_2$ | $v_{d,\nu_6}$ | $-v_c\!\cdot\!\ell_3 - 2v_c\!\cdot\!\ell_4$ | $g_{\nu_6\mu_3}$ 项 |
| $B_3$ | $v_{c,\nu_6}$ | $2v_d\!\cdot\!\ell_3 + v_d\!\cdot\!\ell_4$ | $g_{\mu_4\nu_6}$ 项 |

### Step 4：$g^{\nu_5\nu_6}$ 缩并 → 9 项乘积

$V_{2\text{TGV}}$ 的 Lorentz 部分 = $\text{TGV A}_{\nu_5} \cdot \text{TGV B}^{\nu_5} = \sum_{i,j=1}^3 A_i \cdot B_j$。

| $A_i$ | $B_j$ | $A_i \cdot B_j$ | $v_{ij}$ 前因子 | 归类 |
|------|------|----------------|:---:|:---:|
| $A_1$ | $B_1$ | $v_{ab}v_{cd}\,(\ell_2-\ell_1)\!\cdot\!(\ell_4-\ell_3)$ | $v_{ab}v_{cd}$ | **G1** |
| $A_1$ | $B_2$ | $v_{ab}\,(\ell_2-\ell_1)\!\cdot\!v_d\;(-v_c\!\cdot\!\ell_3 - 2v_c\!\cdot\!\ell_4)$ | $v_{ab}$ | G2 |
| $A_1$ | $B_3$ | $v_{ab}\,(\ell_2-\ell_1)\!\cdot\!v_c\;(2v_d\!\cdot\!\ell_3 + v_d\!\cdot\!\ell_4)$ | $v_{ab}$ | G2 |
| $A_2$ | $B_1$ | $v_{cd}\,v_b\!\cdot\!(\ell_4-\ell_3)\;(-v_a\!\cdot\!\ell_1 - 2v_a\!\cdot\!\ell_2)$ | $v_{cd}$ | G3 |
| $A_2$ | $B_2$ | $v_{bd}\,(-v_a\!\cdot\!\ell_1 - 2v_a\!\cdot\!\ell_2)(-v_c\!\cdot\!\ell_3 - 2v_c\!\cdot\!\ell_4)$ | $v_{bd}$ | G4d |
| $A_2$ | $B_3$ | $v_{bc}\,(-v_a\!\cdot\!\ell_1 - 2v_a\!\cdot\!\ell_2)(2v_d\!\cdot\!\ell_3 + v_d\!\cdot\!\ell_4)$ | $v_{bc}$ | G4c |
| $A_3$ | $B_1$ | $v_{cd}\,v_a\!\cdot\!(\ell_4-\ell_3)\;(2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)$ | $v_{cd}$ | G3 |
| $A_3$ | $B_2$ | $v_{ad}\,(2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)(-v_c\!\cdot\!\ell_3 - 2v_c\!\cdot\!\ell_4)$ | $v_{ad}$ | G4b |
| $A_3$ | $B_3$ | $v_{ac}\,(2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)(2v_d\!\cdot\!\ell_3 + v_d\!\cdot\!\ell_4)$ | $v_{ac}$ | G4a |

### Step 5：同类项合并 → 4 组

按 $v_{ij}$ 前因子的类型，9 项合并为 4 组。**以下为 FeynCalc 直接计算结果**（`verify/v01-lorentz-decomposition.m`），与 effective-vertices.md §3b 的 G2/G3/G4a-d 不同。

| 组 | $v_{ij}$ 前因子 | 表达式 | 验证 |
|---|:---:|------|:---:|
| **G1** | $v_{ab}v_{cd}$ | $v_{ab}v_{cd}\,(\ell_2-\ell_1)\!\cdot\!(\ell_4-\ell_3)$ | ✅ |
| **G2** | $v_{ab}$ | $v_{ab}\bigl[v_c\!\cdot\!\ell_3\,v_d\!\cdot\!\ell_1 + 2v_c\!\cdot\!\ell_4\,v_d\!\cdot\!\ell_1 - v_c\!\cdot\!\ell_3\,v_d\!\cdot\!\ell_2 - 2v_c\!\cdot\!\ell_4\,v_d\!\cdot\!\ell_2 - 2v_c\!\cdot\!\ell_1\,v_d\!\cdot\!\ell_3 + 2v_c\!\cdot\!\ell_2\,v_d\!\cdot\!\ell_3 - v_c\!\cdot\!\ell_1\,v_d\!\cdot\!\ell_4 + v_c\!\cdot\!\ell_2\,v_d\!\cdot\!\ell_4\bigr]$ | ✅ |
| **G3** | $v_{cd}$ | $v_{cd}\bigl[2(v_a\!\cdot\!\ell_4-v_a\!\cdot\!\ell_3)v_b\!\cdot\!\ell_1 + (v_a\!\cdot\!\ell_4-v_a\!\cdot\!\ell_3)v_b\!\cdot\!\ell_2 + v_a\!\cdot\!\ell_1(v_b\!\cdot\!\ell_3-v_b\!\cdot\!\ell_4) + 2v_a\!\cdot\!\ell_2(v_b\!\cdot\!\ell_3-v_b\!\cdot\!\ell_4)\bigr]$ | ✅ |
| **G4a** | $v_{ac}$ | $v_{ac}\bigl[4v_b\!\cdot\!\ell_1 v_d\!\cdot\!\ell_3 + 2v_b\!\cdot\!\ell_2 v_d\!\cdot\!\ell_3 + 2v_b\!\cdot\!\ell_1 v_d\!\cdot\!\ell_4 + v_b\!\cdot\!\ell_2 v_d\!\cdot\!\ell_4\bigr]$ | ✅ |
| **G4b** | $v_{ad}$ | $v_{ad}\bigl[-2v_b\!\cdot\!\ell_1 v_c\!\cdot\!\ell_3 - v_b\!\cdot\!\ell_2 v_c\!\cdot\!\ell_3 - 4v_b\!\cdot\!\ell_1 v_c\!\cdot\!\ell_4 - 2v_b\!\cdot\!\ell_2 v_c\!\cdot\!\ell_4\bigr]$ | ✅ |
| **G4c** | $v_{bc}$ | $v_{bc}\bigl[-2v_a\!\cdot\!\ell_1 v_d\!\cdot\!\ell_3 - 4v_a\!\cdot\!\ell_2 v_d\!\cdot\!\ell_3 - v_a\!\cdot\!\ell_1 v_d\!\cdot\!\ell_4 - 2v_a\!\cdot\!\ell_2 v_d\!\cdot\!\ell_4\bigr]$ | ✅ |
| **G4d** | $v_{bd}$ | $v_{bd}\bigl[v_a\!\cdot\!\ell_1 v_c\!\cdot\!\ell_3 + 2v_a\!\cdot\!\ell_2 v_c\!\cdot\!\ell_3 + 2v_a\!\cdot\!\ell_1 v_c\!\cdot\!\ell_4 + 4v_a\!\cdot\!\ell_2 v_c\!\cdot\!\ell_4\bigr]$ | ✅ |

### Step 6：完整 $V_{2\text{TGV}}$ Lorentz 部分

$$
\boxed{
\begin{aligned}
\mathcal{L}(a,b,c,d;\ell_1,\ell_2,\ell_3,\ell_4) &= \text{G1} + \text{G2} + \text{G3} + \text{G4a} + \text{G4b} + \text{G4c} + \text{G4d}.
\end{aligned}}
$$

其中 $v_{ij} \equiv v_i\!\cdot\!v_j$，$v_i\!\cdot\!\ell_j$ 为向量内积。

> **⚠️ 与 effective-vertices.md 的差异**：effective-vertices.md §3b 的 G2/G3/G4a-d 形式与 FeynCalc 直接计算结果不一致。G4a-d 缺少含 $\ell_1,\ell_2$（不含 $\ell_3,\ell_4$）的项。本页以 FeynCalc 直接结果为准。effective-vertices.md 待修正。

### Step 7：退化自动性

| 退化条件 | 消失的组 | 说明 |
|---------|---------|------|
| $a=b$（$v_{ab}=0$） | G1、G2 | $v_{ab}$ 前因子为零 |
| $c=d$（$v_{cd}=0$） | G1、G3 | $v_{cd}$ 前因子为零 |
| $a=c$（$v_{ac}=0$） | G4a | $v_{ac}$ 前因子为零 |
| $a=d$（$v_{ad}=0$） | G4b | $v_{ad}$ 前因子为零 |
| $b=c$（$v_{bc}=0$） | G4c | $v_{bc}$ 前因子为零 |
| $b=d$（$v_{bd}=0$） | G4d | $v_{bd}$ 前因子为零 |

> **原则**：退化不引入新因子、不改变 $C$ 的结构、不要求分支逻辑。$v_{ij}=0$ 时对应项自动消失。

---

## §2a-b: $I(a,b,c,d)$ 动量空间完整积分

> 验证：`verify/v01-lorentz-decomposition.m`（Lorentz 部分）
> 状态：⚠️ 待独立验证（§2b 的整体结构尚未 CAS 验证，仅 Lorentz 部分已验证）

### Feynman 规则表

| 元素 | 表达式 | 说明 |
|------|--------|------|
| Wilson 线顶点 | $ig\,v_{e_i}^{\mu_i}\,t^{A_i}$ | 外腿 $i$ 附着在边 $e_i$ 上，$e_i \in \{a,b,c,d\}$ |
| Wilson 线传播子 | $\frac{-i\,g_{\mu_i\nu_i}\,\delta^{A_i B_i}}{\ell_i^2+i0}$ | 每条外腿的标准标量传播子 |
| 内部传播子 | $\frac{-i\,g_{\nu_5\nu_6}\,\delta^{B_5 B_6}}{\ell_5^2+i0}$ | TGV A 与 TGV B 之间的胶子 |
| 三胶子顶点 | $g\,f^{ABC}\,\Gamma^{\text{kin}}_{\mu\nu\rho}(k_1,k_2,k_3)$ | $g$ 为 QCD 耦合常数 |
| Fourier 相位 | $e^{-i\ell_i\cdot X_{e_i}(s_i)}$ | $X_{e_i}(s_i) = V_{e_i} + s_i v_{e_i}$ |
| 动量守恒 | $\delta^D(\ell_1+\ell_2+\ell_5)\,\delta^D(\ell_3+\ell_4-\ell_5)$ | 强制 $\ell_5 = -(\ell_1+\ell_2) = \ell_3+\ell_4$ |

### 运动学约定（N=4 四边形）

| 量 | 值 | 说明 |
|----|-----|------|
| $v_{14} = v_{23}$ | $v_{\text{ref}}$ | 参考标度 |
| $v_{12} = v_{34}$ | $z\,v_{\text{ref}}$ | 共形比 $z$ |
| $v_{13} = v_{24}$ | $-(1+z)\,v_{\text{ref}}$ | 由 $\sum v_i = 0$ 确定 |
| $v_i^2$ | $0$ | 所有边 lightlike |
| $V_i$ | $\sum_{k=1}^{i-1} v_k$ | $V_1 = \mathbf{0}$ |
| $X_i(s)$ | $V_i + s\,v_i$ | 边 $i$ 上参数 $s \in [0,1]$ 处的位置 |

### 边附着与 Schwinger 参数

四条外腿分别附着在边 $a,b,c,d$ 上（$a,b$ 属 TGV A，$c,d$ 属 TGV B）：

| 外腿 $i$ | 附着边 $e_i$ | 位置 | Schwinger 参数 | 积分域 |
|:---:|:---:|------|------|------|
| 1 | $a$ | $X_a(s_a) = V_a + s_a v_a$ | $s_a$ | 由边 $a$ 上的腿数决定 |
| 2 | $b$ | $X_b(s_b) = V_b + s_b v_b$ | $s_b$ | 由边 $b$ 上的腿数决定 |
| 3 | $c$ | $X_c(s_c) = V_c + s_c v_c$ | $s_c$ | 由边 $c$ 上的腿数决定 |
| 4 | $d$ | $X_d(s_d) = V_d + s_d v_d$ | $s_d$ | 由边 $d$ 上的腿数决定 |

> 积分域：若边 $e$ 上附着 $m$ 条腿，则 Schwinger 参数有序 $0 \le s_{e,1} \le \cdots \le s_{e,m} \le 1$。具体 bound 由 qgraf-pipeline 的 `orbit_definitions.json` 给出。

### 完整积分表达式

将 $V_{2\text{TGV}}$ 有效顶点（§2c Step 6）代入，消去 $\ell_5$ 和 $\delta$ 函数：

$$
\boxed{
\begin{aligned}
I(a,b,c,d) &= \int_0^1 ds_a ds_b ds_c ds_d \;
\prod_{i=1}^4 (ig\,v_{e_i}^{\mu_i} t^{A_i}) \\
&\quad \times \int \frac{d^D\ell_1}{(2\pi)^D}\,\frac{-ig_{\mu_1\nu_1}\,\delta^{A_1 B_1}}{\ell_1^2+i0}\,e^{-i\ell_1\cdot X_a(s_a)} \\
&\quad \times \int \frac{d^D\ell_2}{(2\pi)^D}\,\frac{-ig_{\mu_2\nu_2}\,\delta^{A_2 B_2}}{\ell_2^2+i0}\,e^{-i\ell_2\cdot X_b(s_b)} \\
&\quad \times \int \frac{d^D\ell_3}{(2\pi)^D}\,\frac{-ig_{\mu_3\nu_3}\,\delta^{A_3 B_3}}{\ell_3^2+i0}\,e^{-i\ell_3\cdot X_c(s_c)} \\
&\quad \times \int \frac{d^D\ell_4}{(2\pi)^D}\,\frac{-ig_{\mu_4\nu_4}\,\delta^{A_4 B_4}}{\ell_4^2+i0}\,e^{-i\ell_4\cdot X_d(s_d)} \\
&\quad \times \int \frac{d^D\ell_5}{(2\pi)^D}\,\frac{-ig_{\nu_5\nu_6}\,\delta^{B_5 B_6}}{\ell_5^2+i0} \\
&\quad \times g\,f^{B_1 B_2 B_5}\,\Gamma^{\text{kin}}_{\nu_1\nu_2\nu_5}(-\ell_1,-\ell_2,-\ell_5) \\
&\quad \times g\,f^{B_3 B_4 B_6}\,\Gamma^{\text{kin}}_{\nu_3\nu_4\nu_6}(-\ell_3,-\ell_4,+\ell_5) \\
&\quad \times (2\pi)^D\delta^D(\ell_1+\ell_2+\ell_5)\;(2\pi)^D\delta^D(\ell_3+\ell_4-\ell_5).
\end{aligned}}
$$

### 消去 $\ell_5$ 后的紧凑形式

$\delta$ 函数积分退化为代入 $\ell_5 = -(\ell_1+\ell_2) = \ell_3+\ell_4$。利用 $g_{\mu_i\nu_i}\delta^{A_i B_i}$ 缩并后，颜色与 Lorentz 部分分离：

$$
\boxed{
\begin{aligned}
I(a,b,c,d) &= \int_0^1 ds_a ds_b ds_c ds_d \;
\bigl(t^{A_1}t^{A_2}t^{A_3}t^{A_4}\,f^{A_a A_b B}f^{A_c A_d B}\bigr) \\
&\quad \times (ig)^4 g^2 \int \frac{d^D\ell_1\cdots d^D\ell_4}{(2\pi)^{4D}} \\
&\quad \times \prod_{i=1}^4 \frac{-i\,v_{e_i}^{\nu_i}}{\ell_i^2+i0}\,e^{-i\ell_i\cdot X_{e_i}(s_i)}
   \times \frac{-i}{(\ell_1+\ell_2)^2+i0} \\
&\quad \times \mathcal{L}(a,b,c,d;\ell_1,\ell_2,\ell_3,\ell_4),
\end{aligned}}
$$

其中 $\mathcal{L} = \text{G1} + \text{G2} + \text{G3} + \text{G4a} + \text{G4b} + \text{G4c} + \text{G4d}$（§2c Step 5）。

### 积分结构总览

| 分量 | 内容 | 后续处理 |
|------|------|---------|
| 颜色 | $t^{A_1}t^{A_2}t^{A_3}t^{A_4}\,f^{A_a A_b B}f^{A_c A_d B}$ | §3：FeynCalc `SUNSimplify` 求值 |
| 耦合 | $(ig)^4 g^2 = g^6$（$i^4=1$） | §9b：归入 $C$ 因子 |
| 传播子 | $\prod_{i=1}^4 \frac{-i}{\ell_i^2+i0} \times \frac{-i}{(\ell_1+\ell_2)^2+i0}$ | §6-§7：Fourier 变换 + Symanzik 参数化 |
| $v_{e_i}^{\nu_i}$ | 4 条 Wilson 线向量 | §4-§5：与 $\mathcal{L}$ 缩并 → $\ell \to i\partial_X$ 算子 |
| $\mathcal{L}$ | 7 组 Lorentz 项（G1-G4d） | §4-§5：按 $\ell$ 结构分类 → $\mathcal{O}$ 算子 |
| $e^{-i\ell_i\cdot X_i}$ | Fourier 相位 | §6：FT 产生位置空间传播子 |
| $ds$ 积分 | $\int_0^1 ds_a ds_b ds_c ds_d$ | §9c-1 步骤 9：超立方变换 |

---

## §3: 颜色因子

> 验证：`verify/v02-color-factors.m`
> 状态：✅ S/U 型通过，⚠️ T 型结果依赖轮廓排序

### Step 1：颜色结构

颜色来自三部分：Wilson 线生成元 $t^{A_i}$（外腿顶点）、TGV 结构常数 $f^{A_a A_b B}f^{A_c A_d B}$（有效顶点）、传播子 $\delta^{A_i B_i}$。缩并后：

$$
\boxed{\text{color} = \text{Tr}[T^{A_1} T^{A_2} T^{A_3} T^{A_4}] \; f^{A_a A_b B} f^{A_c A_d B}.}
$$

其中 $T^{A_i}$ 的排序沿 Wilson 圈轮廓（边 $1 \to 2 \to 3 \to 4$），同一边内的腿按 Schwinger 参数 $s_{e,1} < s_{e,2} < \cdots$ 排序。

**配对固定**：TGV A 的腿在边 $(a,b)$（颜色指标 $A_1, A_2$），TGV B 的腿在边 $(c,d)$（颜色指标 $A_3, A_4$）。配对由边分布唯一确定，无 S/T/U 道的选择自由度。

### Step 2：N=4 边 3,4 上的颜色求值

对于 N=4 四边形，外腿仅附着在边 3 和边 4 上（$a,b,c,d \in \{3,4\}$）。

| 边分布 | 配对 | 轮廓排序 | 颜色因子 | 验证 |
|--------|------|---------|---------|:---:|
| $(3,3,4,4)$ | $(33)(44)$ | $A_1,A_2,A_3,A_4$ | $-\frac{1}{4}C_A^2 C_F$ | ✅ |
| $(3,4,3,4)$ | $(34)(34)$ | $A_3,A_1,A_2,A_4$ | $+\frac{1}{4}C_A^2 C_F$ | ✅ |
| $(3,4,4,3)$ | $(34)(43)$ | — | $0$ | ✅ |

> **⚠️ 轮廓排序敏感性**：颜色因子取决于 Schwinger 参数排序。例如 T 型 $(3,4,3,4)$：
>
> | 轮廓排序 | 颜色因子 | 条件 |
> |---------|---------|------|
> | $A_1,A_3,A_2,A_4$ | $0$ | $s_{3,1} < s_{3,3},\; s_{4,2} < s_{4,4}$ |
> | $A_3,A_1,A_2,A_4$ | $+\frac{1}{4}C_A^2 C_F$ | $s_{3,3} < s_{3,1},\; s_{4,2} < s_{4,4}$ |
> | $A_1,A_3,A_4,A_2$ | $+\frac{1}{4}C_A^2 C_F$ | $s_{3,1} < s_{3,3},\; s_{4,4} < s_{4,2}$ |
>
> 实际排序由 qgraf-pipeline 输出的 `bounds` 决定。montopar 的 `build_sunt_indices()` 自动按轮廓排序，`SUNSimplify` 给出正确结果。

### Step 3：Jacobi 恒等式

$$
\text{color}_S + \text{color}_T = -\frac{1}{4}C_A^2 C_F + \frac{1}{4}C_A^2 C_F = 0.
$$

U 型 $(3,4,4,3)$ 的颜色等价于 $\text{color}_S + \text{color}_T$（Jacobi 恒等式），故恒为零——非微扰贡献。

### Step 4：提取颜色后的积分

颜色因子为纯数（不含 $s$、$\ell$、$X$），提到积分外：

$$
\boxed{
\begin{aligned}
I(a,b,c,d) &= \text{color} \times (ig)^4 g^2
   \int_0^1 ds_a ds_b ds_c ds_d \\
   &\quad \times \int \frac{d^D\ell_1\cdots d^D\ell_4}{(2\pi)^{4D}}
   \prod_{i=1}^4 \frac{-i\,v_{e_i}^{\nu_i}}{\ell_i^2+i0}\,e^{-i\ell_i\cdot X_{e_i}(s_i)}
   \times \frac{-i}{(\ell_1+\ell_2)^2+i0} \\
   &\quad \times \mathcal{L}(a,b,c,d;\ell_1,\ell_2,\ell_3,\ell_4).
\end{aligned}}
$$

| 因子 | 值 | 后续处理 |
|------|-----|---------|
| ① color | 依赖边分布 + 轮廓排序 | §9b 归入 $C$ |
| 耦合 | $(ig)^4 g^2 = g^6$ | §9b 归入 $C$ |
| 传播子 + 相位 | $\prod \frac{-i}{\ell_i^2} e^{-i\ell_i\cdot X_i} \times \frac{-i}{(\ell_1+\ell_2)^2}$ | §6-§7 FT + Symanzik |
| $\mathcal{L}$ | G1-G4d（7 组） | §4-§5 $\ell\to i\partial_X$ |

---

## §4: Lorentz 缩并与 $\ell$ 结构分类

> 验证：`verify/v03-operator-table.m`（待建）
> 状态：✅ 分类完成，待 CAS 验证算子表

### Step 1：传播子缩并

$g_{\mu_i\nu_i}$ 已将 Wilson 线向量 $v_{e_i}^{\mu_i}$ 传递到内部 Lorentz 指标 $\nu_i$。$\mathcal{L}$ 中的 $v_a, v_b, v_c, v_d$ 即这些缩并后的向量。$\mathcal{L}$ 的每一项现在是 $(v_i\!\cdot\!\ell_p)(v_j\!\cdot\!\ell_q)$ 或 $\ell_p\!\cdot\!\ell_q$ 的形式。

### Step 2：$\mathcal{L}$ 的因子化形式

将 §2c Step 5 的 G2-G4d 按公共因子整理：

| 组 | $v_{ij}$ 前因子 | 因子化形式 |
|---|:---:|------|
| **G1** | $v_{ab}v_{cd}$ | $v_{ab}v_{cd}\,(\ell_2-\ell_1)\!\cdot\!(\ell_4-\ell_3)$ |
| **G2** | $v_{ab}$ | $v_{ab}\bigl[(v_d\!\cdot\!\ell_1 - v_d\!\cdot\!\ell_2)(v_c\!\cdot\!\ell_3 + 2v_c\!\cdot\!\ell_4) + (v_c\!\cdot\!\ell_2 - v_c\!\cdot\!\ell_1)(2v_d\!\cdot\!\ell_3 + v_d\!\cdot\!\ell_4)\bigr]$ |
| **G3** | $v_{cd}$ | $v_{cd}\bigl[(v_a\!\cdot\!\ell_4 - v_a\!\cdot\!\ell_3)(2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2) + (v_b\!\cdot\!\ell_3 - v_b\!\cdot\!\ell_4)(v_a\!\cdot\!\ell_1 + 2v_a\!\cdot\!\ell_2)\bigr]$ |
| **G4a** | $v_{ac}$ | $v_{ac}\,(2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)(2v_d\!\cdot\!\ell_3 + v_d\!\cdot\!\ell_4)$ |
| **G4b** | $v_{ad}$ | $-v_{ad}\,(2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2)(v_c\!\cdot\!\ell_3 + 2v_c\!\cdot\!\ell_4)$ |
| **G4c** | $v_{bc}$ | $-v_{bc}\,(v_a\!\cdot\!\ell_1 + 2v_a\!\cdot\!\ell_2)(2v_d\!\cdot\!\ell_3 + v_d\!\cdot\!\ell_4)$ |
| **G4d** | $v_{bd}$ | $v_{bd}\,(v_a\!\cdot\!\ell_1 + 2v_a\!\cdot\!\ell_2)(v_c\!\cdot\!\ell_3 + 2v_c\!\cdot\!\ell_4)$ |

引入简记符号使结构更清晰：

| 符号 | 定义 | 含义 |
|------|------|------|
| $A_1$ | $2v_b\!\cdot\!\ell_1 + v_b\!\cdot\!\ell_2$ | TGV A 侧，$v_b$ 方向的 $\ell$ 组合 |
| $A_2$ | $v_a\!\cdot\!\ell_1 + 2v_a\!\cdot\!\ell_2$ | TGV A 侧，$v_a$ 方向的 $\ell$ 组合 |
| $B_1$ | $2v_d\!\cdot\!\ell_3 + v_d\!\cdot\!\ell_4$ | TGV B 侧，$v_d$ 方向的 $\ell$ 组合 |
| $B_2$ | $v_c\!\cdot\!\ell_3 + 2v_c\!\cdot\!\ell_4$ | TGV B 侧，$v_c$ 方向的 $\ell$ 组合 |
| $C_1$ | $v_d\!\cdot\!\ell_1 - v_d\!\cdot\!\ell_2$ | TGV A 侧，$v_d$ 方向的 $\ell$ 差 |
| $C_2$ | $v_c\!\cdot\!\ell_2 - v_c\!\cdot\!\ell_1$ | TGV A 侧，$v_c$ 方向的 $\ell$ 差 |
| $D_1$ | $v_a\!\cdot\!\ell_4 - v_a\!\cdot\!\ell_3$ | TGV B 侧，$v_a$ 方向的 $\ell$ 差 |
| $D_2$ | $v_b\!\cdot\!\ell_3 - v_b\!\cdot\!\ell_4$ | TGV B 侧，$v_b$ 方向的 $\ell$ 差 |

则 $\mathcal{L}$ 的紧凑形式：

$$
\boxed{
\begin{aligned}
\mathcal{L} &= \underbrace{v_{ab}v_{cd}\,(\ell_2-\ell_1)\!\cdot\!(\ell_4-\ell_3)}_{\text{G1}} \\
&\quad + \underbrace{v_{ab}\,(C_1 B_2 + C_2 B_1)}_{\text{G2}}
   + \underbrace{v_{cd}\,(D_1 A_1 + D_2 A_2)}_{\text{G3}} \\
&\quad + \underbrace{v_{ac}\,A_1 B_1}_{\text{G4a}}
   \underbrace{-\,v_{ad}\,A_1 B_2}_{\text{G4b}}
   \underbrace{-\,v_{bc}\,A_2 B_1}_{\text{G4c}}
   \underbrace{+\,v_{bd}\,A_2 B_2}_{\text{G4d}}.
\end{aligned}}
$$

### Step 3：按 $\ell$ 结构分类

| 类型 | $\ell$ 结构 | 组 | 子项数 | 位置空间算子 |
|------|-------------|-----|:---:|-------------|
| **I** | $\ell_p\!\cdot\!\ell_q$ | G1 | 4 | $-\partial_{X_p}\!\cdot\!\partial_{X_q}$ |
| **II** | $(v_i\!\cdot\!\ell)(v_j\!\cdot\!\ell)$，$v_{ab}$ 或 $v_{cd}$ 前因子 | G2, G3 | 各 4 | $-(v_i\!\cdot\!\partial_X)(v_j\!\cdot\!\partial_X)$ |
| **III** | $(v_i\!\cdot\!\ell)(v_j\!\cdot\!\ell)$，交叉 $v_{ij}$ 前因子 | G4a-d | 各 4 | $-(v_i\!\cdot\!\partial_X)(v_j\!\cdot\!\partial_X)$ |

> 类型 II 和 III 的替换规则相同（都是 $(v\!\cdot\!\ell)(v\!\cdot\!\ell)$ 的二次型），区别仅在前因子和动量下标。总子项数：G1 4 项 + G2 8 项 + G3 8 项 + G4a-d 各 4 项 = **36 项**（算子表 28 行，其中 G2/G3 的 8 行每行含 2 个展开子项）。

---

## §5: $\ell \to i\partial_X$ 替换：微分算子 $\mathcal{O}(a,b,c,d)$

### Step 1：核心恒等式

动量空间中的 $v\!\cdot\!\ell$ 因子作用于相位 $e^{-i\ell\cdot X}$：

$$
\boxed{(v\!\cdot\!\ell)\,e^{-i\ell\cdot X} = i\,v\!\cdot\!\partial_X\bigl[e^{-i\ell\cdot X}\bigr].}
$$

所有 $\ell$ 因子提到积分外，积分内部只剩标量 FT。

### Step 2：三类替换规则

| 类型 | $\ell$ 结构 | 替换为 | $i$ 因子 |
|------|-------------|--------|:---:|
| I | $\ell_p\!\cdot\!\ell_q$ | $-\partial_{X_p}\!\cdot\!\partial_{X_q}$ | $i^2 = -1$ |
| II | $(v_i\!\cdot\!\ell_p)(v_j\!\cdot\!\ell_q)$ | $-(v_i\!\cdot\!\partial_{X_p})(v_j\!\cdot\!\partial_{X_q})$ | $i^2 = -1$ |
| III | $(v_i\!\cdot\!\ell_p)(v_j\!\cdot\!\ell_q)$ | $-(v_i\!\cdot\!\partial_{X_p})(v_j\!\cdot\!\partial_{X_q})$ | $i^2 = -1$ |

### Step 3：定义微分算子 $\mathcal{O}(a,b,c,d)$

将 $\mathcal{L}$ 的所有项中 $\ell \to i\partial_X$，定义：

$$
\boxed{\mathcal{O}(a,b,c,d) \equiv \mathcal{L}\bigl(a,b,c,d;\;\ell_p \to i\partial_{X_p}\bigr).}
$$

**$\mathcal{O}_{\text{I}}$**（类型 I，G1）：

$$
\boxed{\mathcal{O}_{\text{I}} = -v_{ab}v_{cd}\;(\partial_{X_2}-\partial_{X_1})\!\cdot\!(\partial_{X_4}-\partial_{X_3})}
$$

展开为 4 个标量 Laplace 型算子：

| 子项 | $C_t$ | 算子 | $(i,j)$ 对 |
|:---:|------|------|:---:|
| I-1 | $-v_{ab}v_{cd}$ | $-\partial_{X_2}\!\cdot\!\partial_{X_4}$ | $(2,4)$ |
| I-2 | $+v_{ab}v_{cd}$ | $-\partial_{X_2}\!\cdot\!\partial_{X_3}$ | $(2,3)$ |
| I-3 | $+v_{ab}v_{cd}$ | $-\partial_{X_1}\!\cdot\!\partial_{X_4}$ | $(1,4)$ |
| I-4 | $-v_{ab}v_{cd}$ | $-\partial_{X_1}\!\cdot\!\partial_{X_3}$ | $(1,3)$ |

**$\mathcal{O}_{\text{II}}$**（类型 II，G2 + G3）：

G2 的 4 个子项（$v_{ab}$ 前因子）：

| 子项 | $C_t$ | $v_p$ | $D^{(p)}$ | $v_q$ | $D^{(q)}$ |
|:---:|------|------|-----------|------|-----------|
| II-1 | $-v_{ab}$ | $v_d$ | $d_1=+1$ | $v_c$ | $d_3=+1,\;d_4=+2$ |
| II-2 | $+v_{ab}$ | $v_d$ | $d_2=+1$ | $v_c$ | $d_3=+1,\;d_4=+2$ |
| II-3 | $-v_{ab}$ | $v_c$ | $d_1=+1$ | $v_d$ | $d_3=+2,\;d_4=+1$ |
| II-4 | $+v_{ab}$ | $v_c$ | $d_2=+1$ | $v_d$ | $d_3=+2,\;d_4=+1$ |

G3 的 4 个子项（$v_{cd}$ 前因子）：

| 子项 | $C_t$ | $v_p$ | $D^{(p)}$ | $v_q$ | $D^{(q)}$ |
|:---:|------|------|-----------|------|-----------|
| II-5 | $-v_{cd}$ | $v_a$ | $d_3=+1,\;d_4=+2$ | $v_b$ | $d_1=+1$ |
| II-6 | $+v_{cd}$ | $v_a$ | $d_3=+1,\;d_4=+2$ | $v_b$ | $d_2=+1$ |
| II-7 | $-v_{cd}$ | $v_b$ | $d_1=+1,\;d_2=+2$ | $v_a$ | $d_3=+1$ |
| II-8 | $+v_{cd}$ | $v_b$ | $d_1=+1,\;d_2=+2$ | $v_a$ | $d_4=+1$ |

**$\mathcal{O}_{\text{III}}$**（类型 III，G4a-d）：

| 子项 | $C_t$ | $v_p$ | $D^{(p)}$ | $v_q$ | $D^{(q)}$ |
|:---:|------|------|-----------|------|-----------|
| IIIa-1 | $-4v_{ac}$ | $v_b$ | $d_1=+1$ | $v_d$ | $d_3=+1$ |
| IIIa-2 | $-2v_{ac}$ | $v_b$ | $d_2=+1$ | $v_d$ | $d_3=+1$ |
| IIIa-3 | $-2v_{ac}$ | $v_b$ | $d_1=+1$ | $v_d$ | $d_4=+1$ |
| IIIa-4 | $-v_{ac}$ | $v_b$ | $d_2=+1$ | $v_d$ | $d_4=+1$ |
| IIIb-1 | $+2v_{ad}$ | $v_b$ | $d_1=+1$ | $v_c$ | $d_3=+1$ |
| IIIb-2 | $+v_{ad}$ | $v_b$ | $d_2=+1$ | $v_c$ | $d_3=+1$ |
| IIIb-3 | $+4v_{ad}$ | $v_b$ | $d_1=+1$ | $v_c$ | $d_4=+1$ |
| IIIb-4 | $+2v_{ad}$ | $v_b$ | $d_2=+1$ | $v_c$ | $d_4=+1$ |
| IIIc-1 | $+2v_{bc}$ | $v_a$ | $d_1=+1$ | $v_d$ | $d_3=+1$ |
| IIIc-2 | $+4v_{bc}$ | $v_a$ | $d_2=+1$ | $v_d$ | $d_3=+1$ |
| IIIc-3 | $+v_{bc}$ | $v_a$ | $d_1=+1$ | $v_d$ | $d_4=+1$ |
| IIIc-4 | $+2v_{bc}$ | $v_a$ | $d_2=+1$ | $v_d$ | $d_4=+1$ |
| IIId-1 | $-v_{bd}$ | $v_a$ | $d_1=+1$ | $v_c$ | $d_3=+1$ |
| IIId-2 | $-2v_{bd}$ | $v_a$ | $d_2=+1$ | $v_c$ | $d_3=+1$ |
| IIId-3 | $-2v_{bd}$ | $v_a$ | $d_1=+1$ | $v_c$ | $d_4=+1$ |
| IIId-4 | $-4v_{bd}$ | $v_a$ | $d_2=+1$ | $v_c$ | $d_4=+1$ |

> 总计：Type I 4 项 + Type II 8 项 + Type III 16 项 = **28 项**。退化时 $C_t=0$ 的项自动消失。

### Step 4：替换后的 $I(a,b,c,d)$

$\mathcal{O}$ 提到所有 $\ell$ 积分外，内部只剩 5 条标量传播子的 FT：

$$
\boxed{
\begin{aligned}
I(a,b,c,d) &= \text{① color} \times (ig)^4 g^2
   \int_0^1 ds_a ds_b ds_c ds_d \\
   &\quad \times \mathcal{O}(a,b,c,d)\!\left[
   \int \frac{d^D\ell_1\cdots d^D\ell_4}{(2\pi)^{4D}}
   \prod_{i=1}^4 \frac{-i\,e^{-i\ell_i\cdot X_{e_i}(s_i)}}{\ell_i^2+i0}
   \times \frac{-i}{(\ell_1+\ell_2)^2+i0}
   \right].
\end{aligned}}
$$

> **操作顺序**（关键）：先做 FT 将 $\ell$ 积分转为位置空间 $z_A, z_B$ 积分和 Symanzik 分母 $B_5$（§6-§7），再让 $\mathcal{O}$ 作用于 $B_5^n$（§8）。$\mathcal{O}$ 的 $\partial_X$ 通过链式法则产生分子 Pnum。

### Step 5：$\mathcal{O}$ 是二次型

| 特征 | TGV $\mathcal{O}$（线性） | 双 TGV $\mathcal{O}$（二次） |
|------|:---:|:---:|
| $\partial_X$ 阶数 | 1 | 2 |
| 总子项数 | 6 | 36 |
| 链式法则 | $\mathcal{O}(B^n) = n B^{n-1} \mathcal{O}(B)$ | $\mathcal{O}(B^n) = n(n-1)B^{n-2}FG + n B^{n-1}H$ |
| 混合导数 $H$ | $H \equiv 0$ | $H \neq 0$（$v_p\!\cdot\!v_q \neq 0$） |

---

## §6: Fourier 变换：$z_A, z_B$ 解耦

> 状态：✅ 已写入，待 CAS 验证

### Step 1：耦合问题

方括号中的 5 条传播子：

$$
\prod_{i=1}^4 \frac{1}{\ell_i^2+i0} \times \frac{1}{(\ell_1+\ell_2)^2+i0}.
$$

前 4 条各含独立 $\ell_i$，可直接逐个 FT。但第 5 条 $1/(\ell_1+\ell_2)^2$ 耦合 $\ell_1$ 和 $\ell_2$——无法分离。

### Step 2：重新引入 $\ell_5$ 和 $\delta$ 函数

$1/(\ell_1+\ell_2)^2$ 源自内部传播子 $1/\ell_5^2$（$\ell_5 = \ell_1+\ell_2$）。FT 时需显式恢复 $\ell_5$ 以解耦：

$$
\frac{1}{(\ell_1+\ell_2)^2+i0} = \int \frac{d^D\ell_5}{(2\pi)^D}\;
\frac{(2\pi)^D\delta^D(\ell_1+\ell_2-\ell_5)}{\ell_5^2+i0}.
$$

TGV B 的动量守恒给出第二个 $\delta$（在 $V_{2\text{TGV}}$ 中已被隐式使用）：

$$
(2\pi)^D\delta^D(\ell_3+\ell_4+\ell_5).
$$

### Step 3：$\delta$ 函数的 Fourier 表示 → $z_A, z_B$

两个 $\delta$ 函数各引入一个位置空间积分变量：

$$
\boxed{
\begin{aligned}
\delta^D(\ell_1+\ell_2-\ell_5) &= \int \frac{d^Dz_A}{(2\pi)^D}\;e^{i(\ell_1+\ell_2-\ell_5)\cdot z_A}, \\[4pt]
\delta^D(\ell_3+\ell_4+\ell_5) &= \int \frac{d^Dz_B}{(2\pi)^D}\;e^{i(\ell_3+\ell_4+\ell_5)\cdot z_B}.
\end{aligned}}
$$

### Step 4：指数重组：5 条 $\ell$ 完全解耦

$z_A, z_B$ 将每条 $\ell$ 的相位重新组合：

| $\ell$ | 原始相位 | 重组后相位 | 净位移 |
|:---:|---------|---------|------|
| $\ell_1$ | $e^{-i\ell_1\cdot X_1}$ | $e^{-i\ell_1\cdot X_1}\,e^{i\ell_1\cdot z_A}$ | $X_1 - z_A$ |
| $\ell_2$ | $e^{-i\ell_2\cdot X_2}$ | $e^{-i\ell_2\cdot X_2}\,e^{i\ell_2\cdot z_A}$ | $X_2 - z_A$ |
| $\ell_5$ | — | $e^{-i\ell_5\cdot z_A}\,e^{i\ell_5\cdot z_B}$ | $z_A - z_B$ |
| $\ell_3$ | $e^{-i\ell_3\cdot X_3}$ | $e^{-i\ell_3\cdot X_3}\,e^{i\ell_3\cdot z_B}$ | $X_3 - z_B$ |
| $\ell_4$ | $e^{-i\ell_4\cdot X_4}$ | $e^{-i\ell_4\cdot X_4}\,e^{i\ell_4\cdot z_B}$ | $X_4 - z_B$ |

5 条 $\ell$ 完全解耦——每条含一个独立传播子 $1/\ell_i^2$ 和一个独立相位 $e^{-i\ell_i\cdot(\text{位移})}$。

### Step 5：逐个 Fourier 变换

标准标量 FT（$a = D/2-1 = 1-\varepsilon$）：

$$
\boxed{\int\frac{d^D\ell}{(2\pi)^D}\,\frac{-i}{\ell^2+i0}\,e^{-i\ell\cdot x}
= -\frac{\Gamma(a)}{4\pi^{D/2}}\,\frac{1}{(-x^2+i0)^a}.}
$$

5 个积分各给出一条位置空间传播子：

| 传播子 | 位置空间分母 | 物理含义 |
|:---:|------|------|
| $\Delta_1$ | $-(X_1-z_A)^2+i0$ | 边 $a$ 上的点 $X_1$ 到 $z_A$ |
| $\Delta_2$ | $-(X_2-z_A)^2+i0$ | 边 $b$ 上的点 $X_2$ 到 $z_A$ |
| $\Delta_5$ | $-(z_A-z_B)^2+i0$ | TGV A 到 TGV B 的内部传播子 |
| $\Delta_3$ | $-(X_3-z_B)^2+i0$ | 边 $c$ 上的点 $X_3$ 到 $z_B$ |
| $\Delta_4$ | $-(X_4-z_B)^2+i0$ | 边 $d$ 上的点 $X_4$ 到 $z_B$ |

### Step 6：5-点拓扑

```
X_1 ──Δ₁── z_A ──Δ₅── z_B ──Δ₃── X_3
  │                    │
  Δ₂                   Δ₄
  │                    │
 X_2                  X_4
```

双 TGV 的拓扑是两个 3-点星型通过 $z_A$-$z_B$ 桥接。$z_A$ 连接 TGV A 的三条腿（$X_1, X_2, \ell_5$），$z_B$ 连接 TGV B 的三条腿（$X_3, X_4, \ell_5$），$\ell_5$（即 $\Delta_5$）桥接两个子图。

### Step 7：FT 后的 $I(a,b,c,d)$

FT prefactor：5 个传播子各贡献 $-\Gamma(a)/4\pi^{D/2}$，$\delta$ 函数的 $1/(2\pi)^{2D}$ 被 $z$ 积分的分子消去：

$$
\boxed{
\begin{aligned}
I(a,b,c,d) &= \text{① color} \times (ig)^4 g^2
   \left(-\frac{\Gamma(a)}{4\pi^{D/2}}\right)^{\!5} \\
   &\quad \times \int_0^1 ds_a ds_b ds_c ds_d
   \int d^Dz_A\,d^Dz_B \\
   &\quad \times \mathcal{O}(a,b,c,d)\!\left[
   \frac{1}{\Delta_1^a\,\Delta_2^a\,\Delta_5^a\,\Delta_3^a\,\Delta_4^a}
   \right].
\end{aligned}}
$$

> **目标**：方括号内是 5-点标量积分。下一步（§7）用 Symanzik 参数化将其转为 Feynman 参数积分 $B_5^{-(1-3\varepsilon)}$，然后 $\mathcal{O}$ 作用于 $B_5$（§8），产生分子 Pnum 和 $C$ 的 derivative 因子。

---

## §7: Symanzik 参数化：5-点恒等式

> 状态：✅ 已写入，待 CAS 验证

### Step 1：一步 Feynman 参数化

对 5 条传播子做一次 Feynman 参数化（引入 $\alpha_1,...,\alpha_5$，$\delta(1-\sum\alpha)$）：

$$
\frac{1}{\Delta_1^a\Delta_2^a\Delta_5^a\Delta_3^a\Delta_4^a}
= \frac{\Gamma(5a)}{\Gamma(a)^5}\int d\alpha\;\frac{\delta(1-\textstyle\sum\alpha)\,\prod\alpha_i^{a-1}}{Q^{5a}},
$$

其中 $a = D/2-1 = 1-\varepsilon$，$Q = \alpha_1\Delta_1 + \alpha_2\Delta_2 + \alpha_5\Delta_5 + \alpha_3\Delta_3 + \alpha_4\Delta_4$。

### Step 2：矩阵形式

$Q$ 是 $z_A, z_B$ 的二次型。代入 $\Delta_i = -(X_i-z)^2$（$\Delta_5 = -(z_A-z_B)^2$）：

$$
\boxed{Q = -\bigl[\vec{z}^T(\Sigma\otimes I_D)\vec{z} - 2\vec{J}\cdot\vec{z} + C\bigr],}
$$

其中 $\vec{z} = (z_A, z_B)^T$，

$$
\boxed{
\Sigma = \begin{pmatrix}\alpha_1+\alpha_2+\alpha_5 & -\alpha_5 \\[2pt] -\alpha_5 & \alpha_3+\alpha_4+\alpha_5\end{pmatrix},\qquad
\vec{J} = \begin{pmatrix}\alpha_1 X_1+\alpha_2 X_2 \\[2pt] \alpha_3 X_3+\alpha_4 X_4\end{pmatrix},\qquad
C = \sum_{i=1}^4 \alpha_i X_i^2.
}
$$

### Step 3：第一 Symanzik 多项式 $\mathcal{U}$

$$
\boxed{\mathcal{U} = \det\Sigma = (\alpha_1+\alpha_2+\alpha_5)(\alpha_3+\alpha_4+\alpha_5) - \alpha_5^2.}
$$

记 $\Sigma_A \equiv \alpha_1+\alpha_2+\alpha_5$，$\Sigma_B \equiv \alpha_3+\alpha_4+\alpha_5$，则 $\mathcal{U} = \Sigma_A\Sigma_B - \alpha_5^2$。

> **与 TGV（3-点）的对比**：3-点是 $L=1$ 圈，$\Sigma$ 是 $1\times 1$ 矩阵，$\mathcal{U}=1$ 平凡。5-点是 $L=2$ 圈，$\mathcal{U}\neq 1$ 非平凡。

### Step 4：Feynman 分母 $B_5$

配平方后 $Q = (\vec{z}-\vec{z}_0)^T(\Sigma\otimes I_D)(\vec{z}-\vec{z}_0) + \mathcal{F}/\mathcal{U}$。第二 Symanzik 多项式 $\mathcal{F} = C\cdot\mathcal{U} - \vec{J}^T\text{adj}(\Sigma)\vec{J}$。展开后所有 $X_i^2$ 项消去。定义 $B_5 = -\mathcal{F}$（物理正定）：

$$
\boxed{
\begin{aligned}
B_5 &= \alpha_1\alpha_2\,\Sigma_B\,\Delta_{12}
    + \alpha_3\alpha_4\,\Sigma_A\,\Delta_{34}
    + \alpha_5\!\!\sum_{\substack{i\in\{1,2\}\\j\in\{3,4\}}}\!\!\alpha_i\alpha_j\,\Delta_{ij}.
\end{aligned}}
$$

其中 $\Delta_{ij} = -(X_i-X_j)^2+i0$。展开为 6 项：

$$
\boxed{
\begin{aligned}
B_5 &= \alpha_1\alpha_2(\alpha_3+\alpha_4+\alpha_5)\,\Delta_{12}
    + \alpha_3\alpha_4(\alpha_1+\alpha_2+\alpha_5)\,\Delta_{34} \\
    &\quad + \alpha_5\bigl(\alpha_1\alpha_3\Delta_{13} + \alpha_1\alpha_4\Delta_{14}
                      + \alpha_2\alpha_3\Delta_{23} + \alpha_2\alpha_4\Delta_{24}\bigr).
\end{aligned}}
$$

**结构解读**：

| $B_5$ 项 | $\Delta_{ij}$ | 系数 | 物理含义 |
|---------|:---:|------|------|
| 1 | $\Delta_{12}$ | $\alpha_1\alpha_2\Sigma_B$ | $z_A$ 侧同组内积 × 对侧 Feynman 参数和 |
| 2 | $\Delta_{34}$ | $\alpha_3\alpha_4\Sigma_A$ | $z_B$ 侧同组内积 × 对侧 Feynman 参数和 |
| 3-6 | $\Delta_{13},\Delta_{14},\Delta_{23},\Delta_{24}$ | $\alpha_i\alpha_j\alpha_5$ | 跨组积 × $\alpha_5$ 桥接因子 |

### Step 5：高斯积分

配平方后，$2D$ 维高斯积分（$\det(\Sigma\otimes I_D) = \mathcal{U}^D$）：

$$
\boxed{
\int d^{2D}z\;Q^{-5a}
= \frac{-\pi^D\,\Gamma(5a-D)}{\Gamma(5a)}\;\frac{1}{\mathcal{U}^{1+2\varepsilon}\;B_5^{\,5a-D}}.
}
$$

代入 $a = 1-\varepsilon$，$5a-D = 5(1-\varepsilon)-(4-2\varepsilon) = 1-3\varepsilon$：

$$
\int d^{2D}z\;Q^{-5a}
= \frac{-\pi^D\,\Gamma(1-3\varepsilon)}{\Gamma(5-5\varepsilon)}\;\frac{1}{\mathcal{U}^{1+2\varepsilon}\;B_5^{\,1-3\varepsilon}}.
$$

### Step 6：5-点恒等式

将 Feynman 参数化与高斯积分合并：

$$
\boxed{
\begin{aligned}
\int d^Dz_A\,d^Dz_B\;\prod_{i=1}^5\frac{1}{\Delta_i^a}
&= \frac{-\pi^D\,\Gamma(1-3\varepsilon)}{\Gamma(1-\varepsilon)^5}
\int d\alpha\;\frac{\delta(1-\textstyle\sum\alpha)\,\prod\alpha_i^{a-1}}
{\mathcal{U}^{1+2\varepsilon}\;B_5^{\,1-3\varepsilon}}.
\end{aligned}}
$$

其中 $\alpha_i^{a-1} = \alpha_i^{-\varepsilon}$（$a=1-\varepsilon$）。

### Step 7：组装到 $I(a,b,c,d)$

FT prefactor 与 Symanzik prefactor 合并：

$$
\left(-\frac{\Gamma(a)}{4\pi^{D/2}}\right)^{\!5} \times \frac{-\pi^D\,\Gamma(1-3\varepsilon)}{\Gamma(1-\varepsilon)^5}
= \frac{\Gamma(1-3\varepsilon)}{1024\,\pi^{6-3\varepsilon}}.
$$

$$
\boxed{
\begin{aligned}
I(a,b,c,d) &= \text{① color} \times (ig)^4 g^2\,
   \frac{\Gamma(1-3\varepsilon)}{1024\,\pi^{6-3\varepsilon}} \\
   &\quad \times \int_0^1 ds_a ds_b ds_c ds_d
   \int d\alpha\;\delta(1-\textstyle\sum\alpha)\,\prod\alpha_i^{-\varepsilon} \\
   &\quad \times \mathcal{O}(a,b,c,d)\!\left[
   \frac{1}{\mathcal{U}^{1+2\varepsilon}\;B_5^{\,1-3\varepsilon}}
   \right].
\end{aligned}}
$$

> **操作顺序**：$\mathcal{O}$ 的 $\partial_X$ 仅作用于 $B_5$（$\mathcal{U}$ 不含 $X$）。链式法则 $\partial_X B_5^n = n B_5^{n-1}\partial_X B_5$ 产生 Pnum（§8）。

---

## §8: $\mathcal{O}$ 作用于 Symanzik 分母：二次型链式法则

> 状态：✅ 已写入，待 CAS 验证

### Step 1：$\mathcal{U}$ 与 $X$ 无关

$\mathcal{U} = \Sigma_A\Sigma_B - \alpha_5^2$ 仅依赖 Feynman 参数 $\alpha_i$。$\mathcal{O}$ 的 $\partial_X$ 仅作用于 $B_5^{\,n}$（$n = -(1-3\varepsilon)$）：

$$
\mathcal{O}\!\left[\frac{1}{\mathcal{U}^{1+2\varepsilon}\,B_5^{\,1-3\varepsilon}}\right]
= \frac{1}{\mathcal{U}^{1+2\varepsilon}}\;\mathcal{O}\!\left[B_5^{\,n}\right].
$$

### Step 2：系数矩阵 $c_{ij}$

$B_5$ 的 6 项可统一写为 $B_5 = \sum_{i < j} \alpha_i\alpha_j\,c_{ij}\,\Delta_{ij}$，其中：

$$
\boxed{
c_{ij} = \begin{cases}
\Sigma_B \equiv \alpha_3+\alpha_4+\alpha_5, & i,j \in \{1,2\} \\[2pt]
\Sigma_A \equiv \alpha_1+\alpha_2+\alpha_5, & i,j \in \{3,4\} \\[2pt]
\alpha_5, & i \in \{1,2\},\; j \in \{3,4\}.
\end{cases}}
$$

### Step 3：一次导数 $F_i^{(p)}$

对任意向量 $v_p \in \{v_a,v_b,v_c,v_d\}$ 和 $i \in \{1,2,3,4\}$：

$$
\boxed{
\begin{aligned}
F_i^{(p)} &\equiv v_p\!\cdot\!\partial_{X_i} B_5 \\
&= -2\alpha_i\sum_{j \neq i} \alpha_j\,c_{ij}\;v_p\!\cdot\!(X_i-X_j).
\end{aligned}}
$$

$F_i^{(p)}$ 是 $X$ 的线性函数（$\propto A_{\text{ref}}$）。

### Step 4：混合二次导数 $H_{ij}^{(pq)}$

$$
H_{ij}^{(pq)} \equiv v_p\!\cdot\!\partial_{X_i}\bigl(v_q\!\cdot\!\partial_{X_j} B_5\bigr).
$$

利用 $\partial_{X_i}(X_j)^{\mu} = g^{\mu\nu}\delta_{ij}$：

$$
\boxed{
H_{ij}^{(pq)} = 2\alpha_i\alpha_j\,c_{ij}\,(v_p\!\cdot\!v_q) \quad (i \neq j),
}
$$

$$
\boxed{
H_{ii}^{(pq)} = -2\alpha_i\,(v_p\!\cdot\!v_q)\sum_{k \neq i}\alpha_k\,c_{ik}.
}
$$

> **关键**：$H_{ij}^{(pq)}$ 不含 $X$（$\propto A_{\text{ref}}^0$）。$v_p\!\cdot\!v_q \neq 0$ 使 $H \neq 0$——即使 $v_p^2 = 0$。这是双 TGV 与 TGV 的本质区别。

$H_{ij}^{(pq)}$ 满足行和为零：$\sum_{j} H_{ij}^{(pq)} = 0$。

### Step 5：二次型链式法则

$\mathcal{O}$ 的一般项为 $C_t \times (v_p\!\cdot\!D^{(p)})(v_q\!\cdot\!D^{(q)})$，作用在 $B_5^n$ 上：

$$
\boxed{
\begin{aligned}
(v_p\!\cdot\!D^{(p)})(v_q\!\cdot\!D^{(q)})\,B_5^{\,n}
&= n(n-1)\,B_5^{\,n-2}\,F^{(p)}G^{(q)} + n\,B_5^{\,n-1}\,H^{(pq)},
\end{aligned}}
$$

其中：

$$
\boxed{
\begin{aligned}
F^{(p)} &\equiv v_p\!\cdot\!D^{(p)} B_5 = \sum_{i=1}^4 d^{(p)}_i\,F_i^{(p)}, \\[6pt]
G^{(q)} &\equiv v_q\!\cdot\!D^{(q)} B_5 = \sum_{j=1}^4 d^{(q)}_j\,F_j^{(q)}, \\[6pt]
H^{(pq)} &\equiv v_p\!\cdot\!D^{(p)} G^{(q)} = \sum_{i,j=1}^4 d^{(p)}_i d^{(q)}_j\,H_{ij}^{(pq)}.
\end{aligned}}
$$

代入 $n = -(1-3\varepsilon)$：

$$
n(n-1) = (1-3\varepsilon)(2-3\varepsilon),\qquad n = -(1-3\varepsilon).
$$

### Step 6：通用 Pnum 表达式

$\mathcal{O}$ 的完整形式为 $\mathcal{O} = \sum_{t} C_t \times (v_{p_t}\!\cdot\!D^{(p_t)})(v_{q_t}\!\cdot\!D^{(q_t)})$（$t$ 遍历 28 个子项）。定义两个 Pnum 分量：

$$
\boxed{
\begin{aligned}
\text{Pnum}^{(1)} &\equiv (2-3\varepsilon)\sum_{t} C_t\;F^{(p_t)}G^{(q_t)}, \\[6pt]
\text{Pnum}^{(2)} &\equiv -\sum_{t} C_t\;H^{(p_t q_t)}.
\end{aligned}}
$$

则：

$$
\boxed{
\mathcal{O}\!\left[B_5^{\,n}\right]
= (1-3\varepsilon)\left[
   \frac{\text{Pnum}^{(1)}}{B_5^{\,3-3\varepsilon}}
   + \frac{\text{Pnum}^{(2)}}{B_5^{\,2-3\varepsilon}}
\right].
}
$$

### Step 7：$A_{\text{ref}}$ 标度

$B_5 \propto A_{\text{ref}}$，$F_i^{(p)} \propto A_{\text{ref}}$，$H_{ij}^{(pq)} \propto A_{\text{ref}}^0$：

| 分量 | $B_5$ 幂次 | 分子 $A_{\text{ref}}$ | 净 $A_{\text{ref}}$ |
|------|:---:|:---:|:---:|
| $\text{Pnum}^{(1)}/B_5^{3-3\varepsilon}$ | $-(3-3\varepsilon)$ | $F G \propto A_{\text{ref}}^2$ | $A_{\text{ref}}^{-1+3\varepsilon}$ |
| $\text{Pnum}^{(2)}/B_5^{2-3\varepsilon}$ | $-(2-3\varepsilon)$ | $H \propto A_{\text{ref}}^0$ | $A_{\text{ref}}^{-2+3\varepsilon}$ |

---

## §9: $C \times J$ 因子分解

> 状态：✅ 已写入，待 CAS 验证

### Step 1：组装 $I(a,b,c,d)$

将 §8 Step 6 的 $\mathcal{O}(B_5^n)$ 代入 §7 Step 7：

$$
\boxed{
\begin{aligned}
I(a,b,c,d) &= \text{① color} \times g^6\,
   \frac{\Gamma(1-3\varepsilon)}{1024\,\pi^{6-3\varepsilon}}
   \times (1-3\varepsilon)
   \times (-2v_{\text{ref}}+i0)^{3\varepsilon} \\
   &\quad \times \int_0^1 ds_a ds_b ds_c ds_d
   \int d\alpha\;\delta(1-\textstyle\sum\alpha)\,
   \frac{\prod\alpha_i^{-\varepsilon}}{\mathcal{U}^{1+2\varepsilon}} \\
   &\quad \times \left[
   \frac{\text{Pnum}^{(1)}}{B_5^{\,3-3\varepsilon}}
   + \frac{\text{Pnum}^{(2)}}{B_5^{\,2-3\varepsilon}}
   \right].
\end{aligned}}
$$

### Step 2：$C$ 因子分解

提取不依赖 Schwinger 参数和 Feynman 参数的前因子：

$$
\boxed{
\begin{aligned}
C &= \underbrace{\text{Tr}[T^{A_1}T^{A_2}T^{A_3}T^{A_4}]\,f^{A_a A_b B}f^{A_c A_d B}}_{\text{① color}}
   \times \underbrace{(4\pi\alpha_s)^3}_{\text{③ coupling}} \\
   &\quad \times \underbrace{\frac{\Gamma(1-3\varepsilon)}{1024\,\pi^{6-3\varepsilon}}}_{\text{④ ft\_base}}
   \times \underbrace{(1-3\varepsilon)}_{\text{⑤ derivative}}
   \times \underbrace{(-2v_{\text{ref}}+i0)^{3\varepsilon}}_{\text{⑥ v14\_scale}}.
\end{aligned}}
$$

| 因子 | 含义 | 依赖 |
|------|------|------|
| ① color | 颜色迹 + 结构常数 | 边分布 $(a,b,c,d)$ + 轮廓排序 |
| ③ coupling | $g^6 = (4\pi\alpha_s)^3$ | 无 |
| ④ ft_base | FT + Symanzik + 高斯积分的合并 prefactor | 无 |
| ⑤ derivative | 二次型链式法则的 $(1-3\varepsilon)$ | 无（通用） |
| ⑥ v14_scale | $A_{\text{ref}}^{3\varepsilon}$ | 无 |

> **$i$ 因子记账**：$i^4=1$（Wilson 线顶点），$i^2=-1$（$2D$ 高斯），$i^2=-1$（$\ell\to i\partial_X$），净 $i^4 \times i^2 \times i^2 = 1$，$C$ 为实数。

### Step 3：$J$ 积分

$$
\boxed{
\begin{aligned}
J(a,b,c,d) &= \int_0^1 ds_a ds_b ds_c ds_d
   \int d\alpha\;\delta(1-\textstyle\sum\alpha)\,
   \frac{(\alpha_1\alpha_2\alpha_3\alpha_4\alpha_5)^{-\varepsilon}}{\mathcal{U}^{1+2\varepsilon}} \\
   &\quad \times \left[
   \frac{\text{Pnum}^{(1)}}{B_5^{\,3-3\varepsilon}}
   + \frac{\text{Pnum}^{(2)}}{B_5^{\,2-3\varepsilon}}
   \right].
\end{aligned}}
$$

### Step 4：与 TGV 的 $J$ 对比

| 特征 | TGV（3e-3g） | 双 TGV（4e-3g+3g） |
|------|:---:|:---:|
| Schwinger 参数 | $d^3 s$ | $d^4 s$ |
| Feynman 参数 | $d^3\alpha$ | $d^5\alpha$ |
| $\mathcal{U}$ | $=1$（平凡） | $\neq 1$（2-圈非平凡） |
| 分母 | $S^{2-2\varepsilon}$（单项） | **两项**：$B_5^{3-3\varepsilon} + B_5^{2-3\varepsilon}$ |
| 分子 | $\frac{1}{4}\sum\alpha_i\alpha_j W_{ij}$ | $\text{Pnum}^{(1)} + \text{Pnum}^{(2)}$ |
| $\mathcal{O}$ 阶数 | 线性（6 项） | 二次（28 项） |
| 退化 | $v_{ij}=0$ 自动处理 | $v_{ij}=0$ 自动处理 |

### Step 5：最终 $C \times J$ 形式

$$
\boxed{I(a,b,c,d) = C \times J(a,b,c,d).}
$$

### Step 6：操作配方总览

| 步骤 | 输入 | 操作 | 输出 |
|:---:|------|------|------|
| 1 | 边分布 $(a,b,c,d)$ | 计算 $v_{ij}$ 值 → 筛选存活算子项 | 存活项列表 |
| 2 | 存活项 + $X_i$ 位置 | 计算 $v_p\!\cdot\!(X_i-X_j)$ 和 $\Delta_{ij}$ | $\Delta_{ij}$ 表达式 |
| 3 | 存活项 | 逐项计算 $F^{(p)}, G^{(q)}, H^{(pq)}$ | 每项的 $F,G,H$ |
| 4 | 所有 $F,G,H$ | 组装 $\text{Pnum}^{(1)}, \text{Pnum}^{(2)}$ | Pnum 表达式 |
| 5 | Pnum + 运动学 | 组装 $C$（§9 Step 2）和 $J$（§9 Step 3） | $C \times J$ |
| 6 | $s_i$ 变量 | 超立方变换 $s_i \to x_i$（乘 Jacobian） | $x$ 变量 $J$ |

---

## §10: 具体 orbit 示例

> 验证：`verify/v06-explicit-orbit-CJ.m`
> 状态：✅ S 道完成，⚠️ (1,2,3,4) 道需 montopar 代码生成

### orbit 04: $(3,3,4,4)$ S 道（双退化）

**运动学**：$v_a=v_b=v_3$, $v_c=v_d=v_4$, $v_{ab}=v_{cd}=0$, $v_{ac}=v_{ad}=v_{bc}=v_{bd}=z$

**存活算子**：G1=G2=G3=0（$v_{ab}=v_{cd}=0$），仅 G4a-d 存活。合并后：

$$
\mathcal{L}_S = z\,v_3\!\cdot\!(\ell_1-\ell_2)\;v_4\!\cdot\!(\ell_3-\ell_4),
\qquad
\mathcal{O}_S = -z\,(v_3\!\cdot\!(\partial_1-\partial_2))(v_4\!\cdot\!(\partial_3-\partial_4)).
$$

**几何**：

| $X_i$ | 表达式 |
|------|--------|
| $X_1$ | $v_1+v_2+s_{3,1}v_3$ |
| $X_2$ | $v_1+v_2+s_{3,2}v_3$ |
| $X_3$ | $v_1+v_2+v_3+s_{4,1}v_4$ |
| $X_4$ | $v_1+v_2+v_3+s_{4,2}v_4$ |

| $\Delta_{ij}$ | 值 |
|:---:|------|
| $\Delta_{12}$ | $0$（同边，$v_3^2=0$） |
| $\Delta_{34}$ | $0$（同边，$v_4^2=0$） |
| $\Delta_{13}$ | $2(1-s_{3,1})s_{4,1}z$ |
| $\Delta_{14}$ | $2(1-s_{3,1})s_{4,2}z$ |
| $\Delta_{23}$ | $2(1-s_{3,2})s_{4,1}z$ |
| $\Delta_{24}$ | $2(1-s_{3,2})s_{4,2}z$ |

**$B_5$ 和 $\mathcal{U}$**：

$$
\boxed{B_5 = 2\alpha_5 z\,[\alpha_1(1-s_{3,1})+\alpha_2(1-s_{3,2})]\,[\alpha_3 s_{4,1}+\alpha_4 s_{4,2}]},
\qquad
\mathcal{U} = (\alpha_1+\alpha_2+\alpha_5)(\alpha_3+\alpha_4+\alpha_5) - \alpha_5^2.
$$

**Building blocks**：

| 量 | 表达式 |
|----|--------|
| $F^{(p)}$ | $2\alpha_5 z(\alpha_1-\alpha_2)(\alpha_3 s_{4,1}+\alpha_4 s_{4,2})$ |
| $F^{(q)}$ | $-2\alpha_5 z(\alpha_3-\alpha_4)(\alpha_1(1-s_{3,1})+\alpha_2(1-s_{3,2}))$ |
| $H^{(pq)}$ | $2\alpha_5 z(\alpha_1-\alpha_2)(\alpha_3-\alpha_4)$ |

**关键简化**：$F^{(p)}G^{(q)} = -z(\alpha_1-\alpha_2)(\alpha_3-\alpha_4)\,B_5$，即 $\text{Pnum}^{(1)} \propto B_5$！

**Pnum**：

$$
\boxed{
\begin{aligned}
\text{Pnum}^{(1)} &= z^2(2-3\varepsilon)(\alpha_1-\alpha_2)(\alpha_3-\alpha_4)\,B_5, \\[4pt]
\text{Pnum}^{(2)} &= 2z^2\alpha_5(\alpha_1-\alpha_2)(\alpha_3-\alpha_4).
\end{aligned}}
$$

**两项合并**：$\text{Pnum}^{(1)}/B_5^{3-3\varepsilon} + \text{Pnum}^{(2)}/B_5^{2-3\varepsilon} = z^2(\alpha_1-\alpha_2)(\alpha_3-\alpha_4)[z(2-3\varepsilon) + 2\alpha_5]/B_5^{2-3\varepsilon}$。

**完整 $C \times J$**：

$$
\boxed{
\begin{aligned}
C_S &= -\frac{1}{4}C_A^2 C_F \times (4\pi\alpha_s)^3 \times \frac{\Gamma(1-3\varepsilon)}{1024\pi^{6-3\varepsilon}} \times (1-3\varepsilon) \times (-2v_{14}+i0)^{3\varepsilon}, \\[16pt]
J_S &= z^2(\alpha_1-\alpha_2)(\alpha_3-\alpha_4) \int_0^1 ds_{3,1}\int_0^{s_{3,2}} ds_{3,2}\int_0^1 ds_{4,1}\int_0^{s_{4,2}} ds_{4,2} \\
    &\quad \times \int d^5\alpha\;\delta(1-{\textstyle\sum}\alpha)\,
       \frac{(\alpha_1\alpha_2\alpha_3\alpha_4\alpha_5)^{-\varepsilon}}{\mathcal{U}^{1+2\varepsilon}}
       \frac{z(2-3\varepsilon) + 2\alpha_5}{B_5^{\,2-3\varepsilon}}.
\end{aligned}}
$$

### orbit 25: $(1,2,3,4)$ 非退化

**运动学**：

| $v_{ij}$ | 值 |
|:---:|:---:|
| $v_{ab}=v_{12}$ | $z$ |
| $v_{cd}=v_{34}$ | $z$ |
| $v_{ac}=v_{13}$ | $-(1+z)$ |
| $v_{ad}=v_{14}$ | $1$ |
| $v_{bc}=v_{23}$ | $1$ |
| $v_{bd}=v_{24}$ | $-(1+z)$ |

所有 7 组存活，28 个算子项全部贡献。

**几何**：

| $X_i$ | 表达式 |
|------|--------|
| $X_1$ | $s_1 v_1$ |
| $X_2$ | $v_1 + s_2 v_2$ |
| $X_3$ | $v_1+v_2 + s_3 v_3$ |
| $X_4$ | $v_1+v_2+v_3 + s_4 v_4$ |

| $\Delta_{ij}$ | 值 |
|:---:|------|
| $\Delta_{12}$ | $(1-s_1)s_2 z$ |
| $\Delta_{23}$ | $(1-s_2)s_3$ |
| $\Delta_{34}$ | $(1-s_3)s_4 z$ |
| $\Delta_{14}$ | $(1-s_4)s_1$ |
| $\Delta_{13}$ | $z(1-s_1)(1-s_3) + s_1 s_3$ |
| $\Delta_{24}$ | $(1-s_2)(1-s_4) + z s_2 s_4$ |

**$B_5$ 和 $\mathcal{U}$**：

$$
\boxed{
\begin{aligned}
B_5 &= \alpha_1\alpha_2(\alpha_3+\alpha_4+\alpha_5)\,(1-s_1)s_2 z
    + \alpha_3\alpha_4(\alpha_1+\alpha_2+\alpha_5)\,(1-s_3)s_4 z \\
    &\quad + \alpha_5\alpha_1\alpha_3\,[z(1-s_1)(1-s_3)+s_1 s_3]
    + \alpha_5\alpha_1\alpha_4\,(1-s_4)s_1 \\
    &\quad + \alpha_5\alpha_2\alpha_3\,(1-s_2)s_3
    + \alpha_5\alpha_2\alpha_4\,[(1-s_2)(1-s_4)+z s_2 s_4], \\[8pt]
\mathcal{U} &= (\alpha_1+\alpha_2+\alpha_5)(\alpha_3+\alpha_4+\alpha_5) - \alpha_5^2.
\end{aligned}}
$$

**$C$ 因子**：

$$
C = \text{① color} \times (4\pi\alpha_s)^3 \times \frac{\Gamma(1-3\varepsilon)}{1024\pi^{6-3\varepsilon}} \times (1-3\varepsilon) \times (-2v_{14}+i0)^{3\varepsilon}.
$$

① color 由 FeynCalc 求值（轮廓排序 $A_1,A_2,A_3,A_4$，结果待定）。

**$J$ 积分**：

$$
J = \int_0^1 d^4s \int d^5\alpha\;\delta(1-{\textstyle\sum}\alpha)\,
    \frac{(\alpha_1\alpha_2\alpha_3\alpha_4\alpha_5)^{-\varepsilon}}{\mathcal{U}^{1+2\varepsilon}}
    \left[\frac{\text{Pnum}^{(1)}}{B_5^{3-3\varepsilon}} + \frac{\text{Pnum}^{(2)}}{B_5^{2-3\varepsilon}}\right].
$$

$\text{Pnum}^{(1)},\text{Pnum}^{(2)}$ 由 28 个算子项贡献组装（§5 Step 3 算子表 + §8 Step 6 Pnum 公式），由 montopar 代码自动生成。

| 因子 | S 道 $(3,3,4,4)$ | $(1,2,3,4)$ 道 |
|------|:---:|:---:|
| ① color | $-\frac{1}{4}C_A^2 C_F$ | 需 FeynCalc 求值 |
| ③ coupling | $(4\pi\alpha_s)^3$ | 同 |
| ④ ft_base | $\frac{\Gamma(1-3\varepsilon)}{1024\pi^{6-3\varepsilon}}$ | 同 |
| ⑤ derivative | $(1-3\varepsilon)$ | 同 |
| ⑥ v14_scale | $(-2v_{14}+i0)^{3\varepsilon}$ | 同 |
| 存活算子项 | 4（仅 G4a-d） | 28（全部） |
| $B_5$ 项数 | 4 | 6 |
| $J$ 分母幂次 | $B_5^{2-3\varepsilon}$（合并后） | $B_5^{3-3\varepsilon} + B_5^{2-3\varepsilon}$ |

---

<!-- 推导完成 -->