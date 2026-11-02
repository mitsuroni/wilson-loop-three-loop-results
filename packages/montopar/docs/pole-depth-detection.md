# 极点深度探测算法：从 $J(z,\varepsilon)$ 到 `poledepth`

> 输入 montopar 输出的 $J_{\rm raw}$（canonical 参数积分形式），输出一个整数 `poledepth`（Laurent 展开中 $1/\varepsilon$ 的最高幂）。
> 全程不做积分，仅靠"边界退化轨迹 + 横截性"判定。
> 日期：2026-08-12

---

## 1. 动机与管线位置

montopar 产出的 $J$ 积分形如

$$
I_{\rm canonical} = \int_{[0,1]^N} d^N\vec{x}\;\Big(\prod_{i=1}^N x_i^{\alpha_i}\Big)\Big(\prod_{k=1}^K [P_k(\vec{x})]^{-\beta_k}\Big)\cdot\mathcal{F}(\vec{x})
$$

其中 $\alpha_i = a_i + \varepsilon b_i$，$P_k$ 是 $\vec{x}$ 的多项式，$\mathcal{F}$ 是无极点的有理/超越因子。

dimreduce 在消去积分变量前，若能**先知道**该积分的极点深度 $p$，就能：

1. 预判 dimreduce 是否会发散（$p>0$）→ 是否需要打回 hyperint；
2. 给 hyperint 的 $\varepsilon$-展开设定截断阶（展开到 $\varepsilon^p$ 即可）；
3. 对 montopar 产物做一致性校验（如 2-loop pure 交换应到 $1/\varepsilon^4$）。

本算法是 dimreduce 的**前置诊断器**：读 `montopar.json` 的 `J_raw.expanded`，写 `poledepth` 整数到 `dimreduce.json`。

---

## 2. 算法

### 输入 / 输出

- **输入**：`J_raw.expanded` 字符串，即 $I(z,\varepsilon)$ 的 canonical 形。
- **输出**：整数 `poledepth`（$=0$ 表示有限，无极点）。

### 步骤

```
poledepth(I):
  1. 化简到 canonical 形
     - 对每个 P_k 因式分解，吸收纯 z_1 幂到外因子
     - 合并同类项：得到每个 x_i 的净指数 α_i = a_i + ε b_i
                   以及"非平凡多项式因子" {Q_m(x), β_m}（Q_m 不是单变量幂）

  2. 枚举边界退化轨迹（divergence loci）
     A. 端点轨迹：对每个 i，若 a_i < 0
        → 轨迹 L_i = {x_i = 0}，初始深度 d_i = |a_i|
     B. 因子退化轨迹：对每个 Q_m，求 Q_m 在 [0,1]^N 边界上的零点集
        - 若 Q_m 仅在单变量端点退化（如 Q = 1 - x_j 退化于 x_j = 1）
          → 轨迹 L_m，深度 = ⌈Re(β_m)⌉ 的整数部分
        - 若 Q_m 在 codim ≥ 2 角点退化（如 Q = 1 - x_i x_j 退化于角点 (1,1)）
          → 标记为"角点重叠"，待 step 3 处理，初始深度待定

  3. 重叠判定（同一角点有 >1 个因子退化）
     - 若多个 Q_m / x_i^a 在同一角点同时为零 → 非正常交叉
       做 1 次 primary sector 分解（x_< = x_> · t），用 jacobian 消幂
       重算该角点净深度 = 角点处同时退化的独立因子数
       （这是 sector decomposition step 3+4 的最小版本）
     - 【jacobian 陷阱】单个线性因子在 codim-2 角点
       （如 Q ≈ y_1 + c·y_2）经 jacobian（y_1 dy_1 dt）后
       变成 ∫ y_1^ε dy_1 = 有限 → 无极点！
       所以"因子退化" ≠ "一定有极点"，必须过 sector + jacobian 这一步

  4. 横截性判定
     - 两轨迹涉及不相交变量集 → 横截 → 深度相乘
     - 两轨迹共享变量且共角点 → 已在 step 3 合并
     - 两轨迹共变量但不共角点（如 x_1 = 0 与 x_1 = 1）→ 取 max，不相乘

  5. poledepth = 所有"横截组"深度之积
     （把相互横截的极点源分组，组内取 max，组间相乘）
```

> **推广**：上述 5 步适用于 2e/4e（纯 n-gluon，$a=-1$，全因子可分解为纯幂）。3e-3g（TGV，$a=-2$，双层积分，混合不可约因子）需追加规则 **R1（Simplex）/R2（$a=-m$ 降幂）/R3（混合因子）**，见 §6.1。

### 判定依据（对应知识库）

| 步骤 | 对应方法 | 出处 |
|------|---------|------|
| 2A / 3 | Taylor 端点极点提取 | Weinzierl Ch.10, Eqs.10.15–10.18 |
| 3 角点 sector | 迭代子扇区分解（Hironaka blowup） | Weinzierl Ch.10, Eqs.10.12–10.14 |
| 4 横截相乘 | Weinberg 子图定理推论 | [[concepts/superficial-degree-of-divergence]], [[concepts/weinberg-theorem]] |
| 整体上界 | contour 散度指标 $\omega=1-N_{\rm ex}$ | [[papers/hep-th/dotsenko-vergeles-1980]] Eq.2.13 |
| 零扇区剪枝 | LiteRed 1.4 参数零扇区判据 | [[papers/hep-ph/1310.1145]] §3 |
| 2-loop 物理上限 | 振幅极点到 $1/\varepsilon^4$ | [[concepts/1-epsilon-poles]] |

---

## 3. 例 1：2e-01（1-loop）→ `poledepth = 2`

### 3.1 输入

`workspaces/1loop/2e/2e-01/montopar.json` 的 `J_raw.expanded`：

$$
J = z_1\int_0^1 dx_1\,dx_2\;\big[x_2\,z_1 - x_1 x_2\,z_1\big]^{-1+\varepsilon}
= z_1\int_0^1 dx_1\,dx_2\;\big[x_2\,(1-x_1)\,z_1\big]^{-1+\varepsilon}
$$

### 3.2 手算

**Step 1 化简**：吸收 $z_1$，因式分解 $P = x_2(1-x_1)z_1$：

$$
J = z_1^{\varepsilon}\;\underbrace{\int_0^1 dx_1\,(1-x_1)^{-1+\varepsilon}}_{I_1}\;\underbrace{\int_0^1 dx_2\,x_2^{-1+\varepsilon}}_{I_2}
$$

- 净指数：$\alpha_1 = 0$（$x_1$ 只出现在 $1-x_1$），$\alpha_2 = -1+\varepsilon$ → $a_2=-1,\,b_2=1$
- 非平凡因子：$Q_1 = 1-x_1$，$\beta_1 = 1-\varepsilon$

**Step 2 轨迹**：

| 轨迹 | 源 | 涉及变量 | 深度 |
|------|----|---------|------|
| $L_2$ | $a_2=-1<0$ | $x_2$ | $\|a_2\| = 1$ |
| $L_{Q_1}$ | $Q_1=1-x_1$ 退化于 $x_1=1$ | $x_1$ | $1$ |

**Step 3 重叠**：$L_2$ 与 $L_{Q_1}$ 在角点 $(x_1,x_2)=(1,0)$ 相遇，但两者是**两个不同单变量**的端点退化（$x_1\to1$ 与 $x_2\to0$），各自已是正常交叉，无需额外 blowup。

**Step 4 横截**：$L_2$ 涉及 $\{x_2\}$，$L_{Q_1}$ 涉及 $\{x_1\}$ → 不相交 → **横截 → 相乘**。

**Step 5**：`poledepth = 1 × 1 = 2`。

### 3.3 验证

$$
\int_0^1 y^{-1+\varepsilon}\,dy = \frac{1}{\varepsilon}
\quad\Rightarrow\quad
J \sim \frac{1}{\varepsilon^2}
$$

两个独立端点相乘 → $1/\varepsilon^2$。✓ 经典 1-loop cusp 发散。

---

## 4. 例 2：4e-01（2-loop）→ `poledepth = 4`

### 4.1 输入

`workspaces/2loop/4e/4e-01/montopar.json` 的 `J_raw.expanded`：

$$
J = z_1^2\int_{[0,1]^4} dx_{1..4}\;(x_1 x_2)\,\big[x_2 x_4 z_1(1-x_1 x_3)\big]^{-1+\varepsilon}\,\big[x_2 z_1(1-x_1)\big]^{-1+\varepsilon}
$$

### 4.2 手算

**Step 1 化简**：$P_1 = x_2 x_4 z_1(1-x_1 x_3)$，$P_2 = x_2 z_1(1-x_1)$，分子 $x_1 x_2$，吸收 $z_1$：

$$
J \propto \int_{[0,1]^4} dx_{1..4}\;
\underbrace{x_1}_{\alpha_1=+1}\,
\underbrace{x_2^{-1+2\varepsilon}}_{\alpha_2:\,a_2=-1,\,b_2=2}\,
\underbrace{x_3^{0}}_{\alpha_3=0}\,
\underbrace{x_4^{-1+\varepsilon}}_{\alpha_4:\,a_4=-1,\,b_4=1}\,
\underbrace{(1-x_1 x_3)^{-1+\varepsilon}}_{Q_1}\,
\underbrace{(1-x_1)^{-1+\varepsilon}}_{Q_2}
$$

**Step 2 轨迹**：

| 轨迹 | 源 | 涉及变量 | 深度 |
|------|----|---------|------|
| $L_2$ | $a_2=-1<0$ | $x_2$ | $1$ |
| $L_4$ | $a_4=-1<0$ | $x_4$ | $1$ |
| $L_{Q_2}$ | $Q_2=1-x_1$ 退化于 $x_1=1$ | $x_1$ | $1$ |
| $L_{Q_1}$ | $Q_1=1-x_1 x_3$ 退化于角点 $(x_1,x_3)=(1,1)$ | $\{x_1,x_3\}$ | 待定 |

**Step 3 重叠（关键）**：$L_{Q_1}$ 与 $L_{Q_2}$ 在角点 $(x_1,x_3)=(1,1)$ **同时退化** → 非正常交叉。做 1 次 primary sector：

令 $y_1 = 1-x_1,\;y_3 = 1-x_3$。角点附近：

$$
Q_2 = y_1,\qquad Q_1 = 1-(1-y_1)(1-y_3) = y_1+y_3-y_1 y_3 \approx y_1+y_3
$$

取 sector $y_3 = y_1 t_3$（jacobian $y_1$）：

$$
Q_1^{-1+\varepsilon} Q_2^{-1+\varepsilon}\cdot y_1 = [y_1(1+t_3)]^{-1+\varepsilon}\,y_1^{-1+\varepsilon}\cdot y_1 = y_1^{-1+2\varepsilon}\,(1+t_3)^{-1+\varepsilon}
$$

- 净 $y_1$ 指数 $-1+2\varepsilon$ → $a=-1,b=2$ → **深度 1**（极点 $1/(2\varepsilon)$），$t_3$ 积分有限。

另一 sector $y_1 = y_3 t_1$（jacobian $y_3$）：

$$
Q_1 \approx y_3(1+t_1),\quad Q_2 = y_3 t_1
\;\Rightarrow\;
y_3^{-1+2\varepsilon}\,t_1^{-1+\varepsilon}\,(1+t_1)^{-1+\varepsilon}
$$

- $y_3^{-1+2\varepsilon}$ → $1/(2\varepsilon)$，$t_1^{-1+\varepsilon}$ → $1/\varepsilon$ → **深度 2**。

取 max → 角点 $(1,1)$ 净深度 = **2**（两个独立因子同时退化的重叠贡献）。

> **为什么 jacobian 陷阱没在这里咬人**：4e-01 在该角点有**两个独立因子** $Q_1,Q_2$ 同时退化，jacobian 只消 1 个幂，仍剩一个 $-1+\varepsilon$ 端点。对比 §5 的 2e-02：单个线性因子在 codim-2 角点会被 jacobian 消成 $y^\varepsilon$ → 无极点。

**Step 4 横截**：

- $L_2$：$\{x_2\}$
- $L_4$：$\{x_4\}$
- $L_{Q_1 Q_2}$（已合并的角点）：$\{x_1,x_3\}$

三组变量集两两不相交 → **三重横截 → 相乘**。

**Step 5**：`poledepth = 1 × 1 × 2 = 4`。

### 4.3 验证

4e 是 2-loop 纯传播子交换。知识库 [[concepts/1-epsilon-poles]] 指出 2-loop 振幅极点到 $1/\varepsilon^4$；每多一条独立 Wilson-line 交换加一阶横截极点。1-loop（2e-01）深度 2，2-loop（4e-01）深度 4，符合 $2L$ 规律。✓

---

## 5. 补充例：2e-02 与 jacobian 陷阱

`workspaces/1loop/2e/2e-02/montopar.json`：

$$
J = -(1+z_1)\int_0^1 dx_1\,dx_2\;\big[(1-x_1)(1-x_2)+x_1 x_2 z_1\big]^{-1+\varepsilon}
$$

记 $P = (1-x_1)(1-x_2)+x_1 x_2 z_1$。对物理 $z_1>0$，$P\geq0$ 且仅在角点 $(1,0)$、$(0,1)$ 为零——均为 codim-2 单因子退化。

**角点 $(1,0)$ 附近**：令 $y_1=1-x_1\to0,\,y_2=x_2\to0$：

$$
P = y_1(1-y_2)+(1-y_1)y_2 z_1 \approx y_1 + z_1 y_2
$$

做 sector $y_2 = y_1 t$（jacobian $y_1$）：

$$
\int_0^\epsilon dy_1\int_0^1 y_1\,dt\;[y_1(1+z_1 t)]^{-1+\varepsilon}
= \int_0^\epsilon dy_1\,y_1^{\varepsilon}\;\underbrace{\int_0^1(1+z_1 t)^{-1+\varepsilon}dt}_{\text{有限}}
= \frac{\epsilon^{\varepsilon+1}}{\varepsilon+1}\cdot(\text{有限})
\;\to\;\text{有限}
$$

jacobian 把 $y_1^{-1+\varearepsilon}$ 消成 $y_1^{\varepsilon}$ → **无极点**。另一 sector 对称，同样有限。角点 $(0,1)$ 同理。

因此 2e-02 的 `poledepth = 0`（有限）。这与 2e-01 的深度 2 形成对照，说明算法能区分"真端点发散"与"被 jacobian 调节的 codim-2 角点"——正是 sector decomposition 存在的理由。

> ⚠️ 2e-02 的深度 0 是**算法预测**，建议用数值（pySecDec / FIESTA）交叉验证后再定论；本算法的价值在于**不做积分**就给出预测与判定路径。

---

## 6. 3e-3g 系列：双层积分与幂发散

3e-3g（TGV，三胶子顶点，2-loop）暴露了 2e/4e（pure n-gluon）不曾出现的三个新特征，迫使算法推广。先给统一结构，再逐 orbit 分析。

### 6.0 统一结构

$$
J = \int_0^1 d^3x\;\int_\triangle d^3\alpha\,\delta(1-{\textstyle\sum}\alpha)\;
\frac{\text{Jac}\times\text{Pnum}}{(\alpha_1\alpha_2\alpha_3)^{\varepsilon}\;S^{\,2-2\varearepsilon}}
$$

- **双层积分域**：$x\in[0,1]^3$ 超立方 × $\alpha$ 在 2-simplex 上（$\delta(1-\sum\alpha)$）。
- **$S$ 的指数是 $-(2-2\varearepsilon)$**，即 $a=-2$（2-loop 标志，pure n-gluon 是 $a=-1$）。
- $S$ 与 Pnum 均为 $\alpha$ 与 $x$ 的**耦合多项式**，一般不分解为纯幂。
- 四 orbit 共享同一 $C$（退化不改变 $C$），$z$ 依赖全在 $J$。

### 6.1 三条推广规则（3e-3g 引入）

> **规则 R1（Simplex 端点）**：$\alpha_i$ 在 simplex 边界 $\{\alpha_i=0\}$ 退化时，端点判据 $a_i<0$ 与超立方坐标相同（simplex 经 primary sector 映到超立方，但端点发散在映射前即可读出）。注意 simplex 顶点 $\alpha_i=1$（即其余 $\alpha_{j\neq i}=0$）是另一类边界，单独判。

> **规则 R2（$a=-m$ 幂发散，$m>1$）**：因子指数为 $-(m-m\varearepsilon)$ 时（$L$-loop 的 $S$ 给 $m=L$），**不能**用 $|a|=m$ 作深度。在其退化轨迹上做 sector 分解，jacobian（codim $k$）消 $k$ 阶幂，有效指数 $a_{\rm eff}=a+k$：
> - $a_{\rm eff}=-1$ → 该因子贡献深度 1；
> - $a_{\rm eff}>-1$ → jacobian 过度消幂 → 无极点（2e-02 情形）；
> - $a_{\rm eff}<-1$ → 递归再做 sector。

> **规则 R3（混合多项式因子）**：不分解为纯幂的因子（如 $(\alpha_3+\alpha_2 x_3)$），其零轨迹是超曲面而非坐标端点。需先求退化 codim $k$，再走 R2 读深度。这与 2e/4e"全部因子可分解为纯幂"的前提不同。

### 6.2 例 3：3e-3g-01（$(a,b,c)=(3,4,4)$，退化 $b=c$）→ `poledepth = 4`

**输入**：`workspaces/2loop/3e-3g/3e-3g-01/montopar.json`

$$
S = x_2 z\,\alpha_1(1-x_1)(\alpha_3+\alpha_2 x_3),\quad
\text{Pnum} = \tfrac14 z^2(1-x_1)\alpha_1(\alpha_3-\alpha_2),\quad \text{Jac}=x_2
$$

**Step 1 化简合并净指数**（吸收 $z$，$S^{-(2-2\varearepsilon)}$ 拆开，乘分子与 $(\alpha_1\alpha_2\alpha_3)^{-\varearepsilon}$）：

| 变量 | 净指数 $\alpha_i=a_i+\varepsilon b_i$ | 端点轨迹 | 深度 |
|------|------|------|------|
| $x_2$ | $+1-2+2\varearepsilon=-1+2\varearepsilon$ | $\{x_2=0\}$ | 1 |
| $x_1$（$1-x_1$） | $+1-2+2\varearepsilon=-1+2\varearepsilon$ | $\{x_1=1\}$ | 1 |
| $\alpha_1$ | $+1-\varearepsilon-2+2\varearepsilon=-1+\varearepsilon$ | $\{\alpha_1=0\}$ | 1 |
| $\alpha_2,\alpha_3,x_3$ | 混合因子（见下） | — | — |

三个纯幂端点 $a_i=-1<0$，变量集 $\{x_2\},\{x_1\},\{\alpha_1\}$ 两两横截 → 已得 $1\times1\times1=3$。

**Step 2 混合因子 $(\alpha_3+\alpha_2 x_3)^{-2+2\varearepsilon}$（R3）**：

零轨迹（$\alpha_2,\alpha_3\geq0,\,x_3\in[0,1]$）：$\alpha_3+\alpha_2 x_3=0$ iff $\alpha_3=0$ 且 $\alpha_2 x_3=0$ → codim-2 角点 $\{\alpha_3=0,\,x_3=0\}$（$\alpha_2$ 自由）。

sector $x_3\leq\alpha_3$（$x_3=\alpha_3 t$，jacobian $\alpha_3$，R2 降幂）：

$$
\int d\alpha_3\,dx_3\,(\alpha_3+\alpha_2 x_3)^{-2+2\varearepsilon}
\;\to\;
\int d\alpha_3\,\underbrace{\alpha_3\cdot\alpha_3^{-2+2\varearepsilon}}_{\alpha_3^{-1+2\varearepsilon}}\,(1+\alpha_2 t)^{-2+2\varearepsilon}
$$

codim-2 的 jacobian 把 $a=-2$ 降到 $a_{\rm eff}=-1$ → $1/(2\varearepsilon)$ → **深度 1**；$(1+\alpha_2 t)\geq1$ 无发散。另一 sector 对称同深。

**Step 3 横截聚合**：四组变量集 $\{x_2\},\{x_1\},\{\alpha_1\},\{\alpha_3,x_3\}$ 两两不相交 → 四重横截 → `poledepth = 1×1×1×1 = 4`。

与 4e-01 一致，符合 2-loop $1/\varepsilon^4$。✓（待续：02/03/04）

### 6.3 例 4：3e-3g-02（$(a,b,c)=(3,3,4)$，退化 $a=b$）→ `poledepth = 3`

**输入**：`workspaces/2loop/3e-3g/3e-3g-02/montopar.json`

$$
S = x_2 z\,\alpha_3\,R,\quad R = \alpha_1(1-x_1 x_3)+\alpha_2(1-x_1),\quad
\text{Pnum} = \tfrac14 x_2 z^2\alpha_3(\alpha_2-\alpha_1),\quad \text{Jac}=x_1
$$

**关键差别**：$S$ 只因式分解成 **3 个横截因子**（$x_2,\alpha_3,R$），而 3e-3g-01 分解成 4 个——因为 02 里 $(1-x_1)$ **缠在 $R$ 内**（$\alpha_2(1-x_1)$ 项），不能单独析出。

**Step 1 纯幂端点**：

| 变量 | 净指数 | 轨迹 | 深度 |
|------|------|------|------|
| $x_2$ | $+1-2+2\varepsilon=-1+2\varepsilon$ | $\{x_2=0\}$ | 1 |
| $\alpha_3$ | $+1-\varepsilon-2+2\varepsilon=-1+\varepsilon$ | $\{\alpha_3=0\}$ | 1 |

（$\alpha_1,\alpha_2$ 纯幂仅 measure 的 $-\varepsilon$ → $a=0$；$x_1$ 纯幂 $+1$ → 均无极点，皆缠在 $R$ 内。）

**Step 2 混合因子 $R^{-2+2\varearepsilon}$（R2/R3）**：

$R=\alpha_1(1-x_1 x_3)+\alpha_2(1-x_1)\geq0$，$R=0$ iff 两非负项同时为零。退化子轨迹：

| 子轨迹 | 条件 | 变量 | codim |
|--------|------|------|-------|
| 顶点 | $\alpha_1=\alpha_2=0$（$\alpha_3=1$） | $\{\alpha_1,\alpha_2\}$ | 2 |
| $x_1x_3$ 角 | $x_1=1,x_3=1$ | $\{x_1,x_3\}$ | 2 |
| 交叉 | $\alpha_1=0,x_1=1$ | $\{\alpha_1,x_1\}$ | 2 |

各 codim-2，jacobian 降 $a=-2\to a_{\rm eff}=-1$ → **各深度 1**。$R$ 是**单因子**，多退化轨迹择一（max）→ $R$ 总贡献 = 1。

> **$x_1x_3$ 角演算**：$y_1=1-x_1,y_3=1-x_3$，$R\approx y_1(\alpha_1+\alpha_2)+\alpha_1 y_3$。sector $y_3=y_1 t$（jac $y_1$）→ $R\approx y_1(\alpha_1+\alpha_2+\alpha_1 t)$，净 $y_1^{-1+2\varearepsilon}$ → $1/(2\varearepsilon)$，深度 1。✓

> **顶点与 $x_1x_3$ 角的交（codim-4）**：$R\approx\alpha_1 y_1+\alpha_1 y_3+\alpha_2 y_1$（4 变量齐二次）。标度 $\vec u=\lambda\hat u$：$R\to\lambda^2,\;d^4u\to\lambda^3d\lambda$ → $\int\lambda^{-1+4\varearepsilon}d\lambda=1/(4\varearepsilon)$ → **深度 1，不是 2**。同一因子的多个退化点共享"奇点预算"，交点不叠加。

**Step 3 横截聚合**：$\{x_2\},\{\alpha_3\},\{\alpha_1,\alpha_2,x_1,x_3\}_{R}$ 两两不相交 → `poledepth = 1×1×1 = 3`。

**与 3e-3g-01 对照**：

| | 3e-3g-01 | 3e-3g-02 |
|--|---------|---------|
| $S$ 横截因子数 | 4（$x_2,(1-x_1),\alpha_1,(\alpha_3+\alpha_2 x_3)$） | 3（$x_2,\alpha_3,R$） |
| 每因子深度 | 1 | 1 |
| `poledepth` | 4 | 3 |

**洞见**：3e-3g 系列的 `poledepth` = $S$ 因式分解出的**横截因子数**（纯幂因子经分子偏置到 $a=-1$ 给 1；不可约混合因子 $a=-2$ 经 jacobian 给 1）。退化结构改变 $S$ 的可分解性 → 改变深度。**单个 $a=-2$ 因子只给深度 1 而非 2**——这是 TGV（单 $S$）与 4e（两个独立 $a=-1$ 传播子）的本质差别：TGV 把两圈揉进一个 $S$，丢一阶极点。（待续：03/04）

### 6.4 例 5：3e-3g-03（$(a,b,c)=(2,4,4)$，退化 $b=c$）→ `poledepth = 2`

**输入**：`workspaces/2loop/3e-3g/3e-3g-03/montopar.json`

$$
S = \alpha_1\,R',\quad R'=\alpha_2 A+\alpha_3 B,\quad
\begin{cases}A=(1-x_1)(1-x_2 x_3)+zx_1 x_2 x_3\\ B=(1-x_1)(1-x_2)+zx_1 x_2\end{cases}
$$
$\text{Pnum}=\tfrac14(1+z)(1-x_1(1+z))\alpha_1(\alpha_2-\alpha_3)$，$\text{Jac}=x_2$。

**关键**：$S$ 只分解成 **2 个横截因子**（$\alpha_1,R'$）——$x_2$ 缠在 $A,B$ 内（不再是 $S$ 的因子），$(1-x_1)$ 也缠在 $A,B$ 里。

**Step 1 纯幂端点**：

| 变量 | 净指数 | 轨迹 | 深度 |
|------|------|------|------|
| $\alpha_1$ | $+1-\varepsilon-2+2\varepsilon=-1+\varepsilon$ | $\{\alpha_1=0\}$ | 1 |
| $x_2$ | $+1$（仅 Jac，$S$ 无纯 $x_2$ 幂） | — | 0 |
| $x_1,x_3,\alpha_2,\alpha_3$ | 缠在 $R'$（或 Pnum 线性因子） | — | 0 |

> 注：$x_2=0$ 时 $A|_{x_2=0}=B|_{x_2=0}=1-x_1$，$R'|_{x_2=0}=\alpha_1(\alpha_2+\alpha_3)(1-x_1)\neq0$ → $x_2=0$ **非**退化轨迹。
> Pnum 的 $(1-x_1(1+z))$、$(\alpha_2-\alpha_3)$ 是**分子零点**（域内部），不产生极点。

**Step 2 混合因子 $R'^{-2+2\varearepsilon}$（R2/R3）**：

$R'=\alpha_2 A+\alpha_3 B\geq0$（$A,B\geq0$ 对物理 $z>0$），$R'=0$ iff $\alpha_2 A=0$ 且 $\alpha_3 B=0$。主导退化轨迹：

| 子轨迹 | 条件 | 变量 | codim | 深度 |
|--------|------|------|-------|------|
| 顶点 | $\alpha_2=\alpha_3=0$（$\alpha_1=1$） | $\{\alpha_2,\alpha_3\}$ | 2 | 1 |
| $x$-线 | $A=B=0$，含直线 $\{x_1=1,x_2=0\}$（任意 $x_3,\alpha$） | $\{x_1,x_2\}$ | 2 | 1 |

> **$x$-线演算**：近 $\{x_1=1,x_2=0\}$，$y_1=1-x_1,y_2=x_2$，$A\approx y_1+zy_2 x_3,\,B\approx y_1+zy_2$ → $R'\approx y_1(\alpha_2+\alpha_3)+y_2 z(\alpha_2 x_3+\alpha_3)$。sector $y_2=y_1 t$（jac $y_1$）→ $y_1^{-1+2\varearex}$ → 深度 1。✓
> **顶点演算**：$R'\approx\alpha_2 A+\alpha_3 B$（$A,B>0$），sector $\alpha_3=\alpha_2 t$ → $\alpha_2^{-1+2\varearex}$ → 深度 1。✓
> 两轨迹交（codim-4）：$R'$ 双线性，标度 $\lambda$ → $\int\lambda^{-1+4\varearepsilon}d\lambda=1/(4\varearepsilon)$ → 深度 1（同因子奇点预算不叠加）。

$R'$ 单因子 → 多轨迹取 max → 深度 1。

**Step 3 横截聚合**：$\alpha_1$-因子（$\{\alpha_1\}$）与 $R'$-因子（$\{\alpha_2,\alpha_3,x_1,x_2\}$）变量集不相交 → `poledepth = 1 × 1 = 2`。

**模式确认**：

| Orbit | $(a,b,c)$ | $S$ 横截因子数 | `poledepth` |
|-------|-----------|---------------|------------|
| 01 | (3,4,4) | 4 | 4 |
| 02 | (3,3,4) | 3 | 3 |
| 03 | (2,4,4) | 2 | 2 |

**模式**：3e-3g 的 `poledepth` = $S$ 因式分解出的**横截因子数**。退化越深（同边附着越多），$S$ 越不可分，深度越低——R3（混合因子）的累积效应。待 04（非退化）验证。（待续：04）

### 6.5 例 6：3e-3g-04（$(a,b,c)=(2,3,4)$，非退化）→ `poledepth = 1`

**输入**：`workspaces/2loop/3e-3g/3e-3g-04/montopar.json`

$$
S = \alpha_1\alpha_3\,C + \alpha_1\alpha_2\,x_2(1-x_1) + \alpha_2\alpha_3\,x_3 z(1-x_2),\quad C=(1-x_1)(1-x_3)+zx_1 x_3,\;\text{Jac}=1
$$

**Step 1 纯幂端点 — 无**：$S$ 三项分别含 $\alpha_1\alpha_3$、$\alpha_1\alpha_2$、$\alpha_2\alpha_3$，**无任何 $\alpha_i$ 或 $x_i$ 是 $S$ 的公因子**（每项缺一个 $\alpha$）。Pnum 同样纠缠。纯幂端点全部 $a_i\geq0$。

**Step 2 不可约混合因子 $S^{-2+2\varearepsilon}$（R2/R3）**：

$S\geq0$（$z>0$，三项非负），$S=0$ iff 三项同时为零。**$S=0$ 处处 codim ≥ 2**（无 codim-1 分量）——验证：候选超曲面 $\{\alpha_1=0\}$ 代入第三项 $\alpha_2\alpha_3 x_3 z(1-x_2)\neq0$（一般）→ $S\neq0$。

codim-2 退化轨迹：

| 轨迹 | 变量 |
|------|------|
| $\{\alpha_1=0,x_3=0\}$ | $\{\alpha_1,x_3\}$ |
| $\{\alpha_1=0,x_2=1\}$ | $\{\alpha_1,x_2\}$ |
| $\{\alpha_3=0,x_2=0\}$ | $\{\alpha_3,x_2\}$ |
| $\{\alpha_3=0,x_1=1\}$ | $\{\alpha_3,x_1\}$ |
| simplex 顶点 | $\{\alpha_i,\alpha_j\}$ |

> **典型演算**（$\{\alpha_3=0,x_2=0\}$）：近该轨迹 $\alpha_3,x_2\to0$，
> $$S \approx \underbrace{(\alpha_1 C+\alpha_2 x_3 z)}_{c>0}\alpha_3+\underbrace{\alpha_1\alpha_2(1-x_1)}_{d\geq0}x_2\approx c\,\alpha_3+d\,x_2$$
> 2 变量线性组合，$a=-2$，codim-2 → sector $x_2=\alpha_3 t$（jac $\alpha_3$）→ $\alpha_3^{-1+2\varearex}$ → **深度 1**。✓
> **交点不叠加**：两轨迹交（4 变量小）$S$ 双线性，标度 $\lambda$ → $\int\lambda^{-1+4\varearex}d\lambda=1/(4\varearex)$ → 深度 1。单因子奇点预算不叠加。

$S$ 单因子 → 多轨迹 max → 深度 1。

**Step 3 横截聚合**：仅一个不可约因子 $S$ → `poledepth = 1`。

**物理解释**：非退化 = 三胶子分附三**不同** Wilson-line 边 → 无两胶子共边的 cusp 增强 → 最低深度。退化轨道（共边越多）$S$ 越可分 → 深度越高。

### 6.6 3e-3g 系列总结

| Orbit | $(a,b,c)$ | 退化 | $S$ 横截因子数 $N_{\rm fac}$ | `poledepth` |
|-------|-----------|------|----------------------------|------------|
| 01 | (3,4,4) | $b=c$ | 4 | **4** |
| 02 | (3,3,4) | $a=b$ | 3 | **3** |
| 03 | (2,4,4) | $b=c$ | 2 | **2** |
| 04 | (2,3,4) | 无 | 1（不可约） | **1** |

**统一规律**：3e-3g 的 `poledepth` = $S$ 因式分解出的**横截因子数** $N_{\rm fac}$。每因子贡献 1（纯幂经分子偏置到 $a=-1$；不可约混合经 jacobian 把 $a=-2$ 降到 $a_{\rm eff}=-1$）。退化（共边）使 $S$ 可分性增强 → $N_{\rm fac}$ 增大 → 深度增大。

**与 4e 对照**：4e-01 有两独立 $a=-1$ 传播子（横截 → 深度 2）加两纯幂端点 → 4；3e-3g 把两圈揉进单 $S$（$a=-2$），单因子只给 1。**极点深度由 $S$ 的可因式分解结构决定，非单纯圈数**——这是本算法在 2e/4e 之上、由 3e-3g 揭示的核心推广。

---

## 7. 算法要点总结

1. **深度来自三类源**：纯幂端点 $a_i<0$（Step 2A）、单变量因子退化（Step 2B）、混合因子退化角点（Step 3 / R3，需 sector+jacobian 重算）。
2. **组合规则**：**横截相乘**（不同因子的退化轨迹涉及不相交变量集）、**同因子多轨迹取 max**（单因子的奇点预算不叠加，交点经标度分析仍深度 1）、**共变量不共角点取 max**。全程不做积分。
3. **jacobian 陷阱 / R2 降幂**：因子在 codim-$k$ 轨迹退化时，jacobian 消 $k$ 阶幂，$a\to a_{\rm eff}=a+k$。$a_{\rm eff}=-1$→深度 1；$>-1$→无极点（2e-02）；$<-1$→递归 sector。
4. **$a=-m$ 幂发散**（$m>1$，3e-3g 的 $S$ 给 $m=2$）：**不可**用 $|a|=m$ 作深度。单因子经 R2 降到 $a_{\rm eff}=-1$ → 深度 **1 而非 $m$**。这是 TGV（单 $S$，$a=-2$）与 4e（两独立 $a=-1$ 传播子，横截 → 深度 2）的本质差别。
5. **核心规律（3e-3g 揭示）**：`poledepth` = $S$ 因式分解出的**横截因子数** $N_{\rm fac}$（每因子深度 1）。退化（共边附着）使 $S$ 可分性增强 → $N_{\rm fac}$↑ → 深度↑。3e-3g 四例：$N_{\rm fac}=4,3,2,1$ → 深度 $4,3,2,1$。
6. **分子零点不产极点**：Pnum 的内部零点（如 $(1-x_1(1+z))$、$(\alpha_2-\alpha_3)$）不产生极点，只可能对消极点。
7. **预测 vs 验证**：本算法给出**预测深度**（代数/组合，无积分）；真正的 Laurent 系数仍需 hyperint（线性可约时）或数值 sector decomposition（不可约时）验证。3e-3g 四例深度 $4,3,2,1$ 均为**算法预测**，建议数值交叉验证。

---

## 8. 与管线的接口

```
montopar.json (J_raw.expanded)
        │
        ▼
   poledepth 探测器 (本算法)
        │
        ▼  写入 poledepth 字段
   dimreduce.json
        │
        ├─ poledepth = 0 → 期望有限, dimreduce 可直接消变量, 跳过 hyperint
        ├─ poledepth > 0 → dimreduce 需留意发散; 残余 J 交 hyperint 展开到 ε^poledepth
        └─ hyperint 实际展开阶 > poledepth → 标记 kicked_back, 打回 dimreduce
```

**实现计划**（待写）：把 §2 的 `poledepth(I)` 写成 SymPy 原型，输入 `J_raw.expanded` 字符串，自动完成：因式分解 → 净指数 → 轨迹枚举 → 横截分组 → 输出整数。先在 2e-01 / 4e-01 / 2e-02 三个已知例上验证，再推广到全部 4e orbit 与 3e-3g。

---

## 知识库来源

- [[papers/hep-ph/weinzierl-ch-10-sector-decomposition]] — §10.1, Eqs.10.3–10.19：通用积分形式、6 步算法、Taylor 极点提取、终止性
- [[papers/hep-ph/weinzierl-ch-08-multiple-polylogarithms]] — §8.9, Eqs.8.165–8.182：线性可约 ε-展开与 HyperInt 主公式（构造性验证路线）
- [[papers/hep-th/dotsenko-vergeles-1980]] — §2, Eq.2.13：contour 图散度指标 $\omega=1-N_{\rm ex}$（整体上界）
- [[papers/hep-ph/1310.1145]] — §3, Eqs.13–16：LiteRed 1.4 参数零扇区判据（剪枝）
- [[papers/hep-ph/analytic-tools-feynman-integrals]] — Smirnov 2012, Ch.9：expansion-by-regions 与 Newton 多面体（候选轨迹枚举的几何视角）
- [[concepts/1-epsilon-poles]] — 2-loop $1/\varepsilon^4$ 物理上限（sanity check）
- [[concepts/superficial-degree-of-divergence]] / [[concepts/weinberg-theorem]] — 横截相乘的子图定理依据
- [[concepts/primary-sector-decomposition]] — jacobian 消幂的代数来源（Eq.10.7）

> 注：知识库中无名为 "pole-depth detection" 的现成算法页；本算法是把 sector decomposition step 4、primary sector 的 jacobian 消幂、Weinberg 横截相乘、contour 散度指标**组合**成的分层探测器。各组件均有出处，组合方式是设计建议，需在具体 $J_{\rm raw}$ 上验证。
