<!-- double-opposite-current:begin -->
# 6e-32：[1,1] 族 I 的 Laurent 解析解

当前输入为 prepared 0.3。完整 epsilon 表达及 epsilon^-6 至有限项已完成计算和辅助检查；本批没有新的 Lean/MCP 整链收据。

## 1. 输入、域与目标

- 输入：[prepared 6e](../../prepared/6e.json)，orbit_id=6e-32。
- 包 SHA-256：9699532d1ee3a4027c2eb93c961a5b467a9bceda229c0204c13775cfa2451406。
- n=3，a1=1/4，sym=4，kappa=-1/16。
- z=z1>0，J 在 epsilon>0 收敛；I 先取 0<epsilon<1，再作固定正 z 的局部亚纯展开。
- 目标：全部极点至 epsilon^0，余项 O(epsilon)；更高正幂和物理复 z/i0 未纳入完成声明。
- 原外到内顺序：s11 → s21 → s32 → s42 → s31 → s41。

| 变量 | 下界 | 上界 | 端点类型 |
|---|---|---|---|
| s11 | 0 | 1 | O |
| s21 | 0 | 1 | O |
| s32 | 0 | 1 | B |
| s42 | 0 | 1 | A |
| s31 | 0 | s32 | O |
| s41 | 0 | s42 | O |

## 2. 完整结果

H0、H1 及 H1bar=H1(epsilon,1/z) 沿用 [上一批定义](../../families/single-opposite/basis.md)，G 是已补偿 Gamma 因子。

$$
J=\frac{z^{2 \varepsilon} \left(H_{0} H_{1} - H_{1}^{2}\right)}{\varepsilon^{2}}
$$

$$
\boxed{I=- \frac{G^{3} \left(\left(H_{0} H_{1} - H_{1}^{2}\right) e^{\frac{\ell \varepsilon}{2}} + \left(H_{0} \bar{H1} - \bar{H1}^{2}\right) e^{- \frac{\ell \varepsilon}{2}}\right)}{32 \varepsilon^{2} \cosh{\left(\frac{3 \ell \varepsilon}{2} \right)}}}
$$

[结果 r001](../../results/6e-32/r001.json) 包含完全展开 H0/H1 定义的标准 Meijer G 乘积，没有未求值参数积分。

## 3. 原始 J 与分组依据

原被积函数：

$$
z \left(s_{42} z \left(1 - s_{32}\right)\right)^{\varepsilon - 1} \left(- z - 1\right)^{2} \left(s_{11} s_{31} + z \left(1 - s_{11}\right) \left(1 - s_{31}\right)\right)^{\varepsilon - 1} \left(s_{21} s_{41} z + \left(1 - s_{21}\right) \left(1 - s_{41}\right)\right)^{\varepsilon - 1}
$$

两个对边核的变量互不相同；相邻端点的先后顺序由 bounds 保留。放宽到完整立方体后绝对可积，因此可按 [公共方法 §2](../../families/double-opposite/method.md) 使用 Tonelli/Fubini 分组。这里给出数学依据，不把它称为已有的形式换序收据。

| 有序端点 | 字 | 相邻端点正规化权重 |
|---|---|---|
| s11 | O | 1 |
| s21 | O | 1 |
| s31 < s32 | OB | b |
| s41 < s42 | OA | 1 - a |

每个相邻端点的 epsilon 分母单独提出；权重中的 a=s^epsilon、b=(1-s)^epsilon。

## 4. 反向核与 C/I

| 对边核 | 原变量 | 是否 tilde D | 反射后已知因子 |
|---|---|---|---|
| f1 | s11, s31 | True | H1 |
| f2 | s21, s41 | False | H0 - H1 |

tilde D=z D_(1/z)；将其归一化并作两端反射，权重也必须反射。等价的直接检查是 tilde D_z(1-u,1-v)=D_z(u,v)，Jacobian 绝对值为 1。两个核各给 z^(epsilon/2)，邻边再给 z^(r epsilon)。

$$
J=z^{(r+1)\varepsilon}\varepsilon^{-P}S(z),\qquad P=2n-4,\quad S=H_{0} H_{1} - H_{1}^{2}.
$$

C 已含且仅含一次 Euler 补偿。分别计算 J(z) 与 J(1/z)，代入

$$
I=\kappa G^n\frac{J(z)+z^{n\varepsilon}J(1/z)}{1+z^{n\varepsilon}},\qquad\kappa=\frac{\mathrm{sym}\,a_1}{2(-2)^n}.
$$

没有先假设 J 自镜像，也没有从另一个图移用颜色或重数。

## 5. 全部 Laurent 系数（M 主表示）

本批至有限项只需 M_(0,1)、M_(1,1)、M_(2,1) 与 zeta 常数的乘积；这些 M 保持原定义。

| epsilon 幂次 | M 解析系数 |
|---:|---|
| -6 | $0$ |
| -5 | $0$ |
| -4 | $0$ |
| -3 | $0$ |
| -2 | $0$ |
| -1 | $\frac{M_{0,1} M_{1,1}}{32}$ |
| 0 | $\frac{5 M_{0,1}^{3}}{1152} + \frac{7 \pi^{2} M_{0,1}^{2}}{1152} + \frac{3 M_{0,1} M_{2,1}}{128} + \frac{M_{1,1}^{2}}{32}$ |

$$I=\sum_{p=-6}^0c_p(z)\varepsilon^p+O(\varepsilon).$$

显式零项均已保存。复用 H0/H1 至二阶足以得到三圈 epsilon^-2 前因子后的有限项；最终极点阶由分子抵消决定。

## 6. tau 对照

$t_0=\tau_{0,0}=M_{0,1}$，$t_1=\tau_{0,1}=M_{1,1}/2$，$t_2=\tau_{0,2}=M_{2,1}/8-M_{0,1}^2/24$。

| epsilon 幂次 | tau 表达 |
|---:|---|
| -1 | $\frac{\tau_{0,0} \tau_{0,1}}{16}$ |
| 0 | $\frac{14 \tau_{0,0}^{3} + 7 \pi^{2} \tau_{0,0}^{2} + 216 \tau_{0,0} \tau_{0,2} + 144 \tau_{0,1}^{2}}{1152}$ |

M/tau 两式均精确反代到同一 Hlog/period 坐标；没有用数值拟合生成解析答案。完整 Hlog 记录保存在结果 JSON。

## 7. 检查、复用和剩余范围

[本图检查](../../../../workflow/checks/double-opposite/runs/20260916T110141611811Z/6e-32-checks.json)，复用及源哈希见结果 JSON。

- 原顺序 J 在 epsilon=1,2,3、z=2/3 的精确多项式比较残差为零。
- 每个对边核的反射权重另做精确积分比较，避免遗漏反向核的 z 因子。
- 两个非整数参数点的原域 Sobol 检查通过所列统计容差。
- 两个正 z 点，epsilon 逐次减半，闭式与有限项截断的差按 O(epsilon) 缩小。
- H0/H1 与 M 的已有求值按精确版本、规范和 q 范围复用，未重新运行 HyperInt。
- 形式化、物理复分支以及更高正幂未执行；不从本批有限阶纯 M 表达推断任意阶闭合。

<!-- double-opposite-current:end -->
