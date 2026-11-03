<!-- single-opposite-current:begin -->
# 6e-26：[1] 族 I 的解析结果

当前结果绑定 prepared 0.3 原参数输入。解析计算和辅助检查完成，形式验证未执行。

## 1. 输入、条件与目标

- 来源：[prepared 6e](../../prepared/6e.json)，orbit_id=6e-26。
- 包 SHA-256：9699532d1ee3a4027c2eb93c961a5b467a9bceda229c0204c13775cfa2451406。
- n=3，a1=1/2，sym=4，r=0。
- z=z1>0；完整 I 先取 0<epsilon<1，固定正 z 后在 epsilon=0 作局部亚纯展开。
- 目标为完整 epsilon 的 I 及所有极点至有限项 O(epsilon)。物理复 z/i0 不在本批范围。
- 原外到内顺序：s32 → s23 → s41 → s31 → s22 → s21。

| 变量 | 下界 | 上界 | 类型 |
|---|---|---|---|
| s32 | 0 | 1 | A |
| s23 | 0 | 1 | O |
| s41 | 0 | 1 | O |
| s31 | 0 | s32 | A |
| s22 | 0 | s23 | B |
| s21 | 0 | s22 | B |

## 2. 完整 I

H0、H1、H2、HW、HD、HX 均已求值为标准 Meijer G；定义、T 映射和参数范围见 [公共表示](../../families/single-opposite/basis.md)。bar 表示 z→1/z；G 是已补偿 Gamma 因子，等同该文的 Z。

$$
\boxed{I=- \frac{G^{3} \left(\frac{H_{0} e^{\ell \varepsilon}}{8} + \frac{H_{0} e^{- \ell \varepsilon}}{8} - \frac{H_{1} e^{\ell \varepsilon}}{4} - \frac{\bar{H1} e^{- \ell \varepsilon}}{4} + \frac{H_{2} e^{\ell \varepsilon}}{8} + \frac{\bar{H2} e^{- \ell \varepsilon}}{8}\right)}{8 \varepsilon^{4} \cosh{\left(\frac{3 \ell \varepsilon}{2} \right)}}}
$$

[结果 r001](../../results/6e-26/r001.json) 同时给出完全展开 H 定义的 Meijer G 表达，不含未求值参数积分。
完整 epsilon 层的函数数量是本批共同求值表示，不能当作任意高圈的最小基底维数。

## 3. 原始 J 到有序权重

原被积函数为

$$
\left(s_{31} \left(1 - s_{21}\right)\right)^{\varepsilon - 1} \left(s_{32} \left(1 - s_{22}\right)\right)^{\varepsilon - 1} \left(- z - 1\right) \left(s_{23} s_{41} z + \left(1 - s_{23}\right) \left(1 - s_{41}\right)\right)^{\varepsilon - 1}
$$

传播子和分子逐项重建为唯一对边核 D 与邻边 q(1-u)v；绝对可积性及原迭代次序到分组积分的依据见公共表示 §2。没有套用未验证的超立方变换。

| 原端点（从小到大） | 字 | 正规化权重 |
|---|---|---|
| s21 < s22 < s23 | BBO | (1 - BT)**2/2 |
| s31 < s32 | AA | 1/2 |
| s41 | O | 1 |

令 AT=t^epsilon、BT=(1-t)^epsilon、AS=s^epsilon、BS=(1-s)^epsilon；BAT/BAS 表示相应的正规化 BA 下区间积分，R2=Q。每个相邻端点的 epsilon 已提到外面：

$$
J=-(1+z)z^{0\varepsilon}\varepsilon^{-4}\int_0^1dt\,ds\,D_z^{-1+\varepsilon}W,
$$

$$
W=\frac{BT^{2}}{4} - \frac{BT}{2} + \frac{1}{4}
$$

BA 的反射恒等式、前后区间及有序阶乘见公共表示 §3。r 和权重均从本图输入导出。

## 4. 当前 C、crossing 与求值

当前 C 已含一次 Euler 补偿。使用 $\kappa=\mathrm{sym}\,a_1/[2(-2)^n]$，本图 $\kappa=-1/8$，$\delta=r-(n-1)/2=-1$。
分别计算原权重和反射权重，再作 crossing；没有假设裸 J 是 self-mirror。本批 BA 若出现，其奇部分在此组合中精确消失。剩下的幂权重按同一 Meijer G 公式求值。

## 5. 从 epsilon^-6 到有限项

tau_(m,q) 保留旧定义；tau_W/D/X 和 eta 的定义及完整 Hlog 公式见 [系数函数目录](../../families/single-opposite/basis-functions.md)。外部系数无显式 z 或 log(z)，不引入 M 分母。

| epsilon 幂次 | 精确系数 |
|---:|---|
| -6 | $0$ |
| -5 | $0$ |
| -4 | $0$ |
| -3 | $0$ |
| -2 | $- \frac{\pi^{2} \tau_{0,0} + 24 \tau_{0,2}}{384}$ |
| -1 | $\frac{2 \eta_{1,2} - \eta_{2,2} - \tau_{0,3} + 2 \tau_{1,3} - \tau_{2,3}}{32}$ |
| 0 | $\frac{192 \eta_{1,3} - 96 \eta_{2,3} - 10 \pi^{2} \tau_{0,0}^{2} - 240 \tau_{0,0} \tau_{0,2} - 7 \pi^{4} \tau_{0,0} - 168 \pi^{2} \tau_{0,2} - 96 \tau_{0,4} + 192 \tau_{1,4} - 96 \tau_{2,4}}{3072}$ |

$$I=\sum_{p=-6}^{0}c_p(z)\varepsilon^p+O(\varepsilon).$$

对一/二圈统一列出的更低幂次明确为零。完整 Hlog/period 后端逐系数保存在 JSON；所有紧凑表达已反代到同一 Hlog 坐标检查零残差。

## 6. 与 M/N 的关系

| 幂次 | 纯 M 秩/增广秩 | M 加 log 秩/增广秩 |
|---:|---|---|
| -2 | 4/4 | 4/4 |
| -1 | 7/8 | 8/8 |
| 0 | 11/12 | 12/12 |

候选范围是固定自然权重、crossing 偶、系数为所列 zeta 常数的有限多项式；并未测试任意有理函数或放宽权重的空间。M 加 log 的成功表达可按 ell M=M−N 改写，已保存备选公式；失败的候选列表和见证在检查记录中。

## 7. 实际检查和遗留事项

[本图检查记录](../../../../workflow/checks/single-opposite/runs/20260916T093258061172Z/6e-26-checks.json)：

- 三个整数 epsilon 点的原顺序精确 J 与有序权重降维结果一致。
- 两个非整数 epsilon/z 点的原域 Sobol 检查与全 I 的 Meijer G 值一致于所列统计容差。
- 两个正 z 点逐次减半 epsilon，有限项截断误差按 O(epsilon) 缩小。
- 六种母核已重新由 HyperInt 求值；另用 Euler 数值积分和独立 Meijer G Cauchy 系数提取、Hlog 微分递推级数交叉检查。
- 复运动学和本批 Lean/MCP 整链未执行。统计数值检查不是严格区间包围。

<!-- single-opposite-current:end -->
