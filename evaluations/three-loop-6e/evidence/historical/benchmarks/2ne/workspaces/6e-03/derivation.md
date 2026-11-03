<!-- empty-family-current:begin -->
# 6e-03：空集族 I 的完整解析结果

日期：2026-09-16。当前块对应 prepared 0.3；纸面推导和计算检查完成，本批结果未运行 Lean/MCP 整链验证。

## 1. 输入、目标和约定

- 原输入：[prepared 6e](../../prepared/6e.json)，orbit_id=6e-03。
- 包 SHA-256：9699532d1ee3a4027c2eb93c961a5b467a9bceda229c0204c13775cfa2451406。
- z=z1>0；J 在实 epsilon>0 收敛；I 先取 0<epsilon<1，再作固定正 z 下 epsilon=0 的局部亚纯延拓。
- 目标：完整 epsilon 的 I；全部极点及有限项，余项 O(epsilon)。物理复 z/i0 延拓未验收。
- 圈数 n=3，a1=1/4，sym=4；当前 C 已含一次 Euler 补偿。
- 原始外到内积分顺序：s33 → s43 → s32 → s42 → s31 → s41。

| 变量 | 下界 | 上界 | 端点权重 |
|---|---|---|---|
| s33 | 0 | 1 | B |
| s43 | 0 | 1 | A |
| s32 | 0 | s33 | B |
| s42 | 0 | s43 | A |
| s31 | 0 | s32 | B |
| s41 | 0 | s42 | A |

## 2. 最终结果与当前状态

公共记号 G、R2、R12 及 ell=log(z) 的定义见 [公共方法](../../families/empty/method.md)。

$$
J=\frac{z^{3 \varepsilon}}{36 \varepsilon^{6}}
$$

$$
\boxed{I=- \frac{G^{3}}{576 \varepsilon^{6}}}
$$

所有缩写都是显式 Gamma 函数或 log(z)，最终表达不含未求值积分或待定常数。
机器表达式已展开这些定义：[结果 r002](../../results/6e-03/r002.json)。

## 3. 从原始 J 出发的推导

原始被积函数（上下限和顺序见第 1 节）：

$$
z^{3} \left(s_{41} z \left(1 - s_{32}\right)\right)^{\varepsilon - 1} \left(s_{42} z \left(1 - s_{31}\right)\right)^{\varepsilon - 1} \left(s_{43} z \left(1 - s_{33}\right)\right)^{\varepsilon - 1}
$$

逐传播子拆分如下；只在正底数内部拆幂，端点按可积函数解释。

| 传播子 | q | B 端点 u | A 端点 v | 恒等式 |
|---|---|---|---|---|
| f1 | z | s31 | s42 | q(1-u)v |
| f2 | z | s32 | s41 | q(1-u)v |
| f3 | z | s33 | s43 | q(1-u)v |

故原 integrand 等于

$$
s_{41}^{\varepsilon - 1} s_{42}^{\varepsilon - 1} s_{43}^{\varepsilon - 1} z^{3 \varepsilon} \left(1 - s_{31}\right)^{\varepsilon - 1} \left(1 - s_{32}\right)^{\varepsilon - 1} \left(1 - s_{33}\right)^{\varepsilon - 1}
$$

按 [公共方法 §2](../../families/empty/method.md) 的绝对可积界，原迭代积分可以按以下独立有序链分组。这里使用明确的 Tonelli/Fubini 论证，并非假定根输入已附带换序收据。重组 Jacobian 为 1。

| 从小到大的原端点 | 字 | 已求值字积分 |
|---|---|---|
| s31 < s32 < s33 | BBB | $\frac{1}{6 \varepsilon^{3}}$ |
| s41 < s42 < s43 | AAA | $\frac{1}{6 \varepsilon^{3}}$ |

每个字的积分步骤见 [公共方法 §3](../../families/empty/method.md)，代入本图得到：

$$
K=F_{BBB}\,F_{AAA}=\frac{1}{36 \varepsilon^{6}},\qquad J=z^{3\varepsilon}K.
$$

这里保留了所有排序因子（例如 AA 的 1/2、AAA 的 1/6），没有因传播子配对相似而合并图号。

## 4. C/J 到 I

由 sym=4、a1=1/4、n=3，

$$
\kappa=\frac{\mathrm{sym}\,a_1}{2(-2)^n}=- \frac{1}{16}.
$$

先单独计算 $J(1/z)=z^{-r\varepsilon}K$，再代入当前唯一 C/I 约定：

$$
I=\kappa G^nK\frac{z^{r\varepsilon}+z^{(n-r)\varepsilon}}{1+z^{n\varepsilon}}
=\kappa G^nK\,\frac{\cosh[(r-n/2)\varepsilon\ell]}{\cosh(n\varepsilon\ell/2)}.
$$

本图 r=3，双曲比值为 $1$。C 的完整因子已与 prepared 声明逐项检查；G 只是补偿 Gamma 的缩写，没有再乘一次补偿。

## 5. Laurent 展开到有限项

按显式极点上界，正则部分至少展开到 6 阶；实际首个非零幂次为 -6。
下表含显式零项，ell=log(z)。

| epsilon 幂次 | 精确系数 |
|---:|---|
| -6 | $- \frac{1}{576}$ |
| -5 | $0$ |
| -4 | $- \frac{\pi^{2}}{2304}$ |
| -3 | $- \frac{\zeta\left(3\right)}{576}$ |
| -2 | $- \frac{19 \pi^{4}}{276480}$ |
| -1 | $- \frac{12 \zeta\left(5\right) + 5 \pi^{2} \zeta\left(3\right)}{11520}$ |
| 0 | $- \frac{60480 \zeta^{2}\left(3\right) + 631 \pi^{6}}{69672960}$ |

$$I=\sum_{k=-6}^0c_k(z)\varepsilon^k+O(\varepsilon).$$

余项针对固定正 z；Gamma 对数级数和双曲比值在零点邻域解析，正则部分的高阶余项乘回极点因子后为 O(epsilon)。

## 6. 检查与复现入口

[本图实际检查记录](../../../../workflow/checks/empty-family/runs/20260916T083424975899Z/6e-03-checks.json)：

- 原始 J 在 epsilon=1,2,3、z=2/3 的原顺序精确积分残差为零。
- 本图各有序字在三个非整数 epsilon 点通过 70 位工作精度数值积分检查。
- 原 bounds 的非整数 epsilon 随机化 Sobol 检查通过统计容差；误差不是严格区间包围。
- 完整 Gamma 表达的独立 Cauchy/Fourier 系数提取，与至有限项的符号系数一致。
- 三个 z 点的有限项截断余项除以 epsilon，随 epsilon 减半向下一系数收敛。
- 输入、代码、结果修订和检查输出哈希保存在 JSON 中。

复现命令、依赖版本及检查范围见 [计算说明](../../../../workflow/checks/empty-family/README.md)。

## 7. 遗留问题与历史

当前指定实域的解析计算和有限项已完成。本修订没有新 Lean/MCP 收据；物理复运动学分支和后续形式化仍单独开放。

<!-- empty-family-current:end -->
