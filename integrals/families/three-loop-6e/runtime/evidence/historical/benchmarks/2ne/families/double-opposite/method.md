# [1,1] 族：规划、原参数降维与 I 的 Laurent 系数

日期：2026-09-16。范围为 prepared 0.3 的 4e-06、6e-32、6e-33、6e-34。目标延续前两批：计算完整 epsilon 的 I，并列出 epsilon^-6 至 epsilon^0，余项 O(epsilon)。更高正幂不列为本批已计算阶数。

## 1. 执行路线

1. 从当前原始多项式和 bounds 识别两个不同通道的对边核及相邻端点，保存原外到内顺序。
2. 先完成相邻 A/B 端点的有序积分，说明何时允许分组；不根据 [1,1] 标签直接假定结果因子化。
3. 处理反向对边的尺度及反射，把两核写成已知 H0/H1 的乘积，连接当前唯一 C/I。
4. 复用上一批已经求值并检查的 H0/H1 系数至二阶，卷积得到本批系数；检查 M 与 tau 两种表达的精确等价。
5. 用原顺序整数 epsilon 精确积分、非整数原域数值积分、闭式与截断展开的余项缩放作新连接检查。复用文件按实际哈希冻结，不重算相同 HyperInt 积分。

数学规划以实际范围为限；若发现多个 O 同链或剩余权重跨核耦合，应保留缺口，不能套用下面的乘积式。

## 2. 实际输入的两种核

记

$$
D_z(u,v)=(1-u)(1-v)+zuv,\qquad
\widetilde D_z(u,v)=z(1-u)(1-v)+uv=zD_{1/z}(u,v).
$$

每条输入恰有一个 D、一个 tilde D，两个核的四个端点不同；三圈再有一个相邻传播子 z(1-b)a。分子分别为 (1+z)^2 或 (1+z)^2 z。逐多项式重建残差、原端点对应和全部 bounds 均保存在单图推导及检查记录中。

先取 z>0、epsilon>0。相邻端点权重可积；各对边核只在两个角点线性为零，二维 D^(epsilon-1) 可积。放宽排序域到完整立方体后绝对值积分有限，因为所有不同端点仍分属独立因子。Tonelli/Fubini 因而允许改变分组；这是一条本批明确使用的数学桥，不是 prepared 根自带的形式收据。

## 3. 相邻端点的四种位置

A(x)=x^(epsilon-1)，B(x)=(1-x)^(epsilon-1)。对边端点记为 s：

| 排序 | 相邻端点积分 |
|---|---|
| AO | $s^\varepsilon/\varepsilon$ |
| OA | $(1-s^\varepsilon)/\varepsilon$ |
| BO | $[1-(1-s)^\varepsilon]/\varepsilon$ |
| OB | $(1-s)^\varepsilon/\varepsilon$ |

本批最多一条相邻边，每个标记链长度不超过二。相邻两个端点积分共同给显式 epsilon^-2，但分子可能抵消这个极点。每个上下限直接来自原 bounds，没有乘入历史超立方 Jacobian。

## 4. 反向对边、H0/H1 与因子化

沿用 [上一批固定定义](../single-opposite/basis.md)：

$$
H_0=-(1+z)z^{-\varepsilon/2}\int_0^1du\,dv\,D_z(u,v)^{-1+\varepsilon},
$$

$$
H_1=-(1+z)z^{-\varepsilon/2}\int_0^1du\,dv\,D_z(u,v)^{-1+\varepsilon}v^\varepsilon.
$$

它们已经有标准 Meijer G 求值和实际 Hlog 系数。上横线表示 z→1/z，H0bar=H0。给一个端点加 (1-v)^epsilon 则得到 H1bar。

对 tilde D，必须同时跟踪

$$
-(1+z)\widetilde D_z^{-1+\varepsilon}
=z^\varepsilon[-(1+1/z)D_{1/z}^{-1+\varepsilon}].
$$

因此完成积分后总因子仍为 z^(epsilon/2)，但权重需要在镜像下反射。反射 u→1-u、v→1-v 的 Jacobian 绝对值为 1，将 v^epsilon 与 (1-v)^epsilon 互换。这一步不能只替换核而漏掉尺度或权重。

设 A=H0、B=H1，则从实际排序得到：

| 积分 | 反向对边的因子 | 正向对边的因子 | J |
|---|---|---|---|
| 4e-06 | A | A | $z^\varepsilon A^2$ |
| 6e-32 | B | A−B | $z^{2\varepsilon}B(A-B)/\varepsilon^2$ |
| 6e-33 | A−B | A−B | $z^{2\varepsilon}(A-B)^2/\varepsilon^2$ |
| 6e-34 | A−B | B | $z^{2\varepsilon}B(A-B)/\varepsilon^2$ |

所以 6e-32 与 6e-34 的 J 相同。这一等式由不同输入的端点积分导出，不是按图号预设；还要核对 C/sym 才能比较 I。

4e-06 的因子化与 [旧 4e-06A 参考稿](../../../../references/local/2026-09-16/test0903/4e-06A-generating-function.md) 一致。本批只复用其结构和有限项，不沿用“全局丢弃 Euler 常数”的写法，也不直接采纳该草稿未经本轮复核的更高阶表格。

## 5. 当前 C/I 的完整 epsilon 表达

令

$$
G(\varepsilon)=e^{-\gamma_E\varepsilon}\Gamma(1-\varepsilon),\quad
\ell=\ln z,\quad
\kappa=\frac{\mathrm{sym}\,a_1}{2(-2)^n}.
$$

C 已含且仅含一次 Euler 补偿。当前对称化定义给出

$$
I=\kappa G^n\frac{J(z)+z^{n\varepsilon}J(1/z)}{1+z^{n\varepsilon}}.
$$

实际常数为

| 积分 | n | a1 | sym | kappa |
|---|---:|---:|---:|---:|
| 4e-06 | 2 | −1/2 | 1 | −1/16 |
| 6e-32 | 3 | 1/4 | 4 | −1/16 |
| 6e-33 | 3 | 1/2 | 4 | −1/8 |
| 6e-34 | 3 | 1/4 | 4 | −1/16 |

定义 F=B(A−B)、K=(A−B)^2，仅作两个显式乘积的缩写，则

$$
I_{4e-06}=-\frac{G^2 A^2}{16\cosh(\varepsilon\ell)},
$$

$$
I_{6e-32}=I_{6e-34}
=-\frac{G^3}{32\varepsilon^2\cosh(3\varepsilon\ell/2)}
\left[e^{\varepsilon\ell/2}F(z)+e^{-\varepsilon\ell/2}F(1/z)\right],
$$

$$
I_{6e-33}
=-\frac{G^3}{16\varepsilon^2\cosh(3\varepsilon\ell/2)}
\left[e^{\varepsilon\ell/2}K(z)+e^{-\varepsilon\ell/2}K(1/z)\right].
$$

这些是全 epsilon 表达，H0/H1 的标准函数定义可完全代入。I 先取实 0<epsilon<1；固定正 z 后使用 epsilon=0 的局部亚纯延拓。物理复 z/i0 与 Lean 整链不在本批完成声明内。

## 6. 展开与表示选择

只需三个已有函数

$$t_0=\tau_{0,0},\qquad t_1=\tau_{0,1},\qquad t_2=\tau_{0,2}.$$

上一批 Hlog 系数中的低阶恒等式为

$$
H_0=t_0+\varepsilon t_1+\varepsilon^2t_2+O(\varepsilon^3),
$$

$$
H_1^+=t_0+2\varepsilon t_1+
\varepsilon^2\left(4t_2-\frac{\ell^2t_0}{24}\right)+O(\varepsilon^3),
\qquad
H_1^-=-\frac{\varepsilon\ell t_0}{6}+O(\varepsilon^2).
$$

这些关系会与实际 Hlog 输出精确反代，而不是只靠数值点识别。结合 G=1+O(epsilon^2) 和 crossing，可得到：

$$I_{4e-06}=-\frac{t_0^2}{16}+O(\varepsilon),$$

$$
I_{6e-32}=I_{6e-34}
=\frac{t_0t_1}{16\varepsilon}
+\frac{3t_0t_2}{16}+\frac{t_1^2}{8}
+\frac{7t_0^3}{576}+\frac{7\zeta(2)t_0^2}{192}
+O(\varepsilon),
$$

$$
I_{6e-33}
=-\frac{t_1^2}{8}+\frac{t_0^3}{144}
+\frac{\zeta(2)t_0^2}{48}+O(\varepsilon).
$$

特别地，6e-33 的显式 epsilon^-2 完全抵消；6e-32/34 只留简单极点。以上候选须通过本批原积分、卷积及数值检查后才写入完成态结果。

转换

$$
t_0=M_{0,1},\quad t_1=M_{1,1}/2,\quad
t_2=M_{2,1}/8-M_{0,1}^2/24
$$

后，本批到有限项落在 M_(0,1)、M_(1,1)、M_(2,1) 与 zeta(2) 的有限多项式空间。主表采用这一简洁 M 表达，同时保存 tau 和 Hlog；没有理由因 [1] 的个别三圈图需要扩展，就给这四条新增基函数。也不由本批低阶结果推断任意更高阶都纯 M 闭合。

## 7. 复用与检查

- 复用上一批 H0/H1 的实际 HyperInt 系数及独立 Euler/Cauchy 检查；核对具体文件哈希、规范和 q 范围。
- 对每图从原 expanded J 按原顺序积分，在 epsilon=1,2,3 时与降维式作精确多项式比较；只检查 J，不在 Gamma 极点上评估 I。
- 在非整数 epsilon、多个正 z 点对原 bounds 和原传播子作随机化 Sobol 检查，比较完整 I；统计误差不称为严格区间。
- 用已求值 Meijer G 闭式与 Hlog 展开比较，检查截断误差随 epsilon 减半按 O(epsilon) 缩小。
- M、tau、Hlog 三种系数表达互相精确反代；预期的零极点显式存零。
- 计算源、复用源、结果修订和检查输出均存哈希；不把本轮结果赋予旧 Lean 收据。

复现入口：[计算脚本](../../../../workflow/checks/double-opposite/README.md)。
