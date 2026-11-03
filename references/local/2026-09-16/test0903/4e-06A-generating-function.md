---
title: "4e-06A｜因子化积分与生成函数"
category: "纯交换个例"
status: "个例参考"
tags:
  - wilson-loop
  - pure-exchange
  - generating-function
aliases:
  - "4e-06A 生成函数推导"
catalog_updated: 2026-09-15
---

# 4e-06A｜因子化积分与生成函数

> 从 montopar 原始积分出发，通过精确的代数化简与 $\mathcal{M}_{k,1}$ 对称性证明，将对称化积分严格写为三个独立生成子的乘积。
>
> 前置：`docs/2e-02-generating-function-derivation.md`
>
> **约定**: 全局丢弃 $\gamma_E$。

## 1. 原始积分与因子化

### 1.1 montopar 输出

$$J_{raw} = (1+z)^2 \cdot \underbrace{[x_1x_2 + z(1-x_1)(1-x_2)]^{-1+\varepsilon}}_{f_1(x_1,x_2)} \cdot \underbrace{[(1-x_3)(1-x_4) + x_3x_4z]^{-1+\varepsilon}}_{f_2(x_3,x_4)}$$

$f_1 \equiv$ 2e-02 旋转，$f_2 \equiv$ 2e-02 原始。变量完全分离 → $J_{raw} = (1+z)^2 \cdot K(z)^2$，其中 $K(z)$ 是裸 2e-02 积分核。

### 1.2 降维与 $\alpha$ 因子对消

引入 $\alpha = \frac{1}{1+z}$, $\bar\alpha = \frac{z}{1+z}$。从 2e-02 的降维结果：

$$K(z) = \frac{\alpha}{\varepsilon} \int_0^1 \frac{(\bar\alpha/\alpha)^\varepsilon \tau^\varepsilon - (1-\tau)^\varepsilon}{\tau-\alpha} d\tau$$

作 $\tau = 1-\beta$，利用双曲恒等式：

$$(\bar\alpha\bar\beta)^\varepsilon - (\alpha\beta)^\varepsilon = \varepsilon(\alpha\bar\alpha)^{\frac{\varepsilon}{2}} \cdot e^{\frac{\varepsilon}{2}\ln(\beta\bar\beta)} \cdot \frac{\sinh(\frac{\varepsilon}{2}L)}{\varepsilon/2}, \quad L = \ln\frac{\bar\alpha\bar\beta}{\alpha\beta}$$

定义 $\mathcal{M}_{k,1}$ **生成子**（与 2e-02 完全相同的结构）：

$$\mathcal{M}(\alpha) \equiv \int_0^1 \frac{d\beta}{\beta-\bar\alpha} \cdot e^{\frac{\varepsilon}{2}\ln(\beta\bar\beta)} \cdot \frac{\sinh(\frac{\varepsilon}{2}L)}{\varepsilon/2}$$

则 $K(z) = -\alpha^{1-\varepsilon} (\alpha\bar\alpha)^{\frac{\varepsilon}{2}} \mathcal{M}(\alpha)$。代回 $J_{raw}$：

$$J_{raw} = \alpha^{-2} \cdot \alpha^{2-2\varepsilon} (\alpha\bar\alpha)^\varepsilon \mathcal{M}(\alpha)^2 = \left(\frac{\bar\alpha}{\alpha}\right)^\varepsilon \mathcal{M}(\alpha)^2 = z^\varepsilon \mathcal{M}(\alpha)^2$$

**$\alpha$ 因子完美对消。**

## 2. 交叉对称化

规范化后的分支积分：$J_i(z) = -\frac{1}{8}\Gamma(1-\varepsilon)^2 z^\varepsilon \mathcal{M}(\alpha)^2$

4e-06A sym=1：$I = \frac{1}{2} \cdot \frac{J_i(z) + z^{2\varepsilon} J_i(1/z)}{1+z^{2\varepsilon}}$

$z \to 1/z$ 即 $\alpha \leftrightarrow \bar\alpha$：$J_i(1/z) = -\frac{1}{8}\Gamma(1-\varepsilon)^2 z^{-\varepsilon} \mathcal{M}(\bar\alpha)^2$

$$z^{2\varepsilon}J_i(1/z) = -\frac{1}{8}\Gamma(1-\varepsilon)^2 z^{\varepsilon} \mathcal{M}(\bar\alpha)^2$$

求和并利用 $1+z^{2\varepsilon} = 2z^\varepsilon\cosh(\varepsilon\ln z)$：

$$I = -\frac{1}{32}\Gamma(1-\varepsilon)^2 \operatorname{sech}(\varepsilon L_d) \left[ \mathcal{M}(\alpha)^2 + \mathcal{M}(\bar\alpha)^2 \right], \quad L_d = \ln z$$

## 3. $\mathcal{M}_{k,1}(z) = \mathcal{M}_{k,1}(1/z)$ 的积分证明

### 3.1 定理

对于所有 $k \ge 0$，$\mathcal{M}_{k,1}(z) = \int_0^1 \frac{d\beta}{\beta-\bar\alpha} \ln\frac{\bar\alpha\bar\beta}{\alpha\beta} \ln^k(\beta\bar\beta)$ 在 $z \leftrightarrow 1/z$ 下不变。

### 3.2 证明

从 $\mathcal{M}_{k,1}(1/z)$ 的定义出发：

$$\mathcal{M}_{k,1}(1/z) = \int_0^1 \frac{d\beta}{\beta-\alpha} \ln\frac{\alpha\bar\beta}{\bar\alpha\beta} \ln^k(\beta\bar\beta)$$

作变量代换 $\beta = 1-u$：

| 因子 | 变换前 | 变换后 |
|------|--------|--------|
| 测度 | $d\beta$ | $-du$ |
| 极点 | $\beta-\alpha$ | $(1-u)-\alpha = \bar\alpha-u = -(u-\bar\alpha)$ |
| 对数 | $\ln\frac{\alpha\bar\beta}{\bar\alpha\beta}$ | $\ln\frac{\alpha u}{\bar\alpha(1-u)} = \ln\frac{\alpha u}{\bar\alpha\bar u}$ |
| $\ln^k$ | $\ln^k(\beta\bar\beta)$ | $\ln^k(u\bar u)$（不变） |

积分上下限互换引入负号 $\int_1^0 = -\int_0^1$：

$$\mathcal{M}_{k,1}(1/z) = \int_1^0 \frac{-du}{-(u-\bar\alpha)} \ln\frac{\alpha u}{\bar\alpha\bar u} \ln^k(u\bar u) = -\int_0^1 \frac{du}{u-\bar\alpha} \ln\frac{\alpha u}{\bar\alpha\bar u} \ln^k(u\bar u)$$

提出对数内的负号 $\ln\frac{\alpha u}{\bar\alpha\bar u} = -\ln\frac{\bar\alpha\bar u}{\alpha u}$，前面的负号消去：

$$\mathcal{M}_{k,1}(1/z) = \int_0^1 \frac{du}{u-\bar\alpha} \ln\frac{\bar\alpha\bar u}{\alpha u} \ln^k(u\bar u)$$

将哑变量 $u$ 重新标记为 $\beta$，即得 $\mathcal{M}_{k,1}(z)$ 的定义。$\square$

**三个符号变化恰好抵消**：$d\beta \to -du$（负）、$\beta-\alpha \to -(u-\bar\alpha)$（负）、$\int_1^0 = -\int_0^1$（负）→ 三个负号 = 一个负号，恰好翻转对数参数 $\ln\frac{\alpha u}{\bar\alpha\bar u} \to \ln\frac{\bar\alpha\bar u}{\alpha u}$。

### 3.3 推论

所有 $\mathcal{M}_{k,2j+1}$ 在 $z \leftrightarrow 1/z$ 下均不变。因此：

$$\mathcal{M}(\bar\alpha) = \mathcal{M}(\alpha)$$

$$\Sigma\mathcal{F} \equiv \mathcal{F}(z) + \mathcal{F}(1/z) = 2\mathcal{F}(z)$$

**$\mathcal{M}_{i,1}(1/z)$ 不是独立基底，全部映射到 $\mathcal{M}_{i,1}(z)$。**

## 4. 终极生成函数

利用 $\mathcal{M}(\bar\alpha) = \mathcal{M}(\alpha)$，交叉对称和简化为 $2\mathcal{M}(\alpha)^2$：

$$\boxed{I(\varepsilon) = -\frac{1}{16} \cdot \underbrace{Z(\varepsilon)^2}_{\text{$\zeta$ 生成子}} \cdot \underbrace{\operatorname{sech}(\varepsilon L_d)}_{\text{$\ln z$ 生成子}} \cdot \underbrace{\mathcal{M}(\varepsilon; \alpha)^2}_{\text{$\mathcal{M}_{k,1}$ 乘积生成子}}}$$

其中 $Z(\varepsilon) = \exp\!\left(\sum_{k=2}^\infty \frac{\zeta(k)}{k}\varepsilon^k\right)$。

### 与 2e-02 的对比

| 特征 | 2e-02 | 4e-06A |
|------|-------|--------|
| 结构 | $-\frac{1}{2} Z(\varepsilon) \operatorname{sech}(\frac{\varepsilon}{2}L_d) \mathcal{M}(\alpha)$ | $-\frac{1}{16} Z(\varepsilon)^2 \operatorname{sech}(\varepsilon L_d) \mathcal{M}(\alpha)^2$ |
| $\zeta$ 生成子 | $Z(\varepsilon)$ | $Z(\varepsilon)^2$ |
| $\ln z$ 生成子 | $\operatorname{sech}(\frac{\varepsilon}{2}L_d)$ | $\operatorname{sech}(\varepsilon L_d)$ |
| $\mathcal{M}_{k,1}$ 生成子 | $\mathcal{M}(\alpha)$ | $\mathcal{M}(\alpha)^2$ |
| $\mathcal{M}_{k,1}$ 结构 | 单积分 $\mathcal{M}_{k,1}$ | 乘积 $\mathcal{M}_{j,1} \otimes \mathcal{M}_{k,1}$ |

## 5. 各组分级数展开

### 5.1 $\zeta$ 与 $\ln z$ 生成子

$$Z(\varepsilon)^2 = \exp\!\left(\sum_{k=2}^\infty \frac{2\zeta(k)}{k}\varepsilon^k\right) = 1 + \zeta_2\varepsilon^2 + \frac{2\zeta_3}{3}\varepsilon^3 + \left(\frac{\zeta_4}{2} + \frac{\zeta_2^2}{2}\right)\varepsilon^4 + \mathcal{O}(\varepsilon^5)$$

$$\operatorname{sech}(\varepsilon L_d) = \sum_{m=0}^\infty \frac{E_{2m}}{(2m)!} (\varepsilon L_d)^{2m} = 1 - \frac{L_d^2}{2}\varepsilon^2 + \frac{5L_d^4}{24}\varepsilon^4 + \mathcal{O}(\varepsilon^6)$$

其中 $E_{2m}$ 是 Euler 数。$\operatorname{sech}$ 为偶函数，无 $\varepsilon^1, \varepsilon^3$ 阶的 $\ln z$ 项。

### 5.2 前置组合生成子 $P(\varepsilon) = Z(\varepsilon)^2 \operatorname{sech}(\varepsilon L_d)$

$$P(\varepsilon) = \sum_{N=0}^\infty P^{(N)} \varepsilon^N$$

| $N$ | $P^{(N)}$ |
|-----|----------|
| $0$ | $1$ |
| $1$ | $0$ |
| $2$ | $\zeta_2 - \frac{1}{2}\ln^2 z$ |
| $3$ | $\frac{2}{3}\zeta_3$ |
| $4$ | $\frac{\zeta_4}{2} + \frac{\zeta_2^2}{2} - \frac{\zeta_2}{2}\ln^2 z + \frac{5}{24}\ln^4 z$ |

### 5.3 $\mathcal{M}(\alpha)$ 的展开

$$\mathcal{M}(\alpha) = \sum_{k,j\ge 0} \frac{\varepsilon^{k+2j}}{k!(2j+1)!2^{k+2j}} \cdot \mathcal{M}_{k,2j+1}$$

$$\mathcal{M}(\alpha) = \mathcal{M}_{0,1} + \frac{\varepsilon}{2}\mathcal{M}_{1,1} + \varepsilon^2\left(\frac{\mathcal{M}_{2,1}}{8} + \frac{\mathcal{M}_{0,3}}{24}\right) + \varepsilon^3\left(\frac{\mathcal{M}_{3,1}}{48} + \frac{\mathcal{M}_{1,3}}{96}\right) + \mathcal{O}(\varepsilon^4)$$

其中 $\mathcal{M}_{k,1}$（$k \ge 0$）为标准 M 基底，$\mathcal{M}_{0,3}$ 为 IBP 伴随积分。

### 5.4 $\mathcal{M}(\alpha)^2$ 的展开

$$\mathcal{M}(\alpha)^2 = \sum_{N=0}^\infty S^{(N)} \varepsilon^N$$

| $N$ | $S^{(N)}$ |
|-----|----------|
| $0$ | $\mathcal{M}_{0,1}^2$ |
| $1$ | $\mathcal{M}_{0,1} \mathcal{M}_{1,1}$ |
| $2$ | $\frac{1}{4}\mathcal{M}_{0,1} \mathcal{M}_{2,1} + \frac{1}{12}\mathcal{M}_{0,1} \mathcal{M}_{0,3} + \frac{1}{4}\mathcal{M}_{1,1}^2$ |
| $3$ | $\frac{1}{24}\mathcal{M}_{0,1} \mathcal{M}_{3,1} + \frac{1}{48}\mathcal{M}_{0,1} \mathcal{M}_{1,3} + \frac{1}{8}\mathcal{M}_{1,1} \mathcal{M}_{2,1} + \frac{1}{24}\mathcal{M}_{1,1} \mathcal{M}_{0,3}$ |

## 6. $\ln^2 z$ 在 M 基底内的投影

$\mathcal{M}_{0,1}$ 的解析形式：$\mathcal{M}_{0,1} = -\frac{1}{2}\ln^2 z - 3\zeta_2$

$$\boxed{\ln^2 z = -2\mathcal{M}_{0,1} - 6\zeta_2}$$

因此 $P^{(2)} = \zeta_2 - \frac{1}{2}(-2\mathcal{M}_{0,1} - 6\zeta_2) = \mathcal{M}_{0,1} + 4\zeta_2$。

所有 $\ln^{2m} z$ 均可通过 $\mathcal{M}_{0,1}$ 的幂次展开投影到 M 基底。

## 7. 逐阶 Laurent 系数

$$I^{(N)} = -\frac{1}{16} \sum_{k=0}^N P^{(k)} S^{(N-k)}$$

### 各阶结果

| 阶 | $I_{4e-06A}^{(k)}$ | 备注 |
|----|-------------------|------|
| $\varepsilon^0$ | $\displaystyle -\frac{1}{16} \mathcal{M}_{0,1}^2$ | 与旧框架一致 ✓ |
| $\varepsilon^1$ | $\displaystyle -\frac{1}{16} \mathcal{M}_{0,1} \mathcal{M}_{1,1}$ | 纯 M 基底，无 $\ln z$ 污染 |
| $\varepsilon^2$ | $\displaystyle -\frac{1}{64} \mathcal{M}_{0,1} \mathcal{M}_{2,1} - \frac{1}{192} \mathcal{M}_{0,1} \mathcal{M}_{0,3} - \frac{1}{64} \mathcal{M}_{1,1}^2 - \frac{1}{16} \mathcal{M}_{0,1}^3 - \frac{1}{4}\zeta_2 \mathcal{M}_{0,1}^2$ | $\ln^2 z$ 已投影为 $\mathcal{M}_{0,1}^3 + 4\zeta_2 \mathcal{M}_{0,1}^2$ |
| $\varepsilon^3$ | $\displaystyle -\frac{1}{384} \mathcal{M}_{0,1} \mathcal{M}_{3,1} - \frac{1}{768} \mathcal{M}_{0,1} \mathcal{M}_{1,3} - \frac{1}{128} \mathcal{M}_{1,1} \mathcal{M}_{2,1} - \frac{1}{384} \mathcal{M}_{1,1} \mathcal{M}_{0,3} - \frac{1}{16} \mathcal{M}_{0,1}^2 \mathcal{M}_{1,1} - \frac{1}{4}\zeta_2 \mathcal{M}_{0,1} \mathcal{M}_{1,1} - \frac{1}{24}\zeta_3 \mathcal{M}_{0,1}^2$ | 全部在 M-IBP 闭包内 |

## 8. 验证

| 验证项 | 方法 | 结果 |
|--------|------|------|
| $I^{(0)}$ 与旧框架一致 | `fibrationBasis` 简化 | ✓ |
| $I^{(1)}$ 旧框架扩展 | 交叉对称化到 $\varepsilon^1$ | ✓ |
| $\mathcal{M}_{k,1}(1/z) = \mathcal{M}_{k,1}(z)$ | 积分代换证明 | ✓ |
| 边界条件 $I^{(k)}(-1) = 0$ | $\mathcal{M}_{k,1}(-1)=0$ | ✓ |