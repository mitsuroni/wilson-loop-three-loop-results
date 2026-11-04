# Π₂算符系数的显式 Laurent 展开

令 $e=\varepsilon$、$\zeta_n=\zeta(n)$。下表包含实际用于有限项卷积的全部算符有理函数与 Gamma 比值贡献：

$$
p_c(e)=\sum_{k=-2}^{2}p_{c,k}e^k+O(e^3),\qquad
\phi_c(e)=4e^2p_c(e)=\sum_{j=0}^{4}4p_{c,j-2}e^j+O(e^5).
$$

相邻图需要phi至四阶，对边至三阶。R=G10须先展开至五阶，再与表观三阶极点合并。以下统一提供两图所需的全部p系数。

## AA

| $k$ | $p_{c,k}$ |
|---|---|
| $-2$ | $- \frac{25}{12}$ |
| $-1$ | $\frac{17}{72}$ |
| $0$ | $- \frac{3919}{432} + \zeta\left(3\right)$ |
| $1$ | $\frac{3 \zeta_4}{2} - \frac{48535}{2592} + \frac{103 \zeta\left(3\right)}{6}$ |
| $2$ | $\frac{103 \zeta_4}{4} - \frac{1010863}{15552} - \frac{11 \zeta\left(3\right)}{36} + 7 \zeta\left(5\right)$ |

## Af

| $k$ | $p_{c,k}$ |
|---|---|
| $-2$ | $\frac{5}{3}$ |
| $-1$ | $- \frac{19}{18}$ |
| $0$ | $\frac{257}{108} + 8 \zeta\left(3\right)$ |
| $1$ | $12 \zeta_4 - \frac{74 \zeta\left(3\right)}{3} - \frac{295}{648}$ |
| $2$ | $- 37 \zeta_4 + \frac{30977}{3888} + \frac{149 \zeta\left(3\right)}{9} + 56 \zeta\left(5\right)$ |

## Ff

| $k$ | $p_{c,k}$ |
|---|---|
| $-2$ | $0$ |
| $-1$ | $2$ |
| $0$ | $\frac{31}{3} - 16 \zeta\left(3\right)$ |
| $1$ | $- 24 \zeta_4 + \frac{40 \zeta\left(3\right)}{3} + \frac{535}{18}$ |
| $2$ | $20 \zeta_4 - 112 \zeta\left(5\right) - \frac{244 \zeta\left(3\right)}{9} + \frac{9583}{108}$ |

[精确对象](operator-expansion-r001.json) · [相邻解表](families/adjacent/solutions.md) · [对边解表](families/opposite/solutions.md)。最终解表不留未展开的p、U或Gamma。
