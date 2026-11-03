# 本批沿用的 M、N、τ、E

本页复制既有[统一基底定义](../../../2ne/theory/unified-weight-six.md)以方便阅读，不引入新函数。各图r001还保存这些函数的完整Hlog坐标。

沿用
$$
\ell=\log z,\quad W=\log(\beta(1-\beta)),\quad
R=\log\frac{1-\beta}{\beta},\quad L=\ell+R,
$$
$$
M_{k,r}=\int_0^1\frac{W^kL^r}{\beta-z/(1+z)}\,d\beta,
\qquad
N_{k,r}=\int_0^1\frac{W^kRL^r}{\beta-z/(1+z)}\,d\beta.
$$
本表所用 r≥1，内部零点的可去奇性按原定义处理。两者满足
$$
N_{k,r}=M_{k,r+1}-\ell M_{k,r}.
$$

选择的七个 M 为 $M_{0,1},M_{1,1},M_{2,1},M_{3,1},M_{1,3},M_{4,1},M_{2,3}$；两个 N 为 $N_{1,2},N_{2,2}$。还需原有的 $\tau^+_{D,4}$ 和 E。

令 $D_z(t,s)=(1-t)(1-s)+zts$，
$$
H_D(\varepsilon,z)=-(1+z)z^{-\varepsilon/2}
\int_0^1\int_0^1 D_z(t,s)^{\varepsilon-1}t^\varepsilon s^\varepsilon\,dt\,ds,
$$
$$
\tau^+_{D,4}=[\varepsilon^4]\frac{H_D(\varepsilon,z)+H_D(\varepsilon,1/z)}2.
$$
D 是原 H_D 的标签，不是新函数族；提取系数本身已包含四阶 Taylor 展开的 1/4!。

对 $\Omega_2=\{0<t_1<t_2<1,\;0<s_1<s_2<1\}$，
$$
E=\frac{\ell(1+z)^2}{2}\int_{\Omega_2}
\frac{\log(s_2/(1-s_1))}{D_z(t_1,s_1)D_z(t_2,s_2)}\,d^2t\,d^2s
=\ell[\varepsilon^1]Q^-_{\mathrm{id};0,1}.
$$
两者都为 crossing 偶、权重 6。它们已有完整 Hlog 展开，未以未求值积分冒充新的计算完成。

本批的完整ε结构函数T、K及临时b矩H_j另见[推导](method.md)，与这些既有基底符号分开。
