# [3] 族基函数：完整母核与有限项坐标

日期：2026-09-17。完整函数由配对排列标记的 $R_\pi(\varepsilon,z)$ 表示；有限项采用既有 M 与 [2] 的 E。T3 只是已有求值结果的结构缩写，不是又新增一个独立方向。

## 1. 三标记母核

$$
R_\pi(\varepsilon,z)=-(1+z)^3z^{-3\varepsilon/2}
\int_{0<t_1<t_2<t_3<1\atop0<s_1<s_2<s_3<1}
\prod_{i=1}^3D_z(t_i,s_{\pi(i)})^{\varepsilon-1}\,d^3t\,d^3s.
$$

π=123、132、213 是配对，而非三个不同的积分变量次序。$R_{123}$ 自 crossing，$R_{132}(1/z)=R_{213}(z)$。完整函数的显式收敛三重级数见 [方法 §5](method.md)。

$$T_3(z)=R_{123}(0,z).$$
其完整 Hlog 系数与 E 的坐标保存在 [basis-functions.json](basis-functions.json)，来源为本次新鲜 HyperInt 与此前 [2] 母核输出。

## 2. M 的定义

沿用原记号，不改函数规范：
$$
M_{k,r}(z)=\int_0^1\frac{W(\beta)^kL(\beta,z)^r}
{\beta-z/(1+z)}\,d\beta,
$$
$$
W=\ln[\beta(1-\beta)],\qquad
L=\ln z+\ln\frac{1-\beta}{\beta}.
$$

本批使用 $M_{0,1},M_{1,1},M_{2,1},M_{0,3},M_{0,5}$。r≥1 时分母零点被 L 抵消，应按连续延拓读取，不额外加入主值或虚部处方。正实 z 下这些积分有限，权重为 k+r+1。尤其
$$M_{0,1}=\tau_{0,0}=-\frac{\ln^2z+\pi^2}{2}.$$

## 3. E 的结构定义

E 就是 [2] [表示评审](../same-channel-two/representation-review.md)中已检查的方向
$$
E(z)=\eta^{[2],\mathrm{id}}_{0,1;1}(z)
=\ln z\,[\varepsilon^1]\,
\frac{Q_{\mathrm{id};0,1}(\varepsilon,z)
-Q_{\mathrm{id};0,1}(\varepsilon,1/z)}{2}.
$$
这里
$$
Q_{\mathrm{id};0,1}
=(1+z)^2z^{-\varepsilon}
\int_{\Omega_2}D_z(t_1,s_1)^{\varepsilon-1}
D_z(t_2,s_2)^{\varepsilon-1}s_2^\varepsilon\,d^4s.
$$
所以也有直接的对数矩定义
$$
E(z)=\frac{\ln z}{2}(1+z)^2
\int_{\Omega_2}
\frac{\ln[s_2/(1-s_1)]}{D_z(t_1,s_1)D_z(t_2,s_2)}\,d^4s.
$$

E 是 crossing 偶、权重6的函数，并且 E(1)=0。它由同一个双核的端点插入与镜像差导出，不是为 [3] 的计算结果另起的名字。本批复用其已经算出的 Hlog，不把上面的定义积分当作仍未求值的最终答案。

## 4. 有限项的精确转换

$$
\boxed{
T_3=2E-\frac12M_{0,1}M_{2,1}
+\frac12M_{1,1}^2
+\frac{\pi^2}{9}M_{0,3}-\frac1{12}M_{0,5}.
}
$$

$$R_{132}^+(0,z)=\frac{M_{0,1}^3}{36}-T_3.$$

这两个关系经完整 Hlog/ζ 坐标反代，而不是数值拟合。三条有限项因此只需这一结构组合：
$$
c^{20}_0=-T_3/16,\qquad
c^{21}_0=c^{22}_0=T_3/32-M_{0,1}^3/1152.
$$

本族的新结果采用 M+E；此前 [2] 的 Σ/τ 结果和基底迁移状态未改写。固定 ε 阶的转换不等于完整 ε 三核函数可以化成单核函数。

## 5. Hlog 与数据约定

Hlog(z,[a1,...,ar]) 是字母表 {0,−1} 上、基点为0的迭代积分；单字母0为 ln z，全零字按 ln(z)^r/r!。遵循既有 shuffle 约定和正实 z 分支，周期保存为 ζ2、ζ3、ζ5 的幂积，完整系数为有理数。

[basis-functions.json](basis-functions.json)保存逐项字、周期幂与有理系数；结果修订同时保存紧凑表达和同一 Hlog 坐标。[检查记录](../../../../workflow/checks/same-channel-three/README.md)区分符号等价、数值一致性与尚未运行的形式整链。
