# 3e-Gamma1-02：334 两颜色 I 的完整反射解

## 当前物理积分的入口（2026-10-06）

从共同作用量出发，U′=igAU确定真实Wilson词；FT/Dyson规则在提取C前生效。当前[原生C/J](../../sources/montopar-r002.json)给出下式，J就是本文实际积分的原核，域、分支和IBP全保留：

$$I_{P,c}=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-3\epsilon)}{64}\frac{J_c(z)+z^{3\epsilon}J_c(1/z)}{1+z^{3\epsilon}}.$$

用下文已经求出的J代入此定义，得到[当前物理I](../../results/3e-Gamma1-02/r002.json)；组装系数为+1，Euler只补一次。以下带日期的求值过程与r001公式属于旧根历史，计算J的步骤可按同域积分线性复用，旧形式状态不继承。

2026-10-03，两颜色全部极点至有限项完成；[正式r001](../../results/3e-Gamma1-02/r001.json)。原输入是 [Montopar r001](../../sources/montopar-r001.json) 的02轨道，SHA256 `3b92e23125cb6b654ad63d5c19c17fead9cfeb3bff56739dd9f07c6eda87e65a`。这里复用已完成01的积分，先核对完整核及原有序域，不把几何相似当作相等依据。`formal=not_run`。

三端点是 $(3,1),(3,2),(4,1)$。记 $x'=s_{31},y'=s_{32},w'=s_{41}$，原外到内积分为
$$\int_0^1dy'\int_0^{y'}dx'\int_0^1dw'.$$
与01的344域作
$$
(x,y,w)_{01}=(1-w',1-y',1-x'),\quad
\alpha_1\leftrightarrow\alpha_3,\quad t_1\leftrightarrow t_3.
$$
Jacobian绝对值为1；01中的 $0\le y\le w\le1$ 化为 $0\le x'\le y'\le1$。在共同有限I0与置换对称α/t截断下核有界，允许相应迭代顺序的变换；继而采用同一个Feynman边界值和解析延拓。这不是未经条件核对的任意Fubini换序。

[独立反射记录](../../../../../workflow/checks/gamma1-01/independent/reflection.json)逐项验证q、胶子N0/L0/三个剑鱼、夸克N0/L0的全部分子幂次和权重。共同C、两颜色、Euler一次、C4重数4及外尺度规范均保留，故
$$J_{02,c}(z)=J_{01,c}(z)=z^{3\epsilon}j_c(\epsilon),\qquad I_{02,c}=I_{01,c}.$$
01与02仍是两个不同源轨道，各自按原C4权重计入裸贡献。

每个颜色采用
$$I_{02,c}=\frac{e^{-3\gamma_E\epsilon}\Gamma(1-3\epsilon)}{64}j_c.$$
全部系数与[01完整推导§2–3](../3e-Gamma1-01/derivation.md)相同，已逐项保存到本图r001。胶子最高六阶 $1/(4608\epsilon^6)$，夸克最高五阶 $-1/(864\epsilon^5)$。有限项分别为
$$
I_{02,g,0}=\frac{150911\zeta_2^3}{516096}+\frac{457\zeta_2^2}{27648}
+\frac{\zeta_2\zeta_3}{12}-\frac{18325\zeta_2}{93312}
+\frac{1177\zeta_3^2}{9216}-\frac{347\zeta_3}{5184}
+\frac{1907\zeta_5}{2880}-\frac{3201427}{1679616},
$$
$$
I_{02,q,0}=\frac{1365323}{1679616}+\frac{2921\zeta_2}{93312}
+\frac{305\zeta_3}{2592}-\frac{739\zeta_2^2}{34560}
-\frac{11\zeta_2\zeta_3}{64}-\frac{1331\zeta_5}{1440}.
$$
数值为 $0.0725695427120798869\ldots$ 和 $-0.350329707054862310\ldots$。完成范围是全部极点至ε⁰，余项 $O(\epsilon)$；没有新增全ε标准函数闭式或形式证明。01的解析积分及独立核验按确切输入/结果复用，没有把它们计为02独立重新积分的次数。
