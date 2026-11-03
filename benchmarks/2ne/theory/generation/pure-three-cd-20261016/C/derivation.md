# C：逆序纯对边积分的直接生成结构

2026-10-16。$C=R_{321}$ 的完整颜色为 CF³，在当前 a1 通道为零；本页保留它的未投影几何积分。目标是由自身原排序域求生成结构，不通过已知 A/B 反定义 C。核、正实 z 分支及归一化沿用原纯对边母族。

## 1. 原排序直接给出两个互补矩形

\[
C=\int_{0<t_1<t_2<t_3<1\atop0<s_1<s_2<s_3<1}
k_\epsilon(t_1,s_3)k_\epsilon(t_2,s_2)k_\epsilon(t_3,s_1)\,d^3t\,d^3s,
\]
\[
k_\epsilon(t,s)=-(1+z)z^{-\epsilon/2}
[(1-t)(1-s)+zts]^{\epsilon-1}.
\]

固定中间传播子的 t=t2、s=s2。左边 t1 与上边 s3 独立，形成左上矩形；右边 t3 与下边 s1 独立，形成右下矩形。令

\[
F(t,s)=\int_0^tdx\int_0^sdy\,k_\epsilon(x,y),\qquad f(t)=F(t,1),
\]

则四个外端点已经积去，得到

\[
\boxed{C=\int_0^1dt\int_0^1ds\,k_\epsilon(t,s)
[f(t)-F(t,s)][f(s)-F(t,s)].}
\tag{C1}
\]

这与同序图 A 的左下/右上矩形不同；两种排序没有被互相替换。

## 2. 矩形原函数显式求值

置 $p=(1+z)t-1,q=(1+z)s-1,\rho=z/(1+z)^2$，并记

\[
J_\epsilon(x)=\int_0^x\frac{(1-u)^\epsilon-1}{u}du
=-\epsilon x\,{}_3F_2(1,1,1-\epsilon;2,2;x).
\]

已有单核原函数给出

\[
F=-\frac{\rho^{\epsilon/2}}\epsilon
[J_\epsilon(-pq/z)-J_\epsilon(p/z)-J_\epsilon(q/z)+J_\epsilon(-1/z)].
\]

因此两个互补矩形恰好为

\[
f(t)-F=-\rho^{\epsilon/2}\mathcal U_\epsilon/\epsilon,\qquad
f(s)-F=-\rho^{\epsilon/2}\mathcal V_\epsilon/\epsilon,
\]
\[
\begin{aligned}
\mathcal U_\epsilon&=J_\epsilon(-p)-J_\epsilon(-pq/z)+J_\epsilon(q/z)-J_\epsilon(1),\\
\mathcal V_\epsilon&=J_\epsilon(-q)-J_\epsilon(-pq/z)+J_\epsilon(p/z)-J_\epsilon(1).
\end{aligned}
\]

代回 C1：

\[
\boxed{\displaystyle C=-\frac{\rho^{3\epsilon/2}}{z\epsilon^2}
\int_{-1}^zdp\int_{-1}^zdq\,(1+pq/z)^{\epsilon-1}
\mathcal U_\epsilon\mathcal V_\epsilon.}
\tag{C2}
\]

两个 J 差各从 epsilon 一次开始，因此这里只需显式 epsilon^-2；不引入额外 epsilon^-3 后再相消的一阶分子。C2 中 J 已求值，但 C2 单独仍是二维积分表示；下一步才实际完成固定阶求值。

## 3. 任意阶的有限 GPL 规则

将 p/q 的四个符号区间映射到单位方形，原核的剩余有理分母为 1−cuv。它使 u 积分前添运动字母 1/(cv)。两个 J 差及对数幂全部展开为 GPL，作有限 shuffle；随后用字母变动公式积完 v。所有边界由原零点和共同端点固定。

[完整有限词公式 F2–F10](../finite-word-rule.md)明确列出每个源的字母、符号、指标范围和系数。用其中已求值的词和 $\mathscr C_N$，有

\[
\boxed{\displaystyle C^{(q)}=-\sum_{n=0}^q
\frac{(3\log\rho/2)^{q-n}}{(q-n)!}\mathscr C_{n+2}.}
\tag{C3}
\]

源次数只需到 q+2，各项总权重 q+6。这是按任意 q 生成有限 GPL 的规则，未留下待解的微分方程或新的积分常数；并未声称它已压缩成少数 M/N 母函数。

## 4. 前两阶揭示的函数层级

保持原项目的 M/N 指标和阶乘约定，$N_{k,j}=M_{k,j+1}-\log z\,M_{k,j}$。直接结果的零阶为

\[
\boxed{C^{(0)}=M_{0,1}^3/18=(H^{(0)})^3/18.}
\]

一阶可压成

\[
\begin{aligned}
C^{(1)}={}&\frac16M_{1,5}+\frac{13}{12}M_{0,1}M_{1,3}
-\frac32M_{1,1}M_{0,3}-\frac54M_{0,2}N_{1,1}\\
&-\frac{13}{6}\zeta_2M_{1,3}+\frac52\zeta_2N_{1,2}
-\frac{13}{6}\zeta_3M_{0,3}+40\zeta_5M_{0,1}\\
&-10\zeta_2^2M_{1,1}-6\zeta_2\zeta_3M_{0,1}.
\end{aligned}
\tag{C4}
\]

在声明的有理常系数、所有齐次乘积和登记周期空间中，C1 不能只用 M；即使加入已知 log(rho) 也不能。加入 N 则可以。二阶 C2 已超出声明的 M/N 空间，加入 log(rho) 和已有权重六 E/tau 的乘积仍未闭合。rho 因子剥离后的系数也作了单独比较；这个不足不只是已知归一化造成。每个非成员结论保留精确对偶见证，不能将它外推到所有可能的一维源或允许任意有理 z 系数的空间。

这使 C 成为较适合观察结构逐阶出现的对象：有限项退化到 M，一阶需要 N，二阶需要更丰富的已声明源空间。进一步比较表明，C2 仍落在此前独立构造的 Euler/Beta/Mellin 源及其乘积内；C3 才超出这个较大空间，而加入独立双 J 源后可以表示。它没有证明二阶或三阶出现了最小意义下一个全新主积分。

[全部零至三阶显式系数](coefficients.md)已列出，分别有 4、24、84、234 个登记 Hlog/周期坐标；这些数目描述当前表示，不是主积分数。[共同结构对照](../derivation.md)记录各候选空间的边界。

## 5. 验证与物理连接

C1 的降维由原积分域得到，并以 epsilon=1,2,3、一般符号 z 的原六维多项式积分独立核对；实际 Wolfbook/MCP 另核 epsilon=2。非整数 epsilon=0.23、0.37，在 z=0.7、1.4 的完整 C2 与原有序核积分一致，最大观测差约 3.83e−12，非区间误差证书。

直接词构造不使用旧 A/B。Maple 先求值并写出全部新源和零至三阶系数，再读入预先序列化的旧 A/B/H/Q 比较文件，用 $C=H^3/6-HQ+2A+2B$ 做精确回归；四阶位置全部一致。实际状态和系数文件见[计算记录](../../../../../../workflow/checks/pure-three-cd-20261016/README.md)。

最后可与 D 一起恢复 A、B 和原物理 I；C 本身不乘 a1=0 后当成零函数使用。所有新结果是几何 C 的系数，未重复引入 Gamma、Euler、C4 或颜色因子。formal=not_run。
