# D：循环配对纯对边积分的直接生成结构

2026-10-16。保留辅助几何积分 $D=R_{231}=R_{312}$，其完整颜色为 CF³−CA CF²/2，当前 a1 为零。原端点排序和单核与 [C](../C/derivation.md)相同，下面从 D 自身的配对积去外端点。

## 1. 先保留 R231 的原排序

\[
D=\int_{\Omega_3}k_\epsilon(t_1,s_2)k_\epsilon(t_2,s_3)
k_\epsilon(t_3,s_1)\,d^3t\,d^3s.
\]

固定 t=t2。先积分 t1∈(0,t)、t3∈(t,1)，得到 $F_s(t,s_2)$ 与 $f'(s_1)-F_s(t,s_1)$。再对 s1∈(0,s2) 积分，其结果是 $f(s_2)-F(t,s_2)$。最后把 s3 从 s2 到 1 积去：

\[
\int_{s_2}^1k_\epsilon(t,s_3)ds_3=f'(t)-F_t(t,s_2).
\]

因此四个外端点已实际积去，留下

\[
\boxed{D=\int_0^1dt\int_0^1ds\,
[f'(t)-F_t(t,s)]F_s(t,s)[f(s)-F(t,s)].}
\tag{D1}
\]

这里选用的 s 是原 s2。D1 的三个因子都对应明确有序子域，不能把最后的互补矩形替换成 F。

## 2. 原函数导数变成两条幂差

仍取 $p=(1+z)t-1,q=(1+z)s-1,\rho=z/(1+z)^2$。单核原函数的显式导数给出

\[
f'(t)-F_t=-\frac{(1+z)\rho^{\epsilon/2}}{\epsilon p}
[(1+p)^\epsilon-(1+pq/z)^\epsilon],
\]
\[
F_s=-\frac{(1+z)\rho^{\epsilon/2}}{\epsilon q}
[(1+pq/z)^\epsilon-(1-q/z)^\epsilon],\qquad
f(s)-F=-\frac{\rho^{\epsilon/2}}\epsilon\mathcal V_\epsilon,
\]
\[
\mathcal V_\epsilon=J_\epsilon(-q)-J_\epsilon(-pq/z)+J_\epsilon(p/z)-J_\epsilon(1).
\]

乘回 Jacobian 得到

\[
\boxed{\displaystyle D=-\frac{\rho^{3\epsilon/2}}{\epsilon^3}
\int_{-1}^z\frac{dp}{p}\int_{-1}^z\frac{dq}{q}
[(1+p)^\epsilon-(1+pq/z)^\epsilon]
[(1+pq/z)^\epsilon-(1-q/z)^\epsilon]\mathcal V_\epsilon.}
\tag{D2}
\]

三因子各为 O(epsilon)，源分子从 epsilon 三次开始；p=0、q=0 的分母由原幂差的连续差商解释。D2 和旧 Lambda 源的参数配置不同：两条幂差与 J 边界分别属于互补端点，不能只替换图号就认作同一个积分。

## 3. 任意阶系数由有限词和给出

D2 中每个幂差的 i 阶是两个重复字母 GPL 之差，J 项的 k 阶是含一个前零字母的 GPL。把四个 p/q 扇区映射到单位方形后，两个积分都对应前添零；其余变动字母由有限 Chen 递归处理。

用[有限词公式 F2–F10](../finite-word-rule.md)中的已求值词和 $\mathscr D_N$，得到

\[
\boxed{\displaystyle D^{(q)}=-\sum_{n=0}^q
\frac{(3\log\rho/2)^{q-n}}{(q-n)!}\mathscr D_{n+3}.}
\tag{D3}
\]

F8 明确给出 $\mathscr D_N$ 的指标、shuffle、字母变动、端点和符号。给定任意 q 的运算都是有限词运算，输出自变量 1、字母 0/1/−z/−1/z 的 GPL 乘积，总权重 q+6。它不是由 D 的旧系数拟合出的余量定义。

## 4. 与 C 的不同及已完成检查

D 的有限项已经超出声明的 M/N 齐次多项式空间，C 的有限项则是 M01³/18。D 因此继续携带已有三核有限项中的新结构；恢复颜色零图并没有使该结构消失。更大的独立 Euler/Beta/Mellin 空间可表示 D0、D1，但 D2、D3 均不在其中；加入已有独立双 J 源后，零至三阶均可表示。对已知 log(rho) 及旧低权重源乘积的比较另有记录，未将“非成员”提升为全阶最小性定理。

[全部零至三阶显式系数](coefficients.md)已列出，分别有 15、57、153、357 个登记 Hlog/周期坐标；这些是表示长度，不是主积分数。

D1 的原域积分与二维式在 epsilon=1,2,3 的一般 z 多项式控制中精确一致；另有实际 Wolfbook/MCP epsilon=2 核对。epsilon=0.23、0.37，z=0.7、1.4 的完整 D2 与原排序四维余核求积一致，最大观测差约 1.54e−12。数值是经验精化，不作区间认证。

直接源先求值，再用 $D=(HQ-3A-4B)/2$ 做有限阶精确回归。最后可由

\[
A=HQ-H^3/3+2C+2D,\qquad
B=H^3/4-HQ/2-3C/2-2D
\]

恢复物理代表。所有新 C/D 系数均为裸几何积分，原物理 I 的颜色和 Gamma/cosh 因子随后只乘一次。实际产物及完成阶数见[计算记录](../../../../../../workflow/checks/pure-three-cd-20261016/README.md)，formal=not_run。
