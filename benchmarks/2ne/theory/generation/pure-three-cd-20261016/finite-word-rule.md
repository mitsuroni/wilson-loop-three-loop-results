# C、D 的任意阶有限 GPL 词公式

2026-10-16。本页把原排序积分降维后的被积函数实际展开为有限 GPL 词和。所得公式按目标阶 q 直接构造；不要求先求解微分系统，也不读取旧 A/B 的系数。它是一种显式但未必紧凑的表示，不是原 M/N 的全阶闭合结论。

## 1. 两个直接源及四个矩形

记 $\rho=z/(1+z)^2$，$p=(1+z)t-1,q=(1+z)s-1$。由[原 C 积分](C/derivation.md)和[原 D 积分](D/derivation.md)得到

\[
\begin{aligned}
C&=-\frac{\rho^{3\epsilon/2}}{z\epsilon^2}
\int_{-1}^z dp\int_{-1}^z dq\,
(1+pq/z)^{\epsilon-1}\mathcal U_\epsilon(p,q)\mathcal V_\epsilon(p,q),\\
D&=-\frac{\rho^{3\epsilon/2}}{\epsilon^3}
\int_{-1}^z\frac{dp}{p}\int_{-1}^z\frac{dq}{q}
[(1+p)^\epsilon-(1+pq/z)^\epsilon]
[(1+pq/z)^\epsilon-(1-q/z)^\epsilon]\mathcal V_\epsilon(p,q),
\end{aligned}
\tag{F1}
\]
\[
\begin{aligned}
\mathcal U_\epsilon&=J_\epsilon(-p)-J_\epsilon(-pq/z)+J_\epsilon(q/z)-J_\epsilon(1),\\
\mathcal V_\epsilon&=J_\epsilon(-q)-J_\epsilon(-pq/z)+J_\epsilon(p/z)-J_\epsilon(1).
\end{aligned}
\]

将 p、q 各在零处分开，再置 $p=Pu,q=Qv$，其中 $P,Q\in\{z,-1\}$，$u,v\in[0,1]$。定义

\[
c=-PQ/z,\quad \sigma_{PQ}=\operatorname{sgn}(PQ),\quad
\alpha=-1/P,\ \beta=z/Q,\ \delta=-1/Q,\ \eta=z/P,\quad m(v)=1/(cv).
\tag{F2}
\]

四个扇区的 c 为 −z、1、1、−1/z。以下 G 的核为 $du/(u-a)$，上标表示重复字母；空字 G 为 1。

## 2. 每个 epsilon 因子的系数已经求值

定义以下有限 GPL 组合：

\[
L_i=G(m^i;u),\qquad
X_i=G(\alpha^i;u)-G(m^i;u),\qquad
Y_j=G(m^j;u)-G(\beta^j;v),
\tag{F3}
\]
\[
\begin{aligned}
\mathbb U_j={}&G(0,\alpha^j;u)-G(0,m^j;u)
+G(0,\beta^j;v)-G(0,1^j;1),\\
\mathbb V_k={}&G(0,\delta^k;v)-G(0,m^k;u)
+G(0,\eta^k;u)-G(0,1^k;1).
\end{aligned}
\tag{F4}
\]

于是原核的各项展开为 $\sum_{i\ge0}\epsilon^iL_i$、$\sum_{i\ge1}\epsilon^iX_i$、$\sum_{j\ge1}\epsilon^jY_j$，两个矩形差的 J 组合则为 $\sum_{j\ge1}\epsilon^j\mathbb U_j$ 和 $\sum_{k\ge1}\epsilon^k\mathbb V_k$。重复字母 GPL 已包含参数 Taylor 阶乘，无额外 i!j!k!。

## 3. 把两个剩余积分改成有限词操作

定义 $\mathcal E_{\ell,c}$ 为下述有限代数操作，$\ell=0$ 或 $m$。它并不是一个尚待计算的新积分函数。

1. 将同自变量 u 的 GPL 相乘，用带重数的 shuffle 化为词和，并为每个词前添 $\ell$。
2. 对得到的 $G(w(v);1)$ 使用下面的字母变动递归，化为有限个 $G(a;1)G(b;v)$。
3. 将所有自变量 v 的词 shuffle，在前面添零并取 v=1，乘回原有常数 GPL。

第二步明确如下。令 $a_0=1,a_{r+1}=0$，则

\[
\frac{d}{dv}G(a_1,\ldots,a_r;1)
=\sum_{h=1}^rG(a_1,\ldots,\widehat a_h,\ldots,a_r;1)
\frac{d}{dv}\log\frac{a_{h-1}-a_h}{a_{h+1}-a_h}.
\tag{F5}
\]

唯一运动字母为 m(v)。不涉及 m 的差的 dlog 为零；非零项只需

\[
d\log(m-0)=-dv/v,\qquad
d\log(m-a)=-dv/v+\frac{dv}{v-1/(ca)}\quad(a\ne0).
\tag{F6}
\]

若删除词已经写成 $G(a;1)G(b;v)$，对 F6 求原函数就是在 b 前添相应字母。没有运动字母的词是已知常数 GPL；含运动字母的原零点边界为零。F5 每次减少词长，因此给定任意有限 q 必定终止。相同字母的删除项先合并；相同差的 dlog 在共同极限中处理，不分别赋发散值。

在本页四种物理扇区中，F6 及原固定字母只生成

\[
\{0,1,-z,-1/z\}.
\]

所有最终 GPL 都取自变量 1。它们可能以乘积形式出现；需要的每项系数都是有限 shuffle 重数、符号与有理数。端点为 1 的常数词先保留共同 $\log(1-\xi)$，组合完整源后确认其系数消失，再取 $\xi\to1^-$。不对单项擅自删掉发散部分。

## 4. 完整 C、D 系数的通用有限和

现在定义两个完全由 F3–F6 求值的词和：

\[
\boxed{\displaystyle
\mathscr C_N=\sum_{P,Q\in\{z,-1\}}\sigma_{PQ}
\mathcal E_{m,c}\left[
\sum_{\substack{i+j+k=N\\i\ge0,\ j,k\ge1}}L_i\mathbb U_j\mathbb V_k\right],\quad N\ge2,}
\tag{F7}
\]
\[
\boxed{\displaystyle
\mathscr D_N=\sum_{P,Q\in\{z,-1\}}\sigma_{PQ}
\mathcal E_{0,c}\left[
\sum_{\substack{i+j+k=N\\i,j,k\ge1}}X_iY_j\mathbb V_k\right],\quad N\ge3.}
\tag{F8}
\]

F7 前添运动字母而不是零，来自

\[
\frac1{1-cuv}\,du=-\frac1{cv}\frac{du}{u-m(v)}.
\]

结合原 p/q 扇区的 Jacobian，恰好产生 F7 的 $\sigma_{PQ}$ 和下面的整体负号。F8 保留的是原 $dp/p\,dq/q$ 的两次零字。

完整函数可写为

\[
C=-\rho^{3\epsilon/2}\sum_{n\ge0}\epsilon^n\mathscr C_{n+2},\qquad
D=-\rho^{3\epsilon/2}\sum_{n\ge0}\epsilon^n\mathscr D_{n+3}.
\tag{F9}
\]

因此令 $a=\tfrac32\log\rho$，任意物理几何阶 q 的明确规则为

\[
\boxed{\displaystyle
C^{(q)}=-\sum_{n=0}^q\frac{a^{q-n}}{(q-n)!}\mathscr C_{n+2},\qquad
D^{(q)}=-\sum_{n=0}^q\frac{a^{q-n}}{(q-n)!}\mathscr D_{n+3}.}
\tag{F10}
\]

例如 $C^{(0)}=-\mathscr C_2$，$C^{(1)}=-\mathscr C_3-a\mathscr C_2$，$C^{(2)}=-\mathscr C_4-a\mathscr C_3-a^2\mathscr C_2/2$；D 的源下标整体加一。这里的 $\mathscr C_N,\mathscr D_N$ 已在 F7/F8 由有限词操作定义，不是把未知系数换名。其权重分别为 N+4、N+3，F10 中每项总权重 q+6。

## 5. 合法范围与实际实现

取固定正实 z，先在零附近的小 epsilon 域推导。原域放宽到单位立方体后，三条传播子分别使用不重叠的端点对。每个核在零点角附近是可积的线性和的负一次幂；取 $|\Re\epsilon|<\delta<1$ 并稍加宽幂界，可以同时支配任意给定有限次对数导数。这支持原定义的 Taylor 展开与本页有限系数抽取。p=0/q=0 的差商由完整差分定义，可去点不拆开处理。所有 J 自变量不超过 1，取原正实路径及共同端点极限。

[direct_words.py](../../../../../workflow/checks/pure-three-cd-20261016/direct_words.py)实现 F2–F10。它的输入是目标阶和四个矩形的字，未读取旧 A/B 或 Lambda 系数。生成结果先完整求值并写出，再读取旧结果做回归。[实际计算及验证入口](../../../../../workflow/checks/pure-three-cd-20261016/README.md)。

本表示比“仅给微分系统和边界”更具体，任意 q 的输出是有限 GPL 和；但词数增长较快。寻找少数一维母函数的紧凑全阶公式仍有意义，不能从有限词公式推断 M/N 全阶闭合或最小主积分数。formal=not_run。
