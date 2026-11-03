# 先导出生成函数：三核的一个自然双Beta源

日期：2026-09-28。根据用户纠正，后续顺序改为：从未展开积分的幂与边界导出生成函数，求出其方程、常数和系数规则，最后才与M/N及旧函数比较。补齐旧候选空间和扣除已知项不再作为前置步骤。

本页给出一个实际属于三核表达的子块：两个普通Beta矩自然产生一个三参数5F4生成函数及其源方程。没有读取已有系数或做基底拟合；这尚未求完两个三核目标。formal=not_run，检查见[复验](../../../../workflow/checks/natural-beta-source/README.md)。

## 1. 子块来自哪里

在[三核生成式T11–T14](three-kernel-generators.md)中，行列式剩余积分包含两个端点幂差与乘积自变量原函数。将中心变量分正负区间并缩放到单位区间后，正正扇区包含

$$
\int_0^1\frac{du}{u}\int_0^1\frac{dv}{v}
\big[(1+zu)^a-(1-u)^a\big]
\big[(1+zv)^b-(1-v)^b\big]
J_c(-zuv).
$$

a、b、c先保持独立；还原物理图时取ε指定的参数方向。展开两个差分，最后一项恰好是

$$
\boxed{\displaystyle
\mathfrak B(a,b,c;x)
=\int_0^1\frac{du}{u}\int_0^1\frac{dv}{v}
(1-u)^a(1-v)^bJ_c(xuv),}
\tag{B1}
$$

并在该扇区取x为负z。这里

$$
J_c(y)=\int_0^y\frac{(1-t)^c-1}{t}\,dt.
$$

所以参数和权重来自原积分，不是为了匹配M/N而选择。先取a、b实部大于负1、参数在零附近且x的模小于1，使用从零出发的实支或对应解析支。J在uv为零时线性消失，抵消测度的表观端点极点。

## 2. 一次完成Beta积分与求和

直接展开乘积自变量，交换求和与积分后，两个积分各给一个普通Beta函数：

$$
\begin{aligned}
\mathfrak B
&=\sum_{n\ge1}\frac{(-c)_n x^n}{n!\,n}
\mathrm B(n,1+a)\mathrm B(n,1+b)\\
&=-\frac{cx}{(1+a)(1+b)}
{}_5F_4\!\left(
\begin{matrix}1,1,1,1,1-c\\2,2,2+a,2+b\end{matrix};x\right).
\end{aligned}
\tag{B2}
$$

第一行中每个Beta矩都已经求值。分母中的额外n来自J的原函数；省去它得到另一个函数，不能混用。对负实x的延拓沿原积分的正底数支处理；本轮数值控制在收敛圆内进行，不宣称已验证所有延拓点。

例如c为1时，原J等于负y，两边直接退为

$$
\mathfrak B(a,b,1;x)=-\frac{x}{(1+a)(1+b)}.
$$

这与上参数1减c为零时的超几何终止完全一致。

## 3. 生成方程与边界也由同一求和式确定

令 $\theta=x\partial_x$，并令 $V=\theta^2\mathfrak B$。若 $V=\sum_{n\ge1}v_nx^n$，相邻级数系数之比直接给出

$$
(n+a)(n+b)v_n=(n-1)(n-1-c)v_{n-1}-c\,\delta_{n,1},\qquad v_0=0.
$$

所以在整个收敛圆内有

$$
\boxed{\displaystyle
\big[(1-x)\theta^2+(a+b+cx)\theta+ab\big]V=-cx,
\qquad V=\theta^2\mathfrak B.}
\tag{B3}
$$

选择零点解析解，并要求

$$
V(0)=0,\qquad
\mathfrak B(0)=0,\qquad
\theta\mathfrak B(0)=0.
$$

参数在零附近时，这些条件及方程固定所需幂级数，没有由旧答案拟合的常数。将参数展开写为 $V=\sum_{i,j,k\ge0}V_{i,j,k}(x)a^ib^jc^k$，源系数递推为

$$
\begin{aligned}
(1-x)\theta^2 V_{i,j,k}
={}&-\theta V_{i-1,j,k}-\theta V_{i,j-1,k}
-x\theta V_{i,j,k-1}\\
&-V_{i-1,j-1,k}
-\delta_{i,0}\delta_{j,0}\delta_{k,1}x,
\end{aligned}
\tag{B4}
$$

负指标为零。右边只含更低源次数；再由零边界积分两次得到B的源系数。这才是可执行的生成规则，不是只给一个参数积分换名字。

## 4. 与整个三核任务的距离

这个5F4块是一个已完成的自然来源示例。展开端点差分还有带相关运动学参数的幂权重；矩形耦合还产生乘积自变量幂的插入。这些伴随块必须继续从原式推导，不能因为一个子块属于5F4，就宣布全部三核也属于同一个5F4。

下一步从实际剩余块建立最少的独立指数源，研究其Euler/Beta求值、源微分和边界常数。整图与源函数的解析连接完成之后，才讨论系数族的简洁定义及它们与M/N/τ/E的转换。既有ε⁰…ε²结果作为事后回归，不作为生成函数的输入。

本轮没有新增三核Laurent系数或提高56/329的覆盖，也未运行Lean。旧数学正文与证据保持，主计划只更新当前研究顺序。
