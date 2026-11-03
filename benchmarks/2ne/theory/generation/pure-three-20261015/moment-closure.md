# 纯三对边：完整边界矩源的有限微分闭合

日期：2026-10-15。本文从原端点幂积分推导三核源 \(\mathcal G\) 的一个十分量、四阶微分系统。系统包含全部四个符号区间及全部端点幂差，零点解析边界完全固定；它并非此前单个双 Beta 子块的方程，也不使用 Hlog 答案拟合。

这完成了 \(\mathcal W\) 中剩余源的一个有限微分表示。另一个三核源 \(\Lambda\) 不由本文求解，因此本文单独不宣称两个物理积分已经全部闭合。来源为[三核原函数推导 T11–T14](../three-kernel-generators.md)，物理连接使用[完整 I 链](../natural-beta-source.md)。

## 1. 待求源与直接端点矩

固定 \(z>0\)，首先取小正实 \(\epsilon\)，随后延拓到零附近。令

\[
a_\epsilon(p)=(1+p)^\epsilon-(1-p/z)^\epsilon,
\qquad J_\epsilon(y)=\int_0^y\frac{(1-t)^\epsilon-1}{t}\,dt,
\]
\[
\mathcal G(x)=\int_{-1}^{z}\!\int_{-1}^{z}
\frac{a_\epsilon(p)a_\epsilon(q)}{pq}J_\epsilon(-xpq/z)\,dp\,dq.
\tag{M1}
\]

零分母处由连续差商解释。物理值为从 \(x=0\) 沿正实轴延拓的 \(x\to1^-\) 极限。已有母核推导给出

\[
\mathcal G(x)=\sum_{n\ge1}\frac{(-\epsilon)_n}{n!\,n}
\left(-\frac{x}{z}\right)^n T_n^2,
\qquad T_n=\int_{-1}^{z}p^{n-1}a_\epsilon(p)\,dp.
\tag{M2}
\]

这一展开先在 \(|x|<\min(z,1/z)\) 绝对收敛。新步骤是直接求端点矩的有限递归，而不是逐个计算其中的超几何函数。

定义

\[
U_n=\int_{-1}^{z}p^{n-1}(1+p)^\epsilon dp,
\qquad V_n=\int_{-1}^{z}p^{n-1}(1-p/z)^\epsilon dp.
\]

对 \(p^n(1+p)^{\epsilon+1}\) 与 \(p^n(z-p)^{\epsilon+1}\) 积分求导，两个非零端点分别在 \(p=z\) 与 \(p=-1\)。保留端点符号后得到

\[
\begin{aligned}
(n+\epsilon+1)U_{n+1}&=-nU_n+(1+z)^{\epsilon+1}z^n,\\
(n+\epsilon+1)V_{n+1}&=znV_n+(1+z)^{\epsilon+1}z^{-\epsilon}(-1)^n.
\end{aligned}
\tag{M3}
\]

初值来自普通幂原函数：

\[
U_1=\frac{(1+z)^{\epsilon+1}}{1+\epsilon},\qquad
V_1=\frac{(1+z)^{\epsilon+1}z^{-\epsilon}}{1+\epsilon}.
\tag{M4}
\]

这些关系对 \(\Re\epsilon>-1\) 成立；不把端点幂置零的论证外推到边界幂不消失的参数区域。

## 2. 对角矩平方的十维递归

记 \(A=(1+z)^{\epsilon+1}\)、\(b=z^{-\epsilon}\)，并写
\(U_n=A u_n\)、\(V_n=A b v_n\)。于是

\[
d_nu_{n+1}=-nu_n+z^n,\qquad
d_nv_{n+1}=znv_n+(-1)^n,\qquad d_n=n+\epsilon+1,
\tag{M5}
\]
\[
u_1=v_1=(1+\epsilon)^{-1},\qquad
T_n=A(u_n-bv_n).
\]

将四个量 \((u_n,v_n,Z_n,R_n)=(u_n,v_n,z^n,(-1)^n)\) 的对称平方写为

\[
\mathbf y_n=(u_n^2,u_nv_n,v_n^2,u_nZ_n,u_nR_n,
v_nZ_n,v_nR_n,Z_n^2,Z_nR_n,R_n^2)^T.
\tag{M6}
\]

存在显式二次多项式矩阵 \(S(n)\)，使

\[
d_n^2\mathbf y_{n+1}=S(n)\mathbf y_n.
\tag{M7}
\]

为便于检查，将十行完整列出。右边各分量均在同一 \(n\) 处：

\[
S(n)\mathbf y=\begin{pmatrix}
n^2y_1-2ny_4+y_8\\
-zn^2y_2-ny_5+zny_6+y_9\\
z^2n^2y_3+2zny_7+y_{10}\\
-znd_ny_4+zd_ny_8\\
nd_ny_5-d_ny_9\\
z^2nd_ny_6+zd_ny_9\\
-znd_ny_7-d_ny_{10}\\
z^2d_n^2y_8\\
-zd_n^2y_9\\
d_n^2y_{10}
\end{pmatrix}.
\tag{M8}
\]

矩阵上三角，只有前三行需要乘积矩，后面的量依次由更简单的边界量驱动。这里“十维”描述离散矩状态；下面给出的是十个四阶方程，不能称为十维一阶 Pfaffian 系统。

## 3. 完整源的微分方程及全部初值

定义

\[
h_n=\frac{(-\epsilon)_n}{n!\,n}\left(-\frac1z\right)^n,
\quad \mathbf f_n=h_n\mathbf y_n,
\quad \mathbf F(x)=\sum_{n\ge1}\mathbf f_n x^n,
\quad\theta=x\partial_x.
\]

相邻系数之比给出

\[
z(n+1)^2(n+\epsilon+1)^2\mathbf f_{n+1}
+n(n-\epsilon)S(n)\mathbf f_n=0.
\tag{M9}
\]

因此

\[
\boxed{\left[z\theta^2(\theta+\epsilon)^2\mathbf1
+x\theta(\theta-\epsilon)S(\theta)\right]\mathbf F(x)
=\epsilon(1+\epsilon)^2\mathbf y_1\,x.}
\tag{M10}
\]

式中 \(x\) 在所有 \(\theta\) 算子的左边，不得调换次序。初值为

\[
\mathbf F(0)=0,\qquad
\mathbf y_1=(d^{-2},d^{-2},d^{-2},z/d,-1/d,z/d,-1/d,z^2,-z,1)^T,
\quad d=1+\epsilon.
\tag{M11}
\]

要求在零点解析。对于零附近的 \(\epsilon\)，(M10) 的 \(x^1\) 系数固定
\(\mathbf f_1=(\epsilon/z)\mathbf y_1\)，其余系数依 (M9) 唯一确定。没有留给有限阶拟合的积分常数。

原源由固定投影恢复：

\[
\boxed{\mathcal G(x)=(1+z)^{2\epsilon+2}
\left(F_1(x)-2z^{-\epsilon}F_2(x)+z^{-2\epsilon}F_3(x)\right).}
\tag{M12}
\]

最后三个分量本身已有标准原函数值：

\[
F_8=J_\epsilon(-zx),\qquad F_9=J_\epsilon(x),\qquad
F_{10}=J_\epsilon(-x/z).
\tag{M13}
\]

这进一步说明全系统保留了正正、混合和负负端点区间的三类位置。系统每个分量一般从 \(O(\epsilon)\) 开始，但物理投影满足 \(\mathcal G=O(\epsilon^3)\)：这是原来的两个端点幂差各为 \(O(\epsilon)\)、\(J_\epsilon=O(\epsilon)\) 的结果。不能把投影前的低阶项误当原图极点。

## 4. 任意 epsilon 阶的固定有限递归

将明确矩阵按 \(\epsilon\) 展开为
\(S=S^{[0]}+\epsilon S^{[1]}+\epsilon^2S^{[2]}\)。令 (M10) 左边为
\(\mathscr P_0+\epsilon\mathscr P_1+\epsilon^2\mathscr P_2+\epsilon^3\mathscr P_3\)，其中

\[
\begin{aligned}
\mathscr P_0&=z\theta^4\mathbf1+x\theta^2S^{[0]}(\theta),\\
\mathscr P_1&=2z\theta^3\mathbf1+x\big(\theta^2S^{[1]}(\theta)-\theta S^{[0]}(\theta)\big),\\
\mathscr P_2&=z\theta^2\mathbf1+x\big(\theta^2S^{[2]}(\theta)-\theta S^{[1]}(\theta)\big),\\
\mathscr P_3&=-x\theta S^{[2]}(\theta).
\end{aligned}
\tag{M14}
\]

右边 \(\mathbf R=\epsilon(1+\epsilon)^2\mathbf y_1x\) 是次数三的 \(\epsilon\) 多项式。若 \(\mathbf F=\sum_{q\ge0}\epsilon^q\mathbf F_q\)，则

\[
\boxed{\mathscr P_0\mathbf F_q=\mathbf R_q
-\mathscr P_1\mathbf F_{q-1}-\mathscr P_2\mathbf F_{q-2}
-\mathscr P_3\mathbf F_{q-3}.}
\tag{M15}
\]

负指标为零，每阶使用零点解析且值为零的解。这是一个不随目标阶数增长的有限源系统及其任意阶规则。

这组系统还允许证明任意阶可以逐次积分求出，而不只是递归地指定一个微分方程。令

\[
(r_1,\ldots,r_{10})=(0,0,0,1,1,1,1,2,2,2),
\quad(d_1,\ldots,d_{10})=(1,-z,z^2,-z,1,z^2,-z,z^2,-z,1).
\]

直接代入 (M8)，得到每个对角零阶算子的因式分解

\[
(\mathscr P_0)_{ii}=\theta^{r_i}(z+d_ix)\theta^{4-r_i}.
\tag{M15a}
\]

注意这里 \(d_i\) 是十个常系数，与 (M5) 中带指标 \(n\) 的 \(d_n\) 无关。恒等式使用 \(\theta^r x=x(\theta+1)^r\)，故不能当作普通交换多项式因式分解。

按 epsilon 阶数递增、同阶分量从十到一处理。上三角结构保证右边的同阶分量已经求得，低 epsilon 阶也已知。对第 \(i\) 个方程，依次做 \(r_i\) 次 \(dx/x\) 积分、除以 \(z+d_ix\)、再做 \(4-r_i\) 次 \(dx/x\) 积分即可；每个积分选零点解析且值为零的原函数。源在零点消失，因此不需要新边界常数。

这些操作、求导及有理函数乘法保持有限超对数表达，辅助变量的字母只可能为

\[
\{0,-z,1,-1/z\}.
\tag{M15b}
\]

因此每个有限 epsilon 阶均有可执行的有限有理函数／超对数表达，且字母集不随阶数增长。这是辅助变量 \(x\) 的结论；物理 \(x=1\) 的端点值仍需选定极限。本文不据此宣称变换为运动学变量 \(z\) 后已经得到最小字母表，也不宣称某个有限 M/N 多项式空间全阶封闭。本文没有求规范 epsilon 形式或最小一阶系统。

## 5. 延拓到物理点与返回 I

写 \(S(n)=S_0+nS_1+n^2S_2\)。最高导数矩阵为 \(z\mathbf1+xS_2\)，其行列式是

\[
\det(z\mathbf1+xS_2)=z^7(1-x)^4(x+z)^3(1+zx)^3.
\tag{M16}
\]

所以除零点外，有限候选奇点只在 \(x=1,-z,-1/z\)。对正实 \(z\)，从零到一的开区间没有该系统的奇点；零点解析解因此给出原始小 \(x\) 级数的唯一延拓。\(x=1\) 处使用原积分指定的单侧极限。源本身有限并不保证四阶导数也有限，不在该点额外施加解析条件。

已知桥接式为

\[
\mathcal W=-\frac{\rho^{3\epsilon/2}}{\epsilon^3}
\left[\mathcal G(1)-2T_0T_J+J_\epsilon(-1/z)T_0^2\right],
\qquad \rho=\frac z{(1+z)^2}.
\tag{M17}
\]

其中 \(T_0,T_J\) 已在低圈求值。于是 (M10)–(M12) 给出了完整 \(\mathcal W\) 的有限生成系统，包含此前 \(\mathcal G_{\rm rest}\) 的全部内容。对应物理组合

\[
I_{6e-20}+4I_{6e-21}=\frac{P}{16}\left(\frac{H^3}{6}-\mathcal W\right)
\tag{M18}
\]

获得此生成表示；其中 \(P=Z(\epsilon)^3/\cosh(3\epsilon\log z/2)\)。恢复两个物理图各自的结果，还需要独立的 \(\Lambda\) 源。本文不把闭合一个组合扩大为全部图的完成。

## 6. 核验与范围

[可复现脚本](../../../../../workflow/checks/pure-three-generators-20261015/moment_check.py)使用标准有理数运算核对正整数 epsilon 下的直接多项式端点积分及整个十维传递，再在非整数 epsilon 下用独立端点求积检查 (M3)，并对 (M12) 的二阶 \(x\) 导数作原二维积分对照：

\[
\theta^2\mathcal G(x)=\frac{\epsilon x}{z}\int_{-1}^{z}\!\int_{-1}^{z}
a_\epsilon(p)a_\epsilon(q)(1+xpq/z)^{\epsilon-1}\,dp\,dq.
\tag{M19}
\]

实际结果见 [moment_report.json](../../../../../workflow/checks/pure-three-generators-20261015/moment_report.json)。脚本的有限点控制是推导的复核，不替代上述一般参数证明；其中数值精度是经验误差，没有区间认证。formal=not_run，未修改 Wilson/Lean 或旧结果对象。
