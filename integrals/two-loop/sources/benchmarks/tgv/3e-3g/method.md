# 3e-3g：规范代表、归一化与解析求解

## 输入与唯一 I

本批只保存每图一个规范代表（原 A），不保存独立旋转输入。这里取消的是输入的 A/B 双份制度；颜色生成元和其他数学字母不因此更名。

令 e=ε、z=z1>0、λ=1−2e，定义
$$
G_2(e)=e^{-2\gamma_Ee}\Gamma(1-2e),\quad
C=\frac14 C_AC_F(\alpha_s/\pi)^2
e^{-2\gamma_Ee}\Gamma(2-2e)
[\pi(-2v_{14}+i0)]^{2e}.
$$
补偿仅施加一次。原 C、未补偿 IBP 前因子和原始文件哈希均保留。四个三发射拓扑的 C4 轨道大小均为4；本批按与2ne相同的系数提取约定定义
$$
I(z,e)=\frac{\lambda G_2(e)}2
\frac{J(z,e)+z^{2e}J(1/z,e)}{1+z^{2e}}.
$$
03/04 的原域反射给出 $J(1/z)=z^{-2e}J(z)$，从而简化为 $\lambda G_2J/(1+z^{2e})$。01/02 则具有 $J=z^{2e}j(e)$，不能不经核对便套用同一 self-crossing 简化。

J 使用原 s 参数及依赖边界。单纯形测度显式写为
$$
d\mu_e=d^3\alpha\,\delta(1-\alpha_1-\alpha_2-\alpha_3)
(\alpha_1\alpha_2\alpha_3)^{-e}.
$$
上游空的 delta_constraints 并不表示三个α独立；本次在规范输入中补明单纯形。超立方版本仅作为换元核验资料。

## 01、02：全 Gamma 表达

反射有序域并交换α1、α3可将01严格映至02，C和轨道倍数相同，因此完整I相等。把同边两点记为0<x<y<1，孤立的Feynman参数首先积分得到1/e；另外两个尺度积分各给1/(2e)。剩余原域的一维积分是
$$
\frac{J_{01}}{z^{2e}}=\frac1{16e^3}
\int_0^1 db\,
b^{-e-1}(1-b)^{-e}(2b-1)
\frac{1-(1-b)^{-1+2e}}{-1+2e}.
$$
端点差必须整体保留。按Beta恒等式延拓后，设 $B_e=B(1-e,1-e)$，
$$
J_{01}=J_{02}
=\frac{z^{2e}}{16(1-2e)e^4}
[2\Gamma(1+e)\Gamma(1-e)-B_e],
$$
$$
\boxed{I_{01}=I_{02}=
\frac{G_2(e)}{32e^4}
[2\Gamma(1+e)\Gamma(1-e)-B_e].}
$$
在 e=1/2 等中间表示的可去点取完整组合的极限。三个非整数e点的原一维积分已与该式直接比较。

## 03：同边接触面与单核生成函数

记 $D_z(t,u)=(1-t)(1-u)+ztu$。原同边参数0<x<y<1，星型二次式为
$$
S=\alpha_1[\alpha_2D_z(t,x)+\alpha_3D_z(t,y)].
$$
先用α1尺度积分，再在原有序三角域做IBP。两个端点面和出现两次的接触面对角共同产生
$$
F_e(u)=\int_0^u db\,b^{-e}(1-b)^{-1-e}
+\int_u^1db\,b^{-1-e}(1-b)^{-e}-2B_e
=\frac{[u(1-u)]^{-e}-B_e}{e}.
$$
最后一个等式可由导数相同、积分平均值为零证明；没有丢弃接触面。

沿用2ne的M定义，并设
$$
\mathscr C(e,z)=\frac{\pi}{e\sin(\pi e)}
[\cos(\pi e)-\cosh(e\log z)],
$$
$$
H_0(\delta,z)=
-\frac{(1+z)z^{-\delta/2}}
{\Gamma(1-\delta)\Gamma(1+\delta)^2}
G^{3,3}_{3,3}\left(z\middle|
\begin{matrix}\delta,0,0\\0,\delta,\delta\end{matrix}\right).
$$
于是
$$
\boxed{I_{03}=\frac{G_2(e)}
{8e^2\cosh(e\log z)}
[\mathscr C(e,z)-B_eH_0(2e,z)].}
$$
这里 $\mathscr C=M_{0,1}+e^2M_{0,3}/6+\cdots$，
$H_0(2e)=M_{0,1}+eM_{1,1}+e^2(M_{2,1}/2+M_{0,3}/6)+\cdots$。
两者零阶相消，产生简单极点。直接原有序域四维Sobol检查（只解析消去一个α尺度）与完整式一致，未只比较两个答案之间的展开。

## 04：新会话重算 Hlog

源核为
$$
S=\alpha_1\alpha_2y(1-t)+\alpha_1\alpha_3D_z(t,v)
+z\alpha_2\alpha_3v(1-y).
$$
完整P与本次独立Gram导数一致。原IBP保留六个cusp面和recoil。本轮在新目录重跑了旧资料的输入证书、cusp积分、对边边界、困难边界、优化recoil、Hlog组装和M匹配；没有用已知有限项作为积分输入。

关键正规化后的recoil在正半轴图表中只有显式简单极点；HyperInt的发散检查开启，原生dimregPartial再次复核Euler算子。与两个cusp双极点合并后，I最高为二阶极点。新根的S、P、sym=4以及Euler/Gamma规范均逐项连接到此次重算，不自动继承旧材料的状态。

## 04 的完整ε解析表示：逐项积完的收敛级数

记 $p=2-2e,L=1+z$，将有限多项式写为
$$
P(\boldsymbol s,\boldsymbol\alpha)
\left(1-\frac{S(\boldsymbol s,\boldsymbol\alpha)}L\right)^k
=\sum_{\boldsymbol a,\boldsymbol b}
c^{(k)}_{\boldsymbol a,\boldsymbol b}(z)
\boldsymbol s^{\boldsymbol a}\boldsymbol\alpha^{\boldsymbol b}.
$$
04的三个Wilson参数独立在[0,1]，故定义已积完的有限块
$$
B_k(e,z)=\sum_{\boldsymbol a,\boldsymbol b}
c^{(k)}_{\boldsymbol a,\boldsymbol b}(z)
\prod_{i=1}^3\frac1{a_i+1}
\frac{\prod_{j=1}^3\Gamma(b_j+1-e)}
{\Gamma(b_1+b_2+b_3+3-3e)}.
$$
完整解为
$$
J_{04}(z,e)=L^{-p}\sum_{k=0}^{\infty}
\frac{(p)_k}{k!}B_k(e,z),\qquad
\boxed{I_{04}=\frac{e^{-2\gamma_Ee}\Gamma(2-2e)}
{1+z^{2e}}J_{04}.}
$$
有限块内先做全部有符号加和，不能把它拆成独立的无穷多重求和。

在正z下 $0\le S\le L/3$。先取 $1/2<\Re e<1$，每个分子αiαj项可用对应距离的幂作绝对可积界，Wilson相邻距离的幂小于1，Dirichlet指数也大于−1。二项级数因而可逐项积分；随后以该解析函数的亚纯延拓定义e=0附近。Laurent系数另由上述正规化/HyperInt路线计算，不在收敛域外直接逐项展开这个无穷和。

这给出完整ε依赖的显式、逐项积完级数，不声称有限特殊函数压缩或快速数值收敛。12个独立多项式幂的终止级数检查通过。

## 输出边界

四图均给出 e⁻⁶至e⁰的完整表（不出现的负幂为零），余项O(e)，并保留完整Hlog。TGV系数可以混合不同超越权重，不能机械套用纯交换图的统一权重规则。

这里的解析结果与核验均为计算侧成果，formal=not_run；没有增加Kernel/Lean支持或证明物理复分支。源C保留i0，归一化I的求值范围为正实z。
