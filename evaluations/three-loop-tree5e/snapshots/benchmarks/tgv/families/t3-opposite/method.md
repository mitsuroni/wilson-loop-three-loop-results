# T3[1]：有序双核与可因子化通道

本批为5e-3g-22/23/24/25（13通道）及44（24通道）。五图完整ε表示及ε⁻⁶…ε⁰系数已计算；固定z>0，余项O(ε)，formal=not_run。

[解表](solutions.md) · [T3[]母核](../t3-empty/method.md) · [实际核验](../../../../workflow/checks/tgv-t3-opposite/README.md)

## 1. 从原始J开始

记e=ε、λ=1−2e、ℓ=log z，
\[
F_z(t,u)=z(1-t)(1-u)+tu,\quad A_z(t)=(1+z)t-z.
\]
TGV孤立端点为t，同边端点为x<y；pure端点记p、q。原核为
\[
S=\alpha_1[\alpha_2F_z(t,x)+\alpha_3F_z(t,y)],
\quad P=-\frac{1+z}{4}A_z(t)\alpha_1(\alpha_2-\alpha_3).
\]
α1=a、α2=(1−a)b、α3=(1−a)(1−b)，Jacobian=1−a；尺度积分为1/e。对13通道，pure核是−(1+z)F_z(p,q)^(e−1)，因此
\[
J_g=\frac{(1+z)^2}{4e}
\int_{\Omega_g}dt\,dp\,dx\,dy\,dq\int_0^1db\,
\frac{A_z(t)(2b-1)b^{-e}(1-b)^{-e}}
{F_z(t,bx+(1-b)y)^{2-2e}F_z(p,q)^{1-e}}.
\]

| 图 | 边1顺序 | 边3顺序 |
|---|---|---|
| 22 | t<p | x<y<q |
| 23 | t<p | x<q<y |
| 24 | p<t | x<q<y |
| 25 | p<t | q<x<y |

所有端点在(0,1)。44的pure两个端点分别在边4和边2，与TGV三端点无排序耦合，其核为−(1+z)D_z(p,q)^(e−1)，其中D=(1−p)(1−q)+zpq。F与D不可混同；F_z=zD_(1/z)。

五图均为相同固定颜色投影、sym=4，prepared C已经含一次Euler补偿。唯一I为
\[
I_g=\frac{\lambda G(e)}8
\frac{J_g(z,e)+z^{3e}J_g(1/z,e)}{1+z^{3e}},
\qquad
G=e^{-3\gamma_Ee}\Gamma(1-e)\Gamma(1-2e).
\]

## 2. 两组完整反射

同时反向所有边内端点顺序，交换α2、α3并取z→1/z，把25映到22、24映到23。S、P、pure距离、pure分子分别缩放1/z、1/z²、1/z、1/z；两个绝对Jacobian均为1，C相同。所以
\[
J_{25}(z)=z^{3e}J_{22}(1/z),\quad
J_{24}(z)=z^{3e}J_{23}(1/z),
\]
\[
\boxed{I_{22}=I_{25},\qquad I_{23}=I_{24}.}
\]
这是完整ε等式，不仅是有限阶系数一致。不能省略α交换或同边端点逆序。

## 3. 母核密度与保留接触面的消元

沿用T3[]推导中的推前密度ρ0：在u=bx+(1−b)y下，
\[
\int_0^1db\,b^{-e}(1-b)^{-e}(2b-1)
\int_{x<y}dx\,dy\,f(bx+(1-b)y)
=\int_0^1du\,\rho_0(u)f(u).
\]
令B_e=B(1−e,1−e)、T_e(u)=[u(1−u)]^(−e)，则
\[
\rho_0'(u)=\frac{T_e(u)-B_e}{e}.
\]
ρ0两端为零。这是二圈T3的同边接触面消元，不能仅保留外端点面。

固定pure端点q，有序域x<y<q产生
\[
\rho_{\rm low}(u)=q\,\rho_0(u/q)\,\mathbf1_{0<u<q}.
\]
q<x<y产生
\[
\rho_{\rm high}(u)=(1-q)\rho_0((u-q)/(1-q))
\,\mathbf1_{q<u<1}.
\]
中间矩形x<q<y严格满足
\[
\rho_{\rm mid}=\rho_0-\rho_{\rm low}-\rho_{\rm high}.
\]
分段边界处ρ0为零，不额外遗漏δ接触项。由
\(A_z(t)F_z(t,u)^{-2+2e}
=-\lambda^{-1}\partial_uF_z(t,u)^{-\lambda}\)
完成IBP。

定义权重
\[
W_{\rm low}=T_e(u/q)-B_e,\qquad 0<u<q,
\]
\[
W_{\rm mid}=
\begin{cases}
T_e(u)-T_e(u/q),&u<q,\\
T_e(u)-T_e((u-q)/(1-q)),&u>q.
\end{cases}
\]
两个权重都从e¹开始；其零阶消失并非事后数值抵消。

再定义
\[
\mathcal H_{22}=(1+z)^2
\int_{0<t<p<1}\!dt\,dp\int_{0<u<q<1}\!du\,dq\,
F_z(t,u)^{-\lambda}F_z(p,q)^{e-1}W_{\rm low},
\]
\[
\mathcal H_{23}=(1+z)^2
\int_{0<t<p<1}\!dt\,dp\int_{[0,1]^2}\!du\,dq\,
F_z(t,u)^{-\lambda}F_z(p,q)^{e-1}W_{\rm mid}.
\]
则J22=H22/(4λe²)、J23=H23/(4λe²)。这些结构函数保留几何域和核；不是新增基底生成元。

## 4. 23/24的组合先消除一层排序

原J23+J24的t,p两域恰好拼成整个单位正方形。设
\[
Z_\delta(z,u)=(1+z)\int_0^1dt\,F_z(t,u)^{-1+\delta}
=\frac{u^\delta-[z(1-u)]^\delta}
{\delta[u-z/(1+z)]},
\]
\[
\mathcal K=\int_{[0,1]^2}du\,dq\,
Z_{2e}(z,u)Z_e(z,q)W_{\rm mid}(e,u,q).
\]
于是
\[
J_{23}+J_{24}=\frac{\mathcal K}{4\lambda e^2},
\qquad
\boxed{I_{23}=I_{24}=\frac{G\,\mathcal K}
{32e^2(1+z^{3e})}.}
\]
这是23最方便的求值路线，直接把两个不同排序的J组合成物理I。独立原J数值检查仍分别保留23与24。

22先积p：
\[
(1+z)\int_t^1dp\,F_z(p,q)^{e-1}
=\frac{q^e-F_z(t,q)^e}{e[q-z/(1+z)]}.
\]
再对t→u→q积分。Z与此式的内部表观极点都是可去点，必须保留分子差分。

令
\[
\mathcal E_{22}=\frac{z^{-3e/2}\mathcal H_{22}(z)
+z^{3e/2}\mathcal H_{22}(1/z)}2,\quad
\mathcal E_{23}=z^{-3e/2}\mathcal K,
\]
则
\[
I_{22}=\frac{G\,\mathcal E_{22}}{32e^2\cosh(3e\ell/2)},
\quad
I_{23}=\frac{G\,\mathcal E_{23}}{64e^2\cosh(3e\ell/2)}.
\]
本次新鲜计算H22、K的e¹/e²，并投影到物理偶结构。因为G/cosh=1+O(e²)，
\[
c_{22,-1}=\mathcal E_{22}^{(1)}/32,\quad
c_{22,0}=\mathcal E_{22}^{(2)}/32,
\]
\[
c_{23,-1}=\mathcal E_{23}^{(1)}/64,\quad
c_{23,0}=\mathcal E_{23}^{(2)}/64.
\]
所以22–25均只有简单极点。23/24两个非零系数权重分别为5、6；22/25还含较低权重。

## 5. 44：J因子化后的唯一I

44的域和核严格分离，
\(J_{44}=J_{\rm T3}^{F\text{框架}}J_{\rm pure}^{D\text{框架}}\)。
旧二圈T3输入在D框架，本次精确验证
\(S_F=zS_D(1/z)\)、\(P_F=z^2P_D(1/z)\)，并保留原排序和α测度。母图自交叉关系使J_T3^F=J_T3^D，但不是把两个I直接相乘。

沿用已有单对边母函数
\[
H_0(\delta,z)=
-\frac{(1+z)z^{-\delta/2}}
{\Gamma(1-\delta)\Gamma(1+\delta)^2}
G^{3,3}_{3,3}\!\left(z\middle|
\begin{matrix}\delta,0,0\\0,\delta,\delta\end{matrix}\right),
\]
并令
\[
\mathscr C(e,z)=\frac{\pi}{e\sin(\pi e)}
[\cos(\pi e)-\cosh(e\ell)],\quad
\Delta=\mathscr C-B_eH_0(2e,z).
\]
则
\[
J_{\rm T3}=\frac{z^e\Delta}{4\lambda e^2},\quad
J_{\rm pure}=z^{e/2}H_0(e,z),
\]
\[
\boxed{I_{44}=
\frac{G\,H_0(e,z)\Delta}
{32e^2\cosh(3e\ell/2)}.}
\]
这给出44的完整ε特殊函数表达。令M0=M_(0,1)、M1=M_(1,1)、M2=M_(2,1)，则
\[
H_0(e)=M0+\frac e2M1+O(e^2),
\]
\[
\Delta=-e(2M0+M1)
-e^2[(4-\zeta_2)M0+2M1+M2/2]+O(e^3).
\]
因此44也是简单极点，并有紧凑式
\[
\boxed{c_{44,-1}=-\frac{M0(2M0+M1)}{32},}
\]
\[
\boxed{c_{44,0}=-\frac{(4-\zeta_2)M0^2+3M0M1
+(M0M2+M1^2)/2}{32}.}
\]

## 6. 22–25的完整ε：已积分双核级数

结构积分方便推导；完整交付另给不含剩余积分的显式级数。设
a_T=2−2e、a_P=1−e、Λ=1+z，展开有限多项式
\[
A_z(t)(2b-1)
[1-F_z(t,bx+(1-b)y)/\Lambda]^k
[1-F_z(p,q)/\Lambda]^l
=\sum_{a,h,c,d,j,i}q^{k,l}_{a h c d j i}(z)
t^ap^hx^cy^dq^jb^i.
\]
b单项式积完为B(i+1−e,1−e)。边1矩为
\[
m^{(1)}_{22,23}=\frac1{(a+1)(a+h+2)},\quad
m^{(1)}_{24,25}=\frac1{(h+1)(a+h+2)}.
\]
边3矩为
\[
m^{(3)}_{22}=\frac1{(c+1)(c+d+2)(c+d+j+3)},
\]
\[
m^{(3)}_{23,24}=\frac1{(c+1)(c+j+2)(c+d+j+3)},\quad
m^{(3)}_{25}=\frac1{(j+1)(j+c+2)(j+c+d+3)}.
\]
记有限完整和
\[
B_{g;k,l}=\sum q^{k,l}_{a h c d j i}
B(i+1-e,1-e)m^{(1)}_g m^{(3)}_g.
\]
那么
\[
\boxed{J_g=\frac{(1+z)^2}{4e\Lambda^{3-3e}}
\sum_{N=0}^\infty\sum_{k=0}^N
\frac{(a_T)_k(a_P)_{N-k}}{k!(N-k)!}\,B_{g;k,N-k}.}
\]
每块先作有限有符号求和。结合唯一I定义即可得到四图的完整ε；不声称有限特殊函数压缩。

在z>0、实1/2<e<1，FT、FP都不超过max(1,z)<Λ。分开b<1/2与b>1/2后，FT分别由F(t,y)/2或F(t,x)/2控制；再用独立(t,x或y)、(p,q)单位方块上的角点可积界控制乘积。β端点指数−e>−1。由支配收敛交换级数和积分，复条带紧集同样可控，再亚纯延拓到e=0邻域。未在初始收敛条带外直接展开无穷和。

## 7. 基底、余项和核验边界

35个Laurent位置中10个非零；ε⁻⁶…ε⁻²全零。22/25有限项需要既有
\((\tau^+_{D,4}-6E)/64\)，23/24有限项需要既有
\((7\tau^+_{D,4}+10E)/192\)；其余部分均在M/N及周期中。44仅需M及ζ2。无新生成元，不宣称全ε IBP最小性。

展开后的权重只有端点/对角线对数。取|e|<δ、δ足够小（例如δ<1/16），在u=q r图表中用r^(−δ)(1−r)^(−δ)控制低三角；中间矩形使用相同两图表。角点FT与两个非负局部变量之和同阶，纯核积分后的差商仅增加可积对数；对角q=u的对数亦可积。由共同可积界，H/K的e¹/e²展开余项为O(e³)，最后I余项O(e)，固定z>0。

实际完成：5原核/C连接、2完整反射与4负例、8个密度矩、3个t/p分区探针；新鲜HyperInt四个原结构系数，发散检查开启，Maple CPU约32秒；35项基底/Hlog精确回代。9个完整I余项样本、10组原七维QMC及12个有限级数块检查通过。原七维最大相对差约0.105%，8次scramble×2^17点，经验误差。22/23用独立Gauss–Jacobi低维求积并加密节点检查，44用55/65位一维求积；不把全部样本称为高精度或严格包围。

原输入、分类README/catalog、历史工作区及此前成果原字节保留。未修改Kernel/Lean或外部montopar，TGV未进入Kernel默认根；formal=not_run。没有物理复z/i0延拓或新形式整链证明。
