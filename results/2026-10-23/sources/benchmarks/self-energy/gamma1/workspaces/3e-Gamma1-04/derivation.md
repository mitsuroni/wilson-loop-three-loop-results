# 3e-Gamma1-04：234三边两颜色积分

## 当前物理积分的入口（2026-10-06）

从共同作用量出发，U′=igAU确定真实Wilson词；FT/Dyson规则在提取C前生效。当前[原生C/J](../../sources/montopar-r002.json)给出下式，J就是本文实际积分的原核，域、分支和IBP全保留：

$$I_{P,c}=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-3\epsilon)}{64}\frac{J_c(z)+z^{3\epsilon}J_c(1/z)}{1+z^{3\epsilon}}.$$

用下文已经求出的J代入此定义，得到[当前物理I](../../results/3e-Gamma1-04/r002.json)；组装系数为+1，Euler只补一次。以下带日期的求值过程与r001公式属于旧根历史，计算J的步骤可按同域积分线性复用，旧形式状态不继承。

2026-10-04开始求解。输入为正式 [Montopar raw r001](../../sources/montopar-r001.json)，SHA256 `3b92e23125cb6b654ad63d5c19c17fead9cfeb3bff56739dd9f07c6eda87e65a`。目标是CF CA²与CF CA TF nf两个颜色的全部极点至有限项，最终Hlog使用z为自变量。计算证据在 [gamma1-04](../../../../../workflow/checks/gamma1-04/)，工作区仅保存本数学推导。

2026-10-06完成：两个颜色的全部极点至有限项已解析求出，正式结果为[r001](../../results/3e-Gamma1-04/r001.json)。最终自变量z、字母0/−1，最高三阶；ε⁻⁶…ε⁻⁴为零。原始源和C/J不改，formal=not_run。§3–8保留逐步求解时的检查点，当前完整结果和剩余范围以§9为准。

## 1. 正式初始参数积分

三个点位于(2,1)、(3,1)、(4,1)，令x=s21、y=s31、w=s41。原外到内顺序为
$$\int_0^1dx\int_0^1dy\int_0^1dw.$$
令a=α₁、b=α₂、c=1−a−b，t₃=1−t₁−t₂，两组内部参数都在单纯形上。星核为
$$q=ab\,y(1-x)+bc\,zw(1-y)+ac[(1-x)(1-w)+zxw].$$
这是三个独立路径变量的多重仿射式，不能使用03的同边有效点/间隔域。定义P=abc、m=max(at₂t₃,bt₁t₃,ct₁t₂)。完整源为
$$J_c=\int_{x,y,w}\int_{\Delta_\alpha}da\,db\int_{\Delta_t}dt_1dt_2
\left[-\frac{P^{-\epsilon}m^{-1-\epsilon}}{3\epsilon(1+\epsilon)}K_{c,N}
+\frac{P^{-\epsilon}m^{-\epsilon}}{\epsilon}K_{c,L}\right]+J_{c,M}.$$
K是正式源全D算符作用后的完整核；q通道没有原剑鱼M，g通道有三个原剑鱼。环境导数先在独立坐标中完成，再限制到轮廓和面，度规迹D=4−2ε。

源的可选IBP包含六个端点面和非零体余项；不能只把六面相加当作完整积分。最终规范化沿用
$$I_c=\frac{G_3}{64}\frac{J_c(z)+z^{3\epsilon}J_c(1/z)}{1+z^{3\epsilon}},\qquad
G_3=e^{-3\gamma_E\epsilon}\Gamma(1-3\epsilon).$$
这是C4重数4和一次Euler后的唯一I；原raw C不改。外腿Π₁及顶点反项不在本图中。

## 2. 一般三方向分子的内部pinch

没有采用n₂=n₃假设。读完整未移位N_g、N_q，并逐项匹配正式04的N₀、L₀后，以d₁=(k−p₂)²、d₂=(k+p₁)²、d₃=k²作精确分解：
$$N_g+\frac{D-2}{4}N_q=\sum_{j=1}^3A_j(k)d_j+R(k).$$
三个A_j及R均至多对k一次；一般三方向的三次k项实际完全消去。取三个泡的重心−p₁/2、p₂/2、(p₂−p₁)/2，可先积分掉奇部；再加原剑鱼得到P_j。完整多项式保存在 [decomposition-general.json](../../../../../workflow/checks/gamma1-04/decomposition-general.json)。五项精确源连接和分解检查通过。

因此
$$I_g=-\frac{1-\epsilon}{2}I_q+I_p+I_R,$$
其中Iₚ包含pinch与原剑鱼，Iᵣ来自新的完整一般三方向R。Iᵣ的内部Laplacian为零、移位R仅线性依赖t，但其外动量算符仍含三阶，路径体余项仍在。

以上是开始时的输入与降秩；实际求值进度见下文。完整04尚未求完，不能把已完成子项当作完整I。

## 3. 已求完的完整 pinch／原剑鱼子块

2026-10-05检查点。完整原核反射已逐块核验，故本图两个颜色都满足
$$J(z)=z^{3\epsilon}J(1/z),\qquad I=\frac{G_3J}{32(1+z^{3\epsilon})}.$$
这里的$I_p$包括三个内部 pinch 和三个原剑鱼，原剑鱼只计一次。

令$\beta=3\epsilon-1$、$\kappa=1+\epsilon$，$P_t$为普通树级顶点标量核，$R_t$为已有Π₁的完整 recoil，取
$$A_*=ac(1-x),\qquad B_*=aczw.$$
两个独立 pinch 的精确普通核满足
$$C_3=\frac52P_t-\frac{4\epsilon+13}{8}q_w
 -(1+z)(\epsilon+3)A_*-\frac12R_t,$$
$$C_2=\frac32P_t+\frac{(4\epsilon+11)(1+z)}8q_y
 -\frac{4\epsilon+11}{2}R_t+(\epsilon+3)(1+z)(A_*+B_*).$$
$C_1$通过$a\leftrightarrow c$的完整路径反射与$C_3$相连。三个标记权重分别是$(abc\alpha_j)^{-\epsilon}$，不能只交换分子。

复用的是Π₁10／11的裸几何$J$中的标量、纵向与 recoil，未复用其自能物理前因子或最终$I$。实际源的$q$、$L=-P_t$、两个完整分子、标记权重、原域／原序和旧证据哈希共20项核验见 [源连接](../../../../../workflow/checks/gamma1-04/pinch-reuse-binding.json)。六个新矩是
$$\mathcal A_j=\int_{[0,1]^3\times\Delta_\alpha}
 \frac{ac(1-x)(abc\alpha_j)^{-\epsilon}}{q^{2-3\epsilon}},\quad j=2,3,$$
$$\mathcal B_2=\int_{[0,1]^3\times\Delta_\alpha}
 \frac{aczw(abc b)^{-\epsilon}}{q^{2-3\epsilon}}.$$
三者的ε⁰、ε¹均实际积分完；18个原五维正核Sobol数值点提供采样误差范围内支持，数值核验不是形式证明。

记$S_j$为去掉纵向后的Π₁标量块，$R_{c,j}$为其完整 recoil，$W_3$、$Y_2$为相应两端面差。两个裸 pinch 积分为
$$B_3=-\frac{5\beta}{4\kappa}S_3-\frac{4\epsilon+13}{8}W_3
 -\beta(1+z)(\epsilon+3)\mathcal A_3+\frac{\beta}{4\kappa}R_{c,3},$$
$$B_2=-\frac{3\beta}{4\kappa}S_2+\frac{(4\epsilon+11)(1+z)}8Y_2
 +\frac{\beta(4\epsilon+11)}{4\kappa}R_{c,2}
 +\beta(\epsilon+3)(1+z)(\mathcal A_2+\mathcal B_2).$$
因此
$$I_p=\frac{G_3B(1-\epsilon,1-\epsilon)}{8\epsilon(1+z^{3\epsilon})}
 [B_3(z)+z^{3\epsilon}B_3(1/z)+B_2(z)].$$
已得到ε⁻³至ε⁰的四个解析系数，均使用$G(\mathbf a;z)=\mathrm{Hlog}(z,[\mathbf a])$，字母0、−1。例如
$$I_{p,-3}=-\frac{5\zeta_2}{384},$$
$$I_{p,-2}=\frac{7\zeta_3}{96}-\frac{\zeta_2}{64}
 -\frac{15\zeta_2G(0;z)}{256}+\frac{5G(-1,0,0;z)}{128}
 +\frac{15\zeta_2G(-1;z)}{128}-\frac{5G(0,0,0;z)}{256}.$$
四个完整表达式见 [pinch-coefficients.txt](../../../../../workflow/checks/gamma1-04/pinch-coefficients.txt)，均通过crossing。这里的最高三阶只属于$I_p$，不据此判定整图最高极点。

## 4. 内部三角圈压缩为四个函数

在第$k$个Max扇区令$\alpha_k=1/H$、另两腿$\alpha_i=r/H$、$\alpha_j=s/H$，$H=1+r+s$。对内部参数换元后
$$t_k=1/Z,\quad t_i=(r+P)/Z,\quad t_j=(s+Q)/Z,\quad Z=H+P+Q.$$
需要的内部矩只有
$$\Lambda_{nij}=\int_0^\infty\!dP\int_0^\infty\!dQ\;
 \frac{P^iQ^jZ^{2\epsilon-n}}{(r+P)^{1+\epsilon}(s+Q)^{1+\epsilon}}.$$
夸克用$n=4,i+j\leq3$，余项用$n=2,i+j\leq1$，共13种。内部IBP将它们精确约到
$$\mu=\Lambda_{100},\quad
 a_r=r^{-\epsilon}\int_0^\infty\frac{(H+Q)^{-1+2\epsilon}}{(s+Q)^{1+\epsilon}}dQ,
 \quad a_s=a_r(s,r),\quad c_0=H^{2\epsilon}(rs)^{-\epsilon}.$$
52个内部矩ε系数与32个微分关系检查通过。这里只在固定内部参数的正域核验展开；外端点积分仍必须保留$r^{-\epsilon}$等幂，不能直接先展开再积。

令$\gamma_i=\gamma_j=-\epsilon,\gamma_k=1$。对齐次的一次／三次$t$单项式，$\mu$系数有紧凑投影
$$\langle t_1^{m_1}t_2^{m_2}t_3^{m_3}\rangle_\mu
 =\frac{\prod_{\ell=1}^3(\gamma_\ell)_{m_\ell}}
 {(1-2\epsilon)_{m_1+m_2+m_3}}.$$
13个投影均与实际矩约化相符。这允许先压缩动量算符，再作用环境偏导，减少几何展开。

零阶$\mu$也已从实际积分求成经典二重对数。令$h=1/(1+s)$，先积$P$得到
$$\mu_0=\int_0^h\frac{\log(1+ru)-\log r-\log u}{1-u}\,du,$$
故
$$\mu_0=[\log r-\log(1+r)+\log h]\log(1-h)
 +\operatorname{Li}_2(h)+\operatorname{Li}_2\!\left(\frac{r(1-h)}{1+r}\right)
 -\operatorname{Li}_2\!\left(\frac r{1+r}\right).$$
该式与原二维积分的HyperInt结果精确一致；混合Li₂不能删去。

## 5. 新星核换元及已求完的夸克端点子项

使用
$$a=(1-V)A,\quad c=VC,\quad b=(1-V)(1-A)+V(1-C),$$
$$x=1-X/A,\quad y=V(1-C)/b,\quad w=Y/C.$$
其完整域及正雅可比是
$$0<V<1,\quad 0<X<A<1,\quad0<Y<C<1,\qquad
 |\mathrm{Jac}|=\frac{V(1-V)}{ACb}.$$
精确得到
$$q=V(1-V)\mathcal D_z(X,Y),\qquad
 \mathcal D_z=X+zY-(1+z)XY.$$
星核与完整五维Jacobian分别精确核验；新的32组动量算符在三个扇区先约内部矩，再用带完整D迹的环境导数作用。27个原核／相位精确有理点检查通过。换积分次序仍在共同紧域调节下进行，再取相同亚纯延拓。

最终形式为有限和
$$\sum_{k,f,p}\int_0^1dX\int_0^1dY\int_0^1dV\int_X^1dA\int_Y^1dC\;
 K_{kfp}(A,C,V,X,Y,z;\epsilon)\,
 \mathcal F_f\,\mathcal D_z^{3\epsilon-p}\psi_k^{\epsilon},$$
其中$p=2,3$，$K$为有理式；$X,Y,z$的170种单项式已分出，剩下内部系数不含$z$。对$k=3$，
$$\psi_3=\frac{V(1-V)^2}{AC^2b},\qquad
 \psi_2=\frac{V^2(1-V)^2}{ACb^2}.$$
四个场可取$\mathcal F=(\epsilon^2H\mu,\epsilon(1+r)a_r,
\epsilon(1+s)a_s,c_0)$。所有源算符和非零路径体项仍包含在$K$中。

端点不能只取有限参数展开。定义
$$\mathcal I_b(\epsilon)=\int_b^1t^{-1-\epsilon}(1-t)^{-\epsilon}dt,
 \qquad B_\gamma=B(-\epsilon,1-\epsilon).$$
例如第3扇区$a\to0$时，$\mathcal F_1=e^2H\mu$、$\mathcal F_2=e(1+r)a_r$等满足
$$\mathcal F_2\sim\epsilon\left(\frac{1-b}{a}\right)^{\epsilon}\mathcal I_b,
 \quad\mathcal F_3\sim(ab)^{-\epsilon}+\epsilon B_\gamma b^{-\epsilon},$$
$$\mathcal F_1\sim\frac{\mathcal F_2+\epsilon B_\gamma(b^{-\epsilon}-1)}{1-b}.$$
保留这些完整幂后，夸克通道所有慢支路的$V^{-1}$系数精确相消；留下的端点子项全由$B_\gamma$驱动。用端点减除形状$V$或$1-V$得到完整$B(2\epsilon,2)=1/[2\epsilon(1+2\epsilon)]$质量，三扇区全部端点已实际积完至ε⁰。

记此特定减除下的子项为$I_q^{\mathrm{UV}}$，它不是完整$I_q$。例如
$$I_{q,-2}^{\mathrm{UV}}=\frac{G(-1,0,0;z)}{24}
 +\frac{\zeta_2G(-1;z)}8+\frac{\zeta_3}{24}
 -\frac{\zeta_2G(0;z)}{16}-\frac{G(0,0,0;z)}{48}.$$
其三系数完整保存于 [q-uv-all-coefficients.txt](../../../../../workflow/checks/gamma1-04/q-uv-all-coefficients.txt)，均满足crossing。两个独立扇区零阶矩还由一维正核积分独立验证6点，最大差约$6.9\times10^{-47}$；该核验不自动覆盖余项或更高系数。

### 5.1 四重角点的局部谱已积完

2026-10-06校正：旧稿把此角点的线性距离写成$\xi+z\eta$，它在一般$z$下不正确。由原$\mathcal D_z$直接换元得到$z\xi+\eta-(1+z)\xi\eta$，以下使用正确领先距离。旧角点输出保留为历史，不能复用作原核证据。

第三扇区在$A,C,X,Y\to1$时取$\xi=1-X$、$\eta=1-Y$、$U=1-A$、$W=1-C$，原域给出$0<U<\xi$、$0<W<\eta$。对数发散的候选核是
$$\frac{\mathcal D_z^{3\epsilon-3}}{b^{1+\epsilon}},\qquad
 b\sim(1-V)U+VW,\quad\mathcal D_z\sim z\xi+\eta.$$
完整慢支路系数相消，剩余为$B_\gamma\gamma_c(V,\epsilon,z)[(1-V)^{-\epsilon}-1]$。这里$\gamma_c$来自源的实际角点系数，未加新的物理反项；第二扇区同一角点的三个完整领先谱系数都为零。

选择$\rho=z\xi+\eta<\rho_0=z/(1+z)$，该减除支持域对任意$z>0$都在原立方体内。积$U,W$和径向以外的角度，令$u=(1-V)/z$、$v=V$，得到
$$\mathcal A_\epsilon(V,z)=
 -\frac{u^{-\epsilon}-v^{-\epsilon}}
 {\epsilon(1-\epsilon)(2-\epsilon)z^2(u-v)}.$$
$u=v$处是可去奇点，不能按孤立有理极点选支路；相位仍保留$[V(1-V)^2]^\epsilon$，不能把其中的$1-V$替换成角度变量$u$。径向完整质量为$\rho_0^{2\epsilon}/(2\epsilon)$。所剩一维积分已按正确几何重新积完：夸克的ε⁻¹、ε⁰与$R$的ε⁰见 [校正角点系数](../../../../../workflow/checks/gamma1-04/corner1-chamber3-corrected-coefficients.txt)。它们满足
$$J_{R,\mathrm{corner1},0}=-\frac34J_{q,\mathrm{corner1},-1}.$$
旧9个GPL系数点和3个角度点核验了旧局部模型，不能转移到校正后的原核连接。此次有2个独立精确几何检查、6个正域Euler投影积分点，以及正确径向极点与C投影的零残差连接，见 [新控制](../../../../../workflow/checks/gamma1-04/corrected-corner-numerics.json)。此处仍是规定减除下的第三扇区子项，不是完整$I_R$。

### 5.2 A/X/Y重叠角点的两支路已求到有限项

第三扇区$A,X,Y\to0$时，本节先求$p=2$中不含$X,Y$的部分。夸克的$p=3$领先谱精确为零；$R$的$p=3$领先项非零，补项见下节。$a=A(1-V)$趋零，$b\to1-CV$；夸克的两个慢支路精确相消，$R$的慢支路须保留。用截断$A<\rho_0=Cz/(1+z)$定义局部减除，这始终满足$A<1$。先把$Y$积分延长至∞，再精确减掉$Y>C$尾部，两部分都实际积分，未把扩域当作原积分。

夸克的领先Gamma系数可写为
$$K_q(\epsilon,z)=\frac{-8(3\epsilon-1)(\epsilon-1)z}
 {(1+\epsilon)(2\epsilon-3)(2\epsilon-1)},$$
$$\gamma_{q,b}=2K_q,\qquad
 \gamma_{q,1}=-K_q\left(\frac{1+\epsilon}{b}+1-\epsilon\right).$$
令$w=CV$，先积$w$，再积$V$，整个$C,V$谱闭合为$F_q=K_q\mathcal B_q$，其中
$$\begin{aligned}
\mathcal B_q={}&\frac{2[B(\epsilon,1+2\epsilon)-B(\epsilon,2)]}{1-2\epsilon}
 -B(\epsilon,1+2\epsilon)+B(\epsilon,2+\epsilon)\\
&-\frac{1+\epsilon}{\epsilon}
 [B(\epsilon,1+\epsilon)-B(\epsilon,1+2\epsilon)].
\end{aligned}$$
尾部的一维函数也已明确求成标准Gauss函数，
$$T_q(\epsilon,z)=\int_0^1t^{-1-\epsilon}
 \left[\left(1+\frac{t}{1+z}\right)^{3\epsilon}-1\right]dt
 =-\frac{{}_2F_1(-3\epsilon,-\epsilon;1-\epsilon;-1/(1+z))-1}{\epsilon}.$$
因此**该局部夸克减除子项**的全ε形式是
$$J_{q,\mathrm{A}}=\frac{B_\gamma F_q}{3\epsilon z(1-3\epsilon)}
 \left[\frac{[z/(1+z)]^{2\epsilon}}{2\epsilon}
 -[z^2(1+z)]^{\epsilon}T_q\right].$$
其四个极点至有限项已展开并保存于 [q-corner0-chamber3-coefficients.txt](../../../../../workflow/checks/gamma1-04/q-corner0-chamber3-coefficients.txt)。例如
$$J_{q,\mathrm{A},-3}=\frac49(1-\zeta_2).$$

$R$的慢支路有额外的$A^{-\epsilon}$，径向质量为$1/\epsilon$而不是$1/(2\epsilon)$。其有限谱是
$$F_s(\epsilon)=\int_0^1dC\int_0^1dV\,
 \left[\frac{V^2(1-V)}{1-CV}\right]^{\epsilon}
 \mathcal I_{1-CV}(\epsilon).$$
该谱无ε=0发散；四个所需周期已实际积分：
$$\begin{aligned}
F_{s,0}&=2-\zeta_2,\\
F_{s,1}&=3-\zeta_3-\zeta_2,\\
F_{s,2}&=16+\zeta_3-9\zeta_2-\tfrac45\zeta_2^2,\\
F_{s,3}&=35-3\zeta_3-18\zeta_2-\tfrac15\zeta_2^2
 +\zeta_2\zeta_3-3\zeta_5.
\end{aligned}$$
再加入Gamma支路及两种实际尾部积分，得到$p=2$子项$J_{R,\mathrm{A}}^{(2)}$的ε⁻³…ε⁰，[四系数](../../../../../workflow/checks/gamma1-04/r-corner0-chamber3-coefficients.txt)。例如
$$J_{R,\mathrm{A},-3}^{(2)}=\frac54(2-\zeta_2).$$
2个非整数ε的Beta谱控制、q与R两组各3个完整Gauss／不完整Beta展开预算控制通过。展开与全函数的差随ε减半约减半，支持保留至ε⁰的预算。

此截断的局部子项会引入$G(\cdot;z)$字母−2和$\log2$；它们属于辅助减除，尚须在完整余项及组装中核对消项，不能据此判定最终M/N闭包。所有结果的自变量仍为$z$。

### 5.3 R的三次核补项与范围修正

进一步作$X=A\xi,Y=A\eta,A\to0$的完整缩放检查，发现$R$的$p=3$核也贡献领先项。上节的$R$结果仅是$p=2$子项，不能把它的系数转移给完整角点。夸克的四个$p=3$领先支路确实全为零。

三次核的领先项为$(f_X X+f_Y Y)D_0^{3\epsilon-3}/A$，$D_0=X+zY$，有慢Gauss、纯幂与Gamma三个支路。两种延长域的角度矩为
$$\int_0^A dX\int_0^\infty dY\;X D_0^{3\epsilon-3}
=\frac{A^{3\epsilon}}{3\epsilon z(2-3\epsilon)},$$
$$\int_0^A dX\int_0^\infty dY\;Y D_0^{3\epsilon-3}
=\frac{A^{3\epsilon}}{3\epsilon z^2(2-3\epsilon)(1-3\epsilon)}.$$
超出$Y>C$的尾部精确减掉；除上节$T$外还需要
$$U_\delta=-\frac{{}_2F_1(1-3\epsilon,-\delta\epsilon;1-\delta\epsilon;-1/(1+z))-1}{\delta\epsilon},\quad\delta=1,2.$$
所需$C,V$周期全部实际积完，Gamma混合支路的最低系数与慢／纯幂支路逐点相消后再乘$B_\gamma$。四个独立Beta／一维积分预算控制通过，余项按ε³缩放。

现在$p=2$与$p=3$已正确合并，[完整局部四系数](../../../../../workflow/checks/gamma1-04/r-corner0-chamber3-full-coefficients.txt)。例如
$$J_{R,\mathrm{A},-3}^{\mathrm{full}}=\frac{2\zeta_2+\zeta_3-1}{12}.$$
此系数仍属于规定截断下的第3扇区角点子项，不是完整$J_R$的最高系数。

另一侧$C,X,Y\to0$的完整领先谱也已核对：夸克的$p=2,3$六个系数均为零；$R$的两种Gauss谱均非零，九个周期及局部三系数现已实际积完。[谱](../../../../../workflow/checks/gamma1-04/star-cornerC.txt)。这说明仅处理A角点仍不足以完成$R$。

### 5.4 全部V端点已积完，R有限核心已开始实际降维

$R$的第3扇区两个端点与第2扇区一个独立端点的12个所需核心均实际积完；另一个端点及第1扇区经34项完整原核／相位反射复用。全部V端点子项ε⁻²…ε⁰见 [r-uv-all-coefficients.txt](../../../../../workflow/checks/gamma1-04/r-uv-all-coefficients.txt)，三系数均通过crossing。9个最低核心原四维QMC点在记录的采样预算内支持结果；高阶核心仍以CAS检查为主，未声称QMC给出高精度证书。Gauss二阶核还与经典Li₃表达式精确相等。

全局Euler膨胀后，改在固定$A,C,V$先分出纯$X/Y\to0$的全部领先项：$q^{-2}$的常数项及$q^{-3}$的X、Y一次项都须保留。以$D_0$替代其领先距离，作精确加减；不能只减掉$q^{-2}$常数。

$R$有限点核中，内部主函数、Gauss函数只依赖$A,C,V$。把X/Y/z单项式分开后，19个V谱和10个双参数角度核已实际积完；合并后剩下两个α参数。第3扇区该有限核心也已实际积为GPL，[表达式](../../../../../workflow/checks/gamma1-04/r-bulk-third-chamber.txt)。它对应减去V端点及纯X/Y领先核后的点核；仍需补回完整ε的角点质量，不能把它单独当作$J_R$。这一核心的新独立数值检查尚未完成。

### 5.5 C角点与第二扇区有限核心

第3扇区$C/X/Y\to0$的两个Gauss谱包含$q^{-2}$常数和$q^{-3}$一次项。用$C<A/(1+z)$作局部减除，先延长X至∞、再精确扣除$X>A$尾部。九个所需周期已实际积完，ε⁻²…ε⁰保存于 [C角点三系数](../../../../../workflow/checks/gamma1-04/r-cornerC-chamber3-coefficients.txt)。这只是相应局部子项；尚未完成其新的独立数值预算检查。

第二扇区有限点核心也已按同一路线独立积为GPL，[表达式](../../../../../workflow/checks/gamma1-04/r-bulk-second-chamber.txt)。两个核心最终自变量均为z，字母0、−1；第二扇区的精确反射检查通过。第1扇区可按原核／原域完整反射处理。新有限核心的独立数值比较仍待完成。

## 6. 当前实际剩余任务


完整关系仍是
$$I_g=-\frac{1-\epsilon}{2}I_q+I_p+I_R.$$
已求完$I_p$与$I_q^{\mathrm{UV}}$，并实际求完第3扇区的上述两类局部$q/R$角点子项。R的完整V端点子项已求完。C/X/Y局部子项与两个R有限点核心现也已积完。还需完整ε的X/Y领先核质量、第二扇区相应局部质量、夸克有限核心及全部补偿连接，以及最后两颜色组装与完整交叉核验。不能直接把局部子项相加视作原$I$：定义中的精确补偿余项必须一起求完。

全局Euler虽能生成可积候选核，但会把短式放大为百万项。旧正参数路线的第一项实际积分超时；新星核第三扇区正规化与原表示中参与ε⁰的40个有理点系数相符。新路线已依次实际积分$V,A,C,X,Y$，第三扇区ε⁻⁵候选系数精确为零。这只是一个扇区的一个位置，不能据此宣布整图最高极点、其余系数或有限项完成。

当前优先求完整ε的X/Y领先核质量及第二扇区相应局部质量，再完成夸克有限核心和最终组装／独立检查，保留纯$X/Y\to0$与重叠角点的真实幂；全局Euler的ε⁻⁴等位置只作对照，避免继续放大候选核。第二扇区与$I_R$也须独立完成。截止此检查点，Γ₁整体有限阶仍3/4、全ε0/4、数学缺项11，formal=not_run。生产／冻结源、旧结果、Kernel/Lean均未改，未提交Git。


### 下一步的确切积分对象

设从完整核取出的领先系数为$H_0,H_X,H_Y$：
$$H_0=\sum_f K_{R,3,f,2}(X=Y=0)\mathcal F_f,\quad
H_X=\sum_f\partial_XK_{R,3,f,3}(0,0)\mathcal F_f,\quad
H_Y=\sum_f\partial_YK_{R,3,f,3}(0,0)\mathcal F_f.$$
它们依赖$A,C,V,\epsilon$，包含完整三角／Gauss谱，不能只用固定参数的零阶值补回质量。纯X/Y领先核在原矩形上的完整积分分别为
$$L_2=\int_0^A dX\int_0^C dY\,(X+zY)^{3\epsilon-2},$$
$$L_X=\int_0^A dX\int_0^C dY\,X(X+zY)^{3\epsilon-3},\quad
L_Y=\int_0^A dX\int_0^C dY\,Y(X+zY)^{3\epsilon-3}.$$
其极点是$1/(3\epsilon z),1/(6\epsilon z),1/(6\epsilon z^2)$；有限项令
$\ell=\log A+\log(zC)-\log(A+zC)$，分别为
$$L_{2,0}=(1+\ell)/z,$$
$$L_{X,0}=[\ell+A/(A+zC)+1/2]/(2z),\quad
L_{Y,0}=[\ell+zC/(A+zC)+1/2]/(2z^2).$$
下一步须对$\psi^\epsilon(H_0L_2+H_XL_X+H_YL_Y)$处理A/C重叠角点，求完整极点和有限项，并与已积的有限点核心、V端点及四重角点质量衔接。未经这一补回，现有GPL核心不是原$I_R$。q通道也须完成对应ε预算和有限核心。

## 7. 替代方法评估后的角度主积分组织

2026-10-05按用户指定的参考工具及GitHub进行[方法评估](../../../../../workflow/checks/gamma1-04-methods/README.md)。本轮只实际核验上节的三个角度核，不计新的完整颜色积分。

令$\tau=3\epsilon$、$B=zC$、$S=A+B$。在$A,C,z>0$且$\Re\tau>0$的初始收敛域，原矩形求积给出
$$L_2=\frac{S^\tau-A^\tau-B^\tau}{z\tau(\tau-1)},$$
$$L_X=\frac{A[S^{\tau-1}-A^{\tau-1}]/[z(\tau-1)]-L_2}{\tau-2},$$
$$L_Y=\frac{C[S^{\tau-1}-B^{\tau-1}]/(\tau-1)-L_2}{z(\tau-2)}.$$
随后亚纯延拓至ε=0；τ=1、2的表观奇点按完整组合取极限。两个后一式分别来自X、Y全导数的矩形边界，不删去任何上边界；并有$L_X+zL_Y=L_2$。

因此可先将$H_0L_2+H_XL_X+H_YL_Y$合并成一个$L_2$和两条边界幂支路，再与完整内部谱一起积分。是否出现额外相消须在实际$H$上检查，目前尚未求其外层积分。不能取$H(\epsilon)$的零阶替代完整预算。

2026-10-06实际接续：对q/R、两个独立扇区、四个内部场完成16组完整ε合并检查。q通道的领先三次核为零，因此只需同一个$L_2$幂差；R仍保留三条幂支路。表达式见[实际合并系数](../../../../../workflow/checks/gamma1-04/leading-angular-merged.txt)。下面先求最简单的$q$第四场，不把它当作完整颜色积分。

独立检查把$X=Au,Y=Cv$的单位方形分成$u\ge v$、$v\ge u$，置$(u,v)=(\rho,\rho\eta)$或$(\rho\eta,\rho)$。领先核径向积分为$1/\tau$，角度数值积分与闭式在12点一致；另有15项精确混合导数、边界、齐次及三核关系检查，见[记录](../../../../../workflow/checks/gamma1-04-methods/angular-probe.json)。这些只覆盖角度核，formal=not_run；三个核的完整外层质量、q有限核心和最终两颜色I仍未完成。

## 8. 2026-10-06：实际求完的第四场与夸克第三扇区体块

### 8.1 第四场的完整领先质量

这里的第四场是内部四函数分解中的$\mathcal F_4=c_0$，不是一个额外颜色通道。对第3扇区，$c_0=(ab)^{-\epsilon}$，与$\psi_3^\epsilon$相乘后为
$$\left[\frac{V(1-V)}{A^2C^2b^2}\right]^\epsilon,\qquad b=1-A(1-V)-CV.$$
令$R_4=-K_{q,3,4,2}(X=Y=0)/[z(3\epsilon)(3\epsilon-1)]$。实际源中该场的领先三次核为零，故原矩形角度积分后得到
$$J^{(3)}_{q,4,\mathrm{lead}}=\int_0^1dA\,dC\,dV\;
R_4\left[\frac{V(1-V)}{A^2C^2b^2}\right]^\epsilon
\{A^{3\epsilon}+(zC)^{3\epsilon}-(A+zC)^{3\epsilon}\}.$$
不能把三个幂支路各自展开后丢掉其相消。保持完整组合，固定$C,V>0$时的$A$端点是
$$A^{\epsilon-1}r_4(\epsilon,C,V)
\left[\frac{V(1-V)}{C^2(1-CV)^2}\right]^\epsilon,$$
$$r_4=\frac{8(\epsilon-1)CV}{3\epsilon^2(2\epsilon-1)(2\epsilon-3)(1-CV)}.$$
先精确减去这个完整ε端点，并补回$\int_0^1A^{\epsilon-1}dA=1/\epsilon$的亚纯延拓，余项才展开。四个补回周期和三个正规化三维系数核已全部实际解析积完，约10.71秒；结果见[四系数](../../../../../workflow/checks/gamma1-04/q-leading-field4-coefficients.txt)。例如
$$J^{(3)}_{q,4,\mathrm{lead},-3}=\frac89(1-\zeta_2),\qquad
J^{(3)}_{q,4,\mathrm{lead},-2}=\frac{148}{27}-\frac{88}{27}\zeta_2-\frac83\zeta_3.$$
第一项还直接来自$\int_0^1\!dC\,dV\,CV/(1-CV)=\zeta_2-1$。单极点与有限项是自变量$z$、字母0/−1的显式GPL，无未求常数。4个端点周期及三个z点的9个原域正规化系数积分独立QMC检查通过，见[13项数值](../../../../../workflow/checks/gamma1-04/q-field4-numerics.json)。这些是裸J中的一个场分量，尚未施加唯一I前因子，也不证明整个q通道完成。

### 8.2 夸克体块的单极点与有限项

先从全q谱扣除已求V端点，检查发现其点态单极点不为零。有限点核必须写成
$$Q_p^{(0)}(A,C,V,X,Y,z)+3P_p(A,C,V,X,Y,z)\log\mathcal D_z,$$
其中$P_p$是点态ε⁻¹系数；$P_3=0$，$P_2\ne0$。6个内部未约化Euler谱的正参数数值控制支持这两个系数及其预算，见[源连接](../../../../../workflow/checks/gamma1-04/q-point-connection-numerics.json)。

12个V谱已实际求完。4个无对数角度核按完全相同定义复用，另3个$\log\mathcal D_z$角度核实际求完；对X/Y领先项同时减去$\log D_0/D_0^p$，不能只减去无对数核。最后两个α积分中的单极点点核已积完，有限点核在$C\to1$暴露出真实对数发散。

核对原几何发现§5.1的旧局部模型错用了$\xi+z\eta$。修正为$z\xi+\eta$后，实际C发散系数与校正角点极点严格相接。HyperInt在无穷远使用倒变量展开，报告的$\log(1/cc)$系数与物理$\log cc$差一个负号；用$1/(1+cc)$单独校准后，连接残差为零，减除后的有限核心在发散检查保持开启时完成。

为使补回质量与C减除同域，固定$W=1-C$并积分其余角点变量。由$0<U<\xi$、$W<\eta$及正确距离$z\xi+\eta$，依次积$\xi,\eta$后得到
$$D_q(\epsilon,z)=\int_0^1dV\;
\frac{[V(1-V)^2]^\epsilon B_\gamma\gamma_q(V,\epsilon,z)[(1-V)^{-\epsilon}-1]}
{z(1-3\epsilon)(2-3\epsilon)}
\int_0^\infty dt\;
\frac{(zt+1)^{3\epsilon-1}}{[(1-V)t+V]^{1+\epsilon}}.$$
相应C端点为$W^{2\epsilon-1}D_q$，完整补回质量为$D_q/(2\epsilon)$。$D_{q,0},D_{q,1}$已实际积完，零阶同时匹配正确径向角点与体块端点；两个系数在$z=1/2,1,2$的6个独立正域Euler数值检查通过。

这得到第3扇区扣除V端点及完整X/Y领先质量后的体块，含正确C投影补回，其ε⁻¹与ε⁰保存于[新体块](../../../../../workflow/checks/gamma1-04/q-bulk-third-chamber-with-C-corner.txt)。这里采用C投影减除；旧$\rho<\rho_0$角点有限项不能再额外加一次。完整新体块的独立原域数值比较尚待完成。

### 8.3 当时的剩余对象（现由§9完成）

仍需其余完整ε领先质量（内部场1/2/3及其V端点衔接）、第二扇区对应夸克体块与局部质量、R领先质量和校正角点的完整补偿连接，之后才组装两颜色唯一I并做独立检查。R已经求出的两个有限体核继续复用。Γ₁整族仍3/4、全ε0/4、formal=not_run；本轮25个数值控制覆盖上述子对象，不覆盖完整04。

## 9. 完整解析解与完成性审计

### 9.1 从同一个原核作精确加减

原始动量源到C/J的连接仍采用已审查的[完整来源链](../../momentum-to-parametric-review.md)。本次没有更改Γ₁表达式、动量规则、颜色、D=4−2ε或源权重。需要解决的是既有七维参数根，而不是重新选定另一个积分。

在两个独立Max扇区，把完整星核写成
$$f_c=psi_k^epsilonsum_{p=2,3}K_{cfp}mathcal F_fmathcal D_z^{3epsilon-p}.$$
先精确加减两个已求V端点；再加减原矩形上的全部X/Y领先核
$$ell_c=psi_k^epsilon{H_0D_0^{3epsilon-2}+(H_X X+H_Y Y)D_0^{3epsilon-3}},quad D_0=X+zY.$$
源中三次核的常数项全部为零；V减除的二次常数／三次一次领先项也逐通道为零。因此没有额外的V／平坦领先重叠补偿遗漏，且
$$J_c=J_c^{
m UV}+J_c^{
m lead}+J_c^{
m body}$$
保持同一个原域、同一个完整ε核。106个源级零项检查见[完成性审计](../../../../../workflow/checks/gamma1-04/check_complete_bookkeeping.log)。完整body含正确C投影角点质量；不再额外加入旧ρ截断版本。

### 9.2 完整内圈依赖和必要展开阶数

令Max腿为δ，另两参数为a、β，a+β+δ=1。四个实际场可以直接写成
$$mathcal F_1=epsilon^2mu_{
m phys}(a,eta;epsilon),quad
mathcal F_2=epsilon[(1-eta)/a]^epsilonmathcal I_eta,quad
mathcal F_3=epsilon[(1-a)/eta]^epsilonmathcal I_a,quad
mathcal F_4=(aeta)^{-epsilon},$$
$$mu_{
m phys}=int_0^infty du,dv;rac{(1+u+v)^{-1+2epsilon}}{(a+u)^{1+epsilon}(eta+v)^{1+epsilon}}.$$
第三扇区(a,β,δ)=(A(1−V),b,CV)；第二扇区为(A(1−V),CV,b)。内部比值r,s只是μ的辅助变量，最终结果没有Hlog(r)。

合并角度式后，q点态领先核最低ε⁻²，R最低ε⁻¹；所以q的有限项需要μ₀，R的有限项还需要μ₁。已有μ₁通过正域Euler积分独立核验，不能用μ₀替代。A端点的高阶质量另需要完整不完全Beta谱至ε³。12个μ两阶检查和24个未合并角度／Euler原式检查均通过。

### 9.3 按自然端点方向补回完整质量

固定其他参数，A→0的完整谱分成
$$A^{epsilon-1}D_{A,S}(epsilon,C,V)+A^{2epsilon-1}D_{A,G}(epsilon,C,V).$$
两个谱先分别保留完整幂，再精确补回D₍A,S₎/ε和D₍A,G₎/(2ε)。q的慢支路相消；R的慢与Gamma支路虽然在点态双极点处相消，但两种径向质量的分母不同，不能先把它们删掉。这正是原先缺失的重叠端点贡献。

第三扇区C→0时，内圈的齐次度为−3，Euler分部积分给出F₁=εF₂/a+εF₃/β；剩下两条不完全Beta支路，C径向质量为z^(3ε)/ε。第二扇区的C谱由完整场置换和V反射得到，仍分别保持ε、2ε指数。每个补回周期都实际积分，没有以未知常数命名代替求值。

第三扇区和第二扇区的全部q领先质量分别见[q3](../../../../../workflow/checks/gamma1-04/q-leading-all-third-coefficients.txt)、[q2](../../../../../workflow/checks/gamma1-04/q-leading-all-second-coefficients.txt)；R对应[r3](../../../../../workflow/checks/gamma1-04/r-leading-all-third-coefficients.txt)、[r2](../../../../../workflow/checks/gamma1-04/r-leading-all-second-coefficients.txt)。这些已包含原先单独求过的第四场，组装时没有再加第四场。

剩余有限点核实际求出GPL。独立数值从积分前的正域Euler、端点谱及三维核重算全部58个领先／补回位置；另从原五维星核重算18个body位置，均在明确的采样与求积加深预算内一致。预算不是严格区间界。body的q四重角点使用正确zξ+η和盒域减除减去全C投影；R在该角点的点态项为零，但完整ε谱产生有限补回，同样只计一次。

### 9.4 三个Max扇区、两颜色和唯一归一化

第1扇区按完整源、场置换、相位和原域反射复用第3扇区；第2扇区已经完整，不能再重复其反射。精确重建全部J后采用
$$I_c=rac{e^{-3gamma_Eepsilon}Gamma(1-3epsilon)}{32(1+z^{3epsilon})}J_c,qquad
I_g=-rac{1-epsilon}{2}I_q+I_p+I_R.$$
原剑鱼已经包含在完整Iₚ中。有限项包含+I₍q,−1₎/2，来源是分子中的ε，而不是另一个图。C₄重数4、内圈和Euler补偿各计一次；Gamma前因子展开到ε³，足以覆盖最高三阶。原始J与最终I逐项精确相接，12个最终／R反射检查通过；省略ε乘q或Gamma的ε²项都会留下非零差。

全核的V质量至多二阶，完整A/C领先质量至多三阶，body点态至多单阶，四重角点至多单阶／有限；其余正规化核在ε=0邻域可展开。因此最高三阶不是仅根据未出现更高系数推测，ε⁻⁶…ε⁻⁴确为零。

### 9.5 最终结果

两个颜色为C_F C_A²和C_F C_A T_F n_f。各自的全部系数都已求出：
$$I_c=rac{I_{c,-3}}{epsilon^3}+rac{I_{c,-2}(z)}{epsilon^2}+rac{I_{c,-1}(z)}{epsilon}+I_{c,0}(z)+O(epsilon).$$
$$I_{g,-3}=rac{zeta_2}{576}+rac{zeta_3}{128},qquad I_{q,-3}=-rac{zeta_2}{72}.$$
所有Hlog的自变量为z、字母仅0/−1。完整I/J、源哈希、颜色与核验证据见[不可变r001](../../results/3e-Gamma1-04/r001.json)。下表只给便于检查的z=1有限项，不替代一般z的解析式。

| 通道 | z=1的有限项 |
|---|---:|
| g | -4.62965495644975738314864127854954781552998790268523426421019 |
| q | 10.5826987706101835219583831807390747335633807919483217064009 |

下面逐项列出正式结果的精确Maple表达式；zeta[k]是Riemannζ(k)，Hlog(z,[...])是标准G函数。没有未求积分或待定常数。

g 通道：

```text
I[-3] = 1/576*zeta[2]+1/128*zeta[3]
I[-2] = 1/128*Hlog(z,[0, -1, 0, 0])+(-1/64*zeta[2]-1/64*zeta[3])*Hlog(z,[-1])-1/64*Hlog(z,[-1, -1, 0, 0])+3/128*zeta[2]*Hlog(z,[-1, 0])-3/256*zeta[2]*Hlog(z,[0, 0])+(1/128*zeta[3]+1/128*zeta[2])*Hlog(z,[0])+1/384*Hlog(z,[0, 0, 0])-1/192*Hlog(z,[-1, 0, 0])+1/128*Hlog(z,[-1, 0, 0, 0])+3/128*zeta[2]*Hlog(z,[0, -1])+5/432*zeta[2]-1/72*zeta[3]-41/2560*zeta[2]^2-3/64*zeta[2]*Hlog(z,[-1, -1])-1/256*Hlog(z,[0, 0, 0, 0])
I[-1] = 83/2592*zeta[2]-1/15360*(171*z*zeta[2]^2-640*z*zeta[2]+960*z*zeta[3]+171*zeta[2]^2-480*zeta[2]+960*zeta[3])/(1+z)*Hlog(z,[0])-1/4608*(45*z*zeta[2]-80*z+45*zeta[2]-32)/(1+z)*Hlog(z,[0, 0, 0])+287/768*zeta[2]^2-83/216*zeta[3]+105/256*zeta[5]+3/256*Hlog(z,[-1, 0, 0, 0, 0])-1/16*Hlog(z,[-1, 0, 0, 0])+7/128*Hlog(z,[0, 0, 0, 0])+1/64*Hlog(z,[0, -1, 0, 0, 0])-1/32*Hlog(z,[0, -1, -1, 0, 0])+1/8*Hlog(z,[-1, -1, 0, 0])-1/16*Hlog(z,[0, -1, 0, 0])-1/32*Hlog(z,[-1, 0, -1, 0, 0])+1/16*Hlog(z,[-1, -1, -1, 0, 0])-1/32*Hlog(z,[-1, -1, 0, 0, 0])-3/512*Hlog(z,[0, 0, 0, 0, 0])+1/64*Hlog(z,[0, 0, -1, 0, 0])-1/32*zeta[2]*zeta[3]+(1/8*zeta[3]+57/2560*zeta[2]^2-7/96*zeta[2])*Hlog(z,[-1])+(-3/16*zeta[2]-1/32*zeta[3])*Hlog(z,[-1, 0])+(7/64*zeta[2]-5/256*zeta[3])*Hlog(z,[0, 0])+(5/256*zeta[2]-7/288)*Hlog(z,[-1, 0, 0])+(-3/16*zeta[2]-1/32*zeta[3])*Hlog(z,[0, -1])+(3/8*zeta[2]+1/16*zeta[3])*Hlog(z,[-1, -1])+3/64*zeta[2]*Hlog(z,[0, 0, -1])+3/16*zeta[2]*Hlog(z,[-1, -1, -1])+3/64*zeta[2]*Hlog(z,[0, -1, 0])-3/32*zeta[2]*Hlog(z,[-1, 0, -1])-3/32*zeta[2]*Hlog(z,[0, -1, -1])-3/32*zeta[2]*Hlog(z,[-1, -1, 0])
I[0] = 1/1244160*(422739*z*zeta[2]^3+938088*z*zeta[2]^2-1562760*z*zeta[2]*zeta[3]+571050*z*zeta[3]^2+422739*zeta[2]^3+56800*z*zeta[2]-1880400*z*zeta[3]-651240*z*zeta[5]+541512*zeta[2]^2-1562760*zeta[2]*zeta[3]+571050*zeta[3]^2+56800*zeta[2]-1880400*zeta[3]-651240*zeta[5])/(1+z)+1/13824*(6804*z*zeta[2]-675*z*zeta[3]+376*z+6804*zeta[2]-675*zeta[3]+208)/(1+z)*Hlog(z,[0, 0, 0])+1/15360*(717*z*zeta[2]^2-12960*z*zeta[2]+5480*z*zeta[3]+717*zeta[2]^2-7200*zeta[2]+5480*zeta[3])/(1+z)*Hlog(z,[0, -1])+1/1536*(333*z*zeta[2]+280*z+333*zeta[2]+184)/(1+z)*Hlog(z,[0, 0, 0, 0])+1/4608*(3567*z*zeta[2]^2+369*z*zeta[2]*zeta[3]+320*z*zeta[2]-432*z*zeta[3]+261*z*zeta[5]+3567*zeta[2]^2+369*zeta[2]*zeta[3]+264*zeta[2]-1584*zeta[3]+261*zeta[5])/(1+z)*Hlog(z,[0])-1/512*(85*z*zeta[2]+144*z+85*zeta[2]+80)/(1+z)*Hlog(z,[0, -1, 0, 0])+1/7680*(513*z*zeta[2]^2+3280*z*zeta[2]-1360*z*zeta[3]+513*zeta[2]^2+1840*zeta[2]-1360*zeta[3])/(1+z)*Hlog(z,[0, 0])-103/512*Hlog(z,[0, -1, 0, 0, 0, 0])+3/128*Hlog(z,[-1, 0, -1, 0, 0, 0])+51/256*Hlog(z,[-1, -1, 0, 0, 0, 0])+21/128*Hlog(z,[-1, -1, 0, -1, 0, 0])-29/256*Hlog(z,[0, 0, -1, 0, 0, 0])-19/256*Hlog(z,[0, 0, 0, -1, 0, 0])-95/512*Hlog(z,[-1, 0, 0, 0, 0, 0])+21/128*Hlog(z,[-1, -1, -1, 0, 0, 0])+21/128*Hlog(z,[0, -1, -1, -1, 0, 0])-21/64*Hlog(z,[-1, -1, -1, -1, 0, 0])+41/256*Hlog(z,[0, 0, 0, 0, 0, 0])-1/64*Hlog(z,[-1, 0, 0, -1, 0, 0])+5/128*Hlog(z,[0, -1, -1, 0, 0, 0])+21/128*Hlog(z,[-1, 0, -1, -1, 0, 0])-125/384*Hlog(z,[-1, 0, 0, 0, 0])-7/32*Hlog(z,[0, -1, 0, 0, 0])+137/384*Hlog(z,[0, -1, -1, 0, 0])+137/384*Hlog(z,[-1, 0, -1, 0, 0])-137/192*Hlog(z,[-1, -1, -1, 0, 0])+137/384*Hlog(z,[-1, -1, 0, 0, 0])+191/768*Hlog(z,[0, 0, 0, 0, 0])-119/384*Hlog(z,[0, 0, -1, 0, 0])+(137/128*zeta[2]+21/128*zeta[3])*Hlog(z,[0, -1, -1])+(-783/2560*zeta[2]^2+21/16*zeta[2]-137/192*zeta[3])*Hlog(z,[-1, -1])+(137/128*zeta[2]+21/128*zeta[3])*Hlog(z,[-1, 0, -1])+(-137/64*zeta[2]-21/64*zeta[3])*Hlog(z,[-1, -1, -1])+(-211/384*zeta[2]+15/256*zeta[3]-73/1728)*Hlog(z,[-1, 0, 0])+(-69/512*zeta[2]-7/32)*Hlog(z,[-1, 0, 0, 0])+(-9/256*zeta[2]+7/16)*Hlog(z,[-1, -1, 0, 0])+(75/1024*zeta[2]^2-21/32*zeta[2]+137/384*zeta[3])*Hlog(z,[-1, 0])+(137/128*zeta[2]+21/128*zeta[3])*Hlog(z,[-1, -1, 0])+(-221/384*zeta[2]-1/8*zeta[3])*Hlog(z,[0, -1, 0])+(-1211/1280*zeta[2]^2-9/256*zeta[2]*zeta[3]-73/576*zeta[2]+7/16*zeta[3]+3/256*zeta[5])*Hlog(z,[-1])+63/128*zeta[2]*Hlog(z,[-1, 0, -1, -1])-1/8*zeta[2]*Hlog(z,[0, -1, -1, 0])+63/128*zeta[2]*Hlog(z,[-1, -1, 0, -1])-63/64*zeta[2]*Hlog(z,[-1, -1, -1, -1])+63/128*zeta[2]*Hlog(z,[0, -1, -1, -1])+63/128*zeta[2]*Hlog(z,[-1, -1, -1, 0])-57/256*zeta[2]*Hlog(z,[0, 0, 0, -1])-29/256*zeta[2]*Hlog(z,[0, 0, -1, 0])-9/64*zeta[2]*Hlog(z,[-1, 0, -1, 0])-3/64*zeta[2]*Hlog(z,[-1, 0, 0, -1])-119/128*zeta[2]*Hlog(z,[0, 0, -1])
```

q 通道：

```text
I[-3] = -1/72*zeta[2]
I[-2] = -1/27*zeta[2]+1/24*Hlog(z,[-1, 0, 0])+5/72*zeta[3]-1/48*Hlog(z,[0, 0, 0])-1/16*Hlog(z,[0])*zeta[2]+1/8*Hlog(z,[-1])*zeta[2]
I[-1] = -1/144*(8*z+5)/(1+z)*Hlog(z,[0, 0, 0])-1/4*zeta[2]*Hlog(z,[0, 0])+3/16*Hlog(z,[0, -1, 0, 0])+13/144*Hlog(z,[-1, 0, 0])-3/8*Hlog(z,[-1, -1, 0, 0])-1/8*Hlog(z,[0, 0, 0, 0])-9/8*zeta[2]*Hlog(z,[-1, -1])+(13/48*zeta[2]-3/8*zeta[3])*Hlog(z,[-1])+9/16*zeta[2]*Hlog(z,[-1, 0])-1/48*(7*z*zeta[2]-9*z*zeta[3]+6*zeta[2]-9*zeta[3])/(1+z)*Hlog(z,[0])-11/15*zeta[2]^2-13/162*zeta[2]+83/432*zeta[3]+9/16*zeta[2]*Hlog(z,[0, -1])+3/16*Hlog(z,[-1, 0, 0, 0])
I[0] = -53/486*zeta[2]-1/216*(108*z*zeta[2]+17*z+108*zeta[2]+14)/(1+z)*Hlog(z,[0, 0, 0])-1/288*(420*z*zeta[2]^2+64*z*zeta[2]-99*z*zeta[3]+420*zeta[2]^2+60*zeta[2]-99*zeta[3])/(1+z)*Hlog(z,[0])-365/288*zeta[2]^2+65/648*zeta[3]+19/6*zeta[5]+1/6*Hlog(z,[-1, 0, 0, 0, 0])+11/32*Hlog(z,[-1, 0, 0, 0])-5/24*Hlog(z,[0, 0, 0, 0])+3/16*Hlog(z,[0, -1, 0, 0, 0])-41/48*Hlog(z,[0, -1, -1, 0, 0])-11/16*Hlog(z,[-1, -1, 0, 0])+11/32*Hlog(z,[0, -1, 0, 0])-41/48*Hlog(z,[-1, 0, -1, 0, 0])+41/24*Hlog(z,[-1, -1, -1, 0, 0])-41/48*Hlog(z,[-1, -1, 0, 0, 0])-1/48*Hlog(z,[0, 0, 0, 0, 0])+13/24*Hlog(z,[0, 0, -1, 0, 0])+307/144*zeta[2]*zeta[3]+13/8*zeta[2]*Hlog(z,[0, 0, -1])+41/8*zeta[2]*Hlog(z,[-1, -1, -1])+25/24*zeta[2]*Hlog(z,[0, -1, 0])-41/16*zeta[2]*Hlog(z,[-1, 0, -1])-41/16*zeta[2]*Hlog(z,[0, -1, -1])-41/16*zeta[2]*Hlog(z,[-1, -1, 0])+(-5/12*zeta[2]+5/12*zeta[3])*Hlog(z,[0, 0])+(43/48*zeta[2]+31/216)*Hlog(z,[-1, 0, 0])+(-33/16*zeta[2]+41/24*zeta[3])*Hlog(z,[-1, -1])+(77/40*zeta[2]^2+31/72*zeta[2]-11/16*zeta[3])*Hlog(z,[-1])+(33/32*zeta[2]-41/48*zeta[3])*Hlog(z,[-1, 0])+(33/32*zeta[2]-41/48*zeta[3])*Hlog(z,[0, -1])
```

### 9.6 完成范围

Γ₁04两个颜色的全部极点至有限项完成，整族有限阶4/4。全ε标准函数闭式仍0/4；物理复区间延拓、严格源侧全域去调节证明与Wilson形式收据没有提升。全三圈w3的其他缺图、颜色组装、反项与RG另有任务；本轮不据此宣布完整w3完成。没有Git提交。
