# T4[1]：有序对边核、完整边界与求解

本批是5e-3g-27/28/29（13通道）和41/42（24通道）。沿用正实z、唯一C/I与一次Euler补偿，e=ε、λ=1−2e、$G_3=e^{-3\gamma_Ee}\Gamma(1-e)\Gamma(1-2e)$。原输入、排序和颜色未改；formal=not_run。

[解表](solutions.md) · [新增函数及原有基底](basis-functions.md) · [计算证据](../../../../workflow/checks/tgv-t4-opposite/README.md)

## 1. 原核与两个通道

在与二圈3e-3g-04一致的D框架中，记

$$D_z(p,q)=(1-p)(1-q)+zpq,$$
$$S=\alpha_1\alpha_2y(1-x)+\alpha_1\alpha_3D_z(x,w)+z\alpha_2\alpha_3w(1-y),$$
$$P=\frac14\{\alpha_1\alpha_2[-(1+z)(1-x)+zy]+\alpha_1\alpha_3[1+3z-(1+z)x+z(1+z)w]+z\alpha_2\alpha_3[1-y-(1+z)w]\}.$$

α测度是 $d\mu_e=d^3\alpha\,\delta(1-\sum\alpha_i)(\alpha_1\alpha_2\alpha_3)^{-e}$。

27/28/29的原域分别为

$$\mathcal D_{27}:x<p,\ w<q;\qquad \mathcal D_{28}:x<p,\ q<w;\qquad \mathcal D_{29}:p<x,\ q<w,$$

所有Wilson变量都在[0,1]，y独立。定义

$$\widehat J_g=-(1+z)\int_{\mathcal D_g}dx\,dy\,dw\,dp\,dq\int d\mu_e\,P S^{-2+2e}D_z(p,q)^{-1+e}.$$

当前prepared用F框架，逐原多项式有 $S_F(z)=zS_D(1/z)$、$P_F(z)=z^2P_D(1/z)$；连同pure分母和Lorentz因子，得到

$$J_g^{\rm prepared}(z,e)=z^{3e}\widehat J_g(1/z,e).$$

该替换不改变crossing组合，故

$$I_g=\frac{\lambda G_3}{8}\frac{\widehat J_g(z,e)+z^{3e}\widehat J_g(1/z,e)}{1+z^{3e}}.$$

反射 $x\mapsto1-w,y\mapsto1-y,w\mapsto1-x$、α1↔α3、z↦1/z，连同pure端点与排序反转，严格给出I27=I29；28的J本身满足self-crossing。两组结论均已逐原核、边界、Jacobian和C核验。

## 2. 41/42：先在J层因子化

41/42仅在边2有TGV点与pure点的两种互补排序。初始绝对收敛域内，两域合成完整方形：

$$J_{41}+J_{42}=J_{\rm T4}\,J_{\rm pure},\qquad J_{42}(z)=z^{3e}J_{41}(1/z).$$

$J_{\rm T4}$的S/P已与实际二圈prepared输入精确比较，$J_{\rm pure}=z^{e/2}H_0(e,z)$。设 $G_2=e^{-2\gamma_Ee}\Gamma(1-2e)$，则

$$\boxed{I_{41}=I_{42}=e^{-\gamma_Ee}\Gamma(1-e)\frac{\cosh(e\log z)}{8\cosh(3e\log z/2)}I_{\rm T4}(z,e)H_0(e,z).}$$

系数1/8来自三圈唯一I定义及互补域；这里不能直接相乘两个物理I。所用完整二圈T4级数见[原求解](../../3e-3g/method.md)。单对边全函数为

$$H_0(\delta,z)=-\frac{(1+z)z^{-\delta/2}}{\Gamma(1-\delta)\Gamma(1+\delta)^2}G_{3,3}^{3,3}\!\left(z\middle|\begin{matrix}\delta,0,0\\0,\delta,\delta\end{matrix}\right).$$

有限项只需既有二圈系数和 $H_0(e)=M_{0,1}+eM_{1,1}/2+e^2(M_{2,1}/8-M_{0,1}^2/24)+O(e^3)$。Gamma比值在e¹没有项。旧二圈结果被实际哈希和Hlog回代检查，不仅按文件名复用。

## 3. 保留两端点的T4母核

先只积分α及y，保留x,w。令

$$H_e(x,w)=\int d\mu_e\int_0^1dy\,S^{-\lambda},\qquad
\mathcal R_e(x,w)=\int d\mu_e\int_0^1dy\,R S^{-2+2e},$$
$$R=\frac{1+z}{2}\alpha_3\{\alpha_1(1-x+zw)-z\alpha_2w\}.$$

精确恒等式为

$$P=\frac{-zS_x-(1+z)S_y+S_w}{4}+R.$$

y端点的α尺度恰给1/e，其余一维函数为

$$h_j(x,w;e)=\int_0^1db\,[b(1-b)]^{-e}Q_j^{-1+2e},$$
$$Q_0=(1-b)(1-x)(1-w)+z\{b+(1-b)x\}w,\qquad
Q_1=(1-x)\{1-(1-b)w\}+zx(1-b)w.$$

引入pure矩形权重（不含−(1+z)）

$$W_{27}(x,w;e)=\int_x^1dp\int_w^1dq\,D_z(p,q)^{-1+e},\qquad
W_{28}(x,w;e)=\int_x^1dp\int_0^wdq\,D_z(p,q)^{-1+e}.$$

它们的e⁰/e¹已经直接积分成Hlog。H₀、R的极点/有限项、h的e⁰/e¹也由新HyperInt运行得到。因此剩下的是两端点的积分，而不是将已完全积分的二圈I再次乘权重。

## 4. 展开前保留全局边界

把x,w导数先转移到W，得到

$$\begin{aligned}
\frac{4\lambda\widehat J_g}{-(1+z)}={}&z\int_0^1dw\,[H_eW_g]_{x=0}^{1}-\int_0^1dx\,[H_eW_g]_{w=0}^{1}\\
&+\int dx\,dw\,H_e(-z\partial_xW_g+\partial_wW_g)\\
&+\frac{1+z}{e}\int dx\,dw\,(h_1-h_0)W_g+4\lambda\int dx\,dw\,\mathcal R_eW_g.
\end{aligned}$$

27的w=0面有cusp。若先在固定x,w展开并对外层直接积分，有限项发生对数发散，且遗漏的边界也改变1/e系数。该尝试已保存为失败路线，没有使用其V27m1作为最终极点。

正确保留的cusp项是

$$B_{27}^{\rm cusp}(e,z)=\frac{F_0(e)}{2e^2}\mathcal B(e,z),$$
$$F_0(e)=\int_0^1db\,du\,[b(1-b)]^{-e}(1-bu)^{-1+2e}=\zeta_2+e\zeta_3+O(e^2),$$
$$\mathcal B(e,z)=\int_0^1dp\,dq\,D_z(p,q)^{-1+e}\{1-(1-p)^{2e}\}=e b_1+e^2b_2+O(e^3).$$

其中，设 $L_0=\log(1-p),L_1=\log z+\log p$，

$$b_1=\int_0^1\frac{-2L_0(L_1-L_0)}{(1+z)p-1}\,dp,$$
$$b_2=\int_0^1\frac{-L_0(L_1^2-L_0^2)-2L_0^2(L_1-L_0)}{(1+z)p-1}\,dp.$$

内点分母的零点可去；两个系数已求为Hlog。cusp最高只有简单极点。x=0面用H₀的普通端点极限；28的w=1面由原核反射得到，并独立计算W的完整q区间。

R使用与裸T4同一联合区域 {r,1/u} 的Euler正规化，缩放次数为e。因为pure核不依赖r,u，该算子可原样保留在有序pure域内；没有丢掉recoil。所有HyperInt发散检查开启。

28的直接端点母核路线可积；它与上述完整边界路线的极点及有限项精确相同。27则必须采用带cusp的后一路线。

## 5. 外层可积性与Laurent阶数

初始 $1/2<\Re e<1$ 下，可用去掉排序约束后的全立方体给出绝对值上界：TGV与pure使用不同的Wilson端点，分别可积，因此初始Fubini/分组成立。

ε=0附近先正规化R并解析提取y面与cusp的显式极点，再展开。R和y面余部可由可积的TGV余部乘独立pure核上界控制，附加有限个log不改变收敛。接触面只保留S⁻¹而不是S⁻²；H₀的有理因子为 $1/(1-x+zw)$。与pure端点核在相关角落相乘时，三个小变量承受两个一次分母，仍有正的可积余量；内部正z处由原积分S>0保证正则。非退化边界余部只有可积对数行为。

这支持对已减除表达式逐项展开；不是对未正规化原J直接展开。实际积分全部通过发散检查。27–29最高简单极点；41/42由因子化继承二圈T4的二阶极点。固定z>0，余项O(e)；没有严格数值余项包络。

## 6. 从裸核系数到物理I

记 $\widehat J_g=\mathcal V_g/(4\lambda)$，$\mathcal V_g=v_{g,-1}/e+v_{g,0}+O(e)$。由于G₃没有e¹项，

$$c_{g,-1}=\frac{v_{g,-1}(z)+v_{g,-1}(1/z)}{64},$$
$$c_{g,0}=\frac{v_{g,0}(z)+v_{g,0}(1/z)}{64}+\frac{3\log z}{128}\{v_{g,-1}(1/z)-v_{g,-1}(z)\}.$$

28的self-crossing也在这两个裸系数上精确核验。所有非零系数权重为6+k。

## 7. 27–29的完整ε已积分级数

设 $L=1+z,p_T=2-2e,p_P=1-e$，有限多项式

$$P(1-S/L)^k(1-D_z(p,q)/L)^l=\sum c^{k,l}_{abc,hj;\boldsymbol n}(z)x^ay^bw^cp^hq^j\boldsymbol\alpha^{\boldsymbol n}.$$

α单项积分为

$$d_{\boldsymbol n}(e)=\frac{\prod_{i=1}^3\Gamma(n_i+1-e)}{\Gamma(n_1+n_2+n_3+3-3e)}.$$

有序Wilson单项积分是

$$m_{27}=\frac1{(a+1)(a+h+2)(b+1)(c+1)(c+j+2)},$$
$$m_{28}=\frac1{(a+1)(a+h+2)(b+1)(j+1)(c+j+2)},\quad
m_{29}=\frac1{(h+1)(a+h+2)(b+1)(j+1)(c+j+2)}.$$

记已积完的有限块 $B_{g;k,l}=\sum c^{k,l}d_{\boldsymbol n}m_g$，则

$$\boxed{\widehat J_g=-(1+z)L^{-3+3e}\sum_{N=0}^{\infty}\sum_{k=0}^{N}\frac{(p_T)_k(p_P)_{N-k}}{k!(N-k)!}B_{g;k,N-k}.}$$

每块内先完成全部有符号相加。正z时 $0\le S\le L/3$、$0\le D\le L$，初始绝对收敛strip及其复紧子条带给出支配二项展开，再亚纯延拓到e=0。这是已逐项积分的明确级数，没有剩余待求积分；尚未给加速求和或严格无限尾界。

41/42使用第2节的T4级数×Meijer G×Gamma形式。所有五图的有限项另由正规化与精确Hlog求值获得，而非在无穷和中未经证明逐项作e=0展开。

## 8. 表示与证据边界

27/29和28的有限项需要一个新增的crossing偶函数Ω₆，其定义、完整Hlog、原函数边界及旧空间外的坐标见[基函数](basis-functions.md)。旧14维权重6坐标加Ω₆为15维，五图新目标相对旧空间的商维数为1；这不是全ε的IBP最小性定理。41/42只需原M函数和ζ乘积。

原根、完整有序IBP、因子化、反射、有限块及Ω的独立数值定义检查分别记录。数值误差为经验估计；没有小ε完整I余项检查、第二条完全独立解析路线或Lean整链。TGV仍未注册到Kernel默认47根，物理复z/i0分支仍未验证。