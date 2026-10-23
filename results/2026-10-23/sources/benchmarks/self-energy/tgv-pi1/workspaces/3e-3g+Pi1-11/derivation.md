# 3e-3g+Pi1-11：中间腿自能的全部极点和有限项

## 当前物理积分的入口（2026-10-15，r003）

从共同作用量出发，U′=igAU确定真实Wilson词；FT/Dyson规则在提取C前生效。当前[原生C/J](../../sources/montopar-r002.json)给出下式，J就是本文实际积分的原核，域、分支和IBP全保留：

$$R_P=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad I_P=2R_P\frac{J(z)+z^{3\epsilon}J(1/z)}{1+z^{3\epsilon}}.$$

用下文已经求出的J代入此定义，得到[当前物理I](../../results/3e-3g+Pi1-11/r003.json)；组装系数为+1，Euler只补一次。以下带日期的求值过程与r001公式属于旧根历史，计算J的步骤可按同域积分线性复用，旧形式状态不继承。

2026-10-15表示修订：主I、原J和主Laurent系数保持r002的值；当前源哈希、配对引用及附属公式统一到新修订。当前两个系数中的公共项为+χM₀₂/1024；下文两参数母源已修正，实际使用的a=0切片保持不变。[修订说明](../../revisions/r003.md)。以下2026-10-01/02的r001系数表保留历史符号，不能直接当作当前结果。

2026-10-02。11的三个腿依次落在边2、3、4，自能插在中间腿2。已从实际两分子核求出唯一I的全部极点至有限项，完整纵向端点与recoil都保留。最高极点为ε⁻³；C−1、C0属于原M/N的Q(z)系数代数，含同一个−χM₀₂/1024，χ=(1−z)/(1+z)。全ε有限标准函数闭式尚未求得，状态laurent_complete，formal=not_run。

[结果r001](../../results/3e-3g+Pi1-11/r001.json) · [实际核验及限制](../../../../../workflow/checks/tgv-pi1-11/README.md) · [不可变raw/IBP源](../../sources/montopar-r001.json)。输入SHA256为cc9a45b9a3aecd1437fea5fbe222680f67001f8137c09cb38f46ba6fff337e90。旧01–10、旧分类和生产/Kernel/Lean保持。

## 1. 从原有序积分到保留导数的分解

令$x=s_{21},y=s_{31},w=s_{41}$，$\alpha_3=1-\alpha_1-\alpha_2$。原外到内顺序是x、y、w、α₁、α₂，各路径独立在[0,1]，α₂上界1−α₁；实际权重为
$$\alpha_1^{-\epsilon}\alpha_2^{-2\epsilon}\alpha_3^{-\epsilon}.$$

几何核
$$q=\alpha_1\alpha_2y(1-x)+\alpha_1\alpha_3[(1-x)(1-w)+zxw]
+z\alpha_2\alpha_3w(1-y),$$
$$\nu=1-3\epsilon,\quad m=2-3\epsilon,\quad
\delta=3\epsilon-1=-\nu,\quad\kappa=1+\epsilon,\quad Z=z^{3\epsilon}.$$

J是上述原序域对$K_1/(q+i\eta)^m+K_2/(q+i\eta)^{m+1}$的积分。K₁、K₂与IBP通量采用实际源完整多项式，不取D=4近似。唯一归一化沿用
$$R=\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}
{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad
I_{11}=\frac{2R}{1+Z}[J_{11}(z)+ZJ_{11}(1/z)].$$
C4重数4、一次Euler补偿均按原约定计入；颜色与$\lambda_b^3(S_{14}+S_{12})$在I之外。

定义
$$\begin{aligned}
P=\frac14\{&
\alpha_1\alpha_2[-(1+z)(1-x)+zy]\\
&+\alpha_1\alpha_3[1+3z-(1+z)x+z(1+z)w]
+z\alpha_2\alpha_3[1-y-(1+z)w]\}.
\end{aligned}$$
源L精确等于−P。完整张量分解为
$$\frac{K_1}{Q^m}+\frac{K_2}{Q^{m+1}}
=-\frac{2\kappa P}{Q^m}+\partial_y B_{\rm long},
\qquad Q=q+i\eta,$$
$$B_{\rm long}=\frac{b_{0L}}{\nu}Q^{-\nu}+b_{1y}Q^{-m},
\qquad b_{0L}=\frac{(3-2\epsilon)(1+z)}4(\alpha_3-\alpha_1).$$

$b_{1y}$读取源J_ibp.flux.b1_numerators的第2项，保留完整D维迹。标量部分取等价IBP
$$P=\frac{-zq_x-(1+z)q_y+q_w}{4}+\mathcal R_{\rm tree},\qquad
\mathcal R_{\rm tree}=\frac{1+z}{2}\alpha_3
[\alpha_1(1-x+zw)-z\alpha_2w].$$
标量通量系数$\kappa(-z,-(1+z),1)/2$，余项$-2\kappa\mathcal R_{\rm tree}$。连同$b_{0L}$与实际三个b₁，逐项验证
$$K_1=\sum_i(\partial_i b_{1i}-b_{0i}^{\rm alt}q_i)
-2\kappa\mathcal R_{\rm tree},\qquad
K_2=-m\sum_i b_{1i}q_i.$$
因此没有把recoil删除，也没有把展开后的分子冒充原IBP。

## 2. 完整反射连接两个纵向端点

取
$$\alpha_1\leftrightarrow\alpha_3,\quad\alpha_2\ {\rm fixed},\qquad
(x,y,w)\mapsto(1-w,1-y,1-x),\qquad z\mapsto1/z.$$
实际q、K₁、K₂分别缩放z、z²、z³，标记权重与原域保持，两个坐标变换的绝对Jacobian均为1。有限η同步变为η/z；先在共同η与交换不变的α内部截断上做变换，再取共同极限，得到
$$J_{11}(\epsilon,z)=ZJ_{11}(\epsilon,1/z),\qquad
I_{11}=\frac{4R}{1+Z}J_{11}(\epsilon,z).$$
纵向b₀、b₁在上述反射中分别变为−z、−z²倍，路径方向也反转，故
$$J_{\rm long,upper}(z)=ZJ_{\rm long,lower}(1/z).$$

这两个纵向积分均须保留；10号的上纵向反对称判零不适用于中间腿2。

## 3. 先积α尺度，再用混合偏导降低分母幂

在y=0取
$$\alpha_3=a,\quad\alpha_1=(1-a)b,\quad
\alpha_2=(1-a)(1-b),$$
$$q=a(1-a)\mathcal Q,\qquad
\mathcal Q=b[(1-x)(1-w)+zxw]+z(1-b)w.$$
Jacobian为1−a。第一通量和含一个$a(1-a)$的第二通量均留下$a^{2\epsilon-1}da$，线性多项式平均为
$$p_0/(2\epsilon)+p_1/(1+2\epsilon).$$

含下端外法向负号的完整纵向面变为
$$\int_0^1db\,dx\,dw\ b^{-\epsilon}(1-b)^{-2\epsilon}
[N_0\mathcal Q^{-\nu}+N_1\mathcal Q^{-m}],$$
$$N_0=\frac{(b-2\epsilon)(3-2\epsilon)(1+z)}
{8\epsilon\nu(1+2\epsilon)},$$
$$N_1=-\frac{(b-2\epsilon)(1+z)\mathcal Q}{8\epsilon(1+2\epsilon)}
+\frac{bz(1+z)w}{8\epsilon}
+\frac{z(b+1)}{4(1+2\epsilon)}.$$

两个偏导恒等式提供了降幂工具：
$$\mathcal Q_x\mathcal Q_w=b[(1+z)\mathcal Q-z],
\qquad\mathcal Q_{xw}=b(1+z),$$
$$z\mathcal Q^{-m}=
\frac{3\epsilon(1+z)\mathcal Q^{-\nu}
-(3\epsilon b)^{-1}\partial_x\partial_w\mathcal Q^{3\epsilon}}{\delta},$$
$$(1+z)w\mathcal Q^{-m}
=\mathcal Q^{-m}+\frac1{\delta b}\partial_x\mathcal Q^{-\nu}.$$

合并完整N₀、N₁后，密度精确成为
$$\begin{aligned}
&(1+z)\left[\frac b{4\epsilon}+\frac1{2\delta}\right]\mathcal Q^{-\nu}
+\frac z{8\epsilon\delta}\partial_x\mathcal Q^{-\nu}\\
&\hspace{1cm}-\frac{(1+4\epsilon)b+2\epsilon}
{24\epsilon^2(1+2\epsilon)b\delta}\,
\partial_x\partial_w\mathcal Q^{3\epsilon}.
\end{aligned}$$

混合导数的四个端点为z、z(1−b)、0、b，产生Gamma/Beta边界；单x导数还产生一个可直接识别为原M母源的一维差。这里先在收敛的正ε条带取η→0，再使用$\mathcal Q$的偏导恒等式；它不是有限η下删掉修正项的等式。四个非整数ε、b、z点从原a平均密度直接积w，与偏导边界式的差≤5.31×10⁻³⁶。

令
$$\begin{aligned}
Y&=\int_0^1db\,dx\,dw\ b^{-\epsilon}(1-b)^{-2\epsilon}\mathcal Q^{-\nu},\\
Y_B&=\int_0^1db\,dx\,dw\ b^{1-\epsilon}(1-b)^{-2\epsilon}\mathcal Q^{-\nu}.
\end{aligned}$$
令$\rho=z/(1+z)$，$\ell=\log z$，
$$\mathcal U(a,c;z)=\int_0^1\frac{
\exp\{a\log[\beta(1-\beta)]\}
\left(\exp\{c[\ell+\log((1-\beta)/\beta)]\}-1\right)}
{\beta-\rho}\,d\beta.$$
此定义与09号的同一母源一致。先取$z>0$及$\operatorname{Re}(1+a-c)>0$、$\operatorname{Re}(1+a+c)>0$保证两端可积，再作参数解析延拓。在$\beta=\rho$，$\ell+\log((1-\rho)/\rho)=0$，故整个分子为零，内部点可去。本题只用$\mathcal U(0,\epsilon)-\mathcal U(0,-2\epsilon)$。

2026-10-15勘误：旧稿分子写成$e^{aW+cL}-1$，一般$a\ne0$时在内部点不为零。这里改为$e^{aW}(e^{cL}-1)$；$a=0$切片及此前实际使用的两个母源差均不变，不修改已求的J或I系数。

完整下纵向贡献为
$$\Lambda(z)=L_A+ZL_B+\mathcal F(z),$$
$$\mathcal F(z)=\frac{1+z}{4\epsilon}Y_B
+\frac{1+z}{2\delta}Y
+\frac{z^{1+2\epsilon}}
{24\epsilon^2\delta(1+z)}
[\mathcal U(0,\epsilon;z)-\mathcal U(0,-2\epsilon;z)],$$
$$L_A=-\frac{(1+4\epsilon)\mathrm B(1+2\epsilon,1-2\epsilon)
+2\epsilon\mathrm B(2\epsilon,1-2\epsilon)}
{24\epsilon^2(1+2\epsilon)\delta},$$
$$\begin{aligned}
L_B={}&\frac{\mathrm B(1-\epsilon,1-2\epsilon)}{24\epsilon^2\delta}\\
&-\frac{(1+4\epsilon)[\mathrm B(1-\epsilon,1-2\epsilon)-\mathrm B(1-\epsilon,1+\epsilon)]}
{24\epsilon^2(1+2\epsilon)\delta}\\
&-\frac{2\epsilon[\mathrm B(-\epsilon,1-2\epsilon)-\mathrm B(-\epsilon,1+\epsilon)]}
{24\epsilon^2(1+2\epsilon)\delta}.
\end{aligned}$$

负Beta参数解释为共同端点差的解析延拓，不能把单项视为收敛Euler积分。三个正ε点的9个收敛Beta积分支持上述边界转换。上下纵向合并给$(1+Z)(L_A+L_B)+\mathcal F(z)+Z\mathcal F(1/z)$。

## 4. 六个标量面和完整余项

两个cusp面相连，标量尺度与路径尺度积分给
$$C_W=C_X=\frac{\mathrm B(-2\epsilon,1-\epsilon)
-\mathrm B(-2\epsilon,1+2\epsilon)}{18\epsilon^3}.$$
原w=0面有$\alpha_1=a$尺度1/(2ε)，$(1-x)$尺度1/(3ε)，再积y给端点差/(3ε)，因此这里的−2ε来自实际中间标记腿权重。

另两个非退化面记h_X、h_W，反射给
$$h_W(z)=z^{-\nu}h_X(1/z).$$
y=0,1面均有α尺度1/(2ε)，且上面为$z^{-\nu}Y(1/z)$。它们在原J中保留为
$$(1+z)[Y(z)-z^{-\nu}Y(1/z)]/(2\epsilon).$$
整个配对在$I$的C4平均中相消，不能分别判零。

正图表取
$$\alpha=(r,s,1)/(1+r+s),\quad
x=t/(1+t),\quad y=u/(1+u),\quad w=v/(1+v).$$
令H=1+r+s、X=(1+t)(1+u)(1+v)，
$$F=rsuv+rsu+rtuvz+rtvz+ru+r+stvz+svz,\qquad
F_X=rsuv+rsu+ru+r+svz.$$
$$h_X=\int_{\mathbb R_+^4}r^{-\epsilon}s^{-2\epsilon}
H^{-1-2\epsilon}[(1+u)(1+v)]^{-1-3\epsilon}F_X^{-\nu},$$
$$f_{\mathcal R}=r^{-\epsilon}s^{-2\epsilon}H^{-1-2\epsilon}
X^{-1-3\epsilon}F^{-m}N_{\mathcal R},$$
$$N_{\mathcal R}=-\kappa(1+z)(1+u)
[r(1+v)+z(r-s)v(1+t)].$$

余项的几何F和N与10相同，但实际基因子多$s^{-\epsilon}$。这不是积分相等。对实际新权重重新扫描，唯一区域仍为{r,1/u}，尺度次数2ε。令$E=r\partial_r-u\partial_u$，则
$$f_{\mathcal R}\mapsto f_{\mathcal R}-Ef_{\mathcal R}/(2\epsilon).$$
E消去额外$s^{-\epsilon}$，因而可在精确几何桥接后复用旧正规化分子，新的对数展开因子为旧因子−log s。全εEuler分子恒等式再次验证，正规化后区域扫描为空。原小正ε积分收敛、总导数端点为零；在同一积分上正规化后再做解析延拓和展开。

有限面h_X的ε⁰、ε¹，Y和Y_B各ε⁰、ε¹、ε²，正规化recoil的ε⁻¹、ε⁰、ε¹，共11个所需系数全部实际积分。采用HyperInt遇发散中止，顺序分别[r,u,s,v]、[x,w,b]、[r,u,s,t,v]。不是仅生成待积表达式。

完整原J为
$$\begin{aligned}
J_{11}={}&(1+Z)\left[L_A+L_B-\frac{\kappa C_W}{2\nu}\right]
+\mathcal F(z)+Z\mathcal F(1/z)\\
&+\frac{\kappa}{2\nu}\left[
zh_X(z)+z^{-\nu}h_X(1/z)
+\frac{1+z}{2\epsilon}(Y(z)-z^{-\nu}Y(1/z))\right]
+\mathcal R_\epsilon(z).
\end{aligned}$$

计算从共同η、α内部截断的紧域开始；降维与变换在正ε收敛条带取极限，之后作亚纯延拓。正规化后在固定z>0的小ε邻域按可积幂/对数界逐项展开。保存原有序域与变换记录，没有把CAS扫描标成形式测度桥。

## 5. 完整I的M/N表达式

沿用原定义
$$M_{k,j}(z)=\int_0^1\frac{\log^k[\beta(1-\beta)]
[\log z+\log((1-\beta)/\beta)]^j}{\beta-z/(1+z)}\,d\beta,\quad j\ge1,$$
$$N_{k,j}=M_{k,j+1}-\log z\,M_{k,j},\qquad
\chi=\frac{1-z}{1+z},\quad\zeta_n=\zeta(n).$$
内部点可去，不取主值。母源差展开为
$$\mathcal U(0,\epsilon)-\mathcal U(0,-2\epsilon)
=\sum_{j\ge1}\frac{\epsilon^j-(-2\epsilon)^j}{j!}M_{0,j}.$$

完整结果
$$I_{11}=\frac{c_{-3}}{\epsilon^3}
+\frac{c_{-2}}{\epsilon^2}+\frac{c_{-1}}{\epsilon}+c_0+O(\epsilon),$$
$$c_{-3}=\frac{1-2\zeta_2}{384},$$
$$c_{-2}=-\frac{3M_{0,1}}{256}+\frac{M_{1,1}}{128}
-\frac{\zeta_2}{96}-\frac{\zeta_3}{192}+\frac1{384},$$
$$\begin{aligned}
c_{-1}={}&-\frac{5M_{0,1}^2}{1536}
+\frac{\zeta_2M_{0,1}}{256}-\frac{3M_{0,1}}{256}
-\frac{\chi M_{0,2}}{1024}
-\frac{M_{1,1}}{128}+\frac{7M_{2,1}}{512}\\
&-\frac{193\zeta_2^2}{1920}+\frac{5\zeta_2}{768}
-\frac{\zeta_3}{96}+\frac3{128},
\end{aligned}$$
$$\boxed{\begin{aligned}
c_0={}&-\frac{5M_{0,1}^2}{192}-\frac{M_{0,1}M_{1,1}}{512}
-\frac{41\zeta_2M_{0,1}}{512}-\frac{5\zeta_3M_{0,1}}{32}
-\frac{9M_{0,1}}{256}-\frac{\chi M_{0,2}}{1024}\\
&-\frac{13\zeta_2M_{1,1}}{512}+\frac{M_{1,1}}{128}
-\frac{17M_{1,3}}{3072}+\frac{5M_{2,1}}{512}
+\frac{9M_{3,1}}{1024}+\frac{5N_{1,2}}{256}\\
&-\frac{193\zeta_2^2}{960}+\frac{139\zeta_2\zeta_3}{384}
-\frac{47\zeta_2}{768}+\frac{11\zeta_3}{384}
-\frac{13\zeta_5}{192}+\frac{17}{384}.
\end{aligned}}$$

ε⁻⁶、ε⁻⁵、ε⁻⁴均为零。非恒定项$\chi M_{0,2}$在z反演下保持，因为两个因子分别变号；它与其余常系数M/N项共同满足crossing。现行常系数偶基底不足以表示这项；补入已验证的原M₀₂后，Q(z)消元与独立Maple全回代全部精确通过，没有定义新矩。

| z | c−2 | c−1 | c0 |
|---:|---:|---:|---:|
| 1 | 0.1233531841236616 | −0.8099448261414941 | 2.992624373227649 |
| 2 | 0.1308190529280316 | −0.8477365584668339 | 3.097308592536944 |

z=1/2与z=2相同，c−3≈−0.005963198264834513。15个原M/N数值积分与独立Hlog求值的最大差约2.92×10⁻³⁵。这里只建立至有限项的表示，不宣称全εM/N闭包。

## 6. 颜色与后续范围

两个颜色通道仍为
$$CF\,CA^2\,\frac{10-6\epsilon}{3-2\epsilon}I_{11},
\qquad CF\,CA\,TF\,n_f\,\frac{-8(1-\epsilon)}{3-2\epsilon}I_{11}.$$
每通道ε⁻⁴位置为零，ε⁻³至ε⁰四个非零系数已保存，并以M/N与Hlog分别核验。最高留数为$5(1-2\zeta_2)/576$和$-(1-2\zeta_2)/144$。这些是C乘I的颜色部分，完整物理尺度、图重数与RG仍需另行组装。

本次29源/反射/降幂检查、14源面/尺度检查、5原J反射、5I crossing、5M/N、10颜色、Euler恒等式及正规化扫描通过；13面/Beta数值点、15原M/N点与4组原/正规化recoil控制提供支持。QMC精确积掉w后比较同一四维余项，每sector2¹⁶、2¹⁸点、3扰动；末均值差z=1约0.0098、经验SE约0.0051，z=2约0.0209、SE约0.0111，约为2个SE以内。它不是严格区间证书，也没有独立完整五维原J或完整小εI截断检验。

分支有限阶11/12，全ε标准函数闭式仍7/12；下一12为同234核心的首腿1插入。完整w³数学口径122已得、12零、15未完整交付（双TGV10、腿Π₁1、中央Γ₁4）。无可调用Wilson数学MCP、无形式收据，formal=not_run；未提交Git。
