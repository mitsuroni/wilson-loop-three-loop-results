# 3e-3g+Pi1-10：三边核心的全部极点与有限项

## 当前物理积分的入口（2026-10-15，r003）

从共同作用量出发，U′=igAU确定真实Wilson词；FT/Dyson规则在提取C前生效。当前[原生C/J](../../sources/montopar-r002.json)给出下式，J就是本文实际积分的原核，域、分支和IBP全保留：

$$R_P=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad I_P=2R_P\frac{J(z)+z^{3\epsilon}J(1/z)}{1+z^{3\epsilon}}.$$

用下文已经求出的J代入此定义，得到[当前物理I](../../results/3e-3g+Pi1-10/r003.json)；组装系数为+1，Euler只补一次。以下带日期的求值过程与r001公式属于旧根历史，计算J的步骤可按同域积分线性复用，旧形式状态不继承。

2026-10-15表示修订：主I、原J和主Laurent系数保持r002的值；当前源哈希、配对引用及附属公式统一到新修订。[修订说明](../../revisions/r003.md)。以下2026-10-01/02的r001系数表保留历史符号，不能直接当作当前结果。

2026-10-02。10为234三边核心，自能在末边4、腿3。本图从完整两项raw核独立求出ε⁻⁴至ε⁰，全部系数在原M/N与ζ多项式中，仅常系数。recoil完整保留并正规化后实际积分。全ε的有限标准函数闭式尚未求得；状态laurent_complete，formal=not_run。

[结果r001](../../results/3e-3g+Pi1-10/r001.json) · [实际核验](../../../../../workflow/checks/tgv-pi1-10/README.md)。不可变输入为[raw/IBP包](../../sources/montopar-r001.json)，SHA256=cc9a45b9a3aecd1437fea5fbe222680f67001f8137c09cb38f46ba6fff337e90。原01–09与生产/Kernel/Lean保持。

## 1. 完整I与两个分子的起点

记 $x=s_{21},y=s_{31},w=s_{41}$，三个路径参数独立在[0,1]，原外到内顺序是x、y、w、α₁、α₂。权重
$\alpha_1^{-\epsilon}\alpha_2^{-\epsilon}\alpha_3^{-2\epsilon}$，
$\alpha_3=1-\alpha_1-\alpha_2$。令
$$q=\alpha_1\alpha_2y(1-x)+\alpha_1\alpha_3[(1-x)(1-w)+zxw]
+z\alpha_2\alpha_3w(1-y),$$
$$\nu=1-3\epsilon,\quad m=2-3\epsilon,\quad
\delta=3\epsilon-1=-\nu,\quad\kappa=1+\epsilon.$$

R与唯一I沿用
$$R=\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}
{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad
I_{10}=\frac{2R}{1+z^{3\epsilon}}
[J_{10}(z)+z^{3\epsilon}J_{10}(1/z)].$$

C4重数4和一次三圈Euler补偿已在同一约定中；颜色及
$\lambda_b^3(S_{14}+S_{12})$ 留在I之外。

## 2. 等价IBP：标量部分与标记路径全导数

设
$$\begin{aligned}
P=\frac14\{&
\alpha_1\alpha_2[-(1+z)(1-x)+zy]\\
&+\alpha_1\alpha_3[1+3z-(1+z)x+z(1+z)w]
+z\alpha_2\alpha_3[1-y-(1+z)w]\}.
\end{aligned}$$

实际源的L严格等于−P。完整自能张量的标量/二阶导数分解给出同一个原核
$$\frac{K_1}{Q^m}+\frac{K_2}{Q^{m+1}}
=-\frac{2\kappa P}{Q^m}+\partial_w B_{\rm long},
\qquad Q=q+i\eta,$$
$$B_{\rm long}=\frac{b_{0L}}{\nu}Q^{-\nu}+b_1Q^{-m},
\qquad b_{0L}=-\frac{3-2\epsilon}{4}(\alpha_1-\alpha_2).$$

$b_1=-H_3$ 读取实际源完整多项式，D=4−2ε的迹已保留。可用Gram形式
$H_3=-\alpha_1T_{13}T_{21}+\alpha_2T_{23}T_{12}
+G_{12}(\alpha_1W_{31}-\alpha_2W_{32})$，
$T_{ij}=\bar u_i\cdot\mathcal R_j$、$W_{ij}=\mathcal R_i\cdot\mathcal R_j$。

为了利用已有三边图表，将标量P取另一种精确IBP：
$$P=\frac{-zq_x-(1+z)q_y+q_w}{4}+\mathcal R_{\rm tree},\qquad
\mathcal R_{\rm tree}=\frac{1+z}{2}\alpha_3
[\alpha_1(1-x+zw)-z\alpha_2w].$$

相应标量通量的b₀向量是
$\kappa(-z,-(1+z),1)/2$，余项为 $-2\kappa\mathcal R_{\rm tree}$。
这与生产输出选择的b₀/recoil分配不同，但完整两分子已逐项回代：
$$K_1=\sum_i(\partial_i b_{1i}-b_{0i}^{\rm alt}\partial_iq)
-2\kappa\mathcal R_{\rm tree},\qquad
K_2=-m\sum_i b_{1i}\partial_iq.$$

原输入、C和生产代码没有修改，也没有删除recoil；这是求解侧的等价分项。

## 3. 纵向上面严格为零，下端闭合Gamma

在w=1的纵向部分，取
$$\alpha_1\leftrightarrow\alpha_2,\qquad
(x,y)\mapsto(1-y,1-x).$$

实际q保持不变，b₀L和b₁均变号。两个Jacobian绝对值为1，α₃和实际权重不变。有限共同η及交换不变的α内部截断下，紧域上可换序，积分等于自己的负数，因此这一完整分项为0。取共同调节极限仍为0。**此处仅纵向上面相消，标量上面仍保留。**

w=0纵向下端取 $\alpha_1=a,\alpha_2=(1-a)b,\alpha_3=(1-a)(1-b)$、
$X=1-x$。实际源核变成
$$q=a(1-a)X(1-b+by),$$
$$b_1=\frac{a(1-a)X}{4}
[a(1+by)+(1-a)b^2(1-y)].$$

a测度为 $a^{2\epsilon-1}da$，X尺度给 $1/(3\epsilon)$。a平均后，第二分子
$2\epsilon(1+by)+b^2(1-y)=(2\epsilon-b)(1-b+by)+b(1+2\epsilon)$；
y积分于是只产生端点幂差。定义
$$D_0=\mathrm B(-\epsilon,1-2\epsilon)-\mathrm B(-\epsilon,1+\epsilon),$$
$$D_1=\mathrm B(1-\epsilon,1-2\epsilon)-\mathrm B(1-\epsilon,1+\epsilon),$$
$$D_2=\mathrm B(1-\epsilon,1-2\epsilon)-\mathrm B(1-\epsilon,\epsilon).$$

包括下端外法向负号的纵向贡献为
$$\Lambda(\epsilon)=
\frac{(2+\epsilon)(2\epsilon D_0-D_1)}
{72\epsilon^3(1+2\epsilon)\nu}
-\frac{D_2}{24\epsilon^2\delta}.$$

负Beta参数只解释为共同端点差分的解析延拓。先加共同辅助端点指数让分项收敛，再对完整差分延拓；不能将分项当作普通正ε收敛Euler积分。三个非零ε点的原收敛Euler差分与该Gamma式吻合，差小于约10⁻⁴⁸。

## 4. 六个标量面与recoil全部计算

设 $B_a=\mathrm B(\epsilon,1+\epsilon)$，
$$F_s=\frac{\mathrm B(1-\epsilon,-\epsilon)-\mathrm B(1+2\epsilon,-\epsilon)}
{3\epsilon},\qquad
C_X=\frac{B_aF_s}{3\epsilon},\qquad C_W=\frac{D_0}{18\epsilon^3}.$$

$x=1$、$w=0$ 两个标量cusp面精确积为 $z^{-1+3\epsilon}C_X$、$C_W$。
$y=0,1$ 的α尺度分别给 $B_a$、$1/(2\epsilon)$，余下有限核为
$$\begin{aligned}
Y_0&=\int_0^1db\,dx\,dw\,[b(1-b)]^{-\epsilon}
[(1-b)(1-x)(1-w)+z(b+(1-b)x)w]^{-1+3\epsilon},\\
Y_1&=\int_0^1db\,dx\,dw\,b^{-\epsilon}(1-b)^{-2\epsilon}
[(1-x)(1-(1-b)w)+zx(1-b)w]^{-1+3\epsilon}.
\end{aligned}$$

它们的ε⁰、ε¹、ε²六个系数均已实际积成Hlog。

另两个非退化面用
$\alpha=(r,s,1)/(1+r+s)$、
$x=t/(1+t),y=u/(1+u),w=v/(1+v)$。
αJacobian为 $(1+r+s)^{-3}$，每个路径Jacobian为 $(1+\text{参数})^{-2}$。
令H=1+r+s，
$$F_X=rsuv+rsu+ru+r+svz,\qquad
F_W=rsu+rtuz+rtz+stz+sz.$$
相应有限四维核为
$$h_X=\int_{\mathbb R_+^4}(rs)^{-\epsilon}
H^{-1-2\epsilon}[(1+u)(1+v)]^{-1-3\epsilon}F_X^{-1+3\epsilon},$$
$$h_W=\int_{\mathbb R_+^4}(rs)^{-\epsilon}
H^{-1-2\epsilon}[(1+t)(1+u)]^{-1-3\epsilon}F_W^{-1+3\epsilon}.$$
四个所需ε⁰、ε¹系数已实际积分。

完整recoil图表有
$$F=rsuv+rsu+rtuvz+rtvz+ru+r+stvz+svz,$$
$$f_{\mathcal R}=(rs)^{-\epsilon}
H^{-1-2\epsilon}[(1+t)(1+u)(1+v)]^{-1-3\epsilon}
F^{-2+3\epsilon}N_{\mathcal R},$$
$$N_{\mathcal R}=-\kappa(1+z)(1+u)
[r(1+v)+z(r-s)v(1+t)].$$

唯一临界联合区域是r→0、u→∞，尺度次数为2ε。若
$E=r\partial_r-u\partial_u$，使用保留ε的Euler分部积分
$$f_{\mathcal R}\longmapsto f_{\mathcal R}-\frac{Ef_{\mathcal R}}{2\epsilon}.$$

原积分在小正ε共同条带收敛，端点总导数边界为零；之后作亚纯延拓。不是在发散的ε=0原积分上取有限部分。新核的余部扫描为空，只有显式简单极点。其ε⁻¹、ε⁰、ε¹三个系数在开启发散检查的HyperInt中实际积分完成。正规化算子的完整分子恒等式再次精确核验；ε=.2的原recoil与正规化四维QMC控制也一致于经验误差内。

最终原J是
$$\begin{aligned}
J_{10}={}&G_0+z^{3\epsilon}G_Z+
\frac{\kappa}{2\nu}
[z h_X+h_W+(1+z)(B_aY_0-Y_1/(2\epsilon))]
+\mathcal R_\epsilon,\\
G_0={}&\Lambda-\kappa C_W/(2\nu),\qquad
G_Z=-\kappa C_X/(2\nu).
\end{aligned}$$

$\mathcal R_\epsilon$ 已含标量系数−2κ。Gamma项按所需精度展开；有限h需到ε¹、Y需到ε²、recoil需到ε¹，才足够乘R的1/ε获得I有限项。

调节和换序始终从共同η、α内部截断下的紧域开始；cusp/下纵向面在0<Reε<1/3的原域可积，有限面与正规化余部通过具体尺度界及开启的发散检查。正规化后在固定z>0的小ε邻域用可积幂/对数界逐项展开。各程序输出和失败探索保存在检查目录；扫描与CAS成功没有被标成Lean证明。

## 5. 五个最终系数：原M/N常系数多项式

沿用
$$M_{k,j}(z)=\int_0^1
\frac{\log^k[\beta(1-\beta)]
[\log z+\log((1-\beta)/\beta)]^j}{\beta-z/(1+z)}\,d\beta,\quad j\ge1,$$
$$N_{k,j}=M_{k,j+1}-\log(z)M_{k,j},\qquad \zeta_n=\zeta(n).$$

内部点可去，无主值积分。完整结果
$$I_{10}=\sum_{q=-4}^0c_q(z)\epsilon^q+O(\epsilon),$$
$$c_{-4}=-\frac1{768},\qquad
c_{-3}=-\frac1{384}-\frac{\zeta_2}{128},$$
$$c_{-2}=\frac{3M_{1,1}}{256}
-\frac{7\zeta_2}{512}-\frac{\zeta_3}{192}-\frac1{192},$$
$$\begin{aligned}
c_{-1}={}&-\frac{7M_{0,1}^2}{3072}
+\frac{\zeta_2M_{0,1}}{512}
+\frac{3M_{1,1}}{128}+\frac{19M_{2,1}}{1024}\\
&-\frac{33\zeta_2^2}{320}-\frac{11\zeta_2}{256}
-\frac{\zeta_3}{256}-\frac7{384},
\end{aligned}$$
$$\boxed{\begin{aligned}
c_0={}&-\frac{7M_{0,1}^2}{1536}
+\frac{13M_{0,1}M_{1,1}}{1024}
+\frac{\zeta_2M_{0,1}}{256}-\frac{\zeta_3M_{0,1}}8\\
&+\frac{37\zeta_2M_{1,1}}{1024}+\frac{3M_{1,1}}{64}
+\frac{19M_{2,1}}{512}
+\frac{27M_{3,1}+5M_{1,3}}{2048}
+\frac{7N_{1,2}}{512}\\
&-\frac{1923\zeta_2^2}{10240}
+\frac{23\zeta_2\zeta_3}{64}
-\frac{3\zeta_2}{64}-\frac{5\zeta_3}{128}
+\frac{\zeta_5}{96}-\frac7{192}.
\end{aligned}}$$

ε⁻⁶、ε⁻⁵为零。五项从实际积分生成后才作精确M/N消元，随后独立新CAS完整回代，差均为零；没有用已知系数作积分输入。

| z | c₋₂ | c₋₁ | c₀ |
|---:|---:|---:|---:|
| 1 | 0.09551374128866734 | −0.7551977291479347 | 3.056639946351306 |
| 2 | 0.10248981292758194 | −0.7824685094664424 | 3.173575390538622 |

z=1/2与z=2相同。五个crossing精确通过；三个运动学点直接M/N积分共15系数点，最大差约1.63×10⁻³⁴。全εM/N闭包没有由有限阶匹配推断。

## 6. 两颜色与后续任务

两颜色仍为
$$\frac{10-6\epsilon}{3-2\epsilon}I_{10},\qquad
-\frac{8(1-\epsilon)}{3-2\epsilon}I_{10}.$$
各五个完整系数已保存，十条颜色等式在原M/N与Hlog中精确核验。
最高留数为−5/1152、1/288。外尺度、物理总和与RG尚未组装。

本次Γ边界、四个非退化面系数、六个Y系数、三个完整recoil系数均计算完毕。上纵向面曾作四区域正规化探索，后经严格反对称证明为零；该探索输出保留但不参与最终组装。原recoil数值检查精确积掉标记路径后做四维随机QMC，末次原/正规化均值差约4.5×10⁻⁴、2.7×10⁻³，各在经验标准误差尺度内；不是严格区间证书或完整五维原J积分。

有限阶10/12，全ε标准函数闭式仍7/12。下一11为234核心中间腿2插入，须重新核验权重、纵向面和recoil。当前无可调用Wilson数学MCP，无形式收据，formal=not_run，未提交Git。
