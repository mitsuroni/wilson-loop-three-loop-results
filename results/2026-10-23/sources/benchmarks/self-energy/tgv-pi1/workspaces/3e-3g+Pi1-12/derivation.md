# 3e-3g+Pi1-12：首腿插入与10号的完整积分等式

## 当前物理积分的入口（2026-10-15，r003）

从共同作用量出发，U′=igAU确定真实Wilson词；FT/Dyson规则在提取C前生效。当前[原生C/J](../../sources/montopar-r002.json)给出下式，J就是本文实际积分的原核，域、分支和IBP全保留：

$$R_P=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad I_P=2R_P\frac{J(z)+z^{3\epsilon}J(1/z)}{1+z^{3\epsilon}}.$$

用下文已经求出的J代入此定义，得到[当前物理I](../../results/3e-3g+Pi1-12/r003.json)；组装系数为+1，Euler只补一次。以下带日期的求值过程与r001公式属于旧根历史，计算J的步骤可按同域积分线性复用，旧形式状态不继承。

2026-10-15表示修订：主I、原J和主Laurent系数保持r002的值；当前源哈希、配对引用及附属公式统一到新修订。[修订说明](../../revisions/r003.md)。以下2026-10-01/02的r001系数表保留历史符号，不能直接当作当前结果。

2026-10-02。12的三个腿依次落在边2、3、4，自能插在首腿1。完整反射给出全ε关系
$$J_{12}(\epsilon,z)=z^{3\epsilon}J_{10}(\epsilon,1/z),\qquad
I_{12}(\epsilon,z)=I_{10}(\epsilon,z).$$
已将原J重新变换并组装到I的ε⁰，所有极点、有限项、原M/N表示和两颜色保存。最高四阶极点−1/768。全ε关系已建立，但复用的10号尚无全ε积完标准函数闭式；12的状态为laurent_complete，formal=not_run。

[结果r001](../../results/3e-3g+Pi1-12/r001.json) · [完整桥接与实际核验](../../../../../workflow/checks/tgv-pi1-12/README.md) · [10号原积分与求解](../3e-3g+Pi1-10/derivation.md)。不可变[raw/IBP源](../../sources/montopar-r001.json)的SHA256为cc9a45b9a3aecd1437fea5fbe222680f67001f8137c09cb38f46ba6fff337e90。原源、01–11、旧分类及生产/Kernel/Lean保持。

## 1. 插入位置改变了权重，必须变换整个核

记$x=s_{21},y=s_{31},w=s_{41}$，$\alpha_3=1-\alpha_1-\alpha_2$。原外到内顺序x、y、w、α₁、α₂；三个路径参数独立在[0,1]，α₂上界1−α₁。12的实际权重为
$$W_{12}=\alpha_1^{-2\epsilon}\alpha_2^{-\epsilon}\alpha_3^{-\epsilon},$$
10的权重为$\alpha_1^{-\epsilon}\alpha_2^{-\epsilon}\alpha_3^{-2\epsilon}$。因此不能只调换“自能腿”标签后复制原J。

两图的几何分母为
$$q=\alpha_1\alpha_2y(1-x)+\alpha_1\alpha_3[(1-x)(1-w)+zxw]
+z\alpha_2\alpha_3w(1-y).$$
完整核是
$$W_r\left[\frac{K_{1,r}}{(q+i\eta)^{2-3\epsilon}}
+\frac{K_{2,r}}{(q+i\eta)^{3-3\epsilon}}\right].$$
两个K均读实际源全D多项式；本次没有取D=4或删除第二分母幂。

令从12到10的拉回为
$$\Phi:\quad
(\alpha_1',\alpha_2',\alpha_3')=(\alpha_3,\alpha_2,\alpha_1),$$
$$(x',y',w')=(1-w,1-y,1-x),\qquad
z'=1/z,\qquad\eta'=\eta/z.$$
原路径域、正α单纯形保持，路径与α坐标的绝对Jacobian均为1，Φ为对合。实际源逐项满足
$$q_{12}=z\,\Phi^*q_{10},\quad
K_{1,12}=z^2\Phi^*K_{1,10},\quad
K_{2,12}=z^3\Phi^*K_{2,10},\quad
W_{12}=\Phi^*W_{10}.$$

有限调节满足
$$q_{12}+i\eta=z(\Phi^*q_{10}+i\eta/z).$$
对z>0，正实缩放保持上半平面支路，故两个分母幂分别给
$z^{2-(2-3\epsilon)}=z^{3\epsilon}$与
$z^{3-(3-3\epsilon)}=z^{3\epsilon}$。
两项得到相同因子，才可提出整体J关系。

先在共同有限η、交换不变的α内部截断上做紧域变量变换，再利用10的正ε共同收敛条带取极限并作ε亚纯延拓；没有暗中把原有序输入当成任意可置换的未调节乘积积分。原全色C、相位账本、裸归一化、C4重数4和幂处方均相符。6个非整数ε/z的原加权核点（分别η=0、0.017，旧侧用η/z）也支持此完整等式；数值点不是形式证明。

## 2. IBP保留：第一通量的差由余项严格补偿

令$\nu=1-3\epsilon$、$m=2-3\epsilon$、$\kappa=1+\epsilon$。生产IBP为
$$\frac{K_1}{Q^m}+\frac{K_2}{Q^{m+1}}
=\sum_i\partial_i\left(\frac{b_{0i}}{\nu}Q^{-\nu}+b_{1i}Q^{-m}\right)
+\mathcal R_{\rm prod}Q^{-m},\qquad Q=q+i\eta.$$

若P是路径Jacobian矩阵，则10的通量变换为
$$\widetilde b_0=zP^{-1}\Phi^*b_{0,10},\qquad
\widetilde b_1=z^2P^{-1}\Phi^*b_{1,10}.$$
实际12的第二通量精确等于$\widetilde b_1$，第一通量却有
$$b_{0,12}=\widetilde b_0+D,\qquad
D=\kappa(-z,1+z,-1).$$
同时完整余项满足
$$\mathcal R_{{\rm prod},12}-z^2\Phi^*\mathcal R_{{\rm prod},10}
=D\cdot\nabla q_{12}.$$
因为$\nabla\cdot D=0$，
$$\sum_i\partial_i\left(\frac{D_i}{\nu}Q^{-\nu}\right)
=-(D\cdot\nabla q_{12})Q^{-m}.$$
所以第一通量差与余项差在原完整核上逐点抵消；积分后是六个面的总修正与体余项修正相消。此恒等式在共同有限η下成立，不需要先令η=0。

| 12的面 | 对应10的面 | 第一通量分子额外项 |
|---|---|---:|
| x=0 | w′=1 | κz |
| x=1 | w′=0 | −κz |
| y=0 | y′=1 | −κ(1+z) |
| y=1 | y′=0 | κ(1+z) |
| w=0 | x′=1 | κ |
| w=1 | x′=0 | −κ |

表中分子额外项是$n\cdot D$，仍乘$Q^{-\nu}/\nu$；面方向已计入。六个q、第一通量加修正、第二通量、面测度和C/共同调节组均按源实际面坐标核验。不能只比较b₀就否定积分等式，也不能因b₁正确便省略余项。

## 3. 纵向边界与非零标量余项如何继承

10的成功分解是普通标量顶点加标记路径全导数。令P₀为该普通标量顶点，其实际L=−P₀，且$P_{0,12}=z^2\Phi^*P_{0,10}$。拉回后
$$\frac{K_{1,12}}{Q^m}+\frac{K_{2,12}}{Q^{m+1}}
=-2\kappa P_0Q^{-m}+\partial_x B_{{\rm long},12},$$
$$B_{{\rm long},12}=\frac{b_{0L,12}}{\nu}Q^{-\nu}+b_{1x,12}Q^{-m},
\qquad b_{0L,12}=\frac{z(3-2\epsilon)}4(\alpha_3-\alpha_2).$$
$b_{1x,12}$是实际源完整第二通量；没有改变D维迹。

等价标量IBP可写成
$$P_0=\frac{-zq_x+(1+z)q_y+q_w}{4}+\mathcal R_{{\rm tree},12},$$
$$\mathcal R_{{\rm tree},12}=\frac{1+z}{2}\alpha_1
[z\alpha_3w+(\alpha_3-\alpha_2)(1-x)]
=z^2\Phi^*\mathcal R_{{\rm tree},10}.$$
相应b₀向量为
$$\left(-\frac{\kappa z}{2}+b_{0L,12},
\frac{\kappa(1+z)}2,\frac\kappa2\right),$$
余项为$-2\kappa\mathcal R_{{\rm tree},12}$。它与生产分配相差$\kappa q_w$，已回代到完整K₁；第二分子不变。实际余项始终非零。

12的纵向下端x=0对应10纵向上端w′=1。直接在12下端验证
$$\alpha_2\leftrightarrow\alpha_3,\qquad
(y,w)\mapsto(1-w,1-y).$$
q不变、完整b₁与b₀L变号、实际标记权重保持，两个绝对Jacobian为1。共同η/对称α截断下它等于自身的负数，故纵向下端为零；**标量下端仍保留**。

12的纵向上端对应10纵向下端，所以其贡献是
$$z^{3\epsilon}\Lambda_{10}(\epsilon),$$
其中$\Lambda_{10}$是[10号§3](../3e-3g+Pi1-10/derivation.md#3-纵向上面严格为零下端闭合gamma)已完整积出的Gamma/Beta式。六个标量面和完整recoil则整体反射继承。余项先用Φ换到10的图表，再应用其$\omega=2\epsilon$的Euler正规化；没有声称12自身的规范图表具有同样的区域编号或算符。

10已经完成的13个边界/recoil系数积分及Gamma、原M/N、原/正规化余项数值控制按精确版本复用；本次没有把旧检查重新计为12的独立重积分。

## 4. 重新变换原J，再组装唯一I

沿用
$$R=\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}
{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad Z=z^{3\epsilon},$$
$$I_r=\frac{2R}{1+Z}[J_r(z)+ZJ_r(1/z)].$$
一次Euler补偿、C4重数及外置颜色/尺度均一致。J关系给
$$\begin{aligned}
I_{12}
&=\frac{2R}{1+Z}[ZJ_{10}(1/z)+Z\,Z^{-1}J_{10}(z)]\\
&=I_{10}.
\end{aligned}$$

原J不同，不能原样复制其表。若$J_{10}=\sum_{k=-3}^1j_{10,k}\epsilon^k+O(\epsilon^2)$，则
$$j_{12,n}(z)=\sum_{k=-3}^n
\frac{(3\log z)^{n-k}}{(n-k)!}j_{10,k}(1/z),\quad -3\le n\le1.$$
这五个原J系数已实际由HyperInt反演生成，回变换也精确恢复旧J；再按上式R/C4组装五个I系数，并独立核验I=10、crossing和原M/N回代。原J到ε¹的精度对应I有限项所需的R简单极点。

## 5. 全部极点与有限项的原M/N表达式

沿用
$$M_{k,j}(z)=\int_0^1
\frac{\log^k[\beta(1-\beta)]
[\log z+\log((1-\beta)/\beta)]^j}{\beta-z/(1+z)}\,d\beta,\quad j\ge1,$$
$$N_{k,j}=M_{k,j+1}-\log z\,M_{k,j},\qquad\zeta_n=\zeta(n).$$
内部点可去，无主值处方。完整结果
$$I_{12}=\sum_{q=-4}^0c_q(z)\epsilon^q+O(\epsilon),$$
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
+\frac{27M_{3,1}+5M_{1,3}}{2048}+\frac{7N_{1,2}}{512}\\
&-\frac{1923\zeta_2^2}{10240}
+\frac{23\zeta_2\zeta_3}{64}
-\frac{3\zeta_2}{64}-\frac{5\zeta_3}{128}
+\frac{\zeta_5}{96}-\frac7{192}.
\end{aligned}}$$
ε⁻⁶、ε⁻⁵均为零；仅原M/N与ζ常系数多项式，不新增特殊函数。五条Maple原M/N回代及独立Python shuffle回代均精确为零。

| z | c−2 | c−1 | c0 |
|---:|---:|---:|---:|
| 1 | 0.09551374128866734 | −0.7551977291479347 | 3.056639946351306 |
| 2 | 0.10248981292758194 | −0.7824685094664424 | 3.173575390538622 |

z=1/2与z=2相同。这些数值通过I₁₂=I₁₀继承旧15个原M/N控制，未做新12独立数值积分；旧最大差约1.63×10⁻³⁴。全εM/N闭包未由有限阶等式推断。

## 6. 两颜色与本分支收口

两个颜色部分为
$$CF\,CA^2\,\frac{10-6\epsilon}{3-2\epsilon}I_{12},\qquad
CF\,CA\,TF\,n_f\,\frac{-8(1-\epsilon)}{3-2\epsilon}I_{12}.$$
各五系数重新卷积、保存，与10同色结果及原M/N独立卷积完全相符，最高留数−5/1152、1/288。完整$\lambda_b^3(S_{14}+S_{12})$仍在外，尚未作全图物理组装/RG。

本次102完整桥接检查（含6调节核点）、82旧证据版本/范围核对、37独立系数/颜色控制通过。新Maple检查含5旧J绑定、5旧I绑定、5J回变换、5I等式、5crossing、5M/N及10同色等式。哈希只绑定版本，真正复用依据是完整原核、域、权重、分支与归一化的数学连接。

本分支**12/12全部极点至有限项已得**。已验证的归一化配对为01=06、02=05、03=04、07=08、10=12；这些是等值组，不是主积分数。12个源轨道仍分别保留C4重数。全ε标准函数闭式仍7/12，07/08/10/11/12尚未全ε压缩，formal=not_run。

完整w³数学口径为123已得、12零、14未完整交付（双TGV10、中央Γ₁4）。下一步把12个腿Π₁结果接入现有三圈全颜色裸项账本，核对源号、外尺度、C4/web计权和低阶反项精度；双TGV未解图与中央顶点来源仍须继续。当前没有Wilson数学MCP，无形式收据；未提交Git。
