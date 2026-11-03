# Π₂两图审核：动量源、C/J与全ε的共同生成结构

2026-10-06。本轮从三个已计算物质颜色通道中选择Π₂类审核：相邻01、对边02。它们贡献CF²TFnf和CFCA TFnf，也保留纯CFCA²。现行图集下，CF²TFnf的本阶裸贡献来自这两图，因此这里是验证该通道单极点／有限项输入的一个具体起点；本轮不把局部审核扩成其他图类已独立验证。

## 1. 从实际输出开始，明确审核到哪一层

实际qgraf-pipeline的3loop/orbit_definitions.json存在两条2e-Pi2记录：01在边3、4，02在边2、4；两侧传播子共享一个ℓ，中心是一个vacuum-pol-2loop有效1PI算符。总圈数3、耦合g₀⁶、两个Wilson参数均在[0,1]，顺序分别为s31→s41及s21→s41。相位是exp(+iℓ·Φ)，Φ=−X。

这个输出没有展开Π₂内部的全部原始两圈图。它把其标量交给文献算符。因此，本轮审核的是“有效Π₂输出＋有明确来源的全D算符”到C/J/I的链；内部两圈图的独立重新积分仍是另一项验证。

本轮实际运行安装的Montopar CLI，选2e-Pi2并输出到独立检查目录，重放两图及族包。结果与冻结的旧r001一致。该外部版本还保留旧模板的外部负号；项目通过2.0导入把C改号，J和标量算符保持，现行r002直接用于物理求和。不能把旧输出和现行I混用或再次加一个负号。

采用D=4−2ε、D₀=−ig/(k²+i0)、Π₂=−iδ(k²g−kk)𝒫₂及正号Dyson。真实动量起点是

$$\mathcal A_{ij}=\int_0^1ds\int_0^1dt\,
\frac1{N_c}\operatorname{Tr}[(ig_0v_i^{\mu_1}t^{a_1})(ig_0v_j^{\mu_2}t^{a_2})]
\int\frac{d^D\ell}{(2\pi)^D}\,
\frac{-ig_{\mu_1\nu_1}\delta^{a_1b_1}}{\ell^2+i0}
\Pi^{(2)\nu_1\nu_2}_{b_1b_2}(\ell)
\frac{-ig_{\nu_2\mu_2}\delta^{b_2a_2}}{\ell^2+i0}
e^{-i\ell X_{ij}}.$$

两点迹的顺序不引入三点fTrace的号差。其颜色收缩为CF，内部颜色仍需全部保留。

## 2. 内部算符来源及两个Gamma基函数

作者[Ruijl–Ueda–Vermaseren–Vogt的裸全D数据](https://arxiv.org/abs/1703.08532)以a²块列出Π₂。此轮直接解析已核对SHA256的官方gluonB.h，取作者Feynman规范xi=1，重新提取六个有理函数，逐个比较实际Montopar的A/B；没有读旧系数表作为预测值。

$$\mathcal P_2(q)=\frac{g_0^4}{(4\pi)^D}
\left[\frac{\eta(\epsilon)}{1-2\epsilon}\right]^2
q^{-2\epsilon}P_2(\epsilon),\qquad q=-\ell^2-i0,$$
$$\eta=\frac{\Gamma(1+\epsilon)\Gamma(1-\epsilon)^2}{\Gamma(1-2\epsilon)},\qquad
P_2=C_A^2p_{AA}+C_AT_Fn_fp_{Af}+C_FT_Fn_fp_{Ff},$$
$$p_c=A_c+B_cR,\qquad
R=\frac{\Gamma(1+2\epsilon)\Gamma(1-2\epsilon)^2}
{\Gamma(1+\epsilon)^2\Gamma(1-\epsilon)\Gamma(1-3\epsilon)},\qquad U=R^{-1}.$$

R又由嵌套一圈泡的Gamma公式重建。没有TF²nf²的1PI通道；双Π₁串联应另算。尤其CF²TFnf来自内部CF，不能因首坐标投影将它删除。

以Ff为例，实际重新提取的是

$$A_{Ff}=-\frac{4(\epsilon-1)(2\epsilon^2-\epsilon+2)}{\epsilon^3(2\epsilon-3)},\qquad
B_{Ff}=\frac{8(\epsilon-1)(2\epsilon-1)^2(\epsilon^2-2\epsilon+2)}
{\epsilon^3(2\epsilon-3)(3\epsilon-2)(3\epsilon-1)}.$$

完整A/B在新audit.json中保存，与[现行理论P2.4](theory.md)一致。本次采用作者的裸数据及其G-scheme到项目测度字典；未执行FORM/FORCER重新计算原始Π₂图。

## 3. 独立Fourier推导得到C和J

先保留完整横向张量。颜色之外的显式i乘积、ℓ²=−q、Fourier相位以及参考尺度比值均分别核对。母式为

$$\int\frac{d^D\ell}{(2\pi)^D}e^{-i\ell X}q^{-a}
=i\frac{\Gamma(D/2-a)}{4^a\pi^{D/2}\Gamma(a)}
(-X^2+i0)^{a-D/2}.$$

动量二次项是负的坐标二阶导数。对a=1+2ε和2+2ε，系数之比为−1/[12ε(1+2ε)]。令A₀=−2v14、z=v12/v14>0、f=X²/(2v14)、d=v_i·v_j/v14、u=v_i·X/v14、w=v_j·X/v14，最终

$$\mathcal K_2=
\frac{(1+4\epsilon)df+(1-3\epsilon)uw}{2(1+2\epsilon)}f^{-2+3\epsilon},\qquad
J_{ij}=\int_0^1ds\int_0^1dt\,\mathcal K_2.$$

本轮从实际Φ的四条边系数和Gram矩阵重建f,d,u,w，完整分子与两个实际输出相同。对每个颜色提出𝒯c=(CFCA²,CFCA TFnf,CF²TFnf)，raw前因子为

$$C_c^{raw}=\lambda_b^3\mathcal T_c(\pi A_0)^{3\epsilon}
\frac{\Gamma(1-\epsilon)^3}{32(1-2\epsilon)^2}
[A_cU+B_c],\qquad \lambda_b=\frac{g_0^2}{4\pi^2}.$$

prepared恰乘一次exp(−3γ_Eε)。raw没有C4；从四边旋转重新计数，m01=4、m02=2。定义

$$\mathcal R_* =\frac{G^3U}{128\epsilon^2(1-2\epsilon)^2},\qquad G=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),$$
$$I_o=\frac{m_o}{2}\mathcal R_*\frac{J_o(z)+z^{3\epsilon}J_o(1/z)}{1+z^{3\epsilon}}.$$

物理颜色系数为4ε²p_c I_o，不能把完整算符再乘一次，也不能遗漏U。C4、Euler和颜色恢复在这个公式中各计一次。

## 4. 任意横向标量插入的共同几何核

更一般地，若交换线上标量带q^{-hε}，令ν=(h+1)ε；Π₂和两个Π₁串联均有h=2。Fourier推导给出

$$\mathcal K_h=
\frac{(1+2h\epsilon)df+(1-\nu)uw}{2(1+h\epsilon)}f^{\nu-2}.$$

因为f_s=u、f_t=−w、f_st=−d，完整核严格分解为

$$\boxed{\mathcal K_h=d f^{\nu-1}
+\frac{\partial_s\partial_t f^\nu}{2\nu(1+h\epsilon)}.}$$

这个恒等式保留全部边界，说明内部算符只改变标量前因子和ν，几何不必为每个颜色重新求解。它不把Π₂等同于Π₁²。

相邻图f=z(1−s)t，故

$$J_{01}=\frac{(2+\epsilon)z^{3\epsilon}}{18\epsilon^2(1+2\epsilon)},\qquad
\boxed{I_{01}=\frac{G^3U(2+\epsilon)}{1152\epsilon^4(1-2\epsilon)^2(1+2\epsilon)}.}$$

对边令F=(1−s)(1−t)+zst，四角值依次为1、0、0、z。初始Reν>0下，混合导数的积分是1+z^ν；它不能因两个零角点而全部删除。定义唯一标量母积分

$$H_\nu(z)=\int_0^1ds\int_0^1dt\,F^{\nu-1},\qquad
H_0(\nu,z)=-(1+z)z^{-\nu/2}H_\nu(z).$$

于是ν=3ε时

$$J_{02}=z^{\nu/2}H_0(\nu,z)+\frac{1+z^\nu}{2\nu(1+2\epsilon)},$$
$$\boxed{I_{02}=\mathcal R_*\left[
\frac{H_0(3\epsilon,z)}{\cosh(3\epsilon\log z/2)}
+\frac1{3\epsilon(1+2\epsilon)}\right].}$$

原对边J满足J(z)=z^νJ(1/z)；m=2和两个crossing半项共同产生上式。忘掉其中一个半项会使I少一倍。

## 5. 原2ne的M矩自然生成全部ε阶

先将线性的s精确积掉：

$$H_\nu(z)=\frac1\nu\int_0^1dt\,
\frac{z^\nu t^\nu-(1-t)^\nu}{(1+z)t-1}.$$

内部零点处两项分子严格相消。令β=1−t、r=z/(1+z)，定义w=log[β(1−β)]、Q=log z+log[(1−β)/β]，得到

$$\boxed{H_0(\nu,z)=\int_0^1d\beta\,
\frac{e^{\nu w/2}\,2\sinh(\nu Q/2)}{\nu(\beta-r)}.}$$

Q(r)=0，所以分母的内部零点可去。这正好接到原2ne的

$$M_{k,j}(z)=\int_0^1\frac{w^k Q^j}{\beta-r}\,d\beta.$$

展开指数和sinh直接给出任意阶有限和

$$\boxed{[\nu^n]H_0(\nu,z)=
\sum_{j=0}^{\lfloor n/2\rfloor}
\frac{M_{n-2j,2j+1}(z)}{2^n(n-2j)!(2j+1)!}.}$$

特别地，H₀=M₀₁+νM₁₁/2+ν²(M₂₁/8+M₀₃/24)+…。这里的H₀(ν,z)是二变量函数，H₀(0,z)=M₀₁，而不是把新矩重新命名。对固定z>0，原标量母积分在Reν>−1收敛且局部解析；对ν=0附近的紧邻域，端点有统一可积幂／对数界，内部点可去。因此上述逐阶生成在该局部域有解析依据，再按同一亚纯函数延拓。

β→1−β、z→1/z使Q和β−r同时翻号，w不变，故H₀(ν,1/z)=H₀(ν,z)。展开只含奇数第二指标的M，直接解释全阶crossing。这比仅比较有限阶解表多提供一个固定原积分的生成结构。

## 6. 全ε闭式：Mellin–Barnes、Meijer G与关联Gauss支路

对F的两项使用Mellin–Barnes恒等式，再把s,t分别积为Beta函数，得到

$$H_\nu(z)=\frac1{\Gamma(1-\nu)\Gamma(1+\nu)^2}
\int_{c-i\infty}^{c+i\infty}\frac{du}{2\pi i}\,
\Gamma(-u)\Gamma(1-\nu+u)\Gamma(1+u)^2\Gamma(\nu-u)^2z^u.$$

围道须满足max(−1,Reν−1)<c<min(0,Reν)。初始共同域取z>0、0<Reν<1，与两个J端点及四角边界一致。Gamma核正是

$$H_\nu(z)=\frac{G^{3,3}_{3,3}\!\left(z\mid {0,0,\nu\atop0,\nu,\nu}\right)}
{\Gamma(1-\nu)\Gamma(1+\nu)^2}.$$

由此独立匹配现行r002的完整Meijer G参数和Gamma归一化，而不是通过已知Laurent值拟合一个函数。

线性积分还给出更直接的Gauss表示

$$H_\nu(z)=-\frac1{\nu(\nu+1)}\left[
z^\nu\,{}_2F_1(1,\nu+1;\nu+2;1+z+i0)
+\frac1z\,{}_2F_1(1,\nu+1;\nu+2;(1+z)/z-i0)\right].$$

两个支路必须关联；原分式只有一个可去的内部点，换元使两项的边界侧相反。不能把它们都放在同一侧。相反的关联选择给同一个完整组合；ν=0等中间式可去点取完整组合的极限。复运动学仍需沿原处方联合延拓，本轮只核验正实z。

## 7. 本次实际核验和结论范围

113项精确控制覆盖实际输出、官方算符A/B、G10泡比、i/Fourier、两图Gram和完整J、C、一次Euler、C4、全εI、原M生成的14个几何Laurent位置以及42个全颜色位置。没有用RG来决定这些系数。

在ε=0.05、0.11、0.20和z=1/2、1、2，45组控制比较原二维J的四三角Jacobi求积、Gauss闭式、一重MB、未展开生成密度、C4/I及crossing。3个同侧Gauss支路负控制都留下虚部。数值精度加深是非认证控制，不是严格区间证明。

结论：在现行有效Π₂源、作者裸全D算符和项目2.0边界值约定下，未发现需修改当前C/J、I或颜色系数的确定错误。已重建与原结果相同的全ε解，并找到它与原2ne的共同M生成结构。该结论不替代内部原始两圈子图的独立再积分、图集完备性、严格全域去调节或一般复分支证明。formal=not_run，无Wilson MCP形式收据；原输入／结果、上游和完成计数未改，未提交Git。

[实际计算及复现](../../../workflow/checks/pi2-chain-audit/README.md) · [01逐图推导](workspaces/2e-Pi2-01/derivation.md) · [02逐图推导](workspaces/2e-Pi2-02/derivation.md)。
