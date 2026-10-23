# 3e-3g+Pi1-09：独立边自能的全ε闭式与M系数

## 当前物理积分的入口（2026-10-15，r003）

从共同作用量出发，U′=igAU确定真实Wilson词；FT/Dyson规则在提取C前生效。当前[原生C/J](../../sources/montopar-r002.json)给出下式，J就是本文实际积分的原核，域、分支和IBP全保留：

$$R_P=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad I_P=2R_P\frac{J(z)+z^{3\epsilon}J(1/z)}{1+z^{3\epsilon}}.$$

用下文已经求出的J代入此定义，得到[当前物理I](../../results/3e-3g+Pi1-09/r003.json)；组装系数为+1，Euler只补一次。以下带日期的求值过程与r001公式属于旧根历史，计算J的步骤可按同域积分线性复用，旧形式状态不继承。

2026-10-15表示修订：主I、原J和主Laurent系数保持r002的值；当前源哈希、配对引用及附属公式统一到新修订。[修订说明](../../revisions/r003.md)。以下2026-10-01/02的r001系数表保留历史符号，不能直接当作当前结果。

2026-10-02。09为244对边核心，Π₁在独立边2、第1腿。本图从完整raw核与五个IBP面独立求得**Gamma与Gauss $\,{}_2F_1$ 的全ε闭式**；极点至有限项均可用原M、ζ和乘积表示，仅常系数。状态exact_closed_form，formal=not_run。

[结果r001](../../results/3e-3g+Pi1-09/r001.json) · [实际核验](../../../../../workflow/checks/tgv-pi1-09/README.md)。输入为不可变 [raw/IBP包](../../sources/montopar-r001.json)，SHA256=cc9a45b9a3aecd1437fea5fbe222680f67001f8137c09cb38f46ba6fff337e90。默认47数学根未扩展，生产/Kernel/Lean未改。

## 1. 独立腿标记与α尺度

记 $s=s_{21}$、$t=s_{41}$、$u=s_{42}$。原顺序为 $s,u,t,\alpha_1,\alpha_2$，$0\le t\le u\le1$，α在单纯形中。实际权重与分母为
$$\alpha_1^{-2\epsilon}\alpha_2^{-\epsilon}\alpha_3^{-\epsilon},
\qquad q=\alpha_1[\alpha_2D_z(s,t)+\alpha_3D_z(s,u)],$$
$$D_z(s,b)=(1-s)(1-b)+zsb.$$

取 $\alpha_1=a,\alpha_2=(1-a)v,\alpha_3=(1-a)(1-v)$，Jacobian为 $1-a$。设 $\nu=1-3\epsilon$、$\delta=3\epsilon-1$、$m=2-3\epsilon$。从通量两幂提出 $a(1-a)$ 后，共同测度为
$$a^{\epsilon-1}(1-a)^\epsilon\,da\;w(v)\,dv,\qquad
w(v)=v^{-\epsilon}(1-v)^{-\epsilon}.$$

尺度积分是 $B_a=\mathrm B(\epsilon,1+\epsilon)$，线性分子的a平均为
$\langle p_0+p_1a\rangle=p_0+\epsilon p_1/(1+2\epsilon)$。这里使用真实α₁标记权重，本图没有套用07的尺度结果。

原raw两分子重新核验为
$K_1=\sum_i(\partial_i b_{1,i}-b_{0,i}\partial_iq)$、
$K_2=-(2-3\epsilon)\sum_i b_{1,i}\partial_iq$，
且 $\sum_i\partial_i b_{0,i}=0$、recoil=0。因此完整核是原有序域的五面通量和。积a后，每面仍保留
$B_a\int w(v)[N_0D^{-\nu}+N_1D^{-m}]$：

| 原面 | D | N₀ | N₁ |
|---|---|---|---|
| $s=0$ | $1-vt-(1-v)u$ | 0 | $d(2v-1)$ |
| $s=1$ | $z[vt+(1-v)u]$ | 0 | $-z^2d(2v-1)$ |
| $t=0$ | $D_z(s,(1-v)u)$ | $-A(1+z)$ | 0 |
| $u=1$ | $D_z(s,1-v(1-t))$ | $-A(1+z)$ | 0 |
| $t=u$ | $D_z(s,u)$ | $2A(1+z)$ | 0 |

$$d=\frac{1+\epsilon}{4(1+2\epsilon)},\qquad
A=\frac{1+\epsilon}{2\delta}.$$

同边接触面始终保留。五面几何、两分子、a次数和源面Jacobian均逐项核验。

先保留共同 $i\eta$ 和α内部截断作IBP。取 $z>0$、$0<\Re\epsilon<1/3$：a端点为 $a^{\Re\epsilon-1}(1-a)^{\Re\epsilon}$；三个N₀面的s积分H有界；两个N₁端点面的同边尺度要求 $3\Re\epsilon>0$，剩余v角点要求 $2\Re\epsilon>0$，另一端要求 $\Re\epsilon<1$。五面在共同条带分别绝对可积，调节幂由正q幂控制。因此可取共同调节极限，并在该条带按实际可积条件换序，再按ε亚纯延拓。

## 2. 两个s端点面全部积完

在 $s=0$ 面取 $t=1-x,u=1-xp$，$0\le x,p\le1$，Jacobian为x，
$D=x[v+(1-v)p]$。积x、p得到
$$Q(\epsilon)=\frac{d}{3\epsilon\delta}
\int_0^1dv\,\frac{w(v)(2v-1)}{1-v}[1-v^\delta].$$

v=1处差分消零，v=0由 $v^{2\epsilon-1}$ 控制。共同辅助端点指数下拆成Beta，再延拓：
$$Q(\epsilon)=\frac{1+\epsilon}{12\epsilon(1+2\epsilon)\delta}
\left[\frac{\mathrm B(1-\epsilon,-\epsilon)}{1-2\epsilon}
-3\mathrm B(2\epsilon,-\epsilon)\right].
\tag{1}$$

负Beta参数代表同一收敛差分的解析延拓。具体先给 $(1-v)$ 加共同指数σ，使两项分别收敛；差分在 $\Re\sigma>\Re\epsilon-1$ 解析，最后令σ→0。不能将各分项当作正ε条带中的独立收敛Euler积分。

经 $t\to1-u,u\to1-t,v\to1-v$，上端点面成为 $z^{3\epsilon}Q$；分子负号与 $(2v-1)$ 的反射共同变化。两面贡献
$J_{\rm end}=B_a(1+z^{3\epsilon})Q$。

## 3. 不完全Beta余项为何消失

令
$$H(b)=\int_0^1ds\,D_z(s,b)^{-\nu}
=\frac{(zb)^{3\epsilon}-(1-b)^{3\epsilon}}
{3\epsilon[(1+z)b-1]},\qquad B_c=\mathrm B(1-\epsilon,1-\epsilon).$$

内部点 $b=1/(1+z)$ 可去。下、上面换 $b=(1-v)u$、$b=1-v(1-t)$，再加接触面，得
$$J_{\rm side}=B_a A(1+z)\int_0^1db\,H(b)\rho(b,\epsilon),$$
$$\rho=2B_c-L(b)-L(1-b),\qquad
L(b)=\int_0^b\frac{v^{-\epsilon}(1-v)^{-\epsilon}}{1-v}\,dv.$$

本图w对v↔1−v对称，由此
$$\boxed{\rho(b,\epsilon)=
\frac{B_c-[b(1-b)]^{-\epsilon}}{\epsilon}.}
\tag{2}$$

证明：设 $W=\log[b(1-b)]$，两边b导数均为 $e^{-\epsilon W}W'$。在 $\Re\epsilon<0$ 取端点确定常数，得
$L(b)+L(1-b)=[b(1-b)]^{-\epsilon}/\epsilon+\mathrm B(1-\epsilon,-\epsilon)$；
用 $\mathrm B(1-\epsilon,-\epsilon)=-(1-2\epsilon)B_c/\epsilon$。固定 $0<b<1$ 时按ε解析延拓到起始条带，ε=0的表观极点可去。因此没有未知不完全Beta余项，权重全部回到原M母源。

## 4. 全εI与Gauss闭式

原归一化为
$$R=\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}
{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad
I=\frac{2R[J(z)+z^{3\epsilon}J(1/z)]}{1+z^{3\epsilon}}.$$

(2)的密度对b↔1−b对称，结合H反射和端点面得
$J_{09}(z)=z^{3\epsilon}J_{09}(1/z)$，
所以 $I_{09}=4RJ_{09}/(1+z^{3\epsilon})$。一次三圈Euler补偿、C4重数4沿用原约定。

以下使用原M的变量β，记 $\rho_z=z/(1+z)$、$\ell=\log z$、
$W=\log[\beta(1-\beta)]$、$L=\ell+\log[(1-\beta)/\beta]$。已存在的同一母源为
$$\mathcal U(a,b;z)=\int_0^1
\frac{e^{aW}(e^{bL}-1)}{\beta-\rho_z}\,d\beta,\qquad
\mathcal V_{p,q}=
\frac{\mathcal U(p\epsilon/2,q\epsilon/2)-\mathcal U(p\epsilon/2,-q\epsilon/2)}
{q\epsilon}.$$

$\mathcal V_{3,3}=\mathcal V_M(3\epsilon,z)$；
$\mathcal V_{1,3}$是同一U在另两个参数上的组合，不是新函数族。

设
$$G_A=e^{-3\gamma_E\epsilon}
\frac{\Gamma(1-\epsilon)^4\Gamma(1-3\epsilon)\Gamma(1+\epsilon)^2}
{\Gamma(1-2\epsilon)^2\Gamma(1+2\epsilon)},$$
$$G_B=e^{-3\gamma_E\epsilon}
\frac{\Gamma(1-\epsilon)^3\Gamma(1-3\epsilon)\Gamma(1+\epsilon)}
{\Gamma(1-2\epsilon)},$$
$$G_S=e^{-3\gamma_E\epsilon}
\frac{\Gamma(1-\epsilon)^2\Gamma(1-3\epsilon)\Gamma(1+\epsilon)^2}
{\Gamma(1-2\epsilon)\Gamma(1+2\epsilon)}.$$

完整解是
$$\boxed{\begin{aligned}
I_{09}={}&
\frac{G_A}{192\epsilon^4(1+2\epsilon)(1-2\epsilon)^2}
-\frac{G_B}{128\epsilon^4(1-4\epsilon^2)}\\
&+\frac{G_S}{64\epsilon^3(1-2\epsilon)\cosh(3\epsilon\ell/2)}
[B_c\mathcal V_{3,3}-\mathcal V_{1,3}].
\end{aligned}}
\tag{3}$$

Gamma递推精确核验：(3)前两项为 $4RB_aQ$，最后一项前因子由原侧面式直接得到。C、I或图的外部符号没有改动。

U也已全部写成标准函数。两Gauss项须取同一边界 $\xi=1/\rho_z+i0$：
$$\boxed{\begin{aligned}
\mathcal U(a,b;z)=-\frac1{\rho_z}\big[&
z^b\mathrm B(1+a-b,1+a+b)
{}_2F_1(1,1+a-b;2+2a;\xi)\\
&-\mathrm B(1+a,1+a)
{}_2F_1(1,1+a;2+2a;\xi)\big].
\end{aligned}}
\tag{4}$$

先在 $|\xi|<1$ 将 $(1-\xi\beta)^{-1}$ 展成几何级数，逐项Beta积分产生Gauss级数，再沿上半平面继续到所需边界，即得(4)。两项分开有内部Cauchy极点，但它们的分子在 $\beta=\rho_z$ 严格相等，因为 $e^{bL(\rho_z)}=1$。所以同一边界上的残数抵消，合式正是原可去点积分；不得分开选不同支路，定义没有用取实部来遮去虚部。共同−i0给出相同合式。

(3)–(4)仅含Gamma、Beta、Gauss函数，无剩余参数积分。12组参数核验了U的两侧共同边界极限，三个非零ε点核验了完整I的Gauss式与原Euler表示。

## 5. 全部极点与有限项

任意阶展开直接由原M定义给出：
$$\mathcal V_{p,q}=\sum_{n\ge0}\epsilon^n
\sum_{j=0}^{\lfloor n/2\rfloor}
\frac{p^{n-2j}q^{2j}}{2^n(n-2j)!(2j+1)!}
M_{n-2j,2j+1}(z),$$
$$M_{k,j}(z)=\int_0^1\frac{W^kL^j}{\beta-\rho_z}\,d\beta,\qquad j\ge1.$$

内点可去；固定z>0时，足够小的源参数邻域由可积端点幂与对数支配，可逐项展开。Gamma按亚纯函数展开，双曲因子在ε=0非零。因此任意固定阶只需有限原M矩与有限卷积，不会逐阶新增未知余项。

$B_c=1+2\epsilon+(4-\zeta_2)\epsilon^2+(8-2\zeta_2-2\zeta_3)\epsilon^3+\cdots$；
差 $B_c\mathcal V_{3,3}-\mathcal V_{1,3}$ 从ε¹开始，有限项需它到ε³。
用 $\ell^2=-2M_{0,1}-6\zeta_2$、$M_{0,3}=-M_{0,1}^2$：
$$I_{09}=\sum_{n=-4}^{0}c_n(z)\epsilon^n+O(\epsilon),$$
$$c_{-4}=-\frac1{384},\qquad c_{-3}=\frac1{96},$$
$$c_{-2}=\frac1{96}-\frac{7\zeta_2}{256}
+\frac{M_{0,1}}{32}+\frac{M_{1,1}}{64},$$
$$c_{-1}=\frac1{12}+\frac{\zeta_2}{64}-\frac{7\zeta_3}{384}
+\frac{(8-\zeta_2)M_{0,1}+5M_{1,1}+M_{2,1}}{64},$$
$$\boxed{\begin{aligned}
c_0={}&\frac18-\frac{5\zeta_2}{64}+\frac{7\zeta_3}{96}
-\frac{547\zeta_2^2}{5120}\\
&+\frac{3M_{0,1}}8+\frac{29\zeta_2M_{0,1}}{128}
-\frac{\zeta_3M_{0,1}}{32}
+\frac{M_{1,1}}4+\frac{31\zeta_2M_{1,1}}{256}\\
&+\frac{17M_{2,1}}{256}
+\frac{15M_{0,1}^2+9M_{0,1}M_{1,1}}{256}
+\frac{13M_{3,1}}{1536}+\frac{3M_{1,3}}{512}.
\end{aligned}}$$

ε⁻⁶、ε⁻⁵为零。所有M在同一z取值，系数仅有理常数及ζ。09独立求解，最高四阶与重复边插入的最高三阶不同。

| z | c₋₂ | c₋₁ | c₀ |
|---:|---:|---:|---:|
| 1 | −0.0161450738514280 | −0.0198760656618505 | 0.243763969545486 |
| 2 | −0.0143507233420138 | −0.0266267850832843 | 0.259591006167148 |

SymPy用Gamma对数递推，Maple另从未简化 $4RB_aQ$ 与原母源规则展开，五项精确相符；五crossing、五原J归一化回代通过。既有M的Hlog坐标按版本复用仅用于核验，不是数值拟合的来源。

## 6. 两颜色与完成范围

源C仍为 $C_FC_A B_1$，
$B_1=[(10-6\epsilon)C_A-8(1-\epsilon)T_Fn_f]/(3-2\epsilon)$。
两通道为
$$g=\frac{10-6\epsilon}{3-2\epsilon}I_{09},\qquad
f=-\frac{8(1-\epsilon)}{3-2\epsilon}I_{09}.$$

各五系数、共十个颜色系数保存并以独立shuffle核验。最高留数 $g_{-4}=-5/576$、$f_{-4}=1/144$。
$\lambda_b^3(S_{14}+S_{12})$ 留在外面，这里的有限项不是完成重整化与尺度组装后的物理有限项。

实际检查：54项源/面/尺度/密度连接、4项Gamma等式、15项Maple系数/crossing/J回代；9密度点、3个Q Euler点、12组U分支点、6个全ε/Laurent余项点，以及颜色、直接M、完整Gauss与原三面控制。三面数值在精确a、s积分后进行，另外两端点面用独立Euler核验的Q；每轴64→512节点，末差约6.1×10⁻⁸与9.2×10⁻⁸。未实施独立五维积分或严格区间证明。

本分支有限阶9/12，全ε标准函数闭式7/12（六Gamma式加本图Gamma/Gauss）。244核心07–09全部完成，下一10进入234三边核心、六边界加recoil。Wilson数学MCP不可用，无形式收据，formal=not_run。原输入、01–08结果保持，完整物理总和/RG仍开放。
