# 从原始动量积分推导单个 Π₁ 的 montopar 参数化

2026-09-23。本文从用户指定的相邻边动量积分逐步推导完整坐标核和 C×J，并说明同一个局部核怎样用于两端点及四端点的单 Π₁ 图。它是当前输入约定下的实现依据。[r001](theory-r001.md)与已发布 JSON 保留原字节；本文补全推导和符号定义，没有重新定义它们。

本文的 I 始终指下式给定的对象。它不是 Wilson Kernel 中已经完成颜色投影、Euler 准备或轨道组装的最终 I。

## 1. 固定动量起点及各个符号

考虑四边形上第3、4条边的一次自能交换：

$$
\boxed{
\begin{aligned}
I_{34}={}&-\int_0^1 ds_{3,1}\int_0^1 ds_{4,1}\,
(igv_3^{\mu_1}t^{a_1})(igv_4^{\mu_2}t^{a_2})\\
&\times\int\frac{d^D\ell}{(2\pi)^D}
\frac{-ig_{\mu_1\nu_1}\delta^{a_1b_1}}{\ell^2+i0}
\Pi_{b_1b_2}^{\nu_1\nu_2}(\ell)
\frac{-ig_{\nu_2\mu_2}\delta^{b_2a_2}}{\ell^2+i0}
e^{-i\ell\cdot X_{34}}.
\end{aligned}}
\tag{1}
$$

最外层的负号属于式(1)的定义，不能因习惯使用另一种自能规则而再次添加或删除。Π内部已经包含胶子、鬼及无质量夸克贡献；闭合鬼圈、费米子圈的符号不能在外面重复乘入。

采用 Minkowski 度规 `(+−−−)`，D=4−2ε，并保持 D 维张量运算。Wilson 四边形满足

$$
v_i^2=0,\qquad \sum_{i=1}^4v_i=0,\qquad v_{ij}:=v_i\cdot v_j.
\tag{2}
$$

这里 v_ij 是内积本身，不是两倍内积。围道参数化及相位为

$$
\begin{aligned}
x_3(s)&=v_1+v_2+s v_3,\\
x_4(t)&=v_1+v_2+v_3+t v_4,\\
X_{34}&=x_3(s)-x_4(t)=-[(1-s)v_3+t v_4],
\end{aligned}
\qquad s=s_{3,1},\quad t=s_{4,1}.
\tag{3}
$$

因此源 JSON 中保存的 Φ=(1−s)v₃+tv₄ 与 X 的关系为 Φ=−X，`exp(+iℓ·Φ)`和式(1)的指数相同。整条线的动量路由反向会同时反转 X；距离和两个切向量与 X 的收缩乘积都不变，这不能产生一个额外整体负号。

式(1)保留了颜色矩阵。下面取其单位矩阵系数，等价于施加归一化迹 Tr/Nc；由 tᵃtᵃ=CF·1 得到 CF，不额外产生 Nc 或 1/Nc。

## 2. Π₁ 的完整定义：显式 i 已包含在算符中

沿用用户此前指定、并已核对来源的裸一圈自能：

$$
\begin{aligned}
\Pi_{b_1b_2}^{\nu_1\nu_2}(\ell)
&=\delta_{b_1b_2}(\ell^2g^{\nu_1\nu_2}-\ell^{\nu_1}\ell^{\nu_2})
\Pi_1(\ell^2),\\
\Pi_1(\ell^2)&=i\mathcal A_1(\varepsilon)(-\ell^2-i0)^{-\varepsilon},\\
\mathcal A_1(\varepsilon)
&=\frac{g^2}{2(4\pi)^{2-\varepsilon}}
\frac{\Gamma(\varepsilon)\Gamma^2(1-\varepsilon)}{\Gamma(2-2\varepsilon)}
\mathcal B_1(\varepsilon),\\
\mathcal B_1(\varepsilon)
&=\frac{(10-6\varepsilon)C_A-8(1-\varepsilon)T_Fn_f}{3-2\varepsilon}.
\end{aligned}
\tag{4}
$$

TF=1/2，nf为无质量夸克味数，故

$$
\mathcal B_1=\frac{(3D-2)C_A-2(D-2)n_f}{D-1}\qquad(T_F=1/2).
\tag{5}
$$

式(4)包含一个 g² 和一个 i。A₁只是去掉显式 i 和动量幂后的系数，不能把它当成不含耦合的量。相邻图连同两个 Wilson 顶点总阶数为 g⁴，即总二圈。

文献来源为 Davydychev、Osland、Tarasov，[hep-ph/9605348](https://arxiv.org/pdf/hep-ph/9605348)附录C的(C.1)、(C.3)，其规范参数ξ_D=0对应Feynman规范。与逆传播子修正文献记号的字典及已完成的跨来源核对见[文献说明](../../../workflow/self-energy-literature.md)。本文直接以式(4)定义插入，不以不同论文中相同的Π字母替代这个字典。

这是裸算符；未加入反项。g²=4παs仅作记号转换。本文没有隐含引入μ或MSbar因子；若改用带显式μ的耦合约定，必须从式(1)、(4)一起改写。

完整 `(ℓ²g−ℓℓ)` 张量是横向的。下面将其中 ℓ²g 项称为度规项，将 −ℓℓ 变换产生的项称为导数项；保留这两项，才是式(1)定义的完整单图。

## 3. 先收缩颜色、Lorentz 指标和显式相位

三个颜色δ桥连接两个Wilson生成元：

$$
t^{a_1}t^{a_2}\delta^{a_1b_1}\delta_{b_1b_2}\delta^{b_2a_2}
=C_F\mathbf1.
\tag{6}
$$

Lorentz 收缩给出

$$
v_3^{\mu_1}g_{\mu_1\nu_1}
(\ell^2g^{\nu_1\nu_2}-\ell^{\nu_1}\ell^{\nu_2})
g_{\nu_2\mu_2}v_4^{\mu_2}
=\ell^2v_{34}-(v_3\cdot\ell)(v_4\cdot\ell).
\tag{7}
$$

在代入Π所含的 i 之前，其余相位已经确定：

$$
(-1)_{\rm outer}(ig)^2(-i)^2=-g^2.
\tag{8}
$$

令 q=−ℓ²−i0，代入式(4)：

$$
I_{34}=-ig^2C_F\mathcal A_1
\int_0^1ds\int_0^1dt\int\frac{d^D\ell}{(2\pi)^D}
e^{-i\ell\cdot X}
\frac{\ell^2v_{34}-(v_3\cdot\ell)(v_4\cdot\ell)}{(\ell^2+i0)^2}
q^{-\varepsilon}.
\tag{9}
$$

以统一的Feynman边界值解释传播子幂，式(9)成为

$$
I_{34}=ig^2C_F\mathcal A_1
\int_0^1ds\int_0^1dt\int\frac{d^D\ell}{(2\pi)^D}
e^{-i\ell\cdot X}
\left[v_{34}q^{-1-\varepsilon}
+(v_3\cdot\ell)(v_4\cdot\ell)q^{-2-\varepsilon}\right].
\tag{10}
$$

度规项的负号来自ℓ²向−q的转换，导数项原来的分子负号一起被提出，故方括号里是加号。这里和下文的幂核恒等式均按共同Feynman边界值及维数解析延拓理解，并非在有限η时把ℓ²+iη与ℓ²当成完全相同的数。

## 4. Fourier 母式与两个幂核之间的相对负号

定义不含传播子分子和额外 i 的 Fourier 积分

$$
\mathcal F_a(X):=\int\frac{d^D\ell}{(2\pi)^D}
e^{-i\ell\cdot X}(-\ell^2-i0)^{-a}.
\tag{11}
$$

在可作Wick转动的起始区域，时间动量测度给出 i；欧氏Schwinger表示给出

$$
\begin{aligned}
\int\frac{d^Dk}{(2\pi)^D}\frac{e^{ik\cdot x_E}}{(k^2)^a}
&=\frac1{(4\pi)^{D/2}\Gamma(a)}
\int_0^\infty d\tau\,\tau^{a-1-D/2}e^{-x_E^2/(4\tau)}\\
&=\frac{\Gamma(D/2-a)}{4^a\pi^{D/2}\Gamma(a)}(x_E^2)^{a-D/2}.
\end{aligned}
\tag{12}
$$

先在Schwinger积分收敛的参数域建立式(12)，再作解析延拓，得到所需母式

$$
\boxed{\mathcal F_a(X)=i h_a Q^{a-D/2},\quad
h_a=\frac{\Gamma(D/2-a)}{4^a\pi^{D/2}\Gamma(a)},\quad Q=-X^2+i0.}
\tag{13}
$$

两个Fourier指数为1+ε和2+ε。记h₁=h_(1+ε)、h₂=h_(2+ε)，则

$$
\begin{aligned}
\mathcal F_{1+\varepsilon}&=ih_1Q^{-1+2\varepsilon},
&h_1&=\frac{\Gamma(1-2\varepsilon)}{4^{1+\varepsilon}\pi^{2-\varepsilon}\Gamma(1+\varepsilon)},\\
\mathcal F_{2+\varepsilon}&=ih_2Q^{2\varepsilon},
&h_2&=\frac{\Gamma(-2\varepsilon)}{4^{2+\varepsilon}\pi^{2-\varepsilon}\Gamma(2+\varepsilon)},\\
\frac{h_2}{h_1}&=-\frac1{8\varepsilon(1+\varepsilon)}.
\end{aligned}
\tag{14}
$$

末行只用了Γ(1−2ε)=(−2ε)Γ(−2ε)及Γ(2+ε)=(1+ε)Γ(1+ε)。特别地，这个负号不能用Gamma绝对值或独立的分支选择消掉。

不能要求两个h同时来自一个未经调节、绝对收敛的Fourier积分条带：这里包括母核的解析延拓。例如 `0<Re ε<1/2` 时相邻边积分收敛，但F_(2+ε)仍按式(13)的亚纯族定义。坐标微分应作用于该同一族。

由于

$$
(v_i\cdot\partial_X)(v_j\cdot\partial_X)e^{-i\ell\cdot X}
=-(v_i\cdot\ell)(v_j\cdot\ell)e^{-i\ell\cdot X},
\tag{15}
$$

式(10)精确变为

$$
I_{34}=ig^2C_F\mathcal A_1\int_0^1ds\int_0^1dt
\left[v_{34}\mathcal F_{1+\varepsilon}
-(v_3\cdot\partial_X)(v_4\cdot\partial_X)\mathcal F_{2+\varepsilon}\right].
\tag{16}
$$

其中第二项的减号来自两个动量变为两个坐标导数，不是从论文结果反推的符号。

## 5. 显式完成张量求导

引入便于剥离Fourier公共相位的坐标张量

$$
\begin{aligned}
\mathcal H_{\mu\nu}(X)
&:=\frac1i\left[g_{\mu\nu}\mathcal F_{1+\varepsilon}
-\partial_\mu\partial_\nu\mathcal F_{2+\varepsilon}\right]\\
&=h_1\left[g_{\mu\nu}Q^{-1+2\varepsilon}
+\frac{\partial_\mu\partial_\nu Q^{2\varepsilon}}
{8\varepsilon(1+\varepsilon)}\right].
\end{aligned}
\tag{17}
$$

注意H不含Π的A₁，也不含Wilson顶点的g²。式(16)中的两个 i 相乘，因此

$$
\boxed{I_{34}=-g^2C_F\mathcal A_1
\int_0^1ds\int_0^1dt\,v_3^\mu v_4^\nu\mathcal H_{\mu\nu}(X).}
\tag{18}
$$

从∂μQ=−2Xμ直接得到

$$
\partial_\mu\partial_\nu Q^{2\varepsilon}
=-4\varepsilon g_{\mu\nu}Q^{-1+2\varepsilon}
+8\varepsilon(2\varepsilon-1)X_\mu X_\nu Q^{-2+2\varepsilon}.
\tag{19}
$$

代入即得可直接实现的完整张量：

$$
\boxed{\frac{\mathcal H_{\mu\nu}}{h_1}
=\frac{1+2\varepsilon}{2(1+\varepsilon)}g_{\mu\nu}Q^{-1+2\varepsilon}
-\frac{1-2\varepsilon}{1+\varepsilon}X_\mu X_\nu Q^{-2+2\varepsilon}.}
\tag{20}
$$

两个快速检验来自同一个母式：□F_(2+ε)=F_(1+ε)，故∂μH_(μν)=0。亦可直接对式(20)求散度；XμXν项必须采用D=4−2ε的迹，才能与度规项相消。这检验相对符号，却不能单独确定整个图的整体符号。

式(19)、(20)在离奇点处按普通微分理解，在光锥及端点按统一边界值族理解。不能先把不同项的i0分别改号，再把它们相加。

## 6. montopar 的统一几何与原 s 核

为与现有四边形代码一致，采用

$$
r=v_{14},\qquad z=\frac{v_{12}}{v_{14}}=\frac{v_{34}}{v_{14}},
\qquad A_{\rm ref}=-2r.
\tag{21}
$$

起始区域取r<0、z>0。对任意一条不同边间的SE线，令X=x_i(s_i)−x_j(s_j)，定义

$$
f=\frac{X^2}{2r},\qquad d=\frac{v_i\cdot v_j}{r},\qquad
u=\frac{v_i\cdot X}{r},\qquad w=\frac{v_j\cdot X}{r}.
\tag{22}
$$

所有几何量都可由四维“边向量系数”及无量纲Gram表计算；这里的四个系数不是把时空维数D设为4。Gram表为

$$
G_{ij}=\frac{v_i\cdot v_j}{r}=
\begin{pmatrix}
0&z&-(1+z)&1\\ z&0&1&-(1+z)\\
-(1+z)&1&0&z\\ 1&-(1+z)&z&0
\end{pmatrix}.
\tag{23}
$$

若X=Σk ck vk，便有f=(cᵀGc)/2、u=(Gc)i、w=(Gc)j、d=Gij。原端点定义决定ck；同边多个端点不能仅凭边号合并。

将式(22)代入式(20)，提取同一个尺度，得到

$$
v_i^\mu v_j^\nu\mathcal H_{\mu\nu}
=-\frac{h_1}{2}(A_{\rm ref}+i0)^{2\varepsilon}
\mathcal S(f,d,u,w;\varepsilon),
\tag{24}
$$

$$
\boxed{\mathcal S=
\frac{(1+2\varepsilon)df+(1-2\varepsilon)uw}{2(1+\varepsilon)}
(f+i0)^{-2+2\varepsilon}.}
\tag{25}
$$

式(24)的`−1/2`来自r/A_ref，不来自颜色。尺度拆分首先在A_ref>0、f>0验证，再共同解析延拓；式中每个i0是同一母式的边界值简写，不是可以独立选择的多个处方。

原度规/导数两项保持为

$$
\begin{aligned}
\mathcal S_g&=d(f+i0)^{-1+2\varepsilon},\\
\mathcal S_{\partial}&=\frac{-df+(1-2\varepsilon)uw}{2(1+\varepsilon)}(f+i0)^{-2+2\varepsilon},\\
\mathcal S&=\mathcal S_g+\mathcal S_{\partial}.
\end{aligned}
\tag{26}
$$

注意：求导之后式(20)的度规张量系数已吸收导数的一部分，因此它与式(26)中“原度规项”的系数不同。代码的`tensor_decomposition.metric`对应式(26)第一行。

## 7. C 的逐因子组装及 Gamma 简化

将式(24)放入式(18)，两个负号相消。两个轮廓端点的统一前因子是

$$
C_2=\frac{g^2C_F\mathcal A_1h_1}{2}(A_{\rm ref}+i0)^{2\varepsilon},
\qquad J_{34}=\int_0^1ds\int_0^1dt\,\mathcal S.
\tag{27}
$$

这里的1/2是尺度提取给出的，不是图的对称因子。展开式(27)而暂不做Gamma简化：

$$
C_2=\frac{g^4C_F\mathcal B_1}{4(4\pi)^{2-\varepsilon}}
\frac{\Gamma(\varepsilon)\Gamma^2(1-\varepsilon)}{\Gamma(2-2\varepsilon)}
\frac{\Gamma(1-2\varepsilon)}{4^{1+\varepsilon}\pi^{2-\varepsilon}\Gamma(1+\varepsilon)}
(A_{\rm ref}+i0)^{2\varepsilon}.
\tag{28}
$$

使用

$$
\frac{\Gamma(\varepsilon)}{\Gamma(1+\varepsilon)}=\frac1\varepsilon,
\qquad\frac{\Gamma(1-2\varepsilon)}{\Gamma(2-2\varepsilon)}=\frac1{1-2\varepsilon},
\qquad g^2=4\pi\alpha_s,
\tag{29}
$$

得到最终的raw分解

$$
\boxed{\begin{aligned}
I_{34}&=C_2J_{34},\\
C_2&=\frac{C_F\mathcal B_1(\varepsilon)}{16\varepsilon(1-2\varepsilon)}
\left(\frac{\alpha_s}{\pi}\right)^2
[\pi^\varepsilon\Gamma(1-\varepsilon)]^2
(-2v_{14}+i0)^{2\varepsilon},\\
J_{34}&=\int_0^1ds\int_0^1dt\,
\frac{z}{1+\varepsilon}[z(1-s)t+i0]^{-1+2\varepsilon}.
\end{aligned}}
\tag{30}
$$

C包含完整CF和A₁中的CA、TF nf通道。没有颜色a₁投影，没有Euler补偿，没有C4轨道重数，没有反项。两个不同边上的参数各出现一次，也没有额外1/2!。

## 8. 相邻边化简及精确控制积分

由式(3)可直接算出

$$
X^2=2v_{34}(1-s)t,\quad
v_3\cdot X=-t v_{34},\quad v_4\cdot X=-(1-s)v_{34}.
\tag{31}
$$

从而f=z(1−s)t、d=z、u=−zt、w=−z(1−s)，满足uw=zf。式(25)的分子发生精确相加：

$$
(1+2\varepsilon)df+(1-2\varepsilon)uw=2zf.
\tag{32}
$$

因此

$$
\mathcal S_{34}=\frac{\mathcal S_{34,g}}{1+\varepsilon},\qquad
\mathcal S_{34,\partial}=-\frac{\varepsilon}{1+\varepsilon}\mathcal S_{34,g}.
\tag{33}
$$

这解释了式(30)的J，且不需要先积分才能证明比值。若把式(16)两个Fourier核之间的减号误写成加号，会得到旧的错误比值(1+2ε)/(1+ε)，并破坏横向性。

在Re ε>0的正实起始域，按式(1)的原次序逐次积分：

$$
\int_0^1dt\,t^{-1+2\varepsilon}=\frac1{2\varepsilon},\qquad
\int_0^1ds\,(1-s)^{-1+2\varepsilon}=\frac1{2\varepsilon}.
\tag{34}
$$

所以

$$
J_{34,g}=\frac{z^{2\varepsilon}}{4\varepsilon^2},\quad
J_{34,\partial}=-\frac{z^{2\varepsilon}}{4\varepsilon(1+\varepsilon)},\quad
J_{34}=\frac{z^{2\varepsilon}}{4\varepsilon^2(1+\varepsilon)}.
\tag{35}
$$

由此得到式(1)在当前Π定义下的完整控制结果：

$$
\boxed{I_{34}=
\frac1{64}\left(\frac{\alpha_s}{\pi}\right)^2 C_F\mathcal B_1(\varepsilon)
\frac{\Gamma^2(1-\varepsilon)}{\varepsilon^3(1-2\varepsilon)(1+\varepsilon)}
(-2\pi v_{34}+i0)^{2\varepsilon}.}
\tag{36}
$$

先由收敛积分得到式(35)，再按亚纯函数理解ε=0附近的表达。该结果是原动量起点的C×J，未进行Wilson项目的进一步I归一化或Laurent展开。

## 9. 同一核用于对边及三圈单插入

### 9.1 对边：保留u、w

对2/4边，以s=s21、t=s41代入相同的端点位置，得到

$$
f=(1-s)(1-t)+zst,\quad d=-(1+z),\quad
u=-1+(1+z)t,\quad w=1-(1+z)s.
\tag{37}
$$

仍使用式(25)、(27)的同一个C₂和完整S。此时通常uw≠df，因此不能把full/g=1/(1+ε)当作任意SE线的替换规则。尤其d=0也不足以删除整个自能核。

### 9.2 额外裸线与公共三圈C

四端点单Π₁图只增加一条裸交换，但四个端点共享原有序域Ω。裸线的Fourier系数为h₀=Γ(1−ε)/(4π^(2−ε))；连同其两个Wilson顶点，剥离颜色后为

$$
g^2h_0(v_a\cdot v_b)(-X_b^2+i0)^{-1+\varepsilon}
=-\frac12\left(\frac{\alpha_s}{\pi}\right)\pi^\varepsilon\Gamma(1-\varepsilon)
(-2v_{14}+i0)^\varepsilon\,d_b(f_b+i0)^{-1+\varepsilon}.
\tag{38}
$$

颜色必须用四个Wilson生成元的实际顺序重新收缩，不能直接把两端点的CF再乘一个CF。令T_W为完整归一化迹：交错配对为CF(CF−CA/2)，非交错配对为CF²。当前十一张4e-Pi1来源都是前者。

合并式(28)和裸线(38)，总圈数L=2或3的统一表达为

$$
\boxed{C_L(T_W)=\frac{(-1)^{L-2}T_W\mathcal B_1(\varepsilon)}
{2^{L+2}\varepsilon(1-2\varepsilon)}
\left(\frac{\alpha_s}{\pi}\right)^L
[\pi^\varepsilon\Gamma(1-\varepsilon)]^L
(-2v_{14}+i0)^{L\varepsilon},\quad L\in\{2,3\}.}
\tag{39}
$$

这里L=2时T_W=CF；L=3时

$$
J=\int_{\Omega}^{\rm ordered}d^4s\,
\mathcal S(f_{\rm SE},d_{\rm SE},u_{\rm SE},w_{\rm SE};\varepsilon)
\,d_b(f_b+i0)^{-1+\varepsilon}.
\tag{40}
$$

例如4e-Pi1-01的SE端点为(s31,s41)、裸线端点为(s32,s42)，积分次序是

$$
\int_0^1ds_{32}\int_0^{s_{32}}ds_{31}
\int_0^1ds_{42}\int_0^{s_{42}}ds_{41}.
\tag{41}
$$

4e-Pi1-02在同一域上交换SE/裸线的位置。不能丢掉这种标记而认为它们是同一个输入；也不能把式(41)拆成两个无序二重积分。式(39)在本文只为L=2、3建立，不因此声称更高阶、多重插入或Π₂已经实现。

### 9.3 六边对导数和的用途

局部几何给出f_si=u、f_sj=−w、f_sisj=−d，因此

$$
\mathcal S_{\partial}=
\frac1{4\varepsilon(1+\varepsilon)}\partial_{s_i}\partial_{s_j}(f+i0)^{2\varepsilon}.
\tag{42}
$$

两个端点都在独立[0,1]域时，它积分成四个角点之差。对闭合类光四边形，采用共同零距离端点值及延拓约定，六个不同边对的这些角点项相消；在Re ε>0的起始解释下零距离幂为零。这个二点扇区恒等式用于检查完整和，不能据此把式(1)的单图导数项删除。四端点有序域会产生与裸线端点相接的碰撞边界，式(42)本身没有证明它们相消。

## 10. 论文比较：保留输入I与比较对象的区别

记用户此前给出的度规/约化论文式为P，v=v34，TF=1/2。D=4−2ε使(D−4)³=−8ε³，因此它对应

$$
P=-\frac1{64}\left(\frac{\alpha_s}{\pi}\right)^2 C_F\mathcal B_1
\frac{\Gamma^2(1-\varepsilon)}{\varepsilon^3(1-2\varepsilon)}
(-2\pi v_{34}+i0)^{2\varepsilon}.
\tag{43}
$$

由式(36)直接得到

$$
\boxed{I_{34}^{\text{式(1)}}=-\frac{P}{1+\varepsilon}.}
\tag{44}
$$

本文没有通过改动式(1)去消掉这个负号。若另行定义与P同整体符号的完整贡献，则定义为−I₃₄，并满足P/(1+ε)。这是一项需要明确记录的对象转换，不能同时把式(1)的I、论文的W和Wilson项目进一步归一化后的I当作同一对象。

因此式(30)、(39)是**给定动量起点**的代码依据；式(44)是相对给定论文式的比较结论。新推导不把历史的`paper_comparison=failed_overall_sign`改写为已通过。原文图间全导数分配和物理W完整约定的核对仍与本文代数推导分开，见[上一轮审计](paper-check-r001.md)。

## 11. 从公式到实现

生产代码只需实现本文已经完成的局部运算。下表中的公式表达式同时给出符号和归一化，不能仅保留其中的距离幂。

| 数学对象 | 当前实现位置 | 必须保存的连接 |
|---|---|---|
| 式(1)、(3)、(4) | `self_energy.validate_source`、`parse_self_energy_orbit`及`operator`/`phase_ledger` | 一次Π₁、两段传播子、内部δ桥、共享动量、完整相位、外负号及Π所含i |
| 式(6)及四生成元迹 | `color_data` | 实际Wilson顺序；完整CA、CF、TF nf；不得预先投影 |
| 式(22)、(23) | `geometry_for_line` | SE线/裸线标签、原s端点、f/d/u/w；不因d=0丢整个SE |
| 式(25)、(26) | `build_result`中的三个分子 | `metric_num=2(1+eps)df`；`derivative_num=-df+(1-2eps)uw`；相加为完整分子 |
| 式(28)至(30)、(39) | `C.components`及`C.combined` | 各因子只乘一次；总圈数2或3；所有分量乘积等于combined |
| 式(40)、(41) | `J_raw`、`J_raw_s` | 原s上下限和外到内顺序；不是自动执行超立方或IBP |
| 式(33)、(36)、(44) | 数学控制检查 | 相邻full/g；完整C×J；与论文之间的整体负号都分别检查 |

原JSON的`bounds`数组不能默认解释为积分符号次序。本批采用按围道边递增、同边槽号递减的原积分解释，并已逐图对照来源Markdown；相关说明必须随输出保存，不能把排序算法当作Fubini定理。

C4是带SE/裸线标签、保留端点槽位的拓扑作用。其重数保存在元数据，未乘入式(39)。导出一个代表，不把旋转副本和相同代表的C一起重复计数。

若以后更换Π定义或动量式外部符号，需要同时修订来源契约、operator/phase_ledger、共同C、控制例和发布版本。仅改01的答案或只改文档的符号会破坏式(1)到代码的连接。

## 12. 推导检查及使用边界

本次补充检查覆盖动量分子的符号转换、Fourier系数比、D维散度、显式Hessian、相邻化简、未简化Gamma前因子到C₂的连接、裸线到C₃的连接，以及当前13个生产输出对式(25)、(26)、(39)的逐项吻合。错误Fourier相对号和错误整体号保留为反例。具体记录见[检查目录](../../../workflow/checks/self-energy-pi1-theory/README.md)。

本文给出原始积分到参数化的解析推导，并以精确代数辅助核对；它不是Lean形式证明。没有改变原生产代码、r001输入、默认47条prepared或旧积分结果，也没有由此授予全部物理复运动学的延拓已验证状态。Π₂、外接TGV和多重自能不属于本文已实现的范围。
