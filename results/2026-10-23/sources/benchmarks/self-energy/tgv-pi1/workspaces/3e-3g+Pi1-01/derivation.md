# 3e-3g+Pi1 01的完整解析解

## 当前物理积分的入口（2026-10-15，r003）

从共同作用量出发，U′=igAU确定真实Wilson词；FT/Dyson规则在提取C前生效。当前[原生C/J](../../sources/montopar-r002.json)给出下式，J就是本文实际积分的原核，域、分支和IBP全保留：

$$R_P=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad I_P=2R_P\frac{J(z)+z^{3\epsilon}J(1/z)}{1+z^{3\epsilon}}.$$

用下文已经求出的J代入此定义，得到[当前物理I](../../results/3e-3g+Pi1-01/r003.json)；组装系数为+1，Euler只补一次。以下带日期的求值过程与r001公式属于旧根历史，计算J的步骤可按同域积分线性复用，旧形式状态不继承。

2026-10-15表示修订：主I、原J和主Laurent系数保持r002的值；当前源哈希、配对引用及附属公式统一到新修订。[修订说明](../../revisions/r003.md)。以下2026-10-01/02的r001系数表保留历史符号，不能直接当作当前结果。

2026-10-01。目标是从[实际Montopar raw输入](../../sources/montopar-r001.json)得到本图的唯一I，包括全部极点至有限项。本图属于相邻边344，Π₁插在边4的后槽，第r=3腿。下面得到全ε Gamma闭式；不修改原raw或生产代码。[结果r001](../../results/3e-3g+Pi1-01/r001.json)。

## 1. 原积分与可分离的两个端点尺度

记s=s₃₁、t=s₄₁、u=s₄₂，原有序积分为∫₀¹ds∫₀¹du∫₀ᵘdt，三个α在非负单纯形上。取D=4−2ε、z>0、A₀=−2v₁₄>0，先在0<Re ε<1/3求积分，再作ε亚纯延拓。

令x=1−s、b=α₂t+α₃u，则实际两项核为

$$J_{01}=\int ds\,du\,dt\int d\mu_\alpha\,
\alpha_1^{-\epsilon}\alpha_2^{-\epsilon}\alpha_3^{-2\epsilon}
\left[\frac{(1+2\epsilon+\alpha_3)L-\alpha_3E}{(z\alpha_1xb+i0)^{2-3\epsilon}}
-\frac{2(2-3\epsilon)\alpha_3F}{(z\alpha_1xb+i0)^{3-3\epsilon}}\right],$$

$$L=\frac{z^2\alpha_1x}{4}(\alpha_3-\alpha_2),$$
$$E=\frac{z^2\alpha_1x}{4}\left[(4-2\epsilon)(\alpha_2-\alpha_1)-(1-\alpha_3+\alpha_2)\right],$$
$$F=\frac{z^3\alpha_1^2x^2}{8}\left[(\alpha_1-\alpha_2)b+\alpha_2u\right].$$

这里dμ_α=d³α δ(1−Σα)，t的上界是u。将t=u p、0≤p≤1，Jacobian=u；x,u,p均在[0,1]。记m=2−3ε、b'=α₂p+α₃，则q₀=zα₁xu b'。两个核连同Jacobian都带x^{-1+3ε}u^{-1+3ε}，因此

$$\boxed{J_{01}(\epsilon,z)=\frac{z^{3\epsilon}}{9\epsilon^2}T(\epsilon),}$$

$$T=\frac14\int d\mu_\alpha\,
\alpha_1^{2\epsilon-1}\alpha_2^{-\epsilon}\alpha_3^{-2\epsilon}
\int_0^1dp\left[\frac{V}{b'^m}-\frac{m\alpha_3\alpha_2}{b'^{m+1}}\right],$$

$$V=(1+2\epsilon+\alpha_3)(\alpha_3-\alpha_2)
-\alpha_3\left[(4-2\epsilon)(\alpha_2-\alpha_1)-(1-\alpha_3+\alpha_2)\right]
-m\alpha_3(\alpha_1-\alpha_2).$$

T是仍待求的积分，不作为新特殊函数命名来宣称完成。1/(9ε²)只来自两个路径尺度；完整极点阶数还取决于α端点与共同Gamma前因子。

## 2. 单纯形变换、收敛及精确比例积分

取α₁=a、α₂=(1−a)v、α₃=(1−a)(1−v)，a,v∈[0,1]，单纯形Jacobian=1−a。原二维α域与这个方形通过上述双射连接。令y=1−v、δ=1−m=−1+3ε。

V含因子1−a，其商是低次多项式：

$$\frac V{1-a}=U_0(v)+aU_1(v),$$
$$U_0=(\epsilon+2)v^2-5(\epsilon+1)v+2(\epsilon+1),\qquad
U_1=(\epsilon+2)(1-v^2).$$

此步骤在进入Beta积分之前显式消除了虚假的α总尺度端点极点；把长分子完全拆开再逐项积分会暂时产生Beta第二参数为零的项，须在同一辅助调节下先组合。

在当前收敛条带，两项原核均可控制。a端点的共同幂为a^{2Re ε−1}，1−a端点在上述提取后至多为有界多项式；v→0的权重至多v^{-Re ε}。v→1时先对p积分，(1−v+vp)^{-m}以及带1−v的(m+1)幂项均至多为常数乘(1−v)^{1−Re m}，乘原权重后得到(1−v)^{-1+Re ε}。这些幂与两个路径尺度的x,u幂全部可积。

所以在0<Re ε<1/3可对本图作上述分组积分与换序；共同有限iη下各幂也有相同可积控制，沿正实核取η↓0。这里给的是本图的具体控制，不推广为其他图或任意发散迭代积分的Fubini桥。

p只出现在线性分母中，精确积分为

$$T=\frac14\int_0^1da\,a^{2\epsilon-1}\int_0^1dv\,v^{-\epsilon}y^{-2\epsilon}
\left[(U_0+aU_1)\frac{1-y^\delta}{v\delta}+y-y^\delta\right].$$

1−y^δ在v=0为O(v)，这一整体差分不可省略。a积分给1/(2ε)和1/(1+2ε)，故

$$\boxed{T=\frac14\int_0^1dv\,v^{-\epsilon}(1-v)^{-2\epsilon}
\left[C(v)\frac{1-(1-v)^{-1+3\epsilon}}{v(-1+3\epsilon)}
+\frac{(1-v)-(1-v)^{-1+3\epsilon}}{2\epsilon}\right],}$$

$$C(v)=\frac{U_0(v)}{2\epsilon}+\frac{U_1(v)}{1+2\epsilon}.$$

这是保留端点相消的收敛一维积分，后续数值检查直接使用这一式。

## 3. Beta组合闭合成Gamma函数

将C(v)=Σ_{k=0}²c_k v^k。共同引入v^σ，在Re σ>Re ε时拆开的v=0端点项都可用Euler Beta公式；完整差分在σ=0的邻域仍收敛，因此先得到Gamma表达式，再将整个组合解析延拓至σ=0：

$$T=\frac14\left\{\frac1{3\epsilon-1}
\sum_{k=0}^2c_k[\mathrm B(k-\epsilon,1-2\epsilon)-\mathrm B(k-\epsilon,\epsilon)]
+\frac{\mathrm B(1-\epsilon,2-2\epsilon)-\mathrm B(1-\epsilon,\epsilon)}{2\epsilon}\right\}.$$

式中的负Beta参数表示上述共同延拓，不是单独发散Euler积分。特别是B(−ε,ε)的零来自Gamma分母1/Γ(0)=0的亚纯极限，不是随意删去一个原积分。

定义两个有理多项式

$$A(\epsilon)=-45\epsilon^4-90\epsilon^3+89\epsilon^2+6\epsilon-8,$$
$$B(\epsilon)=36\epsilon^4+8\epsilon^3-47\epsilon^2+11\epsilon+4.$$

用Γ(x+1)=xΓ(x)化简全部Beta项，最终得到

$$\boxed{T(\epsilon)=\frac{\Gamma(1-\epsilon)}
{8\epsilon^2(1+2\epsilon)(2-3\epsilon)(1-3\epsilon)^2}
\left[A(\epsilon)\Gamma(1+\epsilon)
+B(\epsilon)\frac{\Gamma(1-2\epsilon)}{\Gamma(1-3\epsilon)}\right].}$$

至此所有参数均已积分，不留下待求的嵌套矩或辅助主积分。此Gamma式与直接从JSON长分子、共同辅助α端点指数求得的Beta组合精确相同。

## 4. 原始J到唯一I及全色贡献的连接

沿用本分支的归一化，

$$R(\epsilon)=\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}
{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},$$

$$I_{01}=2R\frac{J_{01}(\epsilon,z)+z^{3\epsilon}J_{01}(\epsilon,z^{-1})}{1+z^{3\epsilon}}
=\frac{2R}{9\epsilon^2}T.$$

这个等式显式核对了z^{3ε}尺度与两个旋转尺度的组合；I与z无关，原raw图和它的尺度因子并非与z无关。C4次数4及一次Euler补偿由该I定义组织，不再额外乘4或再次补Euler。

代入T得到本图完整答案

$$\boxed{I_{01}(\epsilon)=
\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3}
{2304\epsilon^5(1+\epsilon)(1+2\epsilon)(1-2\epsilon)(1-3\epsilon)(2-3\epsilon)}
\left[A(\epsilon)\frac{\Gamma(1-3\epsilon)\Gamma(1+\epsilon)}{\Gamma(1-2\epsilon)}+B(\epsilon)\right].}$$

当ε→0时A+B→−4，分母除ε⁵外趋于4608；故最高极点恰为五阶，留数−1/1152，不是从扇区计数猜出的上界。

全颜色prepared图轨道和可恢复为

$$\mathcal W_{01}^{prepared}=\lambda_b^3 C_AC_F\mathcal B_1(\epsilon)
(S_{14}+S_{12})I_{01},$$

其中λ_b=α_{s,b}/π，B₁=[(10−6ε)C_A−8(1−ε)T_Fn_f]/(3−2ε)，S₁₂=z^{3ε}S₁₄。若恢复未补Euler的裸图和需乘回e^{+3γ_Eε}。这里保留两个物质/胶子颜色方向，不把I当作已经装配完整w₃。

## 5. 全部极点至有限项

记ζ_n=ζ(n)。本图

$$I_{01}(\epsilon)=\sum_{k=-5}^0c_k\epsilon^k+O(\epsilon),$$

| k | c_k |
|---:|---|
| −5 | −1/1152 |
| −4 | 1/1536 |
| −3 | 73/9216−5ζ₂/768 |
| −2 | 527/18432−41ζ₂/3072−13ζ₃/1152 |
| −1 | 1447/12288−7ζ₂/6144−43ζ₃/1536−29ζ₂²/1024 |
| 0 | 26863/73728−91ζ₂/4096−167ζ₃/9216−4519ζ₂²/61440−41ζ₂ζ₃/768−47ζ₅/640 |

Gamma对数的一次Euler项与补偿相消；通过Γ(1+cε)的对数展开到ε⁵得到这些系数，并用未展开Gamma式的高精度导数逐项独立核对。所需余项由Gamma组合在ε=0的解析性及显式ε⁻⁵确定。

## 6. 验证及范围

直接JSON核与本页的分离形式精确核对，比例参数原函数在三个内部点数值复核。完整一维Euler式与Gamma式在ε=0.12、0.20、0.28一致，经验相对差分别约8×10⁻⁶²、3×10⁻⁵⁶、1×10⁻⁵⁰；六个Laurent系数另与65位Gamma导数比较。数值差不是严格误差包围，解析依据是上述保留共同调节的积分与Beta等式。

[完整检查记录](../../../../../workflow/checks/tgv-pi1-adjacent/README.md)。本会话没有可调用的Wilson数学MCP，未生成或伪造形式收据，formal=not_run；未做一般复运动学延拓或原五维独立数值积分。结果是指定现行raw输入和上述归一化下的完整解析解，未自动推广到其他11图。
