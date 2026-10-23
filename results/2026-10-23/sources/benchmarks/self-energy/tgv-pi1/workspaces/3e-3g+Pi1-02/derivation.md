# 3e-3g+Pi1 02的完整解析解

## 当前物理积分的入口（2026-10-15，r003）

从共同作用量出发，U′=igAU确定真实Wilson词；FT/Dyson规则在提取C前生效。当前[原生C/J](../../sources/montopar-r002.json)给出下式，J就是本文实际积分的原核，域、分支和IBP全保留：

$$R_P=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad I_P=2R_P\frac{J(z)+z^{3\epsilon}J(1/z)}{1+z^{3\epsilon}}.$$

用下文已经求出的J代入此定义，得到[当前物理I](../../results/3e-3g+Pi1-02/r003.json)；组装系数为+1，Euler只补一次。以下带日期的求值过程与r001公式属于旧根历史，计算J的步骤可按同域积分线性复用，旧形式状态不继承。

2026-10-15表示修订：主I、原J和主Laurent系数保持r002的值；当前源哈希、配对引用及附属公式统一到新修订。[修订说明](../../revisions/r003.md)。以下2026-10-01/02的r001系数表保留历史符号，不能直接当作当前结果。

2026-10-01。本图属于344相邻边核心，自能插在边4的前槽，第r=2腿。与01相比，原排序域与q₀相同，但权重、E₂及F₂不同。本页从[实际raw包](../../sources/montopar-r001.json)的02号完整核独立积分，得到全ε Gamma闭式及ε⁻⁵至ε⁰；没有假定I₀₂=I₀₁。[结果r001](../../results/3e-3g+Pi1-02/r001.json)。

## 1. 原对象与两个路径尺度

取s=s₃₁、t=s₄₁、u=s₄₂，原域为∫₀¹ds∫₀¹du∫₀ᵘdt。三α在非负单纯形，α₁+α₂+α₃=1；完整标记权重是α₁^{-ε}α₂^{-2ε}α₃^{-ε}。D=4−2ε保留在分子，m=2−3ε。

原几何q₀=zα₁(1−s)(α₂t+α₃u)，两项直接核分别为K₁/(q₀+i0)^m及K₂/(q₀+i0)^{m+1}。K₁、K₂是实际JSON的完整张量分子。下面的显式变换也给出了它们的紧凑表达。

令x=1−s、t=u p，x,u,p∈[0,1]，路径Jacobian=u；再令α₁=a、α₂=(1−a)v、α₃=(1−a)(1−v)，单纯形Jacobian=1−a。记β=1−a、y=1−v、f=aβ(y+vp)。两个分子精确变为

$$K_1=z^2x\,k_1,\qquad K_2=z^3x^2u\,k_2,$$

$$k_1=-\frac{a\beta}{4}P(a,v),$$

$$P=2a(\epsilon-2)(v^2-2v)+(4-2\epsilon)v^2+(6\epsilon-1)v-(2\epsilon+1),$$

$$k_2=-\frac{ma^2v\beta^2}{4}\left[\beta y^2(1-p)-a(p+y)\right].$$

q₀=zxu f，两个核乘Jacobian后共同带x^{-1+3ε}u^{-1+3ε}。因此

$$\boxed{J_{02}(\epsilon,z)=\frac{z^{3\epsilon}}{9\epsilon^2}T_2(\epsilon),}$$

$$T_2=\int_0^1da\int_0^1dv\,a^{-\epsilon}\beta^{1-3\epsilon}v^{-2\epsilon}y^{-\epsilon}
\int_0^1dp\,\frac{N}{f^{m+1}},\qquad N=k_1f+k_2.$$

N对p至多一次。这个T₂是中间积分，下一步将全部积完，不作为未知特殊函数保留。

## 2. 收敛与线性比例参数的精确积分

先取z>0、A₀=−2v₁₄>0及0<Re ε<1/3。提出尺度和α总尺度后，a端点共同权重是a^{2Re ε−1}，β端点有界。v→0时分母y+vp有正下界，权重v^{-2Re ε}可积；v→1时k₂括号在p=y=0消失，可界为常数乘p+y。先积p后，两项均至多给y^{-1+2Re ε}，也可积。x、u的幂则是−1+3Re ε。

这给出本图原两项核的共同绝对可积控制，允许上述分组积分及共同iη边界值η↓0。之后按ε亚纯延拓解释，不把ε直接置零，也不推广为其他图的积分换序证明。

令f=f₀+f₁p，其中f₀=aβy、f₁=aβv；令N=n₀+n₁p。通用精确原函数给

$$\int_0^1\frac{N}{f^{m+1}}dp
=c_1[(a\beta)^{1-m}-(a\beta y)^{1-m}]
+c_2[(a\beta)^{-m}-(a\beta y)^{-m}],$$

$$c_1=\frac{n_1}{f_1^2(1-m)},\qquad
c_2=-\frac{n_0-n_1f_0/f_1}{mf_1}.$$

从02的实际分子化简得到

$$c_1=\frac{W_0(v)+aW_1(v)}{4v(3\epsilon-1)},\qquad
c_2=\frac{a\beta y^2}{4v},$$

$$W_0=-(\epsilon+2)v^2-3v+3-\epsilon,\qquad
W_1=(\epsilon+2)(v^2-2v).$$

这些等式与JSON长分子逐点代数相符。差分在v=0抵消其表面的1/v奇性；不能将端点幂单独当普通收敛积分。

## 3. 收敛一维Euler式及Beta闭合

记δ=3ε−1。β幂在两项中都消去，a积分只需1/(2ε)和1/(1+2ε)。定义

$$C_2(v)=\frac{W_0(v)}{2\epsilon}+\frac{W_1(v)}{1+2\epsilon},$$

则剩余积分为

$$\boxed{T_2=\frac14\int_0^1dv\,v^{-1-2\epsilon}(1-v)^{-\epsilon}
\left[C_2(v)\frac{1-(1-v)^{3\epsilon-1}}{3\epsilon-1}
+\frac{(1-v)^2-(1-v)^{3\epsilon}}{2\epsilon}\right].}$$

方括号在v=0为O(v)，在v=1的最强幂为(1−v)^{-1+2ε}，与上节的可积界一致。

写C₂(v)=Σ_{k=0}²d_kv^k。为拆开差分，先共同乘v^σ，Re σ>2Re ε时应用Euler Beta公式，然后对完整组合延拓至σ=0。得到

$$T_2=\frac14\left\{\frac1{3\epsilon-1}\sum_{k=0}^2d_k
[\mathrm B(k-2\epsilon,1-\epsilon)-\mathrm B(k-2\epsilon,2\epsilon)]
+\frac{\mathrm B(-2\epsilon,3-\epsilon)-\mathrm B(-2\epsilon,1+2\epsilon)}{2\epsilon}\right\}.$$

负Beta参数均按上述共同延拓解释；B(−2ε,2ε)=0来自Gamma分母的亚纯极限，不是删除原积分。另一条直接从N的展开式积分的路线保留共同α端点辅助指数，组合后取极限，独立给出相同Gamma式。

定义

$$A_2(\epsilon)=36\epsilon^4+18\epsilon^3-28\epsilon^2-6\epsilon+4,$$
$$B_2(\epsilon)=18\epsilon^4+13\epsilon^3-9\epsilon^2-2\epsilon-2.$$

全部Beta项用Gamma递推化简后得到

$$\boxed{T_2(\epsilon)=-\frac{\Gamma(1-2\epsilon)}
{8\epsilon^2(1+2\epsilon)(2-3\epsilon)(1-3\epsilon)^2}
\left[A_2(\epsilon)\Gamma(1+2\epsilon)
+B_2(\epsilon)\frac{\Gamma(1-\epsilon)}{\Gamma(1-3\epsilon)}\right].}$$

这里已没有待计算的参数积分或矩。

## 4. 唯一I、尺度与两个颜色通道

沿用本分支定义

$$R(\epsilon)=\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}
{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},$$

$$I_{02}=2R\frac{J_{02}(\epsilon,z)+z^{3\epsilon}J_{02}(\epsilon,z^{-1})}{1+z^{3\epsilon}}
=\frac{2R}{9\epsilon^2}T_2.$$

因此唯一I与z无关，完整答案为

$$\boxed{I_{02}(\epsilon)=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2}
{2304\epsilon^5(1+\epsilon)(1+2\epsilon)(1-2\epsilon)(1-3\epsilon)(2-3\epsilon)}
\left[A_2(\epsilon)\Gamma(1-3\epsilon)\Gamma(1+2\epsilon)+B_2(\epsilon)\Gamma(1-\epsilon)\right].}$$

这包含一次Euler补偿及规定的四个C4代表组装，不再额外乘4。原raw单图的z^{3ε}尺度没有被删除。

全颜色prepared轨道和为

$$\mathcal W_{02}^{prepared}=\lambda_b^3(S_{14}+S_{12})
\left[C_FC_A^2\frac{10-6\epsilon}{3-2\epsilon}
-C_FC_AT_Fn_f\frac{8(1-\epsilon)}{3-2\epsilon}\right]I_{02}.$$

两个颜色通道共享I。机器结果另外保存b_gI与b_fI的全部极点至有限项，完整卷积了有理系数的正ε项。共同外尺度S₁₄+S₁₂保持在外，尚未作完整w₃物理组装或RG减除。恢复未补Euler的裸轨道和时乘回e^{+3γ_Eε}。

## 5. 全部极点与有限项

A₂(0)+B₂(0)=2，所以最高极点恰为五阶，留数−1/2304。记ζ_n=ζ(n)，

$$I_{02}=\sum_{k=-5}^0c_k\epsilon^k+O(\epsilon),$$

| k | c_k |
|---:|---|
| −5 | −1/2304 |
| −4 | 1/4608 |
| −3 | 65/9216−3ζ₂/512 |
| −2 | 55/2048−31ζ₂/3072−13ζ₃/2304 |
| −1 | 3925/36864−21ζ₂/2048−47ζ₃/4608−199ζ₂²/6144 |
| 0 | 25199/73728−17ζ₂/12288−127ζ₃/9216−259ζ₂²/4096−23ζ₂ζ₃/512−47ζ₅/1280 |

与01的留数和下一阶已不同，不能因为同一附着几何就认为I相等。

展开时使用Gamma对数的有限阶递推。若G=exp(Σℓ_nε^n)、G=Σg_nε^n，则g₀=1、g_n=(1/n)Σ_{k=1}^n kℓ_kg_{n-k}；再与有理分母的逆级数精确卷积。两组Gamma乘积的Euler一次项均与补偿相消，所需ℓ_n分别为[2+3^n+(−2)^n]ζ_n/n和3ζ_n/n。计算至n=5足以得到I有限项，且Γ组合在ε=0的解析性给出所列O(ε)余项。

## 6. 实际核验

两条独立Beta整理路线精确一致，02的分子和标记权重直接读取不可变raw包。收敛一维Euler式在ε=0.12、0.20、0.28与Gamma式比较通过；六个几何系数及十二个颜色系数另与未展开Gamma式的80位导数核对，三个比例积分点也通过。

数值比较不是严格误差包围，解析依据是保持共同调节的积分恒等式和Gamma/Beta延拓。原五维积分未做独立数值回归，一般复运动学延拓未核验。当前没有可调用的Wilson数学MCP，未生成形式收据，formal=not_run；Montopar、Kernel、Lean源以及01结果均未修改。[检查记录](../../../../../workflow/checks/tgv-pi1-02/README.md)。
