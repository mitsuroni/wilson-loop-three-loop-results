# 3e-3g+Pi1 03的完整解析解

## 当前物理积分的入口（2026-10-15，r003）

从共同作用量出发，U′=igAU确定真实Wilson词；FT/Dyson规则在提取C前生效。当前[原生C/J](../../sources/montopar-r002.json)给出下式，J就是本文实际积分的原核，域、分支和IBP全保留：

$$R_P=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},\qquad I_P=2R_P\frac{J(z)+z^{3\epsilon}J(1/z)}{1+z^{3\epsilon}}.$$

用下文已经求出的J代入此定义，得到[当前物理I](../../results/3e-3g+Pi1-03/r003.json)；组装系数为+1，Euler只补一次。以下带日期的求值过程与r001公式属于旧根历史，计算J的步骤可按同域积分线性复用，旧形式状态不继承。

2026-10-15表示修订：主I、原J和主Laurent系数保持r002的值；当前源哈希、配对引用及附属公式统一到新修订。本图三角函数表示的旧号已修正，与当前Gamma式一致。[修订说明](../../revisions/r003.md)。以下2026-10-01/02的r001系数表保留历史符号，不能直接当作当前结果。

2026-10-01。03仍在344相邻边核心，但Π₁插在独立的边3、第r=1腿。原有序域与01、02相同，标记权重改为α₁^{-2ε}α₂^{-ε}α₃^{-ε}。本页从[实际raw包](../../sources/montopar-r001.json)独立求解，得到全ε闭式、全部极点至有限项及两颜色通道的卷积。[结果r001](../../results/3e-3g+Pi1-03/r001.json)。

## 1. 完整两项核在本图的简化

记s=s₃₁、t=s₄₁、u=s₄₂，原路径域为∫₀¹ds∫₀¹du∫₀ᵘdt；三α在非负单纯形上。q₀=zα₁(1−s)(α₂t+α₃u)。D=4−2ε；本图分子中显式D迹的系数是同边类光内积G₂₃=0，所以相应项自动消失，没有提前把D设为4。

令x=1−s、t=u p，路径Jacobian=u，x,u,p∈[0,1]；再取α₁=a、α₂=(1−a)v、α₃=(1−a)(1−v)，单纯形Jacobian=1−a。记β=1−a、y=1−v、m=2−3ε、δ=1−m=3ε−1，定义f=aβ(y+vp)。

实际完整分母为q₀=zxu f。两个张量分子K₁、K₂精确变为

$$K_1=z^2x\,k_1,\qquad K_2=z^3x^2u\,k_2,$$

$$k_1=\frac{a\beta(a+\epsilon)(1-2v)}2,\qquad
k_2=\frac{ma^2\beta^3(y+vp)(1-2v)}4.$$

因此第二项分子含f，可与第一项按同一边界值精确合并：

$$\boxed{\frac{k_1}{f^m}+\frac{k_2}{f^{m+1}}
=\frac{a\beta(1-2v)[2-\epsilon+3\epsilon a]}{4f^m}.}$$

这里合并的是本图完整两项核，而不是只保留自能度规部分。原两个raw项仍保存，在共同可积控制下先取正实核的iη边界值，再使用此式；没有从有限η的两项中丢弃iη分子。

两个路径尺度共同为x^{-1+3ε}u^{-1+3ε}，积掉后

$$\boxed{J_{03}(\epsilon,z)=\frac{z^{3\epsilon}}{9\epsilon^2}T_3(\epsilon),}$$

$$T_3=\frac14\int_0^1da\,a^{\epsilon-1}(1-a)^\epsilon(2-\epsilon+3\epsilon a)
\int_0^1dv\,v^{-\epsilon}(1-v)^{-\epsilon}(1-2v)
\int_0^1dp\,(1-v+vp)^{-m}.$$

这一步显式给出全部变量、域和Jacobian，不借用01/02的结果或镜像关系。

## 2. 本图的共同收敛区及两个Beta组合

先在z>0、A₀=−2v₁₄>0、0<Re ε<1/3计算。a端点权重a^{Re ε−1}与(1−a)^{Re ε}可积，两个路径尺度也可积。v→0时1−v+vp有正下界；v→1时先积p，幂界至多是(1−v)^{1−Re m−Re ε}=(1−v)^{-1+2Re ε}。两个原raw项分别受相同控制，故本图可以分组积分并取共同η↓0；然后按ε亚纯延拓解释。

p积分是完整端点差分

$$\int_0^1dp\,(1-v+vp)^{-m}=\frac{1-(1-v)^\delta}{v\delta}.$$

它在v=0的表面1/v被分子抵消。于是a与v分离：

$$T_3=\frac{A_a(\epsilon)}{4(3\epsilon-1)}V_v(\epsilon),$$

$$A_a=(2-\epsilon)\mathrm B(\epsilon,1+\epsilon)
+3\epsilon\mathrm B(1+\epsilon,1+\epsilon),$$

$$V_v=\mathrm B(-\epsilon,1-\epsilon)-2\mathrm B(1-\epsilon,1-\epsilon)
-\mathrm B(-\epsilon,2\epsilon)+2\mathrm B(1-\epsilon,2\epsilon).$$

v式由同一个收敛差分推出：先共同乘v^σ，取Re σ>Re ε使拆开的Euler项收敛，再对完整组合延拓至σ=0。B(−ε,2ε)不为零，不能按02的另一个Beta零参数情形删除它。

Gamma递推分别给出

$$A_a=\frac{(1+\epsilon)(2+\epsilon)\Gamma(1+\epsilon)^2}
{\epsilon(1+2\epsilon)\Gamma(1+2\epsilon)},$$

$$V_v=\frac{\Gamma(1-\epsilon)}\epsilon
\left[\frac32\frac{\Gamma(1+2\epsilon)}{\Gamma(1+\epsilon)}
-\frac{\Gamma(1-\epsilon)}{(1-2\epsilon)\Gamma(1-2\epsilon)}\right].$$

因此所有剩余参数已经积完：

$$\boxed{T_3=\frac{(1+\epsilon)(2+\epsilon)\Gamma(1-\epsilon)\Gamma(1+\epsilon)}
{8\epsilon^2(1+2\epsilon)(3\epsilon-1)}
\left[3-\frac{2\Gamma(1-\epsilon)\Gamma(1+\epsilon)}
{(1-2\epsilon)\Gamma(1-2\epsilon)\Gamma(1+2\epsilon)}\right].}$$

此式与从JSON长分子、共同辅助α指数直接计算的Gamma组合精确相同。

## 3. 唯一I及Gamma反射简式

本分支约定为

$$R=\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^2\Gamma(2-3\epsilon)}
{64\epsilon(1+\epsilon)\Gamma(2-2\epsilon)},$$

$$I_{03}=2R\frac{J_{03}(\epsilon,z)+z^{3\epsilon}J_{03}(\epsilon,z^{-1})}{1+z^{3\epsilon}}
=\frac{2R}{9\epsilon^2}T_3.$$

故指定唯一I与z无关；原raw图的z^{3ε}因子及恢复后的共同外尺度没有被删除。Gamma反射公式进一步给出

$$\frac{\Gamma(1-\epsilon)\Gamma(1+\epsilon)}
{\Gamma(1-2\epsilon)\Gamma(1+2\epsilon)}
=\frac{\sin(2\pi\epsilon)}{2\sin(\pi\epsilon)}=\cos(\pi\epsilon).$$

本图完整结果可紧凑写为

$$\boxed{I_{03}(\epsilon)=
-\frac{(2+\epsilon)e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3\Gamma(1+\epsilon)\Gamma(1-3\epsilon)}
{2304\epsilon^5(1+2\epsilon)(1-2\epsilon)^2\Gamma(1-2\epsilon)}
\left[3(1-2\epsilon)-2\cos(\pi\epsilon)\right].}$$

Gamma比形式与cos形式都保存于结果r001，符号反射化简已核对。特殊ε处的可去组合按整体亚纯极限处理，不分别给各因子指定无穷值。

这里已经包含一次Euler补偿及规定的C4组装；不再乘4或补第二次Euler。全颜色prepared轨道和恢复为

$$\mathcal W_{03}^{prepared}=\lambda_b^3(S_{14}+S_{12})
\left[C_FC_A^2\frac{10-6\epsilon}{3-2\epsilon}
-C_FC_AT_Fn_f\frac{8(1-\epsilon)}{3-2\epsilon}\right]I_{03}.$$

两颜色通道的完整有理系数卷积亦已存入机器结果；共同外尺度未展开，完整w₃/RG尚未组装。恢复未补Euler的裸轨道和时乘回e^{+3γ_Eε}。

## 4. 全部极点至有限项

完整式除ε⁻⁵外在ε=0解析，方括号趋于1、2+ε趋于2，因此最高极点恰为五阶，留数−1/1152。记ζ_n=ζ(n)，

$$I_{03}=\sum_{k=-5}^0c_k\epsilon^k+O(\epsilon),$$

| k | c_k |
|---:|---|
| −5 | −1/1152 |
| −4 | 7/2304 |
| −3 | 1/192−7ζ₂/768 |
| −2 | 17/576+ζ₂/1536−7ζ₃/1152 |
| −1 | 1/18−3ζ₂/128+49ζ₃/2304−547ζ₂²/15360 |
| 0 | 3/16+11ζ₂/384+7ζ₃/192−11ζ₂²/30720−49ζ₂ζ₃/768−71ζ₅/1920 |

虽然最高留数与01相同，下一阶已经不同，不能据此认为两个I相等。

展开采用Gamma对数的精确有限阶递推，再与有理分母逆级数卷积。两组Gamma乘积的对数系数分别为[3+(−1)^n+3^n−2^n]ζ_n/n与[4+2(−1)^n+3^n−2·2^n−(−2)^n]ζ_n/n，n≥2；一次Euler项都与补偿相消。计算至n=5得到I有限项，解析性给出O(ε)余项。

## 5. 核验及限制

实际完成两条精确Beta路线及Gamma反射化简；三个收敛一维Euler点、六个未展开Gamma导数系数、三个比例参数原函数点、十二个颜色导数系数检查通过。Euler式的经验相对差约10⁻⁶⁰量级，不是严格数值包围。原五维独立数值积分及一般复运动学延拓未做。

当前没有可调用的Wilson数学MCP，没有生成或伪造形式收据，formal=not_run。原raw、01/02结果及Montopar/Kernel/Lean源均保持。[实际检查](../../../../../workflow/checks/tgv-pi1-03/README.md)。本轮只完成03，不直接继承到334核心的04。
