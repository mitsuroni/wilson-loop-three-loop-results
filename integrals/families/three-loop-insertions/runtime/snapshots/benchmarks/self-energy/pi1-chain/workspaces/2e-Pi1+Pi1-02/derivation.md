# 2e-Pi1+Pi1-02：正式求解推导

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1-chain/workspaces/2e-Pi1+Pi1-02/derivation.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。

以下保留迁移前的求解/准备过程记录；其中r001的C/I公式不作为现行数值入口。现行理论及解表已经更新，J与母函数的推导仍可复用。

状态：完整I与ε⁻⁶…ε⁰完成；包含三个颜色通道的有限项卷积。

## 冻结任务

- [prepared包](../../prepared/2e-Pi1+Pi1.json)，SHA-256 `f022b0fc4153ff5bd13bbf7a9ba53833df7ac6388a1e01b2e2ac5bc4a18ef0b8`。
- [归一化契约](../../prepared/normalization.json)，SHA-256 `0bf157b9b405a482603759a8af97817bc6505fbc82654bd47f4b8fae956a58b2`。
- 总圈数3；提出完整B₁²后a₁=1；C4轨道大小m=2。
- 外到内：`s21 → s41`；各变量上下限均0、1，Jacobian=1。
- C含一次三圈Euler补偿，原C完整保存；J及其原顺序未改变。
- 三颜色通道共享本图I，颜色重构与有限项卷积见[说明](../../normalization-and-colors.md)。

记 $e=\varepsilon$、$L=\log z$，沿用[冻结归一化](../../prepared/normalization.json)和[唯一I定义](../../theory.md)SC23：

$$
G=e^{-\gamma_Ee}\Gamma(1-e),\qquad
U=\frac{\Gamma(1+e)^2\Gamma(1-e)\Gamma(1-3e)}{\Gamma(1-2e)^2\Gamma(1+2e)},\qquad P=G^3U.
$$
$$
I=\frac m2R\frac{J(z,e)+z^{3e}J(z^{-1},e)}{1+z^{3e}},\qquad
R=-\frac{P}{128e^2(1-2e)^2}.
$$

提出完整 $C_F\mathcal B_1(e)^2$ 后a₁=1，prepared已含一次三圈Euler补偿。本页不重复补偿，不把原图CJ另称为I。正实z、v14<0下先取0<Re(e)<1/3，再向e=0附近亚纯延拓。

## 1. 原有序J与完整边界

本图m=2，原外到内顺序s21→s41，均为单位区间，Jacobian=1。令s=s21、t=s41、$D=(1-s)(1-t)+zst$。prepared实际完整核为

$$
K=\frac{-(2+e)(1+z)D+(1-3e)z}{2(1+2e)}D^{-2+3e}
=-(1+z)D^{-1+3e}+\frac{\partial_s\partial_tD^{3e}}{6e(1+2e)}.
$$

这是一条完整核恒等式，未删去任何算符项。原次序先积t，再积s；可在截断方形上实施，然后用Re(e)>0的可积性取端点极限。四角距离依次为1、0、0、z，因此导数项积分为 $(1+z^{3e})/[6e(1+2e)]$。

复用[已求值的H₀](../../../../2ne/families/single-opposite/basis.md)：

$$
H_0(\delta,z)=-\frac{(1+z)z^{-\delta/2}}{\Gamma(1-\delta)\Gamma(1+\delta)^2}
G_{3,3}^{3,3}\!\left(z\;\middle|\;\begin{matrix}\delta,0,0\\0,\delta,\delta\end{matrix}\right).
$$

它对应完整裸对边核的积分，且 $H_0(\delta,1/z)=H_0(\delta,z)$。代入δ=3e得到

$$
J_{02}=z^{3e/2}H_0(3e,z)+\frac{1+z^{3e}}{6e(1+2e)},\qquad
J_{02}(z)=z^{3e}J_{02}(1/z).
$$

按本图m=2组装，最终

$$
\boxed{\begin{aligned}
I_{02}(e,z)={}&-\frac{P(e)}{384e^3(1-2e)^2(1+2e)}\\
&-\frac{P(e)H_0(3e,z)}{128e^2(1-2e)^2\cosh(3eL/2)}.
\end{aligned}}
$$

两项共同构成唯一I。第一项也等于 $3eI_{01}/(2+e)$，可作独立前因子控制。

## 2. 既有M函数的连接

旧H₀的实际Hlog坐标精确满足

$$
H_0(\delta,z)=M_{0,1}+\frac\delta2M_{1,1}
+\delta^2\left(\frac18M_{2,1}+\frac1{24}M_{0,3}\right)+O(\delta^3).
$$

这里M沿用原单对边定义：令 $b=z/(1+z)$、$W=\log\beta+\log(1-\beta)$、$T=L+\log(1-\beta)-\log\beta$，$M_{k,n}=\int_0^1W^kT^n/(\beta-b)\,d\beta$；β=b处为可去奇点。它们均已有已求值Hlog表示，本次绑定实际文件哈希并逐项核对，不新增基函数。特别地 $M_{0,1}=-(L^2+\pi^2)/2$。

## 3. 至有限项

写 $I_{02}=\sum_{k=-6}^0c_ke^k+O(e)$，记 $\zeta_k=\zeta(k)$：

| $k$ | $c_k$ |
|---|---|
| $-6,-5,-4$ | $0$ |
| $-3$ | $-1/384$ |
| $-2$ | $-(2+3M_{0,1})/384$ |
| $-1$ | $-(8+\tfrac32\zeta_2+12M_{0,1}+\tfrac92M_{1,1})/384$ |

有限项为

$$
\begin{aligned}
c_0=-\frac1{384}\Big[&16+3\zeta_2+7\zeta_3+\frac{27}{8}M_{2,1}+\frac98M_{0,3}+18M_{1,1}\\
&+\left(36+\frac92\zeta_2-\frac{27}{8}L^2\right)M_{0,1}\Big].
\end{aligned}
$$

实际三阶极点。H₀(3e)须保留至e²，边界和R须按三阶极点预算展开；U的e³项影响有限项，已经保留。对固定正实z，去除显式极点后的表达在e=0的小邻域解析；选择邻域避开cosh分母的零点，余项为O(e)。不声称极端z与e联合极限的一致余项。

## 颜色重构、数值与证据边界

令 $a=(10-6e)/(3-2e)$、$b=8(1-e)/(3-2e)$。CA²CF、CA·CF·TFnf、CF(TFnf)²三个通道的系数分别为a²I、−2abI、b²I；外部耦合与共同尺度仍按冻结定义提出。本图JSON保存这三种乘积至有限项的完整卷积表达及Hlog坐标，不能只将上述cₖ乘e=0的颜色值。

本批79项精确检查通过，包括实际源核、完整边界、C4组装、旧Hlog系数匹配、14个I系数以及42个颜色通道系数的卷积。48组数值控制/矩评估包含9个完整对边源核扇区求积、3个相邻端点积分、6个已有M矩求值、6次实际prepared C重构和24组I/颜色通道余项检查。对边nf²通道在最初ε范围未达到预设余项比率容差；保持公式与容差不变，将范围缩小至1/2048后通过，首次诊断保留。

原核求积为非认证双精度控制，闭式、M矩与余项使用65位mpmath。一般复z物理延拓未验证。当前没有可调用Wilson MCP，未修改Kernel/Lean，formal=not_run；计算完成不等于形式整链证明。[实际检查](../../../../../workflow/checks/self-energy-pi1-chain-solve/README.md)。

[结果r001](../../results/2e-Pi1+Pi1-02/r001.json)
