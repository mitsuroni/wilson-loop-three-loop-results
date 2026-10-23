# 4e-Pi1-10：I求解工作区

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1/workspaces/4e-Pi1-10/derivation.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。

以下保留迁移前的求解/准备过程记录；其中r001的C/I公式不作为现行数值入口。现行理论及解表已经更新，J与母函数的推导仍可复用。

所属[P3：对边SE＋相邻裸线](../../families/p3-opposite-se/README.md)。本图只使用[主推导](../../theory.md)定义的主稿式(31)的I。

## 固定输入

- 来源：[不可变JSON包](../../sources/4e-Pi1-r001.json)中的`4e-Pi1-10`。
- 源SHA-256：`3b2cab8babd9b1f4630d77d01dc64dffb9b6868c90ccab950a55c4c30230ee0a`。
- 总圈数：3；Wilson颜色迹：`CF*(CF-CA/2)`。
- 以专用prepared C进行a₁/Euler和C4系数提取；完整CA/NF颜色因子位于外部前因子，I不含颜色和尺度。
- 外到内：`s22 → s21 → s31 → s41`。

| 原变量 | 下限 | 上限 |
|---|---|---|
| s22 | 0 | 1 |
| s21 | 0 | s22 |
| s31 | 0 | 1 |
| s41 | 0 | 1 |

## 本图计算及同组关系

采用[唯一主推导](../../theory.md)式(31)，记 $e=\varepsilon$、$L=\log z$、$G=e^{-\gamma_Ee}\Gamma(1-e)$。两成员均有 $n=3,a_1=-1/2,\mathrm{sym}=4$。完整CA/NF因子在外部共同前因子中，沿用prepared中的一次Euler补偿。

## 1. 原有序积分的连接

令 $D_z(x,y)=(1-x)(1-y)+zxy$，完整自能核为

$$
K_e(x,y)=\frac{-2(1+z)D_z+(1-2e)z}{2(1+e)}D_z^{-2+2e}.
$$

03的原核是 $K_e(s_{21},s_{41})z^e(1-s_{31})^{e-1}s_{42}^{e-1}$，排序 $s_{41}<s_{42}$；原外到内为 $s_{21}\to s_{31}\to s_{42}\to s_{41}$。10的原核是 $K_e(s_{22},s_{41})(1-s_{21})^{e-1}s_{31}^{e-1}$，排序 $s_{21}<s_{22}$；原外到内为 $s_{22}\to s_{21}\to s_{31}\to s_{41}$。其他边界均为0、1。

在正实z、Re(e)>0下，先对绝对值作Tonelli估计：有序裸线端点权重被完整区间矩的乘积控制，剩余SE角点径向幂为 $r^{-1+2\Re e}$，可积。因此可以按Fubini重组次序积分裸线端点，而不改变原有序根。

$$
J_{03}=\frac{z^e}{e^2}\int_0^1dx\,dy\,K_e(x,y)(1-y^e),\qquad
J_{10}=\frac1{e^2}\int_0^1dx\,dy\,K_e(x,y)[1-(1-x)^e].
$$

原域反射将03的 $(s_{21},s_{31},s_{42},s_{41})$ 映到10的 $(s_{41},s_{31},s_{21},s_{22})=(1-s_{21},1-s_{31},1-s_{42},1-s_{41})$；排序保持，绝对Jacobian为1。距离和完整分子随 $z\mapsto1/z$ 分别取得z和z²因子，因此 $J_{10}(z)=z^{3e}J_{03}(1/z)$。

## 2. 完整核及边界贡献

精确恒等式为

$$
K_e=-(1+z)D_z^{-1+2e}+\frac{\partial_x\partial_yD_z^{2e}}{4e(1+e)}.
$$

这仅是计算同一个完整核的恒等变换。带权重后导数项不为零，必须保留全部边界贡献。先积x，再使用端点 $D_z(1,y)=zy,D_z(0,y)=1-y$，得到

$$
\begin{aligned}
\int_0^1dx\,dy\,(1-y^e)\frac{\partial_x\partial_yD_z^{2e}}{4e(1+e)}
&=\frac{1}{2(1+e)}\left[\frac{z^{2e}}{6e}+\frac{1}{2e}-B(1+e,2e)\right]\\
&=\frac{z^{2e}+3[1-R(e)]}{12e(1+e)},
\end{aligned}
$$

其中 $R(e)=\Gamma(1+e)\Gamma(1+2e)/\Gamma(1+3e)$。先在严格截断域实施微分/积分，再用Re(e)>0下的端点可积性取极限。没有丢弃接触边界或改变张量算符。

沿用[广义H函数的标准求值](../../../../2ne/families/single-opposite/basis.md)，定义本页的参数缩写

$$
\begin{aligned}
\mathscr H_c(\delta,z)
&=-\frac{(1+z)z^{-\delta/2}}{\Gamma(1-\delta)\Gamma(1+\delta)\Gamma(1+(c+1)\delta)}
G_{3,3}^{3,3}\!\left(z\;\middle|\;\begin{matrix}\delta,0,-c\delta\\0,\delta,\delta\end{matrix}\right),\\
F(e,z)&=\mathscr H_0(2e,z)-\mathscr H_{1/2}(2e,z),\qquad
F^\pm=[F(e,z)\pm F(e,1/z)]/2.
\end{aligned}
$$

半整数参数的原因是SE核指数为−1+2e，而裸线端点权重为yᵉ。由此

$$
J_{03}=\frac{z^{2e}}{e^2}F(e,z)
+\frac{z^e\{z^{2e}+3[1-R(e)]\}}{12e^3(1+e)}.
$$

## 3. 单一完整I

用本图的 $I=G^3[J(z)+z^{3e}J(1/z)]/[32e(1-2e)(1+z^{3e})]$ 组装，得到

$$
\boxed{\begin{aligned}
I_{03}=I_{10}={}&\frac{G^3}{32e^3(1-2e)}
\frac{\cosh(eL/2)F^++\sinh(eL/2)F^-}{\cosh(3eL/2)}\\
&+\frac{G^3}{384e^4(1-2e)(1+e)}
\left[1+3(1-R(e))\frac{\cosh(eL/2)}{\cosh(3eL/2)}\right].
\end{aligned}}
$$

这是完整标准特殊函数/Gamma表示，两项共同定义同一个I。

## 4. 有限项与已求值M系数

沿用原 $M_{k,n}$ 定义，所有所需权重五及以下的M已有精确Hlog结果。令 $b=z/(1+z)$、$W=\log\beta+\log(1-\beta)$、$T=L+\log(1-\beta)-\log\beta$，则 $M_{k,n}=\int_0^1W^kT^n/(\beta-b)\,d\beta$，中间点为可去奇点。积分一个对边端点后，广义H的第q阶系数可直接由多项式

$$
\frac{(A_c+T/2)^{q+1}-(A_c-T/2)^{q+1}}{(q+1)!},\qquad
A_c=\frac{c+1}{2}W+\frac c2(T-L)
$$

逐单项式 $W^kT^n\mapsto M_{k,n}$ 得到。取c=0、2的结果已与旧τ/η至q=3的完整Hlog坐标逐项吻合；本图c=1/2使用同一解析恒等式。没有新增有限项基函数，也没有把未求值积分改名作为答案。

写 $I=\sum_{k=-6}^0c_ke^k+O(e)$：

| $k$ | $c_k$ |
|---|---|
| $-6,-5$ | $0$ |
| $-4$ | $1/384$ |
| $-3$ | $1/384$ |
| $-2$ | $(12+5\pi^2-24M_{1,1})/1536$ |
| $-1$ | $(6L^2M_{0,1}-6M_{0,3}-48M_{1,1}-30M_{2,1}-68\zeta_3+20+5\pi^2)/1536$ |

有限项为

$$
\begin{aligned}
c_0=\frac1{12288}\Big[&96L^2M_{0,1}+(264L^2-768-48\pi^2)M_{1,1}\\
&-96M_{0,3}-104M_{1,3}-480M_{2,1}-152M_{3,1}\\
&-32\pi^2L^2-544\zeta_3+352+120\pi^2+21\pi^4\Big].
\end{aligned}
$$

实际四阶极点。F的常数项为零；其展开到e³，Gamma边界项的解析分子展开到e⁴，才能完整取得有限项。固定正实z下解析Taylor余项经相应极点相除后为O(e)。不声称极端z联合极限的统一余项。

## 5. 检查与范围

28项精确检查包含完整源核、导数恒等式、Beta边界、C4组装、12个已有τ/η的Hlog比对和7个系数的crossing对称性；域和源哈希另作断言。9组数值控制覆盖两个带权完整SE核的128/256阶扇区求积及由M定义独立得到的有限项余项缩放。55位mpmath用于Meijer G和M矩；浮点求积没有认证区间界。

共同起始域 $z>0,v_{14}<0,0<\Re e<1/2$，随后向e=0附近亚纯延拓；一般复z物理延拓未验证。无可调用Wilson MCP，未修改Kernel/Lean，formal=not_run。

[03结果](../../results/4e-Pi1-03/r001.json) · [10结果](../../results/4e-Pi1-10/r001.json) · [计算证据](../../../../../workflow/checks/self-energy-pi1-p3/README.md)
