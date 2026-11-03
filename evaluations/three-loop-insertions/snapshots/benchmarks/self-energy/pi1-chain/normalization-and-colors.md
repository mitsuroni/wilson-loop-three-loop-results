# 串联积分的I结构与颜色通道

现行数学约定为[自能2.0](../conventions.md)，输入和结果见[当前清单](../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1-chain/normalization-and-colors.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。

以下保留迁移前的求解/准备过程记录；其中r001的C/I公式不作为现行数值入口。现行理论及解表已经更新，J与母函数的推导仍可复用。

本页解释[唯一理论主稿](theory.md)的SC20–SC23，并落实到[专用prepared](prepared/README.md)，不另定义第二个I。原图贡献是 $\mathcal A=C^{\rm src}J$；最终I是经过颜色、耦合、尺度和C4组装后提取的标量。

## 1. 一个I与三个颜色通道

记 $e=\varepsilon$、$z=z_1=v_{12}/v_{14}$，

$$
a(e)=\frac{10-6e}{3-2e},\qquad b(e)=\frac{8(1-e)}{3-2e},\qquad
\mathcal B_1=aC_A-bT_Fn_f.
$$

完整算符颜色是

$$
C_F\mathcal B_1^2=a^2C_A^2C_F-2abC_AC_FT_Fn_f+b^2C_F(T_Fn_f)^2.
$$

按每个味道次数分别使用原精确颜色重构器，得到：

| 味道次数 | 完整有序基底 | 完整坐标 | 提出算符标量后的首坐标 |
|---|---|---|---|
| 0 | $(C_A^2C_F,C_AC_F^2,C_F^3)$ | $(a^2,0,0)$ | $1$ |
| 1 | $T_Fn_f(C_AC_F,C_F^2)$ | $(-2ab,0)$ | $1$ |
| 2 | $(T_Fn_f)^2(C_F)$ | $(b^2)$ | $1$ |

重构使用关系 $C_A^2=1+2C_AC_F$，不是大Nc近似或字面选项。此处没有任何非零颜色方向被丢弃。**只有把整个B₁²提出后，剩余a₁才是1**；未提出时，各味道分区的首坐标分别为a²、−2ab、b²。不能移用单Π₁加裸线的a₁=−1/2。

在e=0时，三系数为100/9、−160/9、64/9（对应含TF的基底）。TF=1/2时，若改用不含TF的nf记号，中间和最后系数分别变为−80/9、16/9。prepared仍保留完整e依赖，不用这些零阶值替代它们。

## 2. C/J到唯一I

定义

$$
\begin{aligned}
G(e)&=e^{-\gamma_Ee}\Gamma(1-e),\\
U(e)&=\frac{\Gamma(1+e)^2\Gamma(1-e)\Gamma(1-3e)}{\Gamma(1-2e)^2\Gamma(1+2e)},\\
R(e)&=-\frac{G(e)^3U(e)}{128e^2(1-2e)^2},\\
\mathcal N(e)&=\left(\frac{\alpha_s}{\pi}\right)^3C_F\mathcal B_1(e)^2,\\
S_{14}&=(-2\pi v_{14}+i0)^{3e},\qquad S_{12}=z^{3e}S_{14}.
\end{aligned}
$$

专用准备只在源C上加入一次 $e^{-3\gamma_Ee}$，得到 $C=\mathcal N S_{14}R$。原C及旧normalization保存在`original_C`和`original_normalization`；J及原域逐字段相同。C中仍未装配轨道重数。

两图的m分别是相邻4、对边2。唯一I满足

$$
\frac m2 C[J(z,e)+z^{3e}J(z^{-1},e)]
=\mathcal N(S_{14}+S_{12})I(z,e),
$$
$$
\boxed{I(z,e)=\frac m2R(e)
\frac{J(z,e)+z^{3e}J(z^{-1},e)}{1+z^{3e}}.}
$$

I不含CA、CF、TF、nf、耦合或参考尺度。两个几何各有自己的I，**不是每个颜色通道分别计算一个参数积分**。完整贡献可重构为

$$
\left(\frac{\alpha_s}{\pi}\right)^3(S_{14}+S_{12})
\left[a^2C_A^2C_F-2abC_AC_FT_Fn_f+b^2C_F(T_Fn_f)^2\right]I.
$$

因而三个单项式的系数分别是a²I、−2abI和b²I。它们的Laurent系数通常并不相同，因为a、b依赖e；若要输出颜色分解的有限项，必须作卷积，而不能只乘e=0的颜色值。

若某一图 $I=\sum_{j=-p}^0c_je^j+O(e)$，某颜色标量 $\phi(e)=\sum_{k\ge0}\phi_ke^k$，则

$$
\operatorname{coeff}_{e^j}[\phi I]=\sum_{k=0}^{j+p}\phi_kc_{j-k}.
$$

求该颜色通道的有限项要保留φ至eᵖ。这只是对同一个I的颜色组装，不增加独立积分任务。

## 3. Gamma比与展开预算

U含有Γ(1−e)等因子，不能仅通过字面数Γ幂来判断Euler补偿。准备器检查原始Gamma表达式与分解后的C等价，并验证

$$
\left.\partial_e\log U\right|_{e=0}=0,\qquad
\left.\partial_e\log[\Gamma(1-e)^3U]\right|_{e=0}=3\gamma_E.
$$

因此补偿指数必须使用总圈数3，而不是插入数2。U的展开为

$$
U=1+6\zeta_3e^3+9\zeta_4e^4+O(e^5).
$$

相邻图的四阶极点会使U的e³项影响简单极点、e⁴项影响有限项；对边图的三阶极点会使U的e³项影响有限项。不能将U简化为1后再做有限项。

共同R从e⁻²开始，因此 **J统一需要计算至e²**，区别于此前单插入的e¹预算。

| 任务 | J极点预算 | I极点预算 | J最高保留项 | $e^2R$最高保留项 | 颜色标量最高保留项 |
|---|---:|---:|---:|---:|---:|
| CH_A：01，相邻 | 2 | 4 | $e^2$ | $e^4$ | $e^4$ |
| CH_O：02，对边 | 1 | 3 | $e^2$ | $e^3$ | $e^3$ |

这些预算由完整核的端点结构以及理论控制得到。相邻闭式已作为SC24控制式给出；对边角点的径向行为是 $r\,dr\,r^{-2+3e}$，给出一个端点极点。对边导数恒等式需保留其矩形边界，不能因二点几何就逐图删去。

## 4. 当前完成边界

prepared原输入、C/J、域和哈希保持。两图现均已完成全ε及有限项，正式结果2/2；14个I系数位置（9非零），并保存三个颜色通道至有限项的42个位置（27非零）。[结果总表](families/current-results.md)。

相邻I四阶极点，对边I三阶极点，和上述展开预算一致。79精确检查、48数值控制/矩评估通过；旧准备阶段文件及检查记录保留原字节，当前结果目录为catalog-r002。默认Kernel根未扩展，formal=not_run，一般复z物理延拓未验证。
