# 双TGV的紧凑与解析降维初始C/J

2026-10-09。`--double-tgv-forms`在26图的原始C/J旁追加两个可选表示：`J_compact`保留显式因式与公共定义，`J_reduced`把部分积分实际积完。根schema仍为1.3，各新块独立版本化；不改变旧raw、可选`J_ibp`、颜色、Euler/C4账本或积分求值状态。

```powershell
D:/py/python.exe -m montopar.cli --source ../qgraf-pipeline/outputs/3loop/orbit_definitions.json --type "4e-3g+3g" --double-tgv-forms --package-initial
```

可以同时添加`--tgv-ibp`。新开关只接受双TGV目标，且不能与`--package-only`一起用。打包已有新格式时只需`--package-only --package-initial`。Python接口为`build_double_tgv_result(orbit, forms_enabled=True)`。

## 紧凑原核

`J_compact`版本为`double-tgv-compact-1.0`。`definitions`是有序的显式标量定义，每行包含`symbol`和`expression`；只能引用此前定义及该块的积分变量/参数。保留六个距离、两侧α和、各外点梯度、两顶点收缩，以及D=4−2ε的Hessian。

`fields`含B5、U、T1_hat、T2_hat、pnum1、pnum2。它们与原`J_raw_s`对应字段完全相等，可用`montopar.double_tgv_forms.polynomial_fields(block)`通过已有多项式解析器展开，不需要CAS。该辅助函数拒绝重复、前向和未定义的内部符号。`dtg_`为新表示保留前缀，源路径变量不得使用该前缀。

原路径顺序、依赖上下界、五α单纯形均直接继承。`integrand`消费上述显式定义；`geometric_integrand_before_color_screening`保留颜色筛选前核。`regulated_geometric_integrand`另外给出共同有限辅助调节形式B+iη_reg U；η_reg是无量纲调节量，不是α5。

生产端仍导出原展开字段作兼容基准；紧凑块缩短求解输入，不代表原式积分已经求出或生产速度必然提高。

## 通用中央参数积分

设λ=1−3ε，换元

$$\alpha^{\rm simplex}=\frac{(ar,a(1-r),(1-a)t,(1-a)(1-t),c)}{1+c},\qquad h=a(1-a).$$

绝对Jacobian为h/(1+c)⁵，原密度齐次度−5消去投影总尺度。定义

$$L=r(1-r)d_{12},\quad R=t(1-t)d_{34},$$
$$D_c=rt d_{13}+r(1-t)d_{14}+(1-r)t d_{23}+(1-r)(1-t)d_{24},$$
$$Q_0=aL+(1-a)R,\quad Q_1=a^2L+(1-a)^2R+hD_c.$$

于是U=h+c、B=hQ₀+cQ₁。原完整分子满足

$$P_1=N_0+cN_1+c^2N_2,\qquad P_2=cM_1.$$

代码从gᵢ=gᵢ₀+cgᵢ₁重建N₀/N₁/N₂，保留向量收缩；M₁来自完整D维Hessian，无需先展开长多项式。

令n为分母相对λ的增加量，m为c幂次。四个可能的(n,m)为(2,0)、(2,1)、(2,2)、(1,1)，对应系数N₀、N₁、N₂、M₁。逐项中央积分为

$$\begin{aligned}
&h^{1-2\epsilon}\int_0^\infty dc\,
c^{m-\epsilon}(h+c)^{-1-2\epsilon}(hQ_0+cQ_1)^{-\lambda-n}\\
&=h^{m-n-2\epsilon}Q_0^{-\lambda-n}
\mathrm B(m+1-\epsilon,n-m+1)\\
&\quad\times{}_2F_1\left(\lambda+n,m+1-\epsilon;n+2-\epsilon;1-Q_1/Q_0\right).
\end{aligned}$$

这是c=hx/(1−x)后的Euler积分。使用普通Gauss₂F₁；[DLMF 15.6.1](https://dlmf.nist.gov/15.6.E1)写的是正则化函数，不能省掉上述Beta系数。

每个输出项把变量无关的Beta放入`relative_C`，`combined_C=原C×relative_C`；其余权重及原四个路径积分、a/r/t保留在J中。维数从8降为7。新变量在原路径之后积分，保留原四路径的外到内顺序与所有依赖边界。

固定其余参数在内部且Q₀,Q₁>0时，各中央积分在Re ε<1成立。外层须在共同收敛域中解释，或对完整组合继承原亚纯延拓；这不声称四个分项分别具有同一外层条带。Q₀=Q₁时₂F₁=1；Q₀→0须取完整Euler/Gauss边界值，不能分别代入发散因子。代码对双同边的Q₀≡0情形直接走下节。

## 双同边的四维边界和

仅当左两腿同边、右两腿同边、两边不同，并且原域正好是两条独立有序链时适用。当前非零族中是04与19，判断不依赖图号；相同占据数的06、21不能走这条支路。

每对原有序参数使用

$$\int_{0<s<t<1}(\partial_s-\partial_t)f
=\int_0^1dr[2f(r,r)-f(0,r)-f(r,1)].$$

若原顶点槽位与低→高顺序相反，相应算符贡献反号，源颜色也必须保留其原值。两个顶点的九面系数是(2,−1,−1)的外积，不丢接触面。

两侧α总量积分给G(ε)=B(ε,ε)/(2ε)。设d=u_L·u_R/v_ref，两个边界算符到原J的系数为d/[2(1−3ε)]。故每项

$$C_{ij}=C_{\rm raw}\,
\frac{\sigma_L\sigma_R c_i c_j\,d\,B(\epsilon,\epsilon)}{4\epsilon(1-3\epsilon)},$$
$$J_{ij}=\int_0^1dr_Ldr_Rdxdy\,
[x(1-x)y(1-y)]^{-\epsilon}
d(\bar s_L,\bar s_R)^{-1+3\epsilon}.$$

这里距离函数d(·,·)与上式切向内积d不同，输出使用显式`distance`避免记号混淆；\(\bar s_L=x s_1+(1-x)s_2\)、右侧类似，随后限制到该项实际边界。每项有独立C、核、四个变量、边界映射及排序方向。04进一步积分可恢复原Gamma结果，但此输出仍是待求的初始积分。

## 消费合同与调节边界

`J_reduced`版本为`double-tgv-reduced-1.0`，含共享有序`definitions`和`terms`。每项保存`relative_C`、`combined_C`、`integrand`、筛选前核、变量、上下界、外到内次序、Jacobian说明和变换来源。求和规则是

$$C J_{\rm raw,s}=\sum_t \mathrm{combined\_C}_t\int J_t.$$

消费时同时读取块级定义和项级域，不要再乘一次原C。所有C均未加Euler、C4、crossing或a1投影；几何降维没有重新给边界图计算颜色。零颜色图保持有效J为零，并保留原几何。

紧凑原核可逐点处理共同有限调节。解析降维的上述公式则以正实起始区/共同边界值为条件：原有限调节在投影后变成B+iη_reg U(1+c)，对c一般二次，不能把Gauss线性c公式声称为任意有限η_reg恒等式。`equivalence.finite_auxiliary_regulator_identity=false`明确保存这一边界。

`validate_double_tgv_result`重建并逐字段验证两种可选块；打包保留新块并运行同一校验。修改系数、域、调节声明、Gauss参数或删掉一块均会拒绝。已有Stage 3若只认识原始积分字符串，可以继续使用原`J_raw_s`；使用新表示需显式支持定义替换、分项C和普通Gauss函数。
