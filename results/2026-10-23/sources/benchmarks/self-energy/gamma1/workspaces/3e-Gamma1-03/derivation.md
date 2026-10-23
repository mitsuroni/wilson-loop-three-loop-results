# 3e-Gamma1-03：244 对边两颜色积分

## 当前物理积分的入口（2026-10-06）

从共同作用量出发，U′=igAU确定真实Wilson词；FT/Dyson规则在提取C前生效。当前[原生C/J](../../sources/montopar-r002.json)给出下式，J就是本文实际积分的原核，域、分支和IBP全保留：

$$I_{P,c}=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-3\epsilon)}{64}\frac{J_c(z)+z^{3\epsilon}J_c(1/z)}{1+z^{3\epsilon}}.$$

用下文已经求出的J代入此定义，得到[当前物理I](../../results/3e-Gamma1-03/r002.json)；组装系数为+1，Euler只补一次。以下带日期的求值过程与r001公式属于旧根历史，计算J的步骤可按同域积分线性复用，旧形式状态不继承。

2026-10-04：**两个颜色的全部极点至有限项已解析求出**：
$$I_{03,c}(\epsilon,z)=\frac{i_{c,-2}(z)}{\epsilon^2}+\frac{i_{c,-1}(z)}{\epsilon}+i_{c,0}(z)+O(\epsilon),\quad c=g,q.$$
最高极点为二阶，$\epsilon^{-6},\ldots,\epsilon^{-3}$ 全为零。六个系数的最终表示统一使用自变量 $z$、字母 $\{0,-1\}$ 的有限GPL组合，完整表达式见附录B及 [z变量结果](../../results/3e-Gamma1-03/z-r001.json)。原辅助r结果保留于 [r001](../../results/3e-Gamma1-03/r001.json)。03保留运动学依赖，满足 $I(z)=I(1/z)$。

起点为正式 [Montopar raw r001](../../sources/montopar-r001.json) 的03轨道，SHA256 `3b92e23125cb6b654ad63d5c19c17fead9cfeb3bff56739dd9f07c6eda87e65a`。动量源→C/J的既有核验见 [独立审查](../../momentum-to-parametric-review.md)，本轮实际计算见 [gamma1-03检查](../../../../../workflow/checks/gamma1-03/README.md)。完成声明针对固定raw与维数正规化约定下的有限阶；全ε闭式及形式证明未完成，`formal=not_run`。

## 1. 初始积分和唯一规范化

三个外接点在 $(2,1),(4,1),(4,2)$。令 $x=s_{21},y=s_{41},w=s_{42}$，原外到内顺序为
$$\int_0^1dx\int_0^1dw\int_0^wdy.$$
内部参数 $a=\alpha_1,b_\alpha=\alpha_2,c_\alpha=1-a-b_\alpha$ 及 $t_3=1-t_1-t_2$ 分别在两个二维单纯形上。定义
$$D_z(x,v)=(1-x)(1-v)+zxv,\quad q=a[b_\alpha D_z(x,y)+c_\alpha D_z(x,w)],$$
$$P_\alpha=ab_\alpha c_\alpha,\quad m=\max(at_2t_3,b_\alpha t_1t_3,c_\alpha t_1t_2),\quad\sigma=3\epsilon.$$
源的夸克七维无量纲积分为
$$J_q=\int_{x,w,y}\int_{\Delta_\alpha}da\,db_\alpha\int_{\Delta_t}dt_1dt_2
\left[-\frac{P_\alpha^{-\epsilon}m^{-1-\epsilon}}{3\epsilon(1+\epsilon)}K_{q,N}
+\frac{P_\alpha^{-\epsilon}m^{-\epsilon}}{\epsilon}K_{q,L}\right].$$
$K$ 为完整源偏导算符实际作用后的核，具体式见下一节及附录A。胶子通道含相同权重的胶子三角核及原剑鱼泡图，不能只计算一个标量三角母积分。原路径IBP有 $y=0,w=1,y=w$ 三面，体余项为零；先在独立坐标中求导再代入面，保留接触面的法向差。

两颜色前因子为
$$C_c=\lambda_b^3\mathrm{color}_c\frac{(\pi A_0)^{3\epsilon}\Gamma(1-3\epsilon)}{128},\quad
A_0=-2v_{14},\quad z=v_{12}/v_{14}>0,$$
$$\mathrm{color}_g=C_FC_A^2,\qquad\mathrm{color}_q=C_FC_AT_Fn_f.$$
raw不含Euler补偿与C4重数。唯一 $I$ 为
$$I_c=\frac{G_3}{64}\frac{J_c(z)+z^{3\epsilon}J_c(1/z)}{1+z^{3\epsilon}},\qquad
G_3=e^{-3\gamma_E\epsilon}\Gamma(1-3\epsilon).$$
反射 $x\mapsto1-x,y\mapsto1-w,w\mapsto1-y$，同时 $\alpha_2\leftrightarrow\alpha_3,t_2\leftrightarrow t_3,z\mapsto1/z$，保持原有序域和绝对Jacobian。原核十项精确检查（含两个剑鱼项互换）给出
$$J_c(z)=z^{3\epsilon}J_c(1/z),\qquad\boxed{I_c=\frac{G_3J_c(z)}{32(1+z^{3\epsilon})}.}$$
反射在共同有限调节和内部对称截断上核对，再作共同亚纯延拓。03不是01的单尺度反射。

## 2. 内部分子分解与有序路径降维

本图 $n_2=n_3,n_2^2=0$，$n_1\cdot n_2=(1+z)/2$。完整算符字典与01相同，因此可复用一般外动量分子恒等式：
$$\boxed{J_g=-\frac{1-\epsilon}{2}J_q+J_p+J_R.}$$
$J_p$ 含pinch与原剑鱼，不能再加一遍剑鱼；$J_R$ 是内部圈分子余项，与已为零的路径IBP体余项不同。分解和完整 $R$ 算符见 [decomposition-r001.json](../../../../../workflow/checks/gamma1-03/decomposition-r001.json)。这里复用的是一般分子等式，不是01的积分值。

令 $a=A,b_\alpha=(1-A)B,c_\alpha=(1-A)(1-B)$，$S=A(1-A)$，并令
$$V=By+(1-B)w,\quad d=w-y,\quad y=V-(1-B)d,\quad w=V+Bd.$$
路径Jacobian为1，原有序域成为 $0<V<1,0<d<f_0(V)$，其中
$$f_0(V)=\min\left(\frac{V}{1-B},\frac{1-V}{B}\right),\qquad q=SD_z(x,V).$$
令 $X=(1+z)x-1$，三组实际核均有精确正常形
$$K=(1+z)XP_2q^{\sigma-2}
+X[(1+z)D_zP_D+zP_z+(1+z)XdP_{\rm sep}]q^{\sigma-3}.$$
$P$ 只依赖 $A,B,t_1,t_2,\epsilon$，见附录A。九个恒等式逐项还原源核。

$d$ 积分给 $f_0$ 及 $f_1=f_0^2/2$。结合 $\partial_VD=X$ 和
$$zD^{\sigma-2}=\frac{\sigma(1+z)D^{\sigma-1}-\sigma^{-1}\partial_x\partial_VD^\sigma}{\sigma-1}$$
可降维为分布导数。$V_c=1-B$ 上
$$f'_0=\begin{cases}(1-B)^{-1}&V<V_c,\\-B^{-1}&V>V_c,\end{cases}\quad
f''_1=\begin{cases}(1-B)^{-2}&V<V_c,\\B^{-2}&V>V_c\end{cases}-\frac{\delta(V-V_c)}{B(1-B)}.$$
最后接触项不能丢。提取 $S^{3\epsilon}$ 后三个内部系数为
$$C_F=-\frac{P_2+P_D/S+\sigma P_z/[S(\sigma-2)]}{S^2(\sigma-1)},\quad
C_T=\frac{P_{\rm sep}}{S^3(\sigma-1)(\sigma-2)},\quad
C_E=\frac{P_z}{\sigma S^3(\sigma-1)(\sigma-2)}.$$
乘回原 $N/L$ 权重和单纯形Jacobian、积内部变量，得到 $F_q,F_R,T_R,E_q,E_R$。在 $B\mapsto1-B,t_2\mapsto t_3$ 下，$F,E$ 为奇，$T$ 为偶；十二个精确检查支持补全两半参数域。

$x$ 直接积完给
$$H_\epsilon(V)=\int_0^1D_z(x,V)^{3\epsilon-1}dx
=\frac{(zV)^{3\epsilon}-(1-V)^{3\epsilon}}{3\epsilon[(1+z)V-1]}.$$
其余路径密度为
$$\Phi_F(V)=\int_0^V\frac{F(B)}{1-B}dB+\int_0^{1-V}\frac{F(B)}{1-B}dB,$$
$$\Phi_T(V)=\int_0^V\frac{T(B)}{(1-B)^2}dB+\int_0^{1-V}\frac{T(B)}{(1-B)^2}dB-\frac{T(V)}{V(1-V)}.$$
令 $\Phi_q=\Phi_{F_q},\Phi_R=\Phi_{F_R}+\Phi_{T_R}$，则两者均关于 $V\leftrightarrow1-V$ 对称。$C_E$ 的端点给
$$\mathcal E_c=\int_0^1E_c(B,\epsilon)\frac{1-B^{3\epsilon-1}}{1-B}dB,$$
$$\boxed{J_c=J_c^{\rm core}+(1+z^{3\epsilon})\mathcal E_c,\quad
J_c^{\rm core}=(1+z)\int_0^1H_\epsilon(V)\Phi_c(V,\epsilon)dV},\quad c=q,R.$$
这是完整普通核的端点降维，不再额外添加源内IBP面。原路径和上述降维的六个独立基核数值检查通过，含接触项及 $\mathcal E$，见 [路径检查](../../../../../workflow/checks/gamma1-03/path-reduction-checks.json)。

## 3. 内部参数积分与真实极点预算

固定 $b=B$，令 $u=t_2/(1-t_1)$，先取 $u\le b$，以奇偶性补另一半。换元
$$u=bw_0,\quad\frac{t_1}{1-t_1}=t,\quad\frac A{1-A}=\frac{t\rho}{w_0},\quad
A=\frac{t\rho}{w_0+t\rho},\ t_2=\frac{bw_0}{1+t}.$$
原 $da\,db_\alpha\,dt_1dt_2$ 固定 $b$ 的Jacobian是 $btw_0^2/[(w_0+t\rho)^3(1+t)^3]$。有 $m_1/m_2=\rho$ 及
$m_2-m_3=(1-A)t_1(1-t_1)b(1-w_0)\ge0$。故 $\rho<1$ 取 $m_2$，$\rho>1$ 取 $m_1$；高扇区取 $\rho=1/k$，包含绝对Jacobian $k^{-2}$。低扇区取 $\rho=k$。最后取 $k=U/(1+U),w_0=W/(1+W)$。

共同旁观因子是 $[b^2(1-b)]^{-\epsilon}$，两个正谱底数为
$$B_{\rm lo}=\frac{tU^2(1+t)^2(1+W)^3}{[W(1+U)+tU(1+W)]^2[1+(1-b)W]},$$
$$B_{\rm hi}=\frac{tU(1+U)(1+t)^2(1+W)^3}{[WU+t(1+W)(1+U)]^2[1+(1-b)W]}.$$
十个精确有理核乘 $B^\epsilon$ 的表达式见 [profiles-r001.json](../../../../../workflow/checks/gamma1-03/profiles-r001.json)。Jacobian、最大值比较和两个谱底数共四项核验通过；拆开的所有对数底数均正。

固定 $b\in(0,1)$ 的真实发散度为：

| 密度 | 低扇区 | 高扇区 | 最终所需展开 |
|---|---|---|---|
| $F_q$ | $U:2\epsilon$ | 无新增变量发散 | −2至0 |
| $F_R,T_R$ | $U:2\epsilon,t:\epsilon,\{U,t\}:3\epsilon$ | $t:\epsilon$ | −2至0 |
| $E_q$ | 无 | 无 | −1至1 |
| $E_R$ | 无 | 无 | 0至1 |

不存在零发散度拆项。先对低扇区 $U$ 用度数 $2\epsilon$ 的Euler恒等式，积 $U$ 而保留外部指数；对 $R$ 再用 $t$ 的 $(\epsilon-1-t\partial_t)/\epsilon$。若前一步为 $e^{\epsilon L_{\rm rest}}\sum A_p\epsilon^p$，则
$$G_p=-A_{p+1}-t\partial_tA_{p+1}+(1-t\partial_tL_{\rm rest})A_p,$$
$$L_{\rm rest}=\ln t+2\ln(1+t)+3\ln(1+W)-\ln[1+(1-b)W].$$
所以 $R$ 在 $U$ 阶段需多算一阶至ε¹；$E_q$ 的边界Beta极点也需 $K_q$ 至ε¹。逐变量实际解析积分 $U,t,W$，恢复旁观指数 $-2\ln b-\ln(1-b)$，按奇偶性补全 $b$，再积 $b$。发散检查全程开启，不过滤奇异项。

六个完整路径密度见 [core_density.txt](../../../../../workflow/checks/gamma1-03/core_density.txt)，首项为
$$\Phi_{q,-2}=-\frac83-\frac43\ln V-\frac43\ln(1-V),$$
$$\Phi_{R,-2}=\frac72-\frac54\zeta_2+\frac14\ln[V(1-V)]
-\frac18\ln^2[V(1-V)]-\ln V\ln(1-V).$$

## 4. 端点常数：先减除再展开

恢复旁观指数前的下半密度记为 $\widehat E(b,\epsilon)$，令 $K(\epsilon)=\widehat E(1,\epsilon)$。它是普通收敛极限：原 $qE$ 的 $W$ 尾为 $W^{-3+3\epsilon}$，$RE$ 为 $W^{-2+3\epsilon}$，在零点邻域可积。

用保持奇性的显式项
$$E_{\rm ct}=K(\epsilon)[b^{1-2\epsilon}(1-b)^{-\epsilon}-(1-b)^{1-2\epsilon}b^{-\epsilon}]$$
抽出全部边界质量；余项两端至少线性消失，允许对数和小ε幂，可普通逐阶积分。减除项积分为 $KC_1$，
$$C_1=B(2-2\epsilon,-\epsilon)-B(1+\epsilon,-\epsilon)
-B(1-\epsilon,1-2\epsilon)+B(2\epsilon,1-2\epsilon)$$
$$=\pi\csc(\pi\epsilon)+\pi\csc(2\pi\epsilon)
-\frac{1-\epsilon}{\epsilon}B(1-\epsilon,1-2\epsilon)
=\frac1{2\epsilon}-2+(5\zeta_2-6)\epsilon+O(\epsilon^2).$$
四个Beta项按共同亚纯延拓解释，不能把单项均说成同一正ε区间上的普通积分。

实际求得
$$K_{q,-1}=\frac89\zeta_2,\ K_{q,0}=-\frac{16}9\zeta_3-\frac8{27}\zeta_2,\quad
K_{q,1}=-\frac4{45}\zeta_2^2+\frac{16}{27}\zeta_3+\frac{272}{81}\zeta_2,$$
$$K_{R,0}=-\frac23\zeta_2,\qquad K_{R,1}=\frac43\zeta_3.$$
首项还用直接最大值积分独立验证：
$$\int_0^1\frac{1-a}{\max(at_2,(1-a)t_1)}da
=\frac{\ln[(t_1+t_2)/t_1]}{t_2},$$
再取 $T=t_1+t_2,v=t_1/T$，余下 $\int_0^1[-\ln v/(1-v)]dv=\zeta_2$。最终完整常数为
$$\mathcal E_q=\frac{4\zeta_2}{9\epsilon^2}
+\frac{4/9-4\zeta_2/27-32\zeta_3/9}{\epsilon}
+2-\frac8{81}\zeta_2+\frac{272}{27}\zeta_3-\frac65\zeta_2^2+O(\epsilon),$$
$$\mathcal E_R=-\frac{\zeta_2}{3\epsilon}+\frac73\zeta_3+\frac19+O(\epsilon).$$
改用 $b^{2-2\epsilon}(1-b)^{-\epsilon}$ 及其镜像减除，$C_2-C_1=\epsilon/2+O(\epsilon^2)$，余项恰抵消变化；四项精确检查确认结果不依赖该辅助选择。

极点预算还需控制非一致边界。原核关于 $b$ 的角点分母仅有一阶 $1-bw_0$，$U,t$ 的Euler不提高此阶；$F$ 至多有小ε幂与对数，$T$ 仍带 $b(1-b)$，所以接触项 $T(V)/[V(1-V)]$ 没有整数 $V^{-1}$。$E-E_{\rm ct}$ 乘端点权重后也可积。取足够小的ε邻域，可用端点小幂乘对数作一致可积控制，剩余核心可以逐阶积分。全部边界Beta质量已在 $KC_1$ 中，不能额外产生三至六阶极点。此范围为参数积分的共同延拓；源侧全Minkowski分布去调节桥仍开放。

## 5. 最后一维实际求完，处理pinch和原剑鱼

令 $r=1/(1+z)\in(0,1)$，$\ell=\ln z=\ln(1-r)-\ln r$，$A_V=\ell+\ln V,B_V=\ln(1-V)$。最后六个核心系数为
$$J^{\rm core}_{c,n}=\int_0^1\frac{dV}{V-r}
\sum_{p=-2}^n\Phi_{c,p}(V)\frac{3^{n-p}}{(n-p+1)!}
[A_V^{n-p+1}-B_V^{n-p+1}],\quad n=-2,-1,0.$$
$V=r$ 时 $A_V=B_V$，假极点可去。六个积分均已实际求完并化为字母0、1的GPL，没有仅把它们命名为新母积分。

pinch加原剑鱼可先完成内部Beta积分。令 $(A_j,B_j)=(1,1),(2,1),(1,2)$，
$$S_j=B(\epsilon,1+\epsilon),1/(2\epsilon),1/(2\epsilon),\quad
c_y=(-3/4,2+\epsilon,-5/4),\ c_w=(3/4,5/4,-2-\epsilon),$$
$$\mathcal K_{A,B}(V)=\int_0^V\beta^{-A\epsilon}(1-\beta)^{-1-B\epsilon}d\beta,$$
$$F_p=\sum_j\epsilon S_j[c_{w,j}\mathcal K_{B_j,A_j}(V)-c_{y,j}\mathcal K_{A_j,B_j}(1-V)
+(c_{y,j}-c_{w,j})B(1-A_j\epsilon,1-B_j\epsilon)].$$
于是
$$J_p=\frac{2B(1-\epsilon,1-\epsilon)}{\epsilon^2}\int_0^1\frac{F_p(V,\epsilon)}{V-r}
\frac{(zV)^{3\epsilon}-(1-V)^{3\epsilon}}{3\epsilon}dV.$$
其内部展开为
$$\mathcal K_{A,B}=-G_1+\epsilon(AG_{1,0}+BG_{1,1})
-\epsilon^2[A^2G_{1,0,0}+AB(G_{1,0,1}+G_{1,1,0})+B^2G_{1,1,1}]+O(\epsilon^3).$$
三个系数矩 $m_0,m_1,m_2$ 全部积完；
$$J_{p,-2}=2m_0,\quad J_{p,-1}=2(m_1+2m_0),\quad
J_{p,0}=2[m_2+2m_1+(4-\zeta_2)m_0].$$
独立首项桥是
$$I_{p,-2}=\frac3{256}(M_{11}+2M_{01}),\quad
I_{q,-2}=\frac1{48}(M_{11}+2M_{01})+\frac{\zeta_2}{72},$$
其中 $M_{11}+2M_{01}$ 在本变量下为
$$-\int_0^1\frac{[2+\ln(V(1-V))][\ell+\ln V-\ln(1-V)]}{V-r}dV.$$
这个组合与既有对边矩的变量反射一致；最终系数已展开为GPL，不依赖该未求值记号。

## 6. 最终组装、验证和范围

最高极点仅二阶，规范化需要
$$G_3=1+\frac92\zeta_2\epsilon^2+O(\epsilon^3),\quad
\frac1{1+z^{3\epsilon}}=\frac12-\frac34\ell\epsilon+O(\epsilon^3),$$
其中后一式的ε²系数为零。对 $c=q,R$，
$$i_{c,n}=\frac{J^{\rm core}_{c,n}}{64}+\frac{\mathcal E_{c,n}}{32}
-\frac{3\ell J^{\rm core}_{c,n-1}}{128}
+\delta_{n0}\left[\frac{9\zeta_2J^{\rm core}_{c,-2}}{128}+\frac{9\zeta_2\mathcal E_{c,-2}}{64}\right].$$
缺失下标置零，$I_p$ 使用同一个 $G_3/[32(1+z^{3\epsilon})]$。最后
$$i_{g,n}=-i_{q,n}/2+i_{q,n-1}/2+i_{p,n}+i_{R,n}.$$
原raw $J_g,J_q$ 系数也单独导出。

| 系数 | 胶子 $z=1$ | 夸克 $z=1$ | 胶子 $z=2$ | 夸克 $z=2$ |
|---|---:|---:|---:|---:|
| ε⁻² | −0.0266402700305864 | 0.0474022069803852 | −0.0281736034518368 | 0.0497946743262708 |
| ε⁻¹ | 0.165829139873195 | −0.276577133957307 | 0.173456163348481 | −0.286230824637187 |
| ε⁰ | −0.719038967657415 | 1.00495910153475 | −0.740477922366513 | 1.01537200690734 |

实际通过：原核反射10、正常形9、谱换元4和奇偶性12项；原路径降维6个数值基核；独立首项6项及pinch的M桥；更换端点减除4项；最终crossing6项。独立二进制GPL级数和尾零shuffle在75位精度下直接积分核心18点、pinch9点，再组装两颜色18点。核心最大相对误差 $1.3\times10^{-39}$，最终组装 $1.5\times10^{-40}$；24/40节点比较支持收敛。数值检查对象是降维后的系数积分，不宣称独立七维原积分或形式证明。

03的两个颜色全部极点至有限项已无未求积分或未知常数。全ε标准函数闭式、一般复运动学延拓、源侧严格全域去调节、Wilson/Lean整链、完整三圈物理组装/反项/RG仍开放。04三边几何不满足重复方向条件，需另行推导。

## 附录A：核的精确多项式

`**` 表示幂，变量均已在正文定义；多项式由正式普通核精确提取。

### raw_q_N0

```text
P2 = 48*A*eps*(A - 1)**2*(eps - 1)*(3*eps - 1)*(B*t1 - B + t2)*(A*B**2*t1**2 + 2*A*B*t1*t2 - A*B*t1 + A*t2**2 - A*t2 - B**2*t1**2 + B*t1**2)
PD = 0
Pz = 48*A**3*eps*(A - 1)**3*(3*eps - 2)*(3*eps - 1)*(B*t1 + t2)*(B*t1 - B + t2)*(B*t1 + t2 - 1)
Psep = 0
```

### raw_q_L0

```text
P2 = 16*A*(A - 1)*(eps - 1)*(3*eps - 1)*(B*t1 - B + t2)
PD = 0
Pz = 0
Psep = 0
```

### R

```text
P2 = 3*A*eps*(A - 1)**2*(3*eps - 1)*(6*A*B**3*eps + 24*A*B**3*t1 - 18*A*B**3 - 9*A*B**2*eps - 24*A*B**2*t1 + 24*A*B**2*t2 + 15*A*B**2 + 8*A*B*eps**2*t1 - 8*A*B*eps**2 - 8*A*B*eps*t1 + 21*A*B*eps - 24*A*B*t1 - 24*A*B*t2 - 3*A*B + 8*A*eps**2*t2 - 8*A*eps*t2 - 5*A*eps - 24*A*t2 + 15*A - 6*B**3*eps - 24*B**3*t1 + 18*B**3 + 4*B**2*eps*t1 + 8*B**2*eps*t2 + 5*B**2*eps + 31*B**2*t1 - 10*B**2*t2 - 22*B**2 - 4*B*eps*t1 - 8*B*eps*t2 + B*eps - 7*B*t1 + 10*B*t2 + 4*B)/2
PD = -3*A**3*eps*(A - 1)**3*(3*eps - 2)*(3*eps - 1)*(4*B**3*eps*t1 + 8*B**3*t1 - 6*B**3 - 4*B**2*eps*t1 + 4*B**2*eps*t2 - 2*B**2*eps - 8*B**2*t1 + 8*B**2*t2 + 5*B**2 - 4*B*eps*t1 - 4*B*eps*t2 + 6*B*eps - 8*B*t1 - 8*B*t2 - B - 4*eps*t2 - 8*t2 + 5)/2
Pz = 3*A**3*eps**2*(A - 1)**3*(3*eps - 2)*(3*eps - 1)*(2*B**3*t1 - 2*B**2*t1 + 2*B**2*t2 - B**2 + 2*B*t1 - 2*B*t2 - B + 2*t2)
Psep = -3*A**2*B*eps*(A - 1)**3*(B - 1)*(3*eps - 2)*(3*eps - 1)*(4*B**2*eps*t1 + 8*B**2*t1 - 6*B**2 - 4*B*eps*t1 - 7*B*t1 + 2*B*t2 + 5*B - 3*t1 - t2 - 1)/2
```

## 附录B：六个完整解析系数——最终自变量z

2026-10-04按用户要求，所有最终Hlog改用运动学变量z。求解中的r=1/(1+z)为辅助坐标；Hlog第一个参数是自变量，方括号内才是字母。此处统一用`Hlog(z,[...])`，字母为0、−1，所有有理系数也写成z。

基本关系为
$$G(0;r)=-G(-1;z),\qquad G(1;r)=G(0;z)-G(-1;z).$$
高权重经精确有理换元及fibration转换，六个逆换元和六个crossing差均为零，独立z-GPL级数18点误差小于3e−44。[换元证据](../../../../../workflow/checks/gamma1-03-z/README.md)与[z-r001](../../results/3e-Gamma1-03/z-r001.json)保存连接；原r001字节保持。

`Hlog(z,[a1,...,an])`为标准G函数，核dt/(t−a1)；空字1，纯零字log(z)^n/n!，G(−1;z)=log(1+z)。`zeta[k]`为ζ(k)，`^`为幂。以下是正实z区域的完整系数。

### i_g,-2

```maple
-3/64*Hlog(z,[-1, -1, 0, 0])
+(1/128*zeta[2]+3/128*zeta[3])*Hlog(z,[0])
+3/128*Hlog(z,[-1, 0, 0, 0])
-9/64*zeta[2]*Hlog(z,[-1, -1])
+9/128*zeta[2]*Hlog(z,[-1, 0])
+9/128*zeta[2]*Hlog(z,[0, -1])
+(5/96-7/256*zeta[2])*Hlog(z,[0, 0])
+1/384*Hlog(z,[0, 0, 0])
-1/256*Hlog(z,[0, 0, 0, 0])
+3/128*Hlog(z,[0, -1, 0, 0])
-1/192*Hlog(z,[-1, 0, 0])
-1/192*zeta[3]
-159/2560*zeta[2]^2
+43/288*zeta[2]
+(-1/64*zeta[2]-3/64*zeta[3])*Hlog(z,[-1])
```

### i_g,-1

```maple
-115/256*zeta[2]*zeta[3]
-5/48*zeta[3]
+67/256*zeta[5]
-53/1280*zeta[2]^2
+773/864*zeta[2]
-1/144
-3/512*Hlog(z,[0, 0, 0, 0, 0])
-3/32*Hlog(z,[0, -1, -1, 0, 0])
+1/32*Hlog(z,[0, -1, 0, 0, 0])
-3/32*Hlog(z,[-1, 0, -1, 0, 0])
-3/32*Hlog(z,[-1, -1, 0, 0, 0])
+3/16*Hlog(z,[-1, -1, -1, 0, 0])
+3/256*Hlog(z,[-1, 0, 0, 0, 0])
+1/16*Hlog(z,[0, 0, -1, 0, 0])
+1/64*Hlog(z,[0, -1, 0, 0])
-1/32*Hlog(z,[-1, -1, 0, 0])
+1/64*Hlog(z,[-1, 0, 0, 0])
+5/384*Hlog(z,[0, 0, 0, 0])
+(3/16*zeta[3]-3/32*zeta[2])*Hlog(z,[-1, -1])
+9/16*zeta[2]*Hlog(z,[-1, -1, -1])
-9/32*zeta[2]*Hlog(z,[0, -1, -1])
+(-3/32*zeta[3]+3/64*zeta[2])*Hlog(z,[-1, 0])
-9/32*zeta[2]*Hlog(z,[-1, 0, -1])
+(-3/32*zeta[3]+3/64*zeta[2])*Hlog(z,[0, -1])
+(-1/96*zeta[2]+11/36-33/256*zeta[3])*Hlog(z,[0, 0])
+3/16*zeta[2]*Hlog(z,[0, 0, -1])
-9/32*zeta[2]*Hlog(z,[-1, -1, 0])
+1/8*zeta[2]*Hlog(z,[0, -1, 0])
+(9/256*zeta[2]-23/144)*Hlog(z,[-1, 0, 0])
+(81/2560*zeta[2]^2-23/48*zeta[2]-1/32*zeta[3])*Hlog(z,[-1])
-1/4608*(153*z*zeta[2]-344*z+153*zeta[2]-392)/(1+z)*Hlog(z,[0, 0, 0])
-1/15360*(1059*z*zeta[2]^2-3600*z*zeta[2]-240*z*zeta[3]+1059*zeta[2]^2-3760*zeta[2]-240*zeta[3])/(1+z)*Hlog(z,[0])
```

### i_g,0

```maple
-57/64*zeta[2]*zeta[3]
-619/576*zeta[3]
+257/384*zeta[5]
+1369/1920*zeta[2]^2
-5637/5120*zeta[2]^3
+191/512*zeta[3]^2
-1/48
+20453/5184*zeta[2]
+(-9/256*zeta[2]+37/96)*Hlog(z,[-1, -1, 0, 0])
-9/256*zeta[2]*Hlog(z,[0, 0, 0, -1])
+(-37/192-97/512*zeta[2])*Hlog(z,[-1, 0, 0, 0])
+(9/64*zeta[3]-1647/2560*zeta[2]^2+37/32*zeta[2])*Hlog(z,[-1, -1])
+(27/64*zeta[2]-39/64*zeta[3])*Hlog(z,[-1, -1, -1])
+(-27/128*zeta[2]+39/128*zeta[3])*Hlog(z,[0, -1, -1])
+(-9/128*zeta[3]+3551/5120*zeta[2]^2-37/64*zeta[2])*Hlog(z,[-1, 0])
+(-27/128*zeta[2]+39/128*zeta[3])*Hlog(z,[-1, 0, -1])
+(-9/128*zeta[3]+1919/5120*zeta[2]^2-37/64*zeta[2])*Hlog(z,[0, -1])
+(35/576*zeta[2]-15/64*zeta[3]+2287/1728-793/2560*zeta[2]^2)*Hlog(z,[0, 0])
-15/32*zeta[2]*Hlog(z,[-1, 0, 0, -1])
+117/128*zeta[2]*Hlog(z,[0, -1, -1, -1])
-3/16*zeta[2]*Hlog(z,[0, -1, 0, -1])
-1/8*zeta[2]*Hlog(z,[0, -1, -1, 0])
-117/64*zeta[2]*Hlog(z,[-1, -1, -1, -1])
+(43/128*zeta[2]-1/16*zeta[3])*Hlog(z,[0, 0, -1])
-1/32*zeta[2]*Hlog(z,[-1, 0, -1, 0])
-73/256*zeta[2]*Hlog(z,[0, 0, -1, 0])
+117/128*zeta[2]*Hlog(z,[-1, -1, 0, -1])
+(-253/576+127/512*zeta[2])*Hlog(z,[0, 0, 0, 0])
+(-27/128*zeta[2]+39/128*zeta[3])*Hlog(z,[-1, -1, 0])
+(-145/512*zeta[2]-37/192)*Hlog(z,[0, -1, 0, 0])
-3/16*zeta[2]*Hlog(z,[0, 0, -1, -1])
+(-5/384*zeta[2]+75/256*zeta[3]-1615/1728)*Hlog(z,[-1, 0, 0])
+117/128*zeta[2]*Hlog(z,[-1, -1, -1, 0])
+117/128*zeta[2]*Hlog(z,[-1, 0, -1, -1])
+(-139/1280*zeta[2]^2+283/256*zeta[2]*zeta[3]-1615/576*zeta[2]+37/96*zeta[3]-397/256*zeta[5])*Hlog(z,[-1])
+39/128*Hlog(z,[-1, 0, -1, -1, 0, 0])
+39/128*Hlog(z,[-1, -1, -1, 0, 0, 0])
+39/128*Hlog(z,[0, -1, -1, -1, 0, 0])
-231/512*Hlog(z,[-1, 0, 0, 0, 0, 0])
-231/512*Hlog(z,[0, -1, 0, 0, 0, 0])
-1/16*Hlog(z,[0, -1, 0, -1, 0, 0])
+39/128*Hlog(z,[-1, -1, 0, -1, 0, 0])
-1/23040*(3849*z*zeta[2]^2+13455*z*zeta[2]*zeta[3]-32160*z*zeta[2]+4440*z*zeta[3]-17145*z*zeta[5]+3849*zeta[2]^2+13455*zeta[2]*zeta[3]-32440*zeta[2]+4440*zeta[3]-17145*zeta[5])/(1+z)*Hlog(z,[0])
-1/13824*(1296*z*zeta[2]+3699*z*zeta[3]-6376*z+1296*zeta[2]+3699*zeta[3]-6544)/(1+z)*Hlog(z,[0, 0, 0])
-1/16*Hlog(z,[0, 0, -1, -1, 0, 0])
-5/32*Hlog(z,[-1, 0, 0, -1, 0, 0])
-39/64*Hlog(z,[-1, -1, -1, -1, 0, 0])
-3/256*Hlog(z,[0, 0, 0, -1, 0, 0])
+23/128*Hlog(z,[0, -1, -1, 0, 0, 0])
-89/256*Hlog(z,[0, 0, -1, 0, 0, 0])
-25/768*Hlog(z,[0, 0, 0, 0, 0])
-9/128*Hlog(z,[0, -1, -1, 0, 0])
+1/192*Hlog(z,[0, -1, 0, 0, 0])
-9/128*Hlog(z,[-1, 0, -1, 0, 0])
-9/128*Hlog(z,[-1, -1, 0, 0, 0])
+9/64*Hlog(z,[-1, -1, -1, 0, 0])
+7/384*Hlog(z,[-1, 0, 0, 0, 0])
+123/256*Hlog(z,[-1, -1, 0, 0, 0, 0])
+35/128*Hlog(z,[-1, 0, -1, 0, 0, 0])
+43/384*Hlog(z,[0, 0, -1, 0, 0])
+29/256*Hlog(z,[0, 0, 0, 0, 0, 0])
+29/384*zeta[2]*Hlog(z,[0, -1, 0])
```

### i_q,-2

```maple
-1/16*zeta[2]*Hlog(z,[0])
-1/24*Hlog(z,[0, 0])
-1/48*Hlog(z,[0, 0, 0])
+1/24*Hlog(z,[-1, 0, 0])
+1/24*zeta[3]
-1/9*zeta[2]
+1/8*zeta[2]*Hlog(z,[-1])
```

### i_q,-1

```maple
-1/8*Hlog(z,[-1, -1, 0, 0])
-1/48*(12*z*zeta[2]-3*z*zeta[3]+13*zeta[2]-3*zeta[3])/(1+z)*Hlog(z,[0])
+1/16*Hlog(z,[-1, 0, 0, 0])
-3/8*zeta[2]*Hlog(z,[-1, -1])
+3/16*zeta[2]*Hlog(z,[-1, 0])
+3/16*zeta[2]*Hlog(z,[0, -1])
+(-1/24*zeta[2]-53/288)*Hlog(z,[0, 0])
-1/144*(14+11*z)/(1+z)*Hlog(z,[0, 0, 0])
-1/24*Hlog(z,[0, 0, 0, 0])
+1/16*Hlog(z,[0, -1, 0, 0])
+25/144*Hlog(z,[-1, 0, 0])
+1/16*zeta[3]
-1/40*zeta[2]^2
-445/864*zeta[2]
+1/72
+(25/48*zeta[2]-1/8*zeta[3])*Hlog(z,[-1])
```

### i_q,0

```maple
23/16*zeta[2]*zeta[3]
+283/288*zeta[3]
-5/12*zeta[5]
-323/480*zeta[2]^2
+1/16
-10609/5184*zeta[2]
+17/48*Hlog(z,[0, 0, 0, 0, 0])
-3/16*Hlog(z,[0, -1, -1, 0, 0])
-5/48*Hlog(z,[0, -1, 0, 0, 0])
-3/16*Hlog(z,[-1, 0, -1, 0, 0])
-3/16*Hlog(z,[-1, -1, 0, 0, 0])
+3/8*Hlog(z,[-1, -1, -1, 0, 0])
-1/3*Hlog(z,[-1, 0, 0, 0, 0])
-1/12*Hlog(z,[0, 0, -1, 0, 0])
+25/96*Hlog(z,[0, -1, 0, 0])
-25/48*Hlog(z,[-1, -1, 0, 0])
+25/96*Hlog(z,[-1, 0, 0, 0])
+11/36*Hlog(z,[0, 0, 0, 0])
+9/8*zeta[2]*Hlog(z,[-1, -1, -1])
-9/16*zeta[2]*Hlog(z,[0, -1, -1])
-9/16*zeta[2]*Hlog(z,[-1, 0, -1])
-1/4*zeta[2]*Hlog(z,[0, 0, -1])
-9/16*zeta[2]*Hlog(z,[-1, -1, 0])
+1/12*zeta[2]*Hlog(z,[0, -1, 0])
-1/2880*(192*z*zeta[2]^2+2865*z*zeta[2]-750*z*zeta[3]+192*zeta[2]^2+2905*zeta[2]-750*zeta[3])/(1+z)*Hlog(z,[0])
+1/1728*(540*z*zeta[2]-565*z+540*zeta[2]-589)/(1+z)*Hlog(z,[0, 0, 0])
+(3/8*zeta[3]-25/16*zeta[2])*Hlog(z,[-1, -1])
+(-3/16*zeta[3]+25/32*zeta[2])*Hlog(z,[-1, 0])
+(-3/16*zeta[3]+25/32*zeta[2])*Hlog(z,[0, -1])
+(-19/144*zeta[2]+3/8*zeta[3]-1193/1728)*Hlog(z,[0, 0])
+(5/48*zeta[2]+577/864)*Hlog(z,[-1, 0, 0])
+(7/10*zeta[2]^2+577/288*zeta[2]-25/48*zeta[3])*Hlog(z,[-1])
```

## 附录C：是否全部落在原M/N内（2026-10-04）

**不能把03的完整解称为全部M/N解。** 按原M、N定义，在同一z、自然权重≤6、允许ζ常数与乘积及有理z系数的范围内，六个非零系数有五个已完整约化，胶子有限项的权重六部分仍在该空间之外。[独立检查](../../../../../workflow/checks/gamma1-03-mn/README.md) · [等价表示/范围记录](../../results/3e-Gamma1-03/mn-r001.json)。

| 通道 | ε⁻² | ε⁻¹ | ε⁰ |
|---|---|---|---|
| 夸克 | 原M/N，常系数 | 原M/N，需有理z系数 | 原M/N，需有理z系数 |
| 胶子 | 原M/N，常系数 | 原M/N，需有理z系数 | 权重0–5可约，权重6有非零余项 |

这里M/N沿用项目的W、L矩，N(k,j)=M(k,j+1)−log(z)M(k,j)。χ=(1−z)/(1+z)=2r−1；夸克ε⁻¹和有限项分别含χM02/192、χM02/288，胶子ε⁻¹含−χM02/384。因此若只允许有理常数/ζ系数，连这五项中的后三项也不能直接宣称已落入同一范围。

五个成功式在全新HyperInt会话代回完整GPL，差均为零。胶子有限项权重六的33候选在纯长度六词坐标存在左零空间见证：
$$
3c_{001000}+3c_{001100}-c_{010000}-c_{011000}-2c_{100000}-3c_{100100}+c_{101000}.
$$
全部候选对此见证为0，原胶子有限项为3/8；另一新会话逐项核验通过。该见证不依赖ζ乘积独立性，并排除在上述固定权重代数内的有理z系数拟合。它不排除任意高权重消项、变换M/N自变量或任意超越系数的更强表示。

下列五式中的M_k_j表示M(k,j;z)，N_k_j表示N(k,j;z)，zeta_k表示ζ(k)，所有M/N与有理系数的自变量均为z，`**`为幂。

### g:-2

```text
-(480*M_0_1 + 3*M_0_3 + 24*M_1_1 - 54*M_2_1 - 27*N_0_2 + 64*zeta_2)/9216
```

### g:-1

```text
-(216*M_0_1*M_1_1*z + 216*M_0_1*M_1_1 - 3888*M_0_1*z*zeta_3 + 8448*M_0_1*z - 3888*M_0_1*zeta_3 + 8448*M_0_1 - 72*M_0_2*z + 72*M_0_2 + 60*M_0_3*z + 60*M_0_3 + 2052*M_1_1*z*zeta_2 + 2208*M_1_1*z + 2052*M_1_1*zeta_2 + 2208*M_1_1 + 153*M_1_3*z + 153*M_1_3 - 108*M_2_1*z - 108*M_2_1 - 108*M_3_1*z - 108*M_3_1 - 108*N_0_2*z - 108*N_0_2 - 432*N_1_2*z - 432*N_1_2 + 608*z*zeta_2 - 1536*z*zeta_3 + 192*z + 608*zeta_2 - 1536*zeta_3 + 192)/(27648*(z + 1))
```

### q:-2

```text
(6*M_0_1 + 3*M_1_1 + 2*zeta_2)/144
```

### q:-1

```text
(636*M_0_1*z + 636*M_0_1 - 18*M_0_2*z + 18*M_0_2 - 3*M_0_3*z - 3*M_0_3 + 300*M_1_1*z + 300*M_1_1 + 54*M_2_1*z + 54*M_2_1 - 27*N_0_2*z - 27*N_0_2 + 128*z*zeta_2 - 384*z*zeta_3 + 48*z + 128*zeta_2 - 384*zeta_3 + 48)/(3456*(z + 1))
```

### q:0

```text
(9180*M_0_1*M_1_1*z + 9180*M_0_1*M_1_1 - 60480*M_0_1*z*zeta_3 + 71580*M_0_1*z - 60480*M_0_1*zeta_3 + 71580*M_0_1 - 360*M_0_2*z + 360*M_0_2 - 4875*M_0_3*z - 4875*M_0_3 + 52380*M_1_1*z*zeta_2 + 34620*M_1_1*z + 52380*M_1_1*zeta_2 + 34620*M_1_1 + 2610*M_1_3*z + 2610*M_1_3 + 6750*M_2_1*z + 6750*M_2_1 + 810*M_3_1*z + 810*M_3_1 + 7965*N_0_2*z + 7965*N_0_2 - 8640*N_1_2*z - 8640*N_1_2 + 2592*z*zeta_2**2 + 2560*z*zeta_2 + 32640*z*zeta_3 + 6480*z + 2592*zeta_2**2 + 2560*zeta_2 + 32640*zeta_3 + 6480)/(103680*(z + 1))
```
