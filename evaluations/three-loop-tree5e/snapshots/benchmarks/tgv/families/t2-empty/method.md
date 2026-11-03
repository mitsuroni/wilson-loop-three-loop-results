# T2[]：由原参数反射得到完整解析解

对象为5e-3g-05/06/07/08/20/21/31。已有二圈母图3e-3g-02保留原结果。本次计算采用原积分的精确换元连接到已求出的T1[]，继承其完整ε表达和ε⁻⁶…ε⁰展开；不进行拟合。计算层检查见[复算记录](../../../../workflow/checks/tgv-t2-empty/README.md)，formal=not_run。

## 1. 反射与原始积分的连接

令边反射为
$$
\rho=(1\ 2)(3\ 4),\qquad
s^{(2)}_{i,j}=1-s^{(1)}_{\rho(i),\,n_i+1-j}.
$$
上标(2)表示待求T2图，上标(1)表示对应T1图；边内参数按从小到大编号，反射必须倒转整条边上的排序，包括pure端点。不是只倒转TGV端点。

TGV腿按原输入顺序编号，同步取
$$
(\alpha^{(2)}_1,\alpha^{(2)}_2,\alpha^{(2)}_3)
=(\alpha^{(1)}_3,\alpha^{(1)}_2,\alpha^{(1)}_1).
$$
两个反射都保持积分域：每条有序链
$0<s_1<\cdots<s_n<1$映为反向参数的同型链；alpha单纯形映到自身。五个s变量的绝对Jacobian为1；在alpha3=1-alpha1-alpha2的两独立变量表示中，alpha换元的行列式为−1，绝对值为1。原迭代积分在共同绝对收敛域 $z>0,\ 1/2<\Re\varepsilon<1$ 转为此域积分后换元，最后作ε的亚纯延拓。没有静默修改根输入的迭代顺序。

反射保持无量纲Gram矩阵：$v_{12}=v_{34}=z$、$v_{14}=v_{23}=1$、$v_{13}=v_{24}=-1-z$。参考尺度$v_{14}$映到相等的$v_{23}$，因此本次关系是 **z不变**，不需要再乘z幂或令z→1/z。

具体看T2的两同边端点$a<b$（边3）和孤立端点$c$（边4）：
$$
S_2=zc\alpha_3[\alpha_1(1-a)+\alpha_2(1-b)],\qquad
P_2=\frac{z^2c}{4}\alpha_3(\alpha_2-\alpha_1).
$$
取 $t=1-c,\ x=1-b,\ y=1-a$，及上述alpha换名，恰好得到T1核
$$
S_1=z(1-t)\alpha_1(\alpha_2x+\alpha_3y),\qquad
P_1=\frac{z^2(1-t)}4\alpha_1(\alpha_2-\alpha_3).
$$
这里右侧alpha均为T1坐标。**P的最终符号为正向相等**；不能只计腿置换的奇偶性再手工补一个负号。pure的距离多项式、Lorentz因子、alpha测度和所有传播子指数也逐图相同。

| T2图 | T1图 | T2的pure边 | 反射后的pure边 |
|---|---|---|---|
| 05 | 04 | 3–4 | 3–4 |
| 06 | 03 | 3–4 | 3–4 |
| 07 | 02 | 3–4 | 3–4 |
| 08 | 01 | 3–4 | 3–4 |
| 20 | 19 | 2–3 | 1–4 |
| 21 | 18 | 2–3 | 1–4 |
| 31 | 17 | 1–4 | 2–3 |

完整变量字典逐图保存在[反射结果对应的运行记录](../../../../workflow/checks/tgv-t2-empty/current-run.json)和每图r001的reflection字段。全部七对满足
$$J_{\mathrm{T2},g}(z,\varepsilon)=J_{\mathrm{T1},f(g)}(z,\varepsilon).$$

两个输入的实际C、完整颜色坐标$(-1/4,1/2,0)$、sym=4及唯一一次Euler补偿也逐图一致。利用同一个定义
$$
I=\frac{e^{-3\gamma_E\varepsilon}\Gamma(1-\varepsilon)\Gamma(2-2\varepsilon)}8
\frac{J(z,\varepsilon)+z^{3\varepsilon}J(1/z,\varepsilon)}
{1+z^{3\varepsilon}},
$$
得到相同配对的 **完整I相等**。这既不是只在有限阶相等，也不是由数值近似推断。

## 2. 自包含的完整ε表达

以下令$e=\varepsilon,\ \ell=\log z$，并定义
$$
\lambda=1-2e,\quad
G=e^{-3\gamma_Ee}\Gamma(1-e)\Gamma(1-2e),\quad
R=\frac{\cosh(e\ell/2)}{\cosh(3e\ell/2)},
$$
$$
B=B(1-e,1-e),\quad A=\Gamma(1+e)\Gamma(1-e),\quad U=2A-B,
$$
$$
V=\frac{\Gamma(1+e)\Gamma(1+2e)}{\Gamma(1+3e)},\qquad
W=\frac{\Gamma(1+e)^2}{\Gamma(1+2e)}.
$$
这里B(·,·)为Euler Beta函数；G已经包含prepared C中的一次Euler补偿，不再乘第二次。

有限母函数及两个辅助矩为
$$
\mathcal F(e,\nu)=\frac{B(2-e,1-e)}{1+\nu}
{}_3F_2(1,1,2-e;2+\nu,3-2e;1),
$$
$$
h_0=\frac{U}{\lambda e},\qquad
h_1=\frac{B}{2e}-2\mathcal F(e,2e).
$$
h0/h1只是本计算的常数矩，不是2ne的M/τ生成元。
$$
Z(e)=-\sum_{k=1}^{\infty}
\frac{(-e)_k}{k!(k+2e)}\mathcal F(e,e+k),\qquad
Q(e)=\frac{B}{2e^2}(1-2W+V)-2Z(e).
$$
$(a)_k$是升阶Pochhammer符号。Z的每一项已解析积完；实$0<e<1$时收敛，此后作亚纯延拓。沿用[T1的可积控制与正尾界](../t1-empty/method.md)，不把Z冒充有限Gamma闭式。

七个答案为
$$
I_{05}=\frac{\lambda Gh_1}{576e^5},\qquad
I_{06}=\frac{\lambda G(h_0-h_1)}{576e^5},\qquad
I_{07}=2I_{06},\qquad
I_{08}=\frac{GU}{576e^6},
$$
$$
I_{20}=\frac{\lambda GRQ}{64e^4},\qquad
I_{31}=\frac{GRU(1-V)}{128e^6},\qquad
I_{21}=I_{31}-I_{20}.
$$
因此
$$
I_{08}=I_{05}+I_{06},\qquad
I_{07}=2I_{06},\qquad
I_{31}=I_{20}+I_{21}.
$$
这些关系使本批可用四个全ε函数表示；本次不声称这四个在任意函数恒等式或有理ε系数下最小独立。

05–08的I与z无关；20/21/31只通过R依赖z。最高极点依次为6、6、6、6、3、4、4。所有系数仅含π、ζ(3)、ζ(5)及log(z)多项式，见[完整49项解表](solutions.md)。

## 3. 展开、证据复用与边界

T1的有限F积分已算到e⁵，Z到e⁴，Gamma到e⁶；这些精确矩在零点小邻域的可积控制给出各I余项O(e)。本次完整函数等式不含额外ε奇异因子，故所有49项和固定z>0的O(e)余项原样传递。中间e=1/2等可去奇点须取完整组合的极限。

新执行的检查针对反射连接本身：七组S/P/pure、实际排序域、alpha单纯形、Jacobian、C和I归一化；错误alpha置换、缺少顺序反转、P错号、边界变化均拒绝。另直接计算T2原七维超立方输入，在两个(e,z)点与完整解析J交叉检查。旧T1 HyperInt矩及余项证据按原交付复验并引用，本次没有重新运行Maple。

现有Kernel未注册TGV根，也不支持这条反射与解析延拓整链；本次没有新增正式收据、改动Kernel/Lean或外部montopar。结果formal=not_run，物理复z及i0延拓不在本次结论内。原分类README、catalog及输入阶段derivation.md保持历史字节；当前进度由新解表和r001记录。
