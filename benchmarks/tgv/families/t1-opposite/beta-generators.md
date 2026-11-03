# 5e-3g-12/36与11/37：由TGV端点权重导出的Beta生成函数

32/16已经化成三个 $\mathcal V_M$ 的组合。本稿继续求12/36，再由完整关系 $I_{11}=I_{32}-I_{12}$ 得到11/37。下文所有简写图号都指 **5e-3g**，不指纯6e。

结果是：这四图共享一个三参数生成函数 $\mathcal V_{\mathrm B}$。它由原TGV的端点积分直接产生，零参数层回到 $-\partial_a\mathcal V_M$；新增矩有固定的对数积分定义。完整公式与任意阶规则已经得到，四图的系数已算至 $\epsilon^2$。这不等于证明新矩最小，或证明它们能统一全部纯2ne额外函数。

[正文内的逐阶显式系数](#explicit-mnb-coefficients)、[独立系数页](beta-coefficients.md)、[机器数据](beta-generators.json)、[实际检查](../../../../workflow/checks/tgv12-beta/README.md)。原r001不改，formal=not_run。

## 1. 完整的I是什么形式

固定 $z>0$、$\ell=\log z$，沿用

$$
\begin{gathered}
r=\frac{z}{1+z},\qquad W(\beta)=\log[\beta(1-\beta)],\qquad
L(\beta,z)=\ell+\log\frac{1-\beta}{\beta},\\
\mathcal V_M(a,b;z)=\int_0^1\frac{e^{aW}}{\beta-r}
\frac{\sinh(bL)}{b}\,d\beta.
\end{gathered}
\tag{B1}
$$

新出现的端点权重是一个标准不完全Beta函数：

$$
\boxed{\Phi_c(w)=
\frac{\Gamma(1-c)^2}{\Gamma(1-2c)}
(1-w)^{2c}B_{1-w}(1-2c,c).}
\tag{B2}
$$

$B_x(p,q)$ 表示不完全Beta函数，区别于下文的矩 $\mathfrak B$。B2先在 $0<w<1$、$\Re c<1/2$ 附近定义，再解析延拓；不能在中间奇点分别评价Gamma因子和Beta因子。其不含此类表观奇点的Euler表示见B7。

将这个权重放入原有M核，并保留原图要求的镜像配对，就得到唯一的新生成源：

$$
\boxed{\mathcal V_{\mathrm B}(a,b,c;z)=
\int_0^1\frac{e^{aW}}{\beta-r}\frac{\sinh(bL)}{b}
\left[z^c\Phi_c(1-\beta)+z^{-c}\Phi_c(\beta)\right]d\beta.}
\tag{B3}
$$

$a,b,c$ 是原点附近的独立参数；$b=0$ 取极限。$\beta=r$ 时 $L=0$，内部表观极点可去，无需增加主值。交换 $z\leftrightarrow1/z$、$\beta\leftrightarrow1-\beta$ 可见B3为镜像偶函数。

四个完整物理量为

$$
\boxed{\begin{aligned}
I_{12}=I_{36}={}&
\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)\Gamma(1-2\epsilon)}
{256\epsilon^4\cosh(3\epsilon\ell/2)}\\
&\times\Bigg\{2\epsilon\mathcal V_{\mathrm B}
(\epsilon/2,\epsilon/2,\epsilon;z)\\
&\quad-B(1-\epsilon,1-\epsilon)\\
&\quad\times\Big[2\cosh(\epsilon\ell)\mathcal V_M(\epsilon/2,\epsilon/2;z)\\
&\qquad-3\mathcal V_M(3\epsilon/2,3\epsilon/2;z)\\
&\qquad+\mathcal V_M(3\epsilon/2,\epsilon/2;z)\Big]\Bigg\},
\end{aligned}}
\tag{B4}
$$

$$
\boxed{\begin{aligned}
I_{11}=I_{37}={}&
\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)\Gamma(1-2\epsilon)}
{128\epsilon^4\cosh(3\epsilon\ell/2)}\\
&\times\Bigg\{\Gamma(1+\epsilon)\Gamma(1-\epsilon)\\
&\quad\times\Big[2\cosh(\epsilon\ell)\mathcal V_M(\epsilon/2,\epsilon/2;z)\\
&\qquad-3\mathcal V_M(3\epsilon/2,3\epsilon/2;z)\\
&\qquad+\mathcal V_M(3\epsilon/2,\epsilon/2;z)\Big]\\
&\qquad-\epsilon\mathcal V_{\mathrm B}(\epsilon/2,\epsilon/2,\epsilon;z)\Bigg\}.
\end{aligned}}
\tag{B5}
$$

这里保留了全部Gamma前因子，Euler补偿恰好一次。两式相加立即恢复[32/16的完整VM式](32-16-generators.md)。因此T1/T2的[1]组六图已有统一的 $\mathcal V_M+\mathcal V_{\mathrm B}$ 表达；其中32/16根本不需要后者。

<a id="explicit-mnb-coefficients"></a>
### 1.1 已经展开的M/N/B系数：极点至epsilon平方

**B4、B5的逐阶展开已经完成到 $\epsilon^2$。** 以下直接给出展开后的多项式，所有M、N、$\mathfrak B$的自变量均为 $z$，$\zeta_n=\zeta(n)$。M/N的定义见B1、B18，B矩的固定积分定义见B15；此处已经用B17—B19消去了能直接回到M/N的低层矩。

生成式与下列系数是同一个答案的两种写法。每一阶的右边都已没有待展开的Gamma函数、$\mathcal V_M$、$\mathcal V_{\mathrm B}$或参数导数。B矩本身作为固定积分族保留，正如M/N作为固定积分族保留。

写作

$$
\begin{aligned}
I_{11}=I_{37}=\sum_{q=-3}^{2}\epsilon^q[\epsilon^q]I_{11}+O(\epsilon^3),
\\
I_{12}=I_{36}=\sum_{q=-2}^{2}\epsilon^q[\epsilon^q]I_{12}+O(\epsilon^3).
\end{aligned}
$$

更低的系数均为零。下面的 $[\epsilon^q]$ 只标记所列多项式属于哪一阶，并不表示右边还留有提取操作。

#### 11与37

$$
\begin{aligned}
[\epsilon^{-3}]I_{11}={}&- \frac{M_{1,1}}{128}
\end{aligned}
\tag{X11m3}
$$

$$
\begin{aligned}
[\epsilon^{-2}]I_{11}={}&- \frac{3 M_{2,1}}{512} - \frac{M_{0,1}^{2}}{512}\\
&- \frac{13 M_{0,1} \zeta_{2}}{256}
\end{aligned}
\tag{X11m2}
$$

$$
\begin{aligned}
[\epsilon^{-1}]I_{11}={}&- \frac{21 M_{1,3}}{1024} - \frac{7 M_{3,1}}{3072}\\
&+ \frac{5 N_{1,2}}{256} - \frac{111 M_{1,1} \zeta_{2}}{512}\\
&- \frac{19 M_{0,1} M_{1,1}}{512} + \frac{M_{0,1} \zeta_{3}}{32}
\end{aligned}
\tag{X11m1}
$$

$$
\begin{aligned}
[\epsilon^{0}]I_{11}={}&- \frac{119 M_{0,1}^{3}}{46080} - \frac{29 M_{4,1}}{12288}\\
&- \frac{3 N_{2,2}}{2048} - \frac{M_{2,3}}{128}\\
&- \frac{\mathfrak B_{1,1,2}}{512} - \frac{\mathfrak B_{0,1,3}}{768}\\
&- \frac{2573 M_{0,1} \zeta_{2}^{2}}{5120} - \frac{13 M_{0,1} M_{2,1}}{1024}\\
&- \frac{9 M_{2,1} \zeta_{2}}{128} - \frac{7 M_{0,1}^{2} \zeta_{2}}{60}\\
&- \frac{3 M_{1,1} \zeta_{3}}{128}
\end{aligned}
\tag{X110}
$$

$$
\begin{aligned}
[\epsilon^{1}]I_{11}={}&- \frac{323 M_{3,3}}{73728} - \frac{311 M_{1,5}}{245760}\\
&- \frac{69 M_{5,1}}{81920} - \frac{\mathfrak B_{1,1,3}}{1536}\\
&- \frac{\mathfrak B_{2,1,2}}{2048} - \frac{\mathfrak B_{0,1,4}}{3072}\\
&- \frac{N_{1,4}}{4096} - \frac{N_{3,2}}{4096}\\
&- \frac{\mathfrak B_{0,3,2}}{6144} - \frac{37391 M_{1,1} \zeta_{2}^{2}}{20480}\\
&- \frac{929 M_{0,1}^{2} M_{1,1}}{12288} - \frac{283 M_{0,1} M_{1,3}}{6144}\\
&- \frac{207 M_{1,3} \zeta_{2}}{1024} - \frac{91 M_{3,1} \zeta_{2}}{3072}\\
&- \frac{31 M_{0,1} M_{3,1}}{6144} - \frac{9 M_{2,1} \zeta_{3}}{512}\\
&+ \frac{33 M_{0,1}^{2} \zeta_{3}}{512} + \frac{45 M_{0,1} N_{1,2}}{1024}\\
&+ \frac{185 N_{1,2} \zeta_{2}}{1024} - \frac{1531 M_{0,1} M_{1,1} \zeta_{2}}{2048}\\
&+ \frac{35 M_{0,1} \zeta_{2} \zeta_{3}}{256}
\end{aligned}
\tag{X111}
$$

$$
\begin{aligned}
[\epsilon^{2}]I_{11}={}&- \frac{44881 M_{0,1}^{4}}{10321920} - \frac{1033 M_{2,5}}{983040}\\
&- \frac{1015 M_{4,3}}{589824} - \frac{677 M_{6,1}}{2949120}\\
&- \frac{\mathfrak B_{1,1,4}}{6144} - \frac{\mathfrak B_{2,1,3}}{6144}\\
&- \frac{\mathfrak B_{1,3,2}}{12288} - \frac{\mathfrak B_{3,1,2}}{12288}\\
&- \frac{\mathfrak B_{0,1,5}}{15360} - \frac{N_{2,4}}{16384}\\
&- \frac{\mathfrak B_{0,3,3}}{18432} - \frac{N_{4,2}}{32768}\\
&- \frac{960409 M_{0,1} \zeta_{2}^{3}}{286720} - \frac{912987 M_{0,1}^{2} \zeta_{2}^{2}}{573440}\\
&- \frac{579437 M_{0,1}^{3} \zeta_{2}}{2580480} - \frac{41987 M_{2,1} \zeta_{2}^{2}}{81920}\\
&- \frac{2467 M_{4,1} \zeta_{2}}{98304} - \frac{1339 M_{2,3} \zeta_{2}}{16384}\\
&- \frac{1145 M_{0,1}^{2} M_{2,1}}{49152} - \frac{863 M_{0,1} M_{2,3}}{49152}\\
&- \frac{521 M_{0,1} M_{4,1}}{98304} - \frac{111 N_{2,2} \zeta_{2}}{8192}\\
&- \frac{63 M_{1,3} \zeta_{3}}{1024} - \frac{37 \mathfrak B_{1,1,2} \zeta_{2}}{2048}\\
&- \frac{37 \mathfrak B_{0,1,3} \zeta_{2}}{3072} - \frac{33 M_{1,1} \zeta_{5}}{640}\\
&- \frac{27 M_{0,1} N_{2,2}}{8192} - \frac{9 \mathfrak B_{1,1,2} M_{0,1}}{2048}\\
&- \frac{7 M_{3,1} \zeta_{3}}{1024} - \frac{3 \mathfrak B_{0,1,3} M_{0,1}}{1024}\\
&+ \frac{3 M_{0,1} \zeta_{3}^{2}}{32} + \frac{15 N_{1,2} \zeta_{3}}{256}\\
&- \frac{1745 M_{0,1} M_{2,1} \zeta_{2}}{8192} - \frac{333 M_{1,1} \zeta_{2} \zeta_{3}}{512}\\
&- \frac{57 M_{0,1} M_{1,1} \zeta_{3}}{512}
\end{aligned}
\tag{X112}
$$

#### 12与36

$$
\begin{aligned}
[\epsilon^{-2}]I_{12}={}&- \frac{M_{2,1}}{512} - \frac{M_{0,1}^{2}}{512}\\
&+ \frac{M_{1,1}}{64} + \frac{7 M_{0,1} \zeta_{2}}{256}
\end{aligned}
\tag{X12m2}
$$

$$
\begin{aligned}
[\epsilon^{-1}]I_{12}={}&- \frac{5 N_{1,2}}{256} - \frac{M_{3,1}}{512}\\
&+ \frac{M_{1,1}}{32} + \frac{M_{2,1}}{64}\\
&+ \frac{M_{0,1}^{2}}{128} + \frac{11 M_{1,3}}{768}\\
&- \frac{M_{0,1} \zeta_{3}}{32} + \frac{M_{0,1} M_{1,1}}{64}\\
&+ \frac{3 M_{0,1} \zeta_{2}}{64} + \frac{7 M_{1,1} \zeta_{2}}{64}
\end{aligned}
\tag{X12m1}
$$

$$
\begin{aligned}
[\epsilon^{0}]I_{12}={}&- \frac{241 M_{0,1}^{3}}{46080} + \frac{M_{1,1}}{16}\\
&+ \frac{M_{2,1}}{32} + \frac{M_{0,1}^{2}}{64}\\
&+ \frac{\mathfrak B_{1,1,2}}{512} + \frac{\mathfrak B_{0,1,3}}{768}\\
&+ \frac{3 N_{2,2}}{2048} + \frac{3 M_{4,1}}{4096}\\
&+ \frac{13 M_{3,1}}{1536} + \frac{19 M_{1,3}}{1536}\\
&+ \frac{19 M_{2,3}}{6144} - \frac{29 M_{2,1} \zeta_{2}}{1024}\\
&- \frac{3 M_{0,1} M_{2,1}}{512} - \frac{M_{1,1} \zeta_{3}}{64}\\
&+ \frac{3 M_{0,1} \zeta_{2}}{32} + \frac{11 M_{0,1} M_{1,1}}{256}\\
&+ \frac{39 M_{1,1} \zeta_{2}}{256} + \frac{221 M_{0,1}^{2} \zeta_{2}}{7680}\\
&+ \frac{1163 M_{0,1} \zeta_{2}^{2}}{5120}
\end{aligned}
\tag{X120}
$$

$$
\begin{aligned}
[\epsilon^{1}]I_{12}={}&\frac{M_{1,1}}{8} + \frac{M_{2,1}}{16}\\
&+ \frac{M_{0,1}^{2}}{32} + \frac{M_{0,1}^{3}}{64}\\
&+ \frac{\mathfrak B_{1,1,3}}{1536} + \frac{\mathfrak B_{2,1,2}}{2048}\\
&+ \frac{\mathfrak B_{0,1,4}}{3072} + \frac{N_{1,4}}{4096}\\
&+ \frac{N_{3,2}}{4096} + \frac{\mathfrak B_{0,3,2}}{6144}\\
&+ \frac{5 M_{4,1}}{1536} + \frac{13 M_{3,1}}{768}\\
&+ \frac{13 M_{1,5}}{24576} + \frac{19 M_{1,3}}{768}\\
&+ \frac{29 M_{2,3}}{3072} + \frac{37 M_{3,3}}{18432}\\
&+ \frac{43 M_{5,1}}{122880} - \frac{185 N_{1,2} \zeta_{2}}{1024}\\
&- \frac{93 M_{3,1} \zeta_{2}}{4096} - \frac{45 M_{0,1} N_{1,2}}{1024}\\
&- \frac{43 M_{0,1}^{2} \zeta_{3}}{512} - \frac{19 M_{0,1} M_{3,1}}{4096}\\
&- \frac{11 M_{2,1} \zeta_{3}}{512} + \frac{M_{1,1} \zeta_{3}}{64}\\
&+ \frac{3 M_{0,1} \zeta_{2}}{16} + \frac{11 M_{0,1} M_{1,1}}{128}\\
&+ \frac{19 M_{0,1} M_{2,1}}{512} + \frac{37 M_{0,1}^{2} \zeta_{2}}{256}\\
&+ \frac{39 M_{1,1} \zeta_{2}}{128} + \frac{53 M_{0,1}^{2} M_{1,1}}{1536}\\
&+ \frac{69 M_{2,1} \zeta_{2}}{512} + \frac{93 M_{0,1} \zeta_{2}^{2}}{256}\\
&+ \frac{131 M_{0,1} M_{1,3}}{4096} + \frac{1547 M_{1,3} \zeta_{2}}{12288}\\
&+ \frac{2527 M_{1,1} \zeta_{2}^{2}}{2560} - \frac{65 M_{0,1} \zeta_{2} \zeta_{3}}{256}\\
&+ \frac{49 M_{0,1} M_{1,1} \zeta_{2}}{128}
\end{aligned}
\tag{X121}
$$

$$
\begin{aligned}
[\epsilon^{2}]I_{12}={}&- \frac{20701 M_{0,1}^{4}}{2064384} + \frac{M_{1,1}}{4}\\
&+ \frac{M_{2,1}}{8} + \frac{M_{0,1}^{2}}{16}\\
&+ \frac{M_{0,1}^{3}}{32} + \frac{\mathfrak B_{1,1,4}}{6144}\\
&+ \frac{\mathfrak B_{2,1,3}}{6144} + \frac{\mathfrak B_{1,3,2}}{12288}\\
&+ \frac{\mathfrak B_{3,1,2}}{12288} + \frac{\mathfrak B_{0,1,5}}{15360}\\
&+ \frac{N_{2,4}}{16384} + \frac{\mathfrak B_{0,3,3}}{18432}\\
&+ \frac{N_{4,2}}{32768} + \frac{5 M_{4,1}}{768}\\
&+ \frac{13 M_{3,1}}{384} + \frac{19 M_{1,3}}{384}\\
&+ \frac{29 M_{2,3}}{1536} + \frac{121 M_{5,1}}{122880}\\
&+ \frac{163 M_{4,3}}{196608} + \frac{163 M_{2,5}}{327680}\\
&+ \frac{175 M_{3,3}}{36864} + \frac{181 M_{1,5}}{122880}\\
&+ \frac{313 M_{6,1}}{2949120} - \frac{18253 M_{2,1} \zeta_{2}^{2}}{81920}\\
&- \frac{575 M_{0,1}^{2} M_{2,1}}{49152} - \frac{15 N_{1,2} \zeta_{3}}{256}\\
&- \frac{11 M_{3,1} \zeta_{3}}{768} - \frac{3 M_{0,1} \zeta_{3}^{2}}{32}\\
&- \frac{3 M_{1,1} \zeta_{5}}{64} + \frac{M_{1,1} \zeta_{3}}{32}\\
&+ \frac{M_{2,1} \zeta_{3}}{64} + \frac{M_{0,1}^{2} \zeta_{3}}{128}\\
&+ \frac{3 M_{0,1} \zeta_{2}}{8} + \frac{3 \mathfrak B_{0,1,3} M_{0,1}}{1024}\\
&+ \frac{9 \mathfrak B_{1,1,2} M_{0,1}}{2048} + \frac{11 M_{0,1} M_{1,1}}{64}\\
&+ \frac{19 M_{0,1} M_{2,1}}{256} + \frac{27 M_{0,1} N_{2,2}}{8192}\\
&+ \frac{37 M_{0,1}^{2} \zeta_{2}}{128} + \frac{37 \mathfrak B_{1,1,2} \zeta_{2}}{2048}\\
&+ \frac{37 \mathfrak B_{0,1,3} \zeta_{2}}{3072} + \frac{39 M_{1,1} \zeta_{2}}{64}\\
&+ \frac{47 M_{1,3} \zeta_{3}}{1536} + \frac{53 M_{0,1} M_{4,1}}{32768}\\
&+ \frac{69 M_{2,1} \zeta_{2}}{256} + \frac{93 M_{0,1} \zeta_{2}^{2}}{128}\\
&+ \frac{111 N_{2,2} \zeta_{2}}{8192} + \frac{113 M_{0,1} M_{2,3}}{16384}\\
&+ \frac{119 M_{0,1} M_{3,1}}{6144} + \frac{145 M_{3,1} \zeta_{2}}{2048}\\
&+ \frac{167 M_{4,1} \zeta_{2}}{32768} + \frac{173 M_{0,1} M_{1,3}}{6144}\\
&+ \frac{211 M_{1,3} \zeta_{2}}{2048} + \frac{505 M_{0,1}^{2} M_{1,1}}{6144}\\
&+ \frac{1169 M_{2,3} \zeta_{2}}{49152} + \frac{9783 M_{1,1} \zeta_{2}^{2}}{10240}\\
&+ \frac{43979 M_{0,1}^{3} \zeta_{2}}{2580480} + \frac{311141 M_{0,1}^{2} \zeta_{2}^{2}}{573440}\\
&+ \frac{397987 M_{0,1} \zeta_{2}^{3}}{286720} - \frac{811 M_{0,1} M_{2,1} \zeta_{2}}{8192}\\
&+ \frac{M_{0,1} M_{1,1} \zeta_{3}}{256} + \frac{3 M_{0,1} \zeta_{2} \zeta_{3}}{64}\\
&+ \frac{45 M_{1,1} \zeta_{2} \zeta_{3}}{256} + \frac{571 M_{0,1} M_{1,1} \zeta_{2}}{1024}
\end{aligned}
\tag{X122}
$$

有限项中，11/37的新矩部分是

$$
-\frac{\mathfrak B_{0,1,3}}{768}
-\frac{\mathfrak B_{1,1,2}}{512},
$$

12/36取相反号。它们连同各自明确列出的M/N多项式就是完整有限项，不再需要另加 $\tau$ 或 $E$。到 $\epsilon^1,\epsilon^2$ 只增加同一固定B族的更高指标成员。

本节直接展示[原系数页](beta-coefficients.md)的既有精确结果，未改数学值。[机器数据](beta-generators.json)包含所有系数及完整Hlog坐标。其余两图32/16的展开仍见[VM主稿](32-16-generators.md)。**尚未完成的是B矩的进一步简化和最小性研究，它们不阻碍这里的M/N/B展开。**

## 2. 为什么原图会给出这个Beta权重

原核和排序已在[原始推导](method.md)核对。12的同边域是 $w<x<y<1$，额外对边距离为 $D_z(s,w)=(1-s)(1-w)+zsw$。先积分TGV的尺度、同边 $x,y$ 及其Beta参数，原图留下

$$
\begin{aligned}
\mathcal N(\epsilon,w)
={}&B(1-\epsilon,1-\epsilon)
\left[w^{2\epsilon}-1\right]+2\epsilon\Phi_\epsilon(w),\\
\Phi_c(w)
={}&B(1-c,1-c)(1-w)
{}_2F_1(1-2c,1-c;2-2c;1-w).
\end{aligned}
\tag{B6}
$$

这一拆分并不猜测新基底：第一项就是32已经求出的端点幂差，第二项就是同边三角边界留下的Gauss函数。把后者写为Euler积分，参数依赖变得更清楚：

$$
\Phi_c(w)=(1-w)\int_0^1\frac{du}{1-u+uw}
\exp\!\left\{c\left[-\log u-\log(1-u)+2\log(1-u+uw)\right]\right\}.
\tag{B7}
$$

它和M的生成方式一致：一个确定的核乘上指数，对参数求导产生确定的对数幂。不同之处是，这里的指数来自TGV内部的Beta权重。

还可以把B7实际求成B2。令

$$
t=\frac{1-u}{1-u+uw},\qquad
u=\frac{1-t}{1-t+wt}.
$$

则Beta密度比满足 $(1-u+uw)^2/[u(1-u)]=w/[t(1-t)]$，从而

$$
\Phi_c(w)=-w^c\int_0^1
\frac{[t(1-t)]^{-c}}{t-1/(1-w)}\,dt.
\tag{B8}
$$

这是单个Beta密度的Cauchy变换，极点 $1/(1-w)>1$ 不在积分区间。对 $t$ 分部积分给出

$$
\Phi_c'(w)+\frac{2c}{1-w}\Phi_c(w)
=-\frac{\Gamma(1-c)^2}{\Gamma(1-2c)}w^{c-1}.
\tag{B9}
$$

乘以积分因子 $(1-w)^{-2c}$，从 $w$ 积到1：

$$
\frac{\Phi_c(w)}{(1-w)^{2c}}
=\frac{\Gamma(1-c)^2}{\Gamma(1-2c)}
\int_w^1 t^{c-1}(1-t)^{-2c}\,dt.
\tag{B10}
$$

令 $x=1-t$ 就是B2。这里使用的边界是
$\Phi_c(w)/(1-w)^{2c}\to0$，先取 $\Re c<1/2$；仅写 $\Phi_c(1)=0$ 不足以排除齐次项。标准不完全Beta与Gauss函数的对应关系见[DLMF 8.17.7](https://dlmf.nist.gov/8.17.E7)；原图中的参数、边界和归一化由上面的推导固定。

## 3. 从二维K积分，一直回到同一个I12

为便于核对旧式，暂记 $B_\epsilon=B(1-\epsilon,1-\epsilon)$，
$G=e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)\Gamma(1-2\epsilon)$。
原有的 $\mathcal T$ 是32的端点差源，剩下的新块记为 $\mathcal A$：

$$
\begin{aligned}
\mathcal T(\epsilon,z)&=(1+z)\int_0^1ds\,dw\,D_z^{\epsilon-1}(1-w^{2\epsilon}),\\
\mathcal A(\epsilon,z)&=(1+z)\int_0^1ds\,dw\,D_z^{\epsilon-1}\Phi_\epsilon(w),\\
\mathcal K(\epsilon,z)&=-B_\epsilon\mathcal T(\epsilon,z)+2\epsilon\mathcal A(\epsilon,z).
\end{aligned}
\tag{B11}
$$

将B11直接代回完整物理量，左边始终保留目标：

$$
\begin{aligned}
I_{12}=-\frac{G}{256\epsilon^4\cosh(3\epsilon\ell/2)}
\Big\{&-B_\epsilon[z^{\epsilon/2}\mathcal T(\epsilon,z)
+z^{-\epsilon/2}\mathcal T(\epsilon,1/z)]\\
&+2\epsilon[z^{\epsilon/2}\mathcal A(\epsilon,z)
+z^{-\epsilon/2}\mathcal A(\epsilon,1/z)]\Big\}.
\end{aligned}
\tag{B12}
$$

第一行已经由32的求解化成三个VM。第二行先积 $s$：

$$
(1+z)\int_0^1D_z(s,w)^{\epsilon-1}ds
=\frac{(zw)^\epsilon-(1-w)^\epsilon}
{\epsilon(w-1/(1+z))}.
$$

再取 $\beta=1-w$，幂差正好产生M核的 $\sinh$：

$$
\mathcal A(\epsilon,z)
=-z^{\epsilon/2}\int_0^1\frac{e^{\epsilon W/2}}{\beta-r}
\frac{\sinh(\epsilon L/2)}{\epsilon/2}\Phi_\epsilon(1-\beta)\,d\beta.
\tag{B13}
$$

镜像项作 $\beta\mapsto1-\beta$ 后，给出同一个核、权重 $\Phi_\epsilon(\beta)$ 和外因子 $z^{-\epsilon}$。所以B12中的第二个方括号恰好是 $-\mathcal V_{\mathrm B}(\epsilon/2,\epsilon/2,\epsilon;z)$。若将已知三项VM方括号简写为 $\Delta_M(\epsilon,z)$，B12就变成

$$
I_{12}=\frac{G}{256\epsilon^4\cosh(3\epsilon\ell/2)}
\left[2\epsilon\mathcal V_{\mathrm B}(\epsilon/2,\epsilon/2,\epsilon;z)
-B_\epsilon\Delta_M(\epsilon,z)\right].
\tag{B14}
$$

展开 $G$ 和 $\Delta_M$ 即B4，再以 $I_{11}=I_{32}-I_{12}$ 得到B5。36/37通过[原T2完整反射](../t2-opposite/method.md)分别等于12/11；反射保持 $z$、颜色和完整前因子，没有额外符号。

原排序域操作沿用 $1/2<\Re\epsilon<1$ 的共同收敛条带；B6—B10的标准函数恒等式经解析延拓接到该条带及 $\epsilon=0$。分项可能有表观极点，应评价完整组合。

## 4. 它生成什么固定积分族，与M/N有什么关系

给新矩独立于图号的定义。令 $d_+=1-u\beta$、$d_-=1-u+u\beta$，定义

$$
\boxed{\begin{aligned}
\mathfrak B_{k,j,n}(z)={}&\int_0^1\frac{W^kL^j}{\beta-r}\,d\beta
\int_0^1du\,\Bigg[
\frac{\beta}{d_+}\left(\ell-\log[u(1-u)]+2\log d_+\right)^n\\
&\hspace{26mm}+\frac{1-\beta}{d_-}
\left(-\ell-\log[u(1-u)]+2\log d_-\right)^n\Bigg],\\
&k,n\ge0,\qquad j\ge1.
\end{aligned}}
\tag{B15}
$$

**矩本身不含 $1/n!$。** 这一定义只是把B7的指数展开后代入B3，因而

$$
\mathcal V_{\mathrm B}(a,b,c;z)
=\sum_{k,m,n\ge0}
\frac{a^k b^{2m}c^n}{k!(2m+1)!n!}
\mathfrak B_{k,2m+1,n}(z).
\tag{B16}
$$

这就是类比M的固定矩语言。$\mathfrak B_{k,j,n}$ 的权重是 $k+j+n+2$；物理生成式只取奇数 $j$。B3是一维不完全Beta加权积分，B15则将其系数完全写成二维初等对数积分，不能把后者说成已全部约成一维初等矩。

最低层可直接解析约化。因为 $\Phi_0(w)=-\log w$，有

$$
\boxed{\mathcal V_{\mathrm B}(a,b,0;z)=-\partial_a\mathcal V_M(a,b;z),\qquad
\mathfrak B_{k,j,0}=-M_{k+1,j}.}
\tag{B17}
$$

再由 $[c]\Phi_c(w)=-\tfrac12\log^2w+2\operatorname{Li}_2(1-w)$，配对两端并使用二对数反射恒等式，得到

$$
\boxed{\mathfrak B_{k,j,1}=
-\frac34M_{k+2,j}-\frac54M_{k,j+2}
+\frac32N_{k,j+1}
-\left(\frac52M_{0,1}+\frac{11}2\zeta_2\right)M_{k,j}.}
\tag{B18}
$$

这里 $M_{k,j}=\int W^kL^j\,d\beta/(\beta-r)$，
$N_{k,j}=M_{k,j+1}-\ell M_{k,j}$，且 $\ell^2=-2M_{0,1}-6\zeta_2$。B17、B18来自积分恒等式，没有先求目标Hlog再拟合。

因此“定义新生成源”不意味着它的每个系数都独立于M/N。进一步计算得到

$$
\begin{aligned}
\mathfrak B_{0,1,2}={}&\frac{11}2M_{0,1}M_{1,1}-8\zeta_3M_{0,1}
+\frac{65}2\zeta_2M_{1,1}\\
&+\frac{41}{12}M_{1,3}-\frac7{12}M_{3,1}-\frac{13}2N_{1,2}.
\end{aligned}
\tag{B19}
$$

B19是新矩定义并独立积分之后，与旧极点系数比较得到的字典，已在完整Hlog向量上精确回代。它不是新矩的定义，也不承担全阶推导。由此四图的所有极点仍可写在M/N中。

有限项才出现组合 $2\mathfrak B_{0,1,3}+3\mathfrak B_{1,1,2}$。与旧语言的关系为

$$
2\mathfrak B_{0,1,3}+3\mathfrak B_{1,1,2}
\equiv 8(\tau_{D,4}+10E)\pmod{\text{已声明的M/N齐次多项式空间}}.
\tag{B20}
$$

完整、非商空间的字典存于[检查数据](../../../../workflow/checks/tgv12-beta/compact-coefficients.json)。B20解释了这两种旧额外函数为何在本批图中总以固定组合出现，却不分别消去 $\tau_{D,4}$ 或 $E$。

## 5. 怎样生成任意Laurent阶

首先，端点权重本身有显式的一阶递推。记普通Taylor系数 $h_n(w)=[c^n]\Phi_c(w)$ 和 $a_n=[c^n]\Gamma(1-c)^2/\Gamma(1-2c)$。B9给出

$$
\begin{gathered}
h_n'(w)=-\frac{2h_{n-1}(w)}{1-w}
-\frac1w\sum_{j=0}^n a_{n-j}\frac{\log^j w}{j!},\\
h_{-1}=0,\qquad h_n(1)=0,\qquad h_0=-\log w,\\
\sum_{n\ge0}a_nc^n=
\exp\!\left[\sum_{m\ge2}\frac{2-2^m}{m}\zeta_m c^m\right].
\end{gathered}
\tag{B21}
$$

逐次只需积分 $dw/w$、$dw/(1-w)$。每个 $h_n$ 都有确定的边界和权重 $n+1$，没有待求的整图函数留在递推右边。代回B15，可按既定积分次序求出任意固定指标的矩；本轮已完成所需权重不超过8的矩。该结构沿这条积分次序只产生 $z$ 字母 $0,-1$，并不因新函数族而自动扩充字母表。

为给完整I提取系数，定义以下**系数多项式**，它们不是新的积分族：

$$
\begin{aligned}
v_q&=\sum_{k+2m+n=q}
\frac{\mathfrak B_{k,2m+1,n}}
{2^{k+2m}k!(2m+1)!n!},\qquad v_{-1}=0,\\
d_q&=[\epsilon^q]\Delta_M(\epsilon,z),\\
p_q&=[\epsilon^q]\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)\Gamma(1-2\epsilon)}
{\cosh(3\epsilon\ell/2)},\\
\alpha_q&=[\epsilon^q]\Gamma(1+\epsilon)\Gamma(1-\epsilon),\qquad
\beta_q=[\epsilon^q]B(1-\epsilon,1-\epsilon).
\end{aligned}
\tag{B22}
$$

$d_q$ 的M有限和已在[32稿§4](32-16-generators.md)给出，另外三列只需Gamma展开。于是任意 $q\ge-4$ 有

$$
\boxed{\begin{aligned}
[\epsilon^q]I_{12}=[\epsilon^q]I_{36}
&=\frac1{256}\sum_{j=0}^{q+4}p_{q+4-j}
\left(2v_{j-1}-\sum_{r=0}^j\beta_r d_{j-r}\right),\\
[\epsilon^q]I_{11}=[\epsilon^q]I_{37}
&=\frac1{128}\sum_{j=0}^{q+4}p_{q+4-j}
\left(\sum_{r=0}^j\alpha_r d_{j-r}-v_{j-1}\right).
\end{aligned}}
\tag{B23}
$$

$d_0=0$、$d_1=-2M_{1,1}$、$v_0=-M_{1,1}$，所以12/36的表观四阶与三阶极点都相消。举例：

$$
\begin{aligned}
[\epsilon^{-3}]I_{11}&=-\frac{M_{1,1}}{128},\\
[\epsilon^{-2}]I_{12}&=-\frac{M_{0,1}^2+M_{2,1}}{512}
+\frac7{256}\zeta_2M_{0,1}+\frac{M_{1,1}}{64}.
\end{aligned}
\tag{B24}
$$

B23是任意阶的计算规则；[§1.1](#explicit-mnb-coefficients)已将它实际执行到 $\epsilon^2$，逐项列出M/N/B多项式。继续求更高阶仍使用同一规则。

还有一个全阶性质：B5中每个 $\epsilon^q$ 系数都具有权重 $6+q$。原因是 $v_n$ 权重 $n+3$、$d_n$ 权重 $n+2$，而相应Gamma与双曲前因子按阶齐次。B4的完整Beta因子另含 $1/(1-2\epsilon)$，所以12/36会混入低权重项；它们正好与32/16中的低权重项在差式中抵消。

## 6. 已完成什么，仍然缺什么

本批六图的“生成式 → M/N/B逐阶展开”已经完成到 $\epsilon^2$，不存在必须先消去B矩才能给出这些系数的障碍。下述未完成事项针对进一步约化与跨族推广，不是本节已经列出的展开。

本轮从积分直接导出了B3及B4/B5，并给出B15、B21、B23的固定全阶规则。加上32/16，六张[1]图的42个旧Laurent位置（22个非零）全部恢复，六图的12个正阶位置也已齐备。74项新系数/连接检查、12项代数检查、47项数值控制通过；数值误差为经验比较，不是严格区间证明。

为检验与纯2ne的联系，本轮还将实际生成的13个 $n\ge2$ 新矩及其齐次乘积加入原M/N空间：

| 权重 | M/N加本轮Beta矩的秩 | 再加入所测纯2ne及TGV目标的联合秩 |
|---|---:|---:|
| 6 | 13 | 14 |
| 7 | 24 | 28 |
| 8 | 42 | 48 |

在此明确范围内，TGV目标都可表示，但纯6e的36、10、11、12、20、21代表目标均仍有余项。权重六尤其只能覆盖B20中的组合，不能分别覆盖 $\tau_{D,4}$ 和 $E$。这些结论有精确非成员证书；它们不排除尚未计算的其他指标、更多源操作、一般有理 $z$ 系数或变换自变量。

下一步适合研究B9的端点输运怎样与纯2ne的双端源连接，寻找多种图共同的生成操作；不应直接把 $\mathfrak B$ 宣布为全部三圈所需的唯一新族。本轮已完成这六张TGV图的固定生成式与至 $\epsilon^2$ 的检验，进一步压缩和最小性保留为研究问题。[]组、纯2ne的旧计数、原输入及r001均保持，未做新Wilson/Lean整链证明或原七维QMC。
