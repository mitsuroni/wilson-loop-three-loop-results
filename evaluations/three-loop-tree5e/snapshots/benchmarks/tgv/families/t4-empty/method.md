# T4[]：带端点权重的三边 TGV 核

本批成员为 13/14/15/33/34/35/40/43（均指 5e-3g）。正实 z，e=ε；唯一 C 已含一次 Euler 补偿。计算支撑放在 [检查目录](../../../../workflow/checks/tgv-t4-empty/README.md)，原输入/旧证据保持原字节；formal=not_run。

## 1. 原 J、反射和端点消去

以规范 TGV 三腿参数为 x,y,w，定义 D 框架

$$
S=\alpha_1\alpha_2y(1-x)+\alpha_1\alpha_3[(1-x)(1-w)+zxw]
 +z\alpha_2\alpha_3w(1-y),
$$
$$
P=\tfrac14\{\alpha_1\alpha_2[-(1+z)(1-x)+zy]
+\alpha_1\alpha_3[1+3z-(1+z)x+z(1+z)w]
+z\alpha_2\alpha_3[1-y-(1+z)w]\}.
$$

当前 prepared 的 F 框架满足 $S_F(z)=zS_D(1/z)$、$P_F(z)=z^2P_D(1/z)$；这两个关系已逐八图精确核验。令

$$d\mu_e=d^3\alpha\,\delta(1-\alpha_1-\alpha_2-\alpha_3)(\alpha_1\alpha_2\alpha_3)^{-e},\qquad
K_g(z,e)=\int_{[0,1]^3}dx\,dy\,dw\int d\mu_e\,P S^{-2+2e}W_g.$$

原排序域中的两个 pure 端点可精确积分，四个代表的权重为

| g | $W_g$ |
|---|---|
| 13 | $(1-y)^e(1-w^e)$ |
| 14 | $\bigl(1-(1-y)^e\bigr)(1-w^e)$ |
| 15 | $[1-(1-y)^e]w^e$ |
| 40 | $1-x^e$ |

外部共同因子为 $e^{-2}$；因此 $J_g^F(z,e)=z^{2e}K_g(1/z,e)/e^2$。这里的 $e^{-2}$ 是 ε 的负二次幂，不是指数函数。

原核反射 $x\mapsto1-w,y\mapsto1-y,w\mapsto1-x$、$\alpha_1\leftrightarrow\alpha_3$、$z\mapsto1/z$，连同所有 pure 端点顺序反转，给出

$$J_{35}(z)=z^{3e}J_{13}(1/z),\quad J_{34}(z)=z^{3e}J_{14}(1/z),\quad
J_{33}(z)=z^{3e}J_{15}(1/z),\quad J_{43}(z)=z^{3e}J_{40}(1/z).$$

故四对完整 I 分别相等。不是只根据图形标签或低阶系数推断相等。

## 2. 唯一 I 的连接

设 $\lambda=1-2e$、$G_3=e^{-3\gamma_E e}\Gamma(1-e)\Gamma(1-2e)$、$\ell=\log z$。sym=4，当前输入的定义是

$$I_g=\frac{\lambda G_3}{8}\frac{J_g(z)+z^{3e}J_g(1/z)}{1+z^{3e}}.$$

令 $Q_g=4\lambda K_g$，以及

$$\mathcal E_g(z,e)=\frac{z^{-e/2}Q_g(z,e)+z^{e/2}Q_g(1/z,e)}2.$$

则

$$\boxed{I_g(z,e)=\frac{G_3(e)}{32e^2\cosh(3e\ell/2)}\mathcal E_g(z,e).}$$

必须先作此连接，再解释各个裸核系数。不能再加入另一份 Euler 补偿。

## 3. 加权 IBP 和 recoil

精确多项式恒等式为

$$P=\frac{-z\partial_xS-(1+z)\partial_yS+\partial_wS}{4}+R,\qquad
R=\frac{1+z}{2}\alpha_3[\alpha_1(1-x+zw)-z\alpha_2w].$$

记 $B_{v,b}[W]=\int_{v=b}d\mu_e\,d^2s\,W S^{-\lambda}$，并将 $\mathcal R[W]=\int d\mu_e\,d^3s\,RW S^{-2+2e}$。完整加权分部积分给出

$$4\lambda K[W]=z(B_{x,1}-B_{x,0})+(1+z)(B_{y,1}-B_{y,0})-(B_{w,1}-B_{w,0})
+\int d\mu_e\,d^3s\,[-zW_x-(1+z)W_y+W_w]S^{-\lambda}+4\lambda\mathcal R[W].$$

在初始正 e 收敛域内取面值：例如 $0^e=0$。若先在面上把 $u^e$ 换成 1，会丢失极点。端点权重也改变联合缩放区域的 e 次数。

对 13/14/15，将 w 面与 $W_w$ 项先精确合并：定义

$$T_w[W]=\int d\mu_e\,d^3s\,W(\partial_wS)S^{-2+2e}.$$

于是 $-(B_{w,1}-B_{w,0})+\int W_w S^{-\lambda}=\lambda T_w[W]$。这是原 IBP 的恒等改写，避免单独导数项产生的额外正规化层数；没有删除 w 面或 recoil。40 的原分解已经顺利求积。

## 4. 边界的低维解析表示

以下函数均在 e=0 附近有限：

$$F_A(e)=\int_0^1db\,du\,[b(1-b)]^{-e}u^{Ae}(1-bu)^{-1+2e},$$
$$G_A(e)=\int_0^1db\,du\,[b(1-b)]^{-e}(1-u)^{Ae}(1-bu)^{-1+2e},\qquad F_0=G_0.$$

对于单项权重 $x^{Ce}(1-y)^{Ae}w^{Be}$，存在的 cusp 面是

$$B_{w,0}=\frac{B(1+Ce,2e)}e F_A(e)\quad(B=0),\qquad
B_{x,1}=\frac{z^{-1+2e}}{(2+B)e^2}G_A(e).$$

$B>0$ 时 w=0 面为零。另有

$$B(1+e,2e)=\frac1{2e}\exp\!\left[\sum_{n\ge2}\frac{(-1)^n\zeta_n}{n}(1+2^n-3^n)e^n\right].$$

y=0、1 面的 α 尺度积分分别给 $1/e$；余部核为

$$Q_0=(1-b)(1-x)(1-w)+z[b+(1-b)x]w,\quad
Q_1=(1-x)[1-(1-b)w]+zx(1-b)w.$$

相应 integrand 为 $[b(1-b)]^{-e}Q_j^{-1+2e}$ 乘面上的 W。F/G 计算至 e⁴，y 面余部至 e³，足以组成目标有限项。

x=0 的非退化面采用 $\alpha=(r,s,1)/(1+r+s)$、$y=u/(1+u)$、$w=v/(1+v)$：

$$F_H=r(1+u)+rsu(1+v)+zsv,$$
$$f_H=(rs)^{-e}(1+r+s)^{-1-e}(1+u)^{-1-2e}(1+v)^{-1-2e}F_H^{-1+2e}.$$

其 e⁰,e¹,e² 的加权矩可按 r,u,s,v 顺序直接用 HyperInt 积分；w=1 的困难面由同一反射导出。

## 5. 正半轴正规化及阶数

对全体三 Wilson 参数再取 $x=t/(1+t)$，公共分母为

$$F=rtuvz+rsuv+rtvz+st vz+rsu+svz+ru+r,$$
$$\mathrm{base}=(rs)^{-e}(1+r+s)^{-1-e}[(1+t)(1+u)(1+v)]^{-1-2e}F^{-2+2e}.$$

这里 $st vz$ 表示 $s\,t\,v\,z$。recoil 分子为

$$R_c=\frac{1+z}{2}(1+u)[r(1+v)+z(r-s)v(1+t)],$$

$T_w$ 分子为

$$N_w=(1+v)[r(zt-1)(1+u)+zs(1+t)].$$

对于单项权重 $(1-y)^{Ae}w^{Be}$，recoil 的联合区域 $\{r,1/u\}$ 缩放次数为 $(1+A)e$；$T_w$ 的两区域 $\{r,v\}$、$\{v,1/t\}$ 分别为 $(1+B)e$、$(2+B)e$。使用保留 e 的 Euler IBP 提取显式极点，再展开已经局部可积的积分。余部扫描和 HyperInt 发散检查开启。

若 $E_{\mathcal S}$ 对趋零变量取 $v\partial_v$、趋无穷变量取 $-v\partial_v$，对应测度次数 d，则使用

$$f\longmapsto\frac{(\omega-d)f-E_{\mathcal S}f}{\omega}.$$

该式在初始收敛域由总导数边界为零得到，然后作亚纯延拓；不是给发散的 e=0 积分直接取有限部分。算法及中间产物实际使用哪些路径，以交付清单为准。

## 6. 完整 ε 的已积分级数

设 p=2−2e、L=1+z，有限多项式

$$P(1-S/L)^n=\sum c^{(n)}_{abc,ijk}(z)x^ay^bw^c\alpha_1^i\alpha_2^j\alpha_3^k.$$

α 单项积分为

$$D_{ijk}(e)=\frac{\Gamma(i+1-e)\Gamma(j+1-e)\Gamma(k+1-e)}{\Gamma(i+j+k+3-3e)}.$$

Wilson 单项积分为

$$m_{13}=\frac{B(b+1,1+e)}{a+1}\left(\frac1{c+1}-\frac1{c+1+e}\right),$$
$$m_{14}=\frac1{a+1}\left(\frac1{b+1}-B(b+1,1+e)\right)\left(\frac1{c+1}-\frac1{c+1+e}\right),$$
$$m_{15}=\frac{1/(b+1)-B(b+1,1+e)}{(a+1)(c+1+e)},\qquad
m_{40}=\frac{1/(a+1)-1/(a+1+e)}{(b+1)(c+1)}.$$

因此

$$\boxed{K_g=L^{-p}\sum_{n=0}^{\infty}\frac{(p)_n}{n!}\sum_{abc,ijk}c^{(n)}_{abc,ijk}(z)D_{ijk}(e)m_g(a,b,c;e).}$$

每个有限块先完成全部有符号相加，再对 n 求和。没有残留未求值积分。由第2节公式给出完整 I，并用反射覆盖八图。

正 z 下 $0\le S\le L/3$。先取实 $1/2<e<1$，逐个分子 $\alpha_i\alpha_j$ 项用对应正距离给出绝对可积界；α 指数为 e−1 或 −e，邻边端点幂 p<1，权重有界。二项展开可支配地逐项积分，并在复 e 的紧子条带局部一致；随后亚纯延拓到 e=0。不能在该无穷级数中未经另证直接逐项作 e=0 展开。

完整 ε 的类别是显式收敛级数，尚无有限特殊函数压缩、快速收敛实现或严格尾界。

## 7. 验证范围

八个原核/C、四组反射、端点积分及完整加权 IBP 已检查。新数值检查从原七维 J 与独立五维权重核出发，并另测完整加权 IBP。有限块用五维 Gauss-Jacobi 与解析单项公式比较。它们是计算证据；经验误差不是严格区间。

本批没有 TGV Kernel 根注册、Lean 整链或复 z/i0 分支验证。正 e 的原积分数值吻合也不等于逐 Laurent 系数独立证明；实际系数核验及未执行的检查会在结果/检查说明中分别列明。