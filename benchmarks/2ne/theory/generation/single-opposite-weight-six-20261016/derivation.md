# 36/39 的四个权重六矩：M/N 分离与唯一新增方向

2026-10-16。接续[独立端点生成源](../single-opposite-free-endpoints-20261016/derivation.md)。本页把其中四个尚未约化的矩求成精确 Hlog，并在同一个 M/N 空间内比较，再恢复完整物理有限项。结论是：四个矩共同只增加一个方向，其中一个矩完全属于旧空间。这里的“一个”限于本页的权重六和四个矩，不是全阶或全部 6e 的主积分数。

## 1. 先固定原源和比较空间

取 z>0、实对数，沿用

\[
D_z=(1-t)(1-s)+zts,\qquad q_z=D_z/\sqrt z,\qquad
k_\lambda=-(1+z)z^{-\lambda/2}D_z^{\lambda-1},
\]
\[
E_\lambda(a,b;z)=\int_0^1dt\,ds\,k_\lambda
\frac{1-t^a}{a}\frac{1-s^b}{b},
\]
\[
\mathcal E(\lambda,a,b;z)=\frac12\left[
z^{(a+b)/2}E_\lambda(a,b;z)
+z^{-(a+b)/2}E_\lambda(a,b;1/z)\right].
\]

三个参数在原积分零点邻域独立展开：

\[
\mathcal E(\lambda,a,b;z)=
\sum_{p,m,n\ge0}\frac{\lambda^pa^mb^n}{p!m!n!}\mathcal E_{p,m,n}(z).
\tag{1}
\]

每个矩已经包含 crossing 偶投影，且 m、n 对称，权重为 p+m+n+4。下标是参数导数次数；例如 λ² 的 Taylor 系数是 \(\mathcal E_{2,0,0}/2\)，不是 \(\mathcal E_{2,0,0}\)。

比较空间 \(\mathcal V_6\) 定义为：总权重六、crossing 偶的全部 ζ/M/N 齐次单项式的有理常数线性包。使用原始 \(M_{k,j},N_{k,j}\)，k≥0、j≥1，并包含它们的乘积；不加入任意有理 z 系数或额外原函数。具体定义见[单对边族 C2–C3](../../../families/single-opposite/three-loop.md#definitions)，特别是

\[
N_{k,j}=M_{k,j+1}-\ell M_{k,j},\qquad
\ell=\log z,\qquad \ell^2=-2M_{0,1}-6\zeta_2.
\]

枚举得到 21 个候选列，精确秩 12；只保留 M/ζ 的对应空间秩为 11。允许整体偶的奇函数乘积，不能只枚举每个因子都偶的单项式。下文所有“超出 M/N”都相对于这一定义。

## 2. 四个矩来自母方程，不由物理答案反解

从上一阶段已经确定的原二维幂矩及其方程出发：

\[
F_\lambda(a,b;z)=\int_0^1dt\,ds\,t^as^bD_z^{\lambda-1},
\]
\[
[\theta(\theta-\lambda)+z(\theta+1+a)(\theta+1+b)]F=z^\lambda,
\qquad \theta=z\partial_z.
\]

令

\[
c=a+b+\lambda,\quad d=(a+\lambda)(b+\lambda),\quad
Y=(1+z)^{1+c}F,\quad U=\theta Y,\quad h=\frac z{1+z}.
\]

方程成为三角递推的形式：

\[
\theta Y=U,\qquad
\theta U=\lambda U+hc\,U-hd\,Y+
\exp[\lambda\ell+c\log(1+z)].
\tag{2}
\]

在参数总次数 q 上，右侧 U、Y 只取低阶；源项是已知对数多项式。积分 \(d\log z\) 对应给 Hlog 词添加首字母 0，积分 \(dz/(1+z)\) 对应添加首字母 −1。因此这条递推直接给出任意固定阶的有限 Hlog 向量。

它仍需要边界常数。原 Mellin 表示在 z→0 的两支给出

\[
F=A_0+\frac{z^\lambda}{\lambda}
\left(\ell-\frac1\lambda+d_0\right)
+O(z)+O(z^{1+\lambda}\ell),
\]
\[
A_0=B(1+a,\lambda)B(1+b,\lambda),\qquad
d_0=\gamma_E-\psi(1-\lambda)+\psi(1+a+\lambda)+\psi(1+b+\lambda).
\]

按每个参数系数取切向正则化常数：

\[
\operatorname{reg}_0U=d_0,\qquad
\operatorname{reg}_0Y=A_0-\lambda^{-2}+d_0/\lambda.
\tag{3}
\]

具体展开使用

\[
\log(\lambda^2A_0)=
\sum_{r\ge2}\frac{(-1)^r\zeta_r}{r}
\left[a^r+b^r+2\lambda^r-(a+\lambda)^r-(b+\lambda)^r\right],
\]
\[
d_0=\sum_{r\ge1}\zeta_{r+1}
\left[\lambda^r+(-1)^{r+1}\bigl((a+\lambda)^r+(b+\lambda)^r\bigr)\right].
\tag{4}
\]

式(3)中所有负 λ 次幂精确抵消。例如零阶是 \(Y_0=\ell^2/2+3\zeta_2\)、\(U_0=\ell\)，所以不能把 Y 的边界常数置零。

为提取端点差，先令 \(H=(1+z)F=e^{-c\log(1+z)}Y\)，再取

\[
\frac{H(a,b)-H(a,0)-H(0,b)+H(0,0)}{ab}.
\tag{5}
\]

在 Taylor 向量中，这等于取 a、b 次数至少为一的系数并分别降一次。随后乘 \(-e^{(a+b-\lambda)\ell/2}\)，最后作 z↔1/z 的偶投影和阶乘恢复。这就是本次生成四个矩的完整路径；计算它们需要 Y 到参数总次数四，Gamma 边界到次数六。

生成结果先写入独立文件，再读取旧 M/N 和物理答案进行比较。七个矩的反演均由新 Maple/HyperInt 会话实际执行，并检查反演两次恢复原式。四个权重六矩的完整表达见[Hlog 系数表](hlog-coefficients.md)；机器可读向量见[reduction.json](../../../../../workflow/checks/single-opposite-weight-six-20261016/reduction.json)。

## 3. 三条明确恒等式与一个不可消去的矩

本页暂以 \(\mathcal E_{2,0,0}\) 作为新增方向的代表，不另起固定权重的新名称。精确结果是

\[
\boxed{\begin{aligned}
\mathcal E_{1,1,0}
={}&\frac12\mathcal E_{2,0,0}
-\frac{3}{16}\zeta_2M_{0,1}^2
+\frac12\zeta_2M_{2,1}
-\frac12\zeta_3M_{1,1}\\
&+\frac{M_{0,1}^3}{144}
+\frac{M_{0,1}M_{2,1}}8
-\frac{N_{2,2}}{16}
+\frac{M_{2,3}}{48},
\end{aligned}}
\tag{6}
\]

\[
\boxed{\begin{aligned}
\mathcal E_{0,2,0}
={}&\frac{13}{20}\zeta_2^2M_{0,1}
+\frac{13}{60}\zeta_2M_{0,1}^2
+\frac14\zeta_2M_{2,1}
-\zeta_3M_{1,1}\\
&+\frac{M_{0,1}^3}{180}
-\frac{N_{2,2}}{16}
+\frac{M_{4,1}}{48},
\end{aligned}}
\tag{7}
\]

\[
\boxed{
\mathcal E_{0,1,1}
=\frac12\mathcal E_{2,0,0}
+\frac{69}{80}\zeta_2^2M_{0,1}
+\frac{3}{40}\zeta_2M_{0,1}^2
+\frac{M_{0,1}^3}{120}
-\frac{M_{4,1}}{64}.
}
\tag{8}
\]

因此，在模去 \(\mathcal V_6\) 的意义下：

| 矩 | 相对于 \([\mathcal E_{2,0,0}]\) 的新部分 | 是否属于旧空间 |
|---|---:|---|
| \(\mathcal E_{2,0,0}\) | 1 | 否 |
| \(\mathcal E_{1,1,0}\) | 1/2 | 否 |
| \(\mathcal E_{0,2,0}\) | 0 | 是 |
| \(\mathcal E_{0,1,1}\) | 1/2 | 否 |

21 个旧列加入四个矩后秩从 12 变成 13。上表并不表示 \(\mathcal E_{2,0,0}\) 本身没有旧空间成分；它是商空间的一个方便代表。式(6)–(8)是完整函数恒等式，而非仅比较最高权重部分。

也不能把单端二阶矩和两端一阶矩混同：

\[
[\mathcal E_{0,2,0}-\mathcal E_{0,1,1}]
=-\frac12[\mathcal E_{2,0,0}]\ne0.
\tag{9}
\]

这正是保留独立 a、b 的价值：物理对角线只提供它们的加权和。

### 非成员结论不依赖未知的 zeta 独立性

设 \(c_w(f)\) 是 f 的标准 shuffle 展开中，纯六字母 Hlog \(H_w(z)\) 的系数，不带任何 zeta 因子。构造

\[
\mathcal L(f)=
\frac13c_{(-1,0,-1,0,0,0)}(f)
-c_{(-1,0,0,-1,0,0)}(f)
-\frac13c_{(0,-1,-1,0,0,0)}(f)
+c_{(0,-1,0,-1,0,0)}(f).
\tag{10}
\]

对全部 21 个旧候选逐列有 \(\mathcal L=0\)，而

\[
\mathcal L(\mathcal E_{2,0,0})=-2,\quad
\mathcal L(\mathcal E_{1,1,0})=-1,\quad
\mathcal L(\mathcal E_{0,2,0})=0,\quad
\mathcal L(\mathcal E_{0,1,1})=-1.
\tag{11}
\]

故至少一个新方向确实必要；式(6)–(8)又证明一个已经足够。本结论不以 ζ₃²、ζ₂³ 等周期之间未证明的有理独立性为前提。

## 4. 回到 36/39 的完整有限项

沿物理切片 λ=a=b=ε，记

\[
\mathcal E(\epsilon,\epsilon,\epsilon)
=S_0+\epsilon S_1+\epsilon^2S_2+\cdots.
\]

正确的组合为

\[
S_2=\frac12\mathcal E_{2,0,0}
+2\mathcal E_{1,1,0}+\mathcal E_{0,2,0}+\mathcal E_{0,1,1},
\qquad [S_2]=2[\mathcal E_{2,0,0}].
\tag{12}
\]

完整物理前因子保持为

\[
I_{36}=I_{39}
=-\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3}
{32\epsilon^2\cosh(3\epsilon\ell/2)}
\mathcal E(\epsilon,\epsilon,\epsilon).
\]

其中

\[
S_0=\frac14M_{2,1}-\frac12\zeta_2M_{0,1},
\qquad
\frac{e^{-3\gamma_E\epsilon}\Gamma(1-\epsilon)^3}
{\cosh(3\epsilon\ell/2)}
=1+\epsilon^2\left(\frac32\zeta_2-\frac98\ell^2\right)+O(\epsilon^3).
\]

将式(6)–(8)连同这一前因子一起回代，得到不显含 \(\ell\) 的结果：

\[
\boxed{\begin{aligned}
[\epsilon^0]I_{36}=[\epsilon^0]I_{39}
={}&-\frac{\mathcal E_{2,0,0}}{16}
+\frac{209}{2560}\zeta_2^2M_{0,1}
+\frac{29}{768}\zeta_2M_{0,1}^2
-\frac{53}{512}\zeta_2M_{2,1}\\
&+\frac{\zeta_3M_{1,1}}{16}
-\frac{M_{0,1}^3}{1152}
-\frac{13M_{0,1}M_{2,1}}{512}
+\frac{3N_{2,2}}{512}\\
&-\frac{M_{2,3}}{768}
-\frac{M_{4,1}}{6144}.
\end{aligned}}
\tag{13}
\]

这条式子完整回代到 36 的 r002 和 39 的 r001，Hlog/周期坐标逐项相等；两图的 ε⁻²、ε⁻¹ 也同时一致。这是已有物理有限项的新生成源分解，并非新发布的 Laurent 阶。

## 5. 对自然生成源的具体启发

新增代表没有脱离原参数积分。直接令 a=b=0：

\[
\mathcal E(\lambda,0,0;z)
=-\frac{1+z}{2}\int_0^1\!\!\int_0^1
\frac{q_z^\lambda}{D_z}
\left[\log t\log s+\log(1-t)\log(1-s)\right]dt\,ds.
\tag{14}
\]

所以

\[
\mathcal E_{p,0,0}
=-\frac{1+z}{2}\int_0^1\!\!\int_0^1
\frac{\log^p q_z}{D_z}
\left[\log t\log s+\log(1-t)\log(1-s)\right]dt\,ds.
\tag{15}
\]

这是事先由原图端点减除得到、可跨权重生成的子族。已知

\[
\mathcal E_{0,0,0}=\frac14M_{2,1}-\frac12\zeta_2M_{0,1},
\]
\[
\begin{aligned}
\mathcal E_{1,0,0}={}&\frac74\zeta_2M_{1,1}
-2\zeta_3M_{0,1}
+\frac14M_{0,1}M_{1,1}
-\frac12N_{1,2}\\
&+\frac1{12}M_{1,3}+\frac18M_{3,1},
\end{aligned}
\tag{16}
\]

而 \(\mathcal E_{2,0,0}\) 首次超出旧空间。其序列具有用户希望的特征：同一个自然源的低阶退化为旧函数，高阶出现必要的新部分。

这次结果还说明，权重六时“同一个端点多取一次对数”和“两个端点共同取对数”的作用不同。前者对应的 \(\mathcal E_{0,2,0}\) 可全部约化，后者的 \(\mathcal E_{0,1,1}\) 携带与传播子二次对数相同的新商方向。这里是已算阶的精确关系；尚未从中推出一个通用于任意 p、m、n 的积分恒等式，也没有证明二维表达不可化为一维。

下一可执行检验应保持母源不变：

1. 权重七求 p+m+n=3 的六个端点对称矩，检验是否都能由旧空间和 \(\mathcal E_{3,0,0}\) 表示。
2. 权重八求 p+m+n=4 的九个矩，允许 \(\mathcal E_{4,0,0}\)、\(\zeta_2\mathcal E_{2,0,0}\)、\(M_{0,1}\mathcal E_{2,0,0}\)，检查是否仍需新方向。
3. 若出现余项，优先保留具体混合指标及其原核；若闭合，再寻找参数积分层面的递推，而不把有限阶相等外推为全阶。

## 6. 实际检查范围

- 原母方程递推、权重、端点对称和边界：160 项精确检查。
- 新 Maple/HyperInt 会话：七次精确反演及七次双反演恢复；发散检查开启。
- 旧裸对角线、低阶 M/N、完整物理归一化和两图旧结果：12 项精确连接。
- 同一 21 列空间的精确约化；独立 SymPy 稠密秩复核 12→13；六项完整回代。
- 三个完整非成员见证和一个无周期见证，分别逐列湮灭；五条商空间关系；14 个错误倍数/结构假设控制。
- z=0.3、0.7 的七个矩共 14 项原二维核对照，并各检查 crossing。160/280 节点精化，最大绝对差约 \(2.82\times10^{-11}\)，不作区间误差证书。

[脚本、来源绑定及检查说明](../../../../../workflow/checks/single-opposite-weight-six-20261016/README.md)。本轮没有 Wilson 数学 MCP 可调用，formal=not_run；不修改原 C/J、已发布结果、历史证据、PDF 或 Wilson/Lean。
