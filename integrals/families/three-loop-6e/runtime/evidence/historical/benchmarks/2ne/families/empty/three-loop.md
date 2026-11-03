# 6e的[]族：求解进度与完整生成式

日期：2026-09-28。本文是三圈纯2ne空集族的专门入口，范围为6e-01、02、03、07、08、09、17、18、19、35，共10图。较低圈结果仍见[统一语言主稿](../../theory/generation/unified-language.md)。

**结论：本族10/10已经达到当前全阶生成目标。** 完整I是Gamma比值与双曲函数的组合，任意阶系数由显式递推产生；不需要新增运动学函数族。极点至有限项的旧结果也已全部求出。

本次按“先了解进度”的要求归集现有成果，核对下列完整I与当前结果的连接；没有重算原六维积分。旧结果、清单与检查保持，formal=not_run。

## 1. 分层完成状态

| 层次 | 进度 | 已有内容 |
|---|---:|---|
| 原有序J到完整I的连接 | 10/10 | 端点有序字、Gamma/Beta求值、当前C/I规范 |
| 标准函数全ε解 | 10/10 | 下文的六组完整公式覆盖全部10图 |
| 固定函数族的任意阶生成 | 10/10 | Gamma对数级数、双曲比和有限卷积 |
| 极点至有限项 | 10/10 | 旧结果均存至ε⁰，余项O(ε) |
| Lean/MCP形式整链 | 未运行 | 解析推导与CAS/数值检查不等于形式证明 |

本族全ε生成的完成数不依赖更高正幂系数是否逐项打印。这里10图按原输入编号计数，不按相同函数或主积分数去重。

## 2. 共同定义与求解来源

令 $z>0$、$\ell=\log z$，

$$
Z=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),\qquad
c_n=\cosh\!\left(\frac{n\epsilon\ell}{2}\right),
$$

$$
R_2=\frac{\Gamma(1+\epsilon)^2}{\Gamma(1+2\epsilon)},\qquad
R_{12}=\frac{\Gamma(1+\epsilon)\Gamma(1+2\epsilon)}
{\Gamma(1+3\epsilon)}.
$$

先从小正ε下的原有序域求值，再在ε=0附近延拓。公共[空集方法](method.md)已完成：将相邻核分解为A/B端点幂，按原排序形成长度不超过3的有序字；这些字全部由Beta/Gamma函数求值，再恢复一次Euler补偿和物理crossing。

因此本族新的工作主要是统一表达与后续复用，原积分求值本身已经完成。

## 3. 全部10图的完整I

### 6e-01、02、03

$$
\boxed{\displaystyle
I_{6e-01}=-\frac{Z^3}{288\epsilon^6},\qquad
I_{6e-02}=I_{6e-03}=-\frac{Z^3}{576\epsilon^6}
=\frac12I_{6e-01}.}
\tag{E1}
$$

这三图全ε不依赖z。

### 6e-07、17

$$
\boxed{\displaystyle
I_{6e-07}=I_{6e-17}
=-\frac{Z^3c_1}{32\epsilon^6c_3}(R_2-R_{12}).}
\tag{E2}
$$

### 6e-08、18

$$
\boxed{\displaystyle
I_{6e-08}=I_{6e-18}
=-\frac{Z^3c_1}{32\epsilon^6c_3}(1-2R_2+R_{12}).}
\tag{E3}
$$

### 6e-09、19

$$
\boxed{\displaystyle
I_{6e-09}=I_{6e-19}
=-\frac{Z^3c_1}{64\epsilon^6c_3}(1-2R_2+R_{12})
=\frac12I_{6e-08}.}
\tag{E4}
$$

### 6e-35

$$
\boxed{\displaystyle
I_{6e-35}
=-\frac{Z^3c_1}{16\epsilon^6c_3}(1-R_2)^2.}
\tag{E5}
$$

以上相等及比例均为完整物理I的恒等式，包含归一化；每张原图的身份和来源仍分别保留。

## 4. 任意阶为何已经可生成

三个Gamma因子的对数为

$$
\begin{aligned}
\log Z&=\sum_{k\ge2}\frac{\zeta_k}{k}\epsilon^k,\\
\log R_2&=\sum_{k\ge2}
\frac{(-1)^k(2-2^k)\zeta_k}{k}\epsilon^k,\\
\log R_{12}&=\sum_{k\ge2}
\frac{(-1)^k(1+2^k-3^k)\zeta_k}{k}\epsilon^k .
\end{aligned}
\tag{E6}
$$

对任意 $F(\epsilon)=\exp(\sum_{k\ge1}a_k\epsilon^k)=\sum_{r\ge0}f_r\epsilon^r$，

$$
f_0=1,\qquad f_r=\frac1r\sum_{k=1}^{r}k\,a_kf_{r-k}.
\tag{E7}
$$

双曲因子的系数是显式有理数乘ℓ的偶次幂。用[统一语言的前因子递推](../../theory/generation/unified-language.md)，乘回显式ε⁻⁶并取有限卷积即可生成任意Laurent阶。最后使用

$$
\ell^2=-2M_{0,1}-6\zeta_2
$$

得到ζ常数与M01的多项式。本族甚至不需要一般M/N高指标矩。

显式ε⁻⁶不是所有图的实际极点阶。例如

$$
R_2-R_{12}=\zeta_2\epsilon^2+O(\epsilon^3),\quad
1-2R_2+R_{12}=2\zeta_3\epsilon^3+O(\epsilon^4),\quad
(1-R_2)^2=\zeta_2^2\epsilon^4+O(\epsilon^5).
$$

这分别解释07/17、08/09/18/19、35的实际首项。

## 5. 逐图现有结果

所有图已列至ε⁰；完整极点与有限项见原[解表](solutions.md)，本页不重复长系数表。

| 图 | 实际首项 | 全阶生成 | 当前结果 |
|---|---:|---|---|
| 6e-01 | ε⁻⁶ | 已完成 | [r002](../../results/6e-01/r002.json) |
| 6e-02 | ε⁻⁶ | 已完成 | [r002](../../results/6e-02/r002.json) |
| 6e-03 | ε⁻⁶ | 已完成 | [r002](../../results/6e-03/r002.json) |
| 6e-07 | ε⁻⁴ | 已完成 | [r002](../../results/6e-07/r002.json) |
| 6e-08 | ε⁻³ | 已完成 | [r002](../../results/6e-08/r002.json) |
| 6e-09 | ε⁻³ | 已完成 | [r002](../../results/6e-09/r002.json) |
| 6e-17 | ε⁻⁴ | 已完成 | [r002](../../results/6e-17/r002.json) |
| 6e-18 | ε⁻³ | 已完成 | [r002](../../results/6e-18/r002.json) |
| 6e-19 | ε⁻³ | 已完成 | [r002](../../results/6e-19/r002.json) |
| 6e-35 | ε⁻² | 已完成 | [r002](../../results/6e-35/r002.json) |

## 6. 下一步与证据

本族作为已完成的标量前因子与有序字控制，供[1]族及更复杂族复用；不再把这10图列为待求全阶积分。需要具体正ε阶时按(E6)–(E7)直接生成，不为填表盲目扩算。

原解析与数值依据见[empty-family检查](../../../../workflow/checks/empty-family/README.md)。本次[两族进度与完整I核对](../../../../workflow/checks/6e-empty-single-status/README.md)只复用已有结果、核对25图来源和新文档公式，不产生新的积分求值。

物理复z/i0与形式整链未在本批完成。本轮全项目仍为32/47达到固定函数族全阶生成目标、56/329位置完成新框架有限阶重建；本族复用不重复增加计数。未提交Git。
