# 自能的统一现行约定

2026-10-06原生生产已统一：同一作用量、Wilson 路径序与 Fourier/Wick 规则在源→C/J前生效，版本 `wilson-qcd-action-path-fourier/1.0`。新物理 I 直接求和，`assembly_sign=1`；raw 不含 Euler/C4，准备前验证 Gamma 权重、只补一次 Euler。旧 r001 源与结果原字节保留，以下带日期的旧号／待实施段落按其历史版本阅读。17个二点/四端点自能原生prepared与现行r002精确相同，今后新源不再进入旧1.0→2.0号桥。[检查与读取规则](../../workflow/checks/shared-sign-production/README.md)。

2026-10-06：本文2.0的含i插入、正Dyson连接及17个现行r002符号由[一圈至三圈共同推导](../theory/montopar/shared-sign-conventions.md)接回同一作用量与FT，值保持。共同推导同时完成三点腿Π₁／Γ₁的物理色序说明；两点插入约定与三槽Wilson词是分别导出的操作，不按“有无费米圈”给整族加号。

日期：2026-09-25（新增TGV腿源边界，既有数学约定不变）。版本：`wilson-self-energy-1pi-insertion/2.0`。本文是本项目单Pi1、单Pi2及双Pi1串联的共同符号契约。现行输入和结果由[current.json](current.json)唯一选择；不再在物理求和时临时补自能负号。

## 1. 一个算符定义，一条 Dyson 规则

传播子和算符定义为

$$
D_{0,\mu\nu}^{ab}(k)=\frac{-i\delta^{ab}g_{\mu\nu}}{k^2+i0},\qquad
\boxed{\Pi_L^{ab,\mu\nu}(k)=-i\delta^{ab}(k^2g^{\mu\nu}-k^\mu k^\nu)\mathcal P_L(k^2).}
$$

ΠL是**已经含i的1PI插入**；PL是文献的逆传播子标量修正。不要把ΠL再称作δ(D_F⁻¹)：在这一约定下，δ(D_F⁻¹)=−ΠL。所有内部胶子、鬼、夸克符号和颜色已经包含在PL中。

连接规则固定为

$$
\boxed{D_F=D_0+D_0\Pi D_0+D_0\Pi D_0\Pi D_0+\cdots.}
$$

每条链的外部系数均为+1。两个Pi1串联是两个插入，不是一个Pi2；单个Pi2虽然内部有两圈，仍只有一个插入。展开横向标量分母直接检验：

$$
\frac1{1+\mathcal P_1+\mathcal P_2+\cdots}
=1-\mathcal P_1-\mathcal P_2+\mathcal P_1^2+\cdots.
$$

Wilson顶点保持ig0 v·A。路径排序、C4重数、颜色投影和Euler准备均不改变。来源及显式i的推导见[两圈符号审计](../renormalization/two-loop/source-audit.md)；其r001输入属于旧约定，不能与现行I混用。

## 2. 三个分支的 C 与唯一 I

设r=v14、z=v12/v14，先取r<0、z>0，再保持共同Feynman处方。定义

$$
G=e^{-\gamma_E\varepsilon}\Gamma(1-\varepsilon),\quad
U=\frac{\Gamma(1+\varepsilon)^2\Gamma(1-\varepsilon)\Gamma(1-3\varepsilon)}
{\Gamma(1-2\varepsilon)^2\Gamma(1+2\varepsilon)},
$$
$$
\mathcal B_1=\frac{(10-6\varepsilon)C_A-8(1-\varepsilon)T_Fn_f}{3-2\varepsilon},
\qquad \mathcal B_2=4\varepsilon^2P_2.
$$

P2及pAA、pAf、pFf继续使用原[文献算符展开](pi2/operator-expansion.md)，不翻转这些文献标量。单Pi1的总圈数n=2或3，另外两个分支总圈数n=3。令S14=(−2πr+i0)^(nε)、S12=z^(nε)S14，所有分支沿用

$$
C=a_1\mathcal N_nS_{14}R_n,\qquad
\boxed{I=\frac{ma_1}{2}R_n\frac{J(z)+z^{n\varepsilon}J(1/z)}{1+z^{n\varepsilon}}},
$$
$$
\frac m2C[J(z)+z^{n\varepsilon}J(1/z)]
=\mathcal N_n(S_{14}+S_{12})I.
$$

现行因子为

| 分支 | 外部因子 Nn | Rn | a1 |
|---|---|---|---|
| 单Pi1，n=2/3 | (αs/π)^n CA^(n−2) CF B1 | (−1)^(n−1) G^n/[2^(n+2) ε(1−2ε)] | 1 / −1/2 |
| 双Pi1串联 | (αs/π)^3 CF B1² | −G³U/[128ε²(1−2ε)²] | 1 |
| 单Pi2 | (αs/π)^3 CF B2 | +G³U/[128ε²(1−2ε)²] | 1（完整B2保留） |

Pi2保留CFCA²、CFCA TFnf、CF²TFnf三个通道；串联另有CF(TFnf)²。相同的Pi2几何I乘4ε²pc给各颜色系数，不能删Ff。

Pi2和串联共用完整几何J，但其现行R相反：

$$
\boxed{J_{\Pi_2}=J_{\Pi_1+\Pi_1},\qquad I_{\Pi_2}=-I_{\Pi_1+\Pi_1}.}
$$

这不表示Π2等于Π1²。旧数据字段`I_geom`只是Pi2唯一I的schema名称，仍保留以兼容读取；没有另定义第二种物理积分。

## 3. 旧数据迁移与当前入口

旧源额外乘了每个插入一个−1，而算符已采用上述含i定义。因此在保持J、算符和外部颜色定义时，迁移因子取决于**插入次数**：

$$
C_{\rm current}=(-1)^{N_{\rm ins}}C_{\rm r001},\qquad
I_{\rm current}=(-1)^{N_{\rm ins}}I_{\rm r001}.
$$

| 分支 | 图数 | 相对r001 | 当前结果 |
|---|---:|---:|---|
| 单Pi1 | 13 | −1 | results/各图/r002.json |
| 单Pi2 | 2 | −1；共享I及三个颜色系数同时变号 | results/各图/r002.json |
| 双Pi1串联 | 2 | +1，数学值保持 | results/各图/r002.json |

raw和prepared分别发布新的`*-r002.json`；旧r001及旧无修订文件名的prepared包、manifest、normalization.json全部原字节保留，仅用于历史重放。现行清单是[current-r002.json](current-r002.json)，不是三个旧prepared/manifest.json。

本项目的规范化位于导入边界：[preprocessing/self_energy.py](../../preprocessing/self_energy.py)。其操作是版本化输入转换，不是再次计算积分；拒绝重复规范化、错误插入次数、错误源号及Euler次数。当前结果应通过`load_current_result(orbit_id)`读取，读取时核对raw、prepared和result的绑定哈希及约定。

现行求和直接加现行I。两圈例子为

$$
\widehat w_2=C_F(X_s^{2\varepsilon}+X_t^{2\varepsilon})
\left\{C_A\left[\sum I_{4e}+\sum I_{3e-3g}\right]
+\mathcal B_1\left(I_{B_A}+I_{B_O}\right)\right\}.
$$

现行BA本身带负号，故这个正号组装与此前旧I乘−1的已验证物理总和完全相同。新的[两圈r002结果](../renormalization/two-loop/result-r002.json)保留全部物理系数；旧r001结果仍可复现。

## 4. 范围和证据

本次统一的是D:/wilson-agent内的现行输入、结果、理论、解表和读取入口。外部wilsonloop仓库的生产源码及输出保持原样，本项目导入时必须经过上述规范化；历史源快照也不修改。默认47个纯交换根和其预处理命令不变。

迁移通过Dyson代数、全部17个C/J连接、119个主Laurent位置、84个颜色系数位置，以及两圈完整裸/有限/RG结果不变性检查。旧数学结果通过线性变换复用，不自动继承新定义下的形式证明。formal=not_run；物理复z、三圈完整图集和外部生产迁移不在本次验收范围。证据见[迁移检查](../../workflow/checks/self-energy-signs/README.md)。

## 5. 三圈覆盖与树级TGV腿插入（2026-09-25）

上次17图迁移已包含三圈15图（11个4e-Pi1、2个串联、2个Pi2），不能因本次转向三圈而再翻号。详细范围见[三圈状态表](three-loop-status.md)。

新增3e-3g+Pi1的12个原始图直接采用同一Pi定义与正号D0PiD0连接；完整动量核中(−i)³=i，符号来自内部张量/传播子，不在图外补负号。[新分支理论](tgv-pi1/theory.md)保留投影张量作用于树TGV及CFCA²/CFCA TFnf两个通道。

新族已完成上游源定义、公式和图册，尚未有C/J/I，因此不扩充当前已求解current.json，也不套用上表无TGV的Rn。其2.0源标签表示直接遵守约定，禁止再应用legacy迁移。本项目旧r002保持，外部旧分支仍通过导入边界适配；新增分支的生产器和输出已经回写。[检查](../../workflow/checks/self-energy-three-loop/README.md)。
