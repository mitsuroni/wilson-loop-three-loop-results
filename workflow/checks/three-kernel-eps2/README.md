# 三核三图ε⁰…ε²：实际求值与函数空间比较

数学主稿：[三核至ε²](../../../benchmarks/2ne/theory/generation/three-kernel-eps2.md)。全部系数：[Hlog表达](../../../benchmarks/2ne/theory/generation/three-kernel-eps2-hlog.md)、[JSON数据](../../../benchmarks/2ne/theory/generation/three-kernel-eps2.json)。

## 实际成果

- 6e-20、6e-21、6e-22完整I的ε⁰、ε¹、ε²均已求值；21与22相同，余项O(ε³)，全部负幂为零。
- 三图零阶、两个独立一阶逐项匹配历史结果；二阶为新求值，I20/I21分别有152/153个Hlog-周期坐标，字母只有0、−1。
- 两个自然源矩先从矩形边界推导，再用字算法实际积分。两次字首补零取端点替代通用积分器的变量变换；所有线性、系数和端点条件由脚本拒绝式检查。
- 成功会话的Maple CPU约91秒完成六个源分子系数。四个低阶组合与本轮通用积分器输出精确相同。独立参数求积24例、区间方差遗漏4个负控制、Gamma/crossing前因子检查通过。
- 权重8的M/N空间秩36，加入两目标为38；加旧基底乘积后40→42，再加入现有部分单核/双核源后45→47。每个不包含结论有精确对偶见证，成功表达和规范变换保留完整坐标。

## 证据入口

| 文件 | 内容 |
|---|---|
| [fast-run.json](fast-run.json) | 成功CAS会话；原始脚本、每块内积分/外变量表达、完成输出、周期表哈希与日志 |
| [assembly-order2.json](assembly-order2.json) | 独立新会话完成ρ因子卷积与crossing；全部输入输出绑定 |
| [algebra.json](algebra.json) | 三个原根的排列/κ/δ连接、三次幂边界恒等式、物理前因子及源预算 |
| [coefficients.json](coefficients.json) | 原源矩、完整物理系数、历史回归及来源哈希 |
| [comparison.json](comparison.json) | 完整候选列表、矩阵秩、逐目标精确不相容见证 |
| [numeric.json](numeric.json) | 三维/四维参数积分与55位Hlog独立求值的24个比较 |
| [validation.json](validation.json) | 命令、退出码、旧文件保护、源绑定与产物哈希 |

## 可复验命令

```powershell
D:/py/python.exe -B -X utf8 workflow/checks/three-kernel-eps2/run_checks.py
```

该命令复用已保存的实际CAS输出，核对来源、公式、数值与空间见证，不重复耗时积分。需要重新积分时运行run_fast.py，再将输出会话路径交给assemble_maple.py；每次创建新目录，不覆盖历史会话。

最初两次脚本在积分前因Maple受保护名称degree、inner失败，已保存。之后通用积分路线完成低阶及部分高阶；在快速路线完成并通过四个精确交叉检查后，已主动停止其对应的cmaple进程，见该会话的cancellation.json，不能把该退出码当成数学失败。旧三核研究中的900/300秒超时会话保持不变。

初始稠密矩阵比较完成M/N层后换成稀疏有理消元。新消元器用SymPy小矩阵独立核对，并对每项检查完整重构或对偶见证；不是浮点秩判定。

## 比较范围与限制

系数域为有理数；周期使用ζ2、ζ3、ζ5、ζ7和标准ζ(5,3)的形式坐标，已施加既有偶ζ和交换恒等式。结果不是周期超越性定理，也不排除任意有理函数系数、无限和或未枚举源族。

最大的已测空间仍未包含全部单核六阶系数、双核四阶偶源矩或任意端点字。故需要两个新权重8代表，不等于已证明需要两个新的三核全阶函数族。

数值阈值是网格加密的经验控制，不是严格区间包络。当前正实z分支；物理复z延拓和Lean/MCP整链未执行，formal=not_run。旧根、r001/catalog、此前推导和证据保持原字节。有限项新框架覆盖56/329；正ε专项不重复增加该计数。未提交Git。
