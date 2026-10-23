# Wilson 项目地图

本页负责定位。当前优先级看 [master-plan](master-plan.md)，工作状态看 [任务账本](backlog/)，数学版本看原有 `current.json` / `families/*/catalog.json`，历史验证只对其固定输入有效。数据处理约定见 [数据流](data-flow.md)。

## 新窗口先做什么

1. 读总计划当前部分，向主工作区整合窗口取得一个任务 ID，执行 `./workflow/backlog.ps1 task <ID> --plain`。
2. 读任务卡列出的输入、结果选择器、推导和检查记录。不要根据文件名最大 revision 自行选择数学版本。
3. 运行 `D:/py/python.exe -m workflow.project_index summary`；缺缓存或过期则 `build`。用 `trace` 找固定来源，用 `impact --compare` 判断修订范围。
4. 核对任务负责人、修改范围、实际 Git 状态和共享输出，再开始工作。新写入任务按 [worktree 与交接规则](project-management.md) 分配；既有窗口不搬迁。

## 长期目标与执行窗口

用户授权的三个独立方向见 [长期目标](long-term-goals.md)：6e 母源、当前物质颜色 RG 和 Kernel 可移植恢复。首轮任务为 WILSON-3.1、WILSON-2.3、WILSON-7.2；具体负责人、worktree、输入和下一步由主 [任务账本](backlog/) 及交接记录定位。双 TGV 继续使用其原任务入口。

第二阶段实际交付见[审阅集成记录](checks/parallel-goals-20261017/phase2.md)。6e 权重九向量已冻结并通过三个点的新阶比较；物质颜色入口必须同时读[日志约定勘误](../benchmarks/renormalization/theory/matter-independent-20261017/derivation-r002.md)，不能仅凭 r001 的算术相等宣布外部物理 Γ 验证完成。

## 软件模块

| 位置 | 职责与入口 | 产物/检查 |
|---|---|---|
| 外部 `qgraf-pipeline` | 图生成、拓扑/轨道、端点排序、颜色与倍数；实际安装及代码哈希见 [upstream.json](project-index/upstream.json) | 原始 `integrals.json`、`orbit_definitions.json`；仓内 [三圈源快照](checks/w3-p0/sources/qgraf/3loop/) |
| 外部 `montopar` | 从图到 raw C/J；[模块地图](../benchmarks/theory/montopar/module-map.md)、[I/O 合同](../benchmarks/theory/montopar/io-contract.md) | 各族 `sources/`；[原生符号生产审计](checks/shared-sign-production/README.md) |
| [preprocessing](../preprocessing/README.md) | `python -m preprocessing`；只从已提交的有序导出生成默认 47 根 | [manifest 0.3](../benchmarks/2ne/prepared/manifest.json)，Python tests |
| [wilson_kernel](../wilson_kernel/README.md) | 数学任务/操作与 stdio MCP | 数学 API、任务恢复；不由管理索引改动 |
| [interfaces](../interfaces/wilson-theorem-boundary.md) | 封闭 AST/schema、受控 Lean 生成与边界 | Schema/负例；不是新的解析公式存储 |
| [proof_store](../proof_store/README.md) | 内部证明生成、重放、证据准入 | 收据与 `.wilson/` 本地状态；收据绑定精确陈述/上下文 |
| [lean](../lean/README.md) | 可复用数学及分离的 benchmark 定义 | 固定 `lake-manifest.json`；只按数学改动运行相应检查 |
| `benchmarks/*/workspaces/` | 每积分 `derivation.md`，求解窗口入口 | 结果在相应 `results/`，不可覆盖旧修订 |
| `workflow/checks/` | 框架测试、CAS 脚本、历史运行及证据 | 先读对应 README；很多旧脚本会覆盖原输出 |
| [project_index](project_index/) | 可重建只读管理索引 | [来源登记](project-index/sources.json)，缓存 `.wilson/project-index/`，管理合同与负例测试 |
| [backlog](backlog/) | 本地任务账本 | 根 [配置](../backlog.config.yml)、[工具锁](tools.lock.json)；主工作区统一领取/集成 |

## 科学计算入口

| 范围 | 当前版本/输入入口 | 推导、检查和任务 |
|---|---|---|
| 任意四腿树图联合核 | [2344联合推导](../benchmarks/theory/four-leg-2344-joint-20261023/derivation.md)、[九组状态](../benchmarks/theory/four-leg-remaining-groups-20261023/derivation.md) | WILSON-4.34 In Progress：[本轮证据](checks/four-leg-2344-joint-20261023/README.md)。2344局部接触/边界/体与重复边差商已明确，严格有限F1三步降为五维Li2核；下一合并实际核并检查第17组非线性候选，2234继承。完整有限项、整体亏格零及最终block归属未证，2334与1234缺口保持。 |
| 多边形一圈 | [五／六边形及一般N计划](polygon-one-loop-plan.md)、[输入与N5结果](../benchmarks/polygons/one-loop/README.md) | WILSON-4.15完成P2五边形全ε级数及ε²；下一WILSON-4.17六边形十五项；[P2交接](checks/pentagon-one-loop-20261018/README.md)；新根未准入，整体积分formal待办 |
| 2e/4e/6e 默认 47 根 | [ordered manifest](../benchmarks/2ne/ordered-source-manifest.json)、[prepared manifest](../benchmarks/2ne/prepared/manifest.json)；各 [families](../benchmarks/2ne/families/) catalog 的 `result/latest_result` | [总览](../benchmarks/2ne/README.md)、[全阶研究](../benchmarks/2ne/theory/generation/)，WILSON-3 |
| 两圈 3e-3g | [prepared](../benchmarks/tgv/3e-3g/prepared/manifest.json)、[results](../benchmarks/tgv/3e-3g/results/) | 04以[连续主推导](../benchmarks/tgv/3e-3g/all-orders-problem/derivation.md)为入口；WILSON-4.4已给有限z闭包，下一非齐次Gauss/初值/ε=0连接。[阶段3证据](checks/tgv04-barnes-reduction-20261019/README.md)。仍缺独立 current 选择器，不猜 latest |
| 单TGV无退化通用结构 | [Q结构讨论](../benchmarks/tgv/3e-3g/method-comparison-20261020/q-connection-discussion.md)、[正则组合全阶二维化](../benchmarks/tgv/3e-3g/method-comparison-20261020/regular-sector.md)、[04的M原域来源](../benchmarks/tgv/3e-3g/method-comparison-20261020/mn-origin.md)、[体项分离](../benchmarks/tgv/3e-3g/method-comparison-20261020/body-insertions.md) | WILSON-4.29 In Progress。尖点/y面全阶M已分离；完整E已积为共同二维核，端点幂与Gauss修正各自受控，修正O(eps)。有限项仅来自端点部分，但全阶M闭包与高阶tau对应未证；下一先比较标准Q部分积分与K0/K1及完整E边界响应；未证end=F或corr=Q，保留二维端点约化。完整I与AC4开放。[检查](checks/tgv-method-comparison-20261020/regular-sector-20261021/README.md)。 |
| 5e-3g | [current](../benchmarks/tgv/5e-3g/current.json) | [有效分支和](../benchmarks/recursion-audit/results/5e-effective-sum/r001.json) |
| 4e-4g | [current](../benchmarks/fgv/4e-4g/current.json) | [有效和推导](../benchmarks/fgv/4e-4g/workspaces/effective-sum/derivation.md) WILSON-4.30：[07原域有理分解与源关系](../benchmarks/fgv/4e-4g/workspaces/4e-4g-07/derivation.md)，[检查](checks/fgv07-rational-20261021/README.md)；全阶纯函数三方向闭合仍开放。 |
| 双 TGV | [current](../benchmarks/double-tgv/4e-3g+3g/current.json)、[进度/逐 ID 预算](double-tgv-laurent-plan.md) | WILSON-1；当前 15、16、24 未完全完成；25 的小基底约化另列 |
| 自能/链/Pi2 | [current-r002](../benchmarks/self-energy/current-r002.json) | raw/prepared/result 与 legacy 符号连接逐项登记 |
| TGV 腿 Pi1 | [current](../benchmarks/self-energy/tgv-pi1/current.json) | 12 图；[原生符号检查](checks/shared-sign-production/native-assembly-check.json) |
| 中央 Gamma1 | [current](../benchmarks/self-energy/gamma1/current.json) | [完整报告](../benchmarks/self-energy/gamma1/report/complete/README.md)、[r003 数值修订](checks/gamma1-numeric-r003-20261016/README.md)，WILSON-5 |
| 一/两圈与三圈反项 | [一圈](../benchmarks/renormalization/one-loop/result-r001.json)、[两圈正阶](../benchmarks/renormalization/three-loop/lower-order-inputs-r001.json)、[反项及 beta](../benchmarks/renormalization/three-loop/counterterms-r001.json) | 三圈有限项需要 eps² / eps¹ / eps⁰ |
| 三圈完整颜色/RG | [current](../benchmarks/renormalization/three-loop/current.json)、[旧组装账本](../benchmarks/renormalization/three-loop/assembly-ledger.md) | [31 图输入](checks/w3-ca-nf/inputs.json)、[历史 RG 检查](checks/w3-ca-nf/README.md)，WILSON-2 |
| 纯交换图共同母积分 | [权重七八生成](../benchmarks/theory/integral-geometry-20261017/adjacent-tower/higher-weights/derivation.md)、[负一边界](../benchmarks/theory/integral-geometry-20261017/adjacent-tower/minus-one/derivation.md)、[47图字典](../benchmarks/theory/integral-geometry-20261017/adjacent-tower/graph-reduction/derivation.md) | WILSON-4.28原范围Done，47图有限操作全阶链与边界已接通；Q10反馈、R高权来源和立方关系缺陷已明确，[本轮检查](checks/adjacent-tower-20261023/higher-weights/README.md)。下一固定边界的Green/Euler与Q11/12复合核约化；最小基、全阶zeta周期和整图formal开放。；[32篇公式渲染检查](checks/adjacent-tower-20261023/katex-fix/README.md)已通过。 |
| M/N 函数族的运算 | [微分与积分规律](../benchmarks/theory/integral-geometry-20261017/m-calculus/derivation.md) | WILSON-4.27 Done；任意矩微分、κ积分原函数、2e-02闭合与N源方程，区分裸J逐阶M/N和循环后的两个H；70控制/4局部MCP及[交接](checks/m-calculus-20261019/README.md)。WILSON-4.28已接续到逐次邻边插入的实际算子 |
| 对边＋邻边的单端权重 | [4e-02/05：从J到I](../benchmarks/theory/integral-geometry-20261017/single-endpoint-differential/derivation.md) | WILSON-4.26 Done；正文依次区分原J、辅助F与完整I；自由b方程、Beta源/基点与反射H₂ε保持，参数及矩阵另列附录；[交接与修订](checks/single-endpoint-differential-20261019/README.md) |
| 低圈微分与系数空间 | [2e-02：从J到M与I](../benchmarks/theory/integral-geometry-20261017/differential-structure/2e-02/derivation.md)；[2e/4e/HD综合推导](../benchmarks/theory/integral-geometry-20261017/differential-structure/derivation.md) | WILSON-4.24 Done；独立讲义明确M结构先于循环求和、完整I及递推边界，附至ε²的常系数M表达；[文稿核对](checks/low-loop-differential-20261018/twoe02-writing-20261019/README.md)。原三/九分量系统及tau证据见[交接](checks/low-loop-differential-20261018/README.md)，最小性与hard全阶映射开放 |
| 其他生成研究 | [Gamma 展开](gamma-expansion-plan.md)、[Montopar 算符](montopar-operator-plan.md)、[s3 blocks](../benchmarks/theory/s3-blocks/derivation.md) | WILSON-4，旧窗口是否继续待核实 |

## 可直接查询的身份

```powershell
D:/py/python.exe -m workflow.project_index trace rg:CF-CA-TF-nf:current
D:/py/python.exe -m workflow.project_index impact i:4e-3g+3g-13:r002 --compare i:4e-3g+3g-13:r001
D:/py/python.exe -m workflow.project_index impact i:3e-Gamma1-03:r003 --compare i:3e-Gamma1-03:r002
D:/py/python.exe -m workflow.project_index check --json
```

`file:<仓库相对路径>` 用于普通产物；`i:<图号>:<revision>` 用于单图；`i:<图号>:current` 仅来自已登记选择器。纯胶子 `assembly:CF-CA2:next` 是尚未完成的组装目标，不能当作已有结果。

旧报告和实验不逐条复制为待办。详细历史保留在原路径，以及 [总计划原字节快照](master-plan.history-20261017.md)。管理实施的实际验证见 [本次报告](checks/project-management-20261017/README.md)。
