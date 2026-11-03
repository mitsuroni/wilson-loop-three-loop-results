# [2] 族计算与检查

7 个当前根均已交付：全 ε 显式收敛级数，以及 ε⁻⁶…ε⁰ 的已求值 Hlog。结果入口：[族目录](../../../benchmarks/2ne/families/same-channel-two/README.md)。所有数学状态为计算侧，formal: not_run。

## 实际运行

- 新鲜 Maple 2025 / HyperInt：[20260916T125747177085Z](runs/20260916T125747177085Z/maple.log)。从原积分降维核计算 V0 和三个类型 U1/U2，启用 period 表、检查发散并遇发散中止；无旧结果拟合或 PSLQ。Maple 日志最后的 CPU time 为 21.890 秒，未据此声明端到端墙钟耗时。
- [独立检查](runs/20260916T125747177085Z/independent-checks.json)：21 个原顺序整数 ε 的 J；14 个原域非整数 Sobol（各8×2¹⁷点）；14 个 Hlog/独立三维对数积分比较；8 个运动学场景、24 个余项样本。
- [全 ε 级数检查](runs/20260916T125747177085Z/series-checks.json)：21 个整数截断恒等式，12 个非整数有限块，4 组逐级部分和（每组16/32/64/128）。级数单调收敛检查不等于高精度误差认证；收敛证明在 method.md。
- 7 项单元/负例通过，包括错误 C、错误核、错误分子、显式多项式系数及端点可去值。
- [M/N/τ 比较](basis-analysis.json)、[新增一维扩充](extended-basis.json)：精确有理 Hlog/ζ 坐标运算；成功式全部反代，不相容式保留左零空间见证。排除范围仅为列出的有限候选空间。
- [当前交付指针](current-delivery.json)、[交付概要](delivery-summary.json)：7 份成果、49 个系数完整；源码/主稿/旧矩复用与所有输出按 SHA256 绑定。
- 前35个结果和 Kernel/Lean 未修改。没有执行新的 Lean/MCP 全链，没有扩大物理复分支范围。

旧 M/τ 函数从已验收的单对边矩文件复用；发布脚本逐个比对旧交付 summary 的哈希，保存具体依赖。原始本地 GF 档案只作阅读背景，不是新 checkout 验证的必要文件。

## 文件与重算

使用本项目已有 D:/py/python.exe、SymPy、mpmath、NumPy、SciPy、jsonschema；Maple/HyperInt 只在重新积分时需要。正实数值 Hlog 复用独立微分方程幂/对数级数求值器。

~~~powershell
D:/py/python.exe -m preprocessing
D:/py/python.exe -m unittest discover -s workflow/checks/same-channel-two -p test_two.py -v
D:/py/python.exe workflow/checks/same-channel-two/verify_results.py
~~~

重新计算流程：

1. run_maple.py：建立新的时间戳 run，执行 compute_coefficients.mpl。
2. coefficient_math.py、extend_basis.py：计算各系数并比较指定基底。
3. check_two.py、full_series.py：在新 run 写入独立检查与级数检查（exclusive create；不得覆盖旧报告）。
4. publish.py：保存源码快照、新结果修订和工作区。已有 derivation.md 必须与已发布前驱的哈希及其冻结快照吻合后才可更新；有用户未归档改动时拒绝覆盖。
5. verify_results.py：检查全部哈希、schema、输入、Markdown-only、紧凑系数与 Hlog 的精确一致性、成功式及不相容见证。

非本机环境需配置 Maple 和 HyperInt 路径，不能通过改动历史 provenance 来伪造相同环境。代码放本目录，不进入逐图求解工作区。

开发中的两次数值检查失败及发布时的键类型/UTF-8问题记录在 [探索记录](exploration-notes.md)。失败的首次发布只生成草稿快照，未生成有效结果修订；其文件保留，不作为完成证据。当前成果使用成功交付指针。

当前 r002 仅修正检查目录的相对导航，并纳入修正后的 JSON 规范化检查器；与 r001 的 input、domain、normalization、reduction、I、J、laurent 逐项完全相同。r001 与各次源码快照保留。
