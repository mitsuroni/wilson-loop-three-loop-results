# 权重 ≤6：47 图统一系数基底

当前运行：[20260917T043810203852Z](runs/20260917T043810203852Z/session-checks.json)。结果采用固定正实 z、prepared 0.3 和既有 I 归一化；旧结果修订、族解表、Kernel/Lean 不变。

数学说明见 [统一基底](../../../benchmarks/2ne/theory/unified-weight-six.md)，可读输出为 [47 图解表](../../../benchmarks/2ne/theory/weight-six/solutions.md) 和 [92 项转换](../../../benchmarks/2ne/theory/weight-six/conversions.md)。

## 实际结果

- 329 个 ε⁻⁶…ε⁰ 系数位置：128 非零、201 零，其中 8 个旧文件未显式存储的结构零均注明来源。
- 92 个登记项的精确转换，包括周期、偶 M/N、τ/η、E、双核 ρ/χ 和 Σ。
- 权重 2/3/4/5/6 的实际目标系数秩与候选空间秩均为 2/2/4/8/14。
- 去除低权重乘积后，各权重新增代数生成元数为 2/2/1/4/5；权重 6 的乘积空间秩为 9。
- 选中 3 个周期、7 个 M、2 个 N、原有 τDq4 和 E。Σ 被转换为这些对象的多项式。
- 6 项单元/负例通过；329+92 次精确回代及来源重读一致；最终验证含 161 项当前文件/快照/输出哈希检查。
- 本轮没有新积分、HyperInt、Lean/MCP 形式整链；formal=not_run。

这些秩只针对声明的 Q-Hlog/形式周期坐标。不是一般周期超越性、Q(z,ε) 主积分最小性或高圈闭包结论。

## 文件与复算

- unify.py：只读当前 47 个不可变结果修订，补结构零；登记候选；按权重先取乘积再取新函数，生成精确转换与新快照。
- verify.py：只读复验输入、329/92 转换、秩、乘积商维数和数学页面哈希。
- publish.py：生成新增展示页与原字节出版快照，不修改旧族结果。
- test_unified.py：错误坐标、权重不匹配、遗漏方向、shuffle 重数/系数扰动及零/单位。
- run_checks.py：顺序执行上述测试、生成、发布、复验，保留真实 stdout/stderr、返回码和耗时。
- inspect_results.py：紧凑诊断输出，不影响结果。
- 当前运行内 basis-and-conversions.json 是完整机器可读表；provenance.json 绑定读取依赖和输出；publication.json 绑定数学页面，source/ 与 publication/ 保存原字节。

只读复验：

~~~powershell
D:/py/python.exe -X utf8 workflow/checks/unified-weight-six/verify.py
~~~

生成新的独立运行、更新展示页和当前指针：

~~~powershell
D:/py/python.exe -X utf8 workflow/checks/unified-weight-six/run_checks.py
~~~

20260917T043008971710Z 是初轮精确消元，已保留；当前运行补齐单元/负例、独立复验和可读发布。早期本地接线时遇到 SymPy 的 Zeta 解析及个别旧结果无 loop_order 字段，已分别显式映射 Zeta 并从绑定 prepared 包读取圈数；这些启动失败不包含在成功检查数内。

## 下一步

先扩充独立版本的权重 7 周期与低核数候选空间，再计算裸双核 q=3、裸三核 q=1。不能把缺少装饰矩导致的不相容误认为三核独有函数。尚未执行的新积分与全 ε 最小性继续在总计划登记。
