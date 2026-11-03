# [1,2] / [2,1] 计算与复算

日期：2026-09-17。当前 6e-37、6e-38 两条完成态结果见 [解表](../../../benchmarks/2ne/families/mixed-opposite/solutions.md)，[数学方法](../../../benchmarks/2ne/families/mixed-opposite/method.md)。

本轮先运行 D:/py/python.exe -m preprocessing，47 条输入通过；没有修改 prepared 内容或 Kernel/Lean。原始根按实际 D/tilde D、bounds 与 C 检查，映射到 H0×Qid、H0×Qswap。两个结果均无 ε 极点，有限项 −M01³/96。

## 文件与命令

~~~powershell
D:/py/python.exe -m preprocessing
D:/py/python.exe -m unittest discover -s workflow/checks/mixed-opposite -p test_*.py
D:/py/python.exe workflow/checks/mixed-opposite/check_mixed.py
D:/py/python.exe workflow/checks/mixed-opposite/publish.py
D:/py/python.exe workflow/checks/mixed-opposite/verify_results.py
~~~

- reduce_mixed.py：输入身份、传播子/排序/C 审计、母核表达与精确系数。
- check_mixed.py：原顺序多项式积分、有限块级数、原域 Sobol 与独立降维求积。
- test_mixed.py：7 个单元/负例，含错误分子、跨通道排序与重复 Euler 因子拒绝。
- publish.py：保存文档/代码/复用输出快照，新增不可变结果修订；不会运行 Maple。
- verify_results.py：只读核对 schema、全部公式、依赖/快照、当前结果与已保存的检查。
- [current-checks.json](current-checks.json)、[current-delivery.json](current-delivery.json)：运行指针。
- [delivery-summary.json](delivery-summary.json)：实际范围与输出清单。

check_mixed.py 每次创建新运行目录；publish.py 每次创建新结果修订。只想检查已有成果时运行 verify_results.py，不需重算或发布。

## 本轮实际检查

数学检查运行：runs/20260917T021108676259Z。
交付运行：runs/20260917T021941919630Z。

- 7 个单元/负例通过。
- 6 次原外到内 J 的精确积分：每图 ε=1,2,3、z=2/3。
- 6 次对应的终止双级数精确比较；没有在 Gamma 极点评估 I。
- 4 个非整数参数案例的原六维域 Sobol，每例 8×2^17 样本，种子/均值/容差全部保存。
- 4 个有限项的 Hlog/log 多项式/独立降维求积对照；14 个 Laurent 系数的精确导出。
- 3 组 H0 Meijer G/一维求积、Qid+Qswap=H0²/2 和两类 Q crossing 对照。
- 12 个全 I 与有限项差值的 ε 减半样本，检验 O(ε) 缩放。
- 非零 ε 下两图值不同的实际对照，避免将相同有限项误报为全函数相同。

复用的 H0、M01、Qid00、Qswap00 都是此前实际 HyperInt 的 q=0 输出，具体路径、哈希与副本保存在交付 run。已有 single-opposite/verify_results.py 和 same-channel-two-structure/verify_structure.py 本轮亦已通过；这只是复核旧成果依赖，没有重新运行旧 HyperInt 积分。

所有数值比较都是一致性检查，未产生严格数值包络。full ε 结果为 Meijer G 与显式收敛双级数；级数收敛依据见 method.md，有限项为显式 log/π 多项式。formal=not_run。

发布器初稿曾因 schema 字典括号不匹配而未能启动；修正后才产生上述交付目录和 r001。该故障未运行积分、未改写旧成果。一次中文文档读取遇到终端 GBK 编码错误，改用 Python -X utf8 后读取成功。
