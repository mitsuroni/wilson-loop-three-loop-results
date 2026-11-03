# T4[] 计算与复验

2026-10-08更正：历史组装的`series(evenQ,eps=0,3)`漏掉必要正阶数，影响13/35的−1/0及14/34、15/33的0阶。历史r001和日志保存；当前结果为新r002，读取[5e current](../../../benchmarks/tgv/5e-3g/current.json)。显式预算、独立八图和及全44有效总和见[递归验收](../recursion-5e-remainder/README.md)。以下旧验收不覆盖本次发现的预算错误。

[八图解表](../../../benchmarks/tgv/families/t4-empty/solutions.md) · [原核推导](../../../benchmarks/tgv/families/t4-empty/method.md) · [当前交付指针](current-delivery.json)

八图13/14/15/33/34/35/40/43的完整ε已积分收敛级数与ε⁻⁶…ε⁰均已完成：56位置、32非零。四组完整I关系13=35、14=34、15=33、40=43；最高极点依次4/3/3/2。既有M/N/τ/E足够，没有增加基底生成元，所有非零系数权重为6+k。

## 实际检查

- check_math.py：八图原核、排序域、C、D/F尺度桥；四组完整反射及四项负例；完整P分解保留非零recoil。
- Maple/HyperInt：新算cusp至e⁴、对边面余部至e³、困难面至e²、加权recoil至e²、权重导数及合并T_w至目标阶数。正半轴Euler正规化后扫描无剩余区域；HyperInt发散检查始终开启。
- 56项M/N/τ/E表达精确回代Hlog，四组28项系数相等；11项不同顺序、内部编码及组装复现的一致性检查通过。
- 16组原七维J与独立端点降维对照，最大相对差约0.874%；16组完整加权IBP对照，最大约1.012%。每组8次扰乱、每次2¹⁷个Sobol点。
- 8组原T_w与双Euler正规化后的非零ε对照，最大相对差约1.435%，8×2¹⁶点；六个精确有理点连接独立SymPy算子与实际Maple正规化输出。
- 12个全ε有限块：精确单项求和与独立五维Gauss-Jacobi一致。
- 复用的统一权重6基底、mixed输入证据及原分类复验通过；旧工作区、分类和数学结果未改字节。

上述数值误差是经验估计；未执行小ε完整I余项对照，也未用第二条独立解析路线重算各有限项。不同积分顺序/编码核验不冒充数学路线独立性。没有形式整链、复z/i0验证、全ε快速求和或严格无限尾界。formal=not_run；未改Kernel/Lean、外部montopar，未提交Git。

## 复算入口

```powershell
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t4-empty/verify_results.py
```

只读复验核对源/快照哈希、原根连接、56项精确回代、统计及数值日志，不启动Maple或QMC。

check_math.py、full_series.py、check_numeric.py、check_regularizer_numeric.py分别保存原核、级数块、原域/IBP、正规化数值算法。assemble.py生成唯一I的组装脚本；formulas.py比较既有基底；publish.py仅首次创建r001，后续数学修改必须新建修订。

launch.py为每次Maple运行建立新目录，不覆盖旧输出；后半程新增maple_worker.py保存真实退出码。Maple2025/HyperInt本机路径在脚本和environment.json中显式记录；必要读入的原生.m文件也随交付保存。

## 运行与探索记录

| 指针 | 使用范围 |
|---|---|
| current-scan.json / current-ibp_scan.json | 原核、recoil和导数项的发散扫描 |
| current-recoil.json | 四个代表的加权recoil系数 |
| current-easy.json / current-hard.json | cusp、y面和困难面的全部所需矩 |
| current-derivatives.json | D_x、两种D_y完成；D_w单算路线被后续合并项替代 |
| current-tw_vfirst.json | T_w正规化和13的低阶矩；q1运行异常结束，原启动器未记录退出码，原因未确立 |
| current-tw_14.json / current-tw_15.json | 14/15完整T_w矩 |
| current-tw_words.json | 直接构造内部词列表，补齐13高阶，并复算14/15最后系数；完整结束且退出码0 |
| current-assemble_13_14_15_40.json | 最终四代表I，含全部零项和crossing检查 |
| current-consistency.json | 11项精确一致性比较 |

较慢的r-first、t-first探测已停止，日志和已有有效结果保留。首轮有限项拟合遇到未约化ζ(4)的解析器限制，随后在新组装会话显式使用ζ(4)=2ζ(2)²/5，旧输出没有改写。额外正规化数值检查曾因巨型表达式触发Python编译递归限制，改为独立紧凑Euler表达式加实际Maple有理点对照后通过。

只以所需产物齐备和复验通过判定完成；不从单一日志标记或进程曾成功启动推断完成。完整依赖列表见交付manifest；导航文件独立维护，输入阶段README/catalog/workspace仍记录历史状态。
首次封存的只读复验还捕获了六个算子点未导出的缺口：导出器使用了Maple保留名value，进程退出码仍为0，空列表曾让数值脚本的循环跳过。已改为pointValue，并显式要求恰好六条记录；六点精确比较及八组数值对照已重新执行通过。首次未通过的封存目录保留，当前交付指向修复后的新快照，八个r001数学结果字节未改。