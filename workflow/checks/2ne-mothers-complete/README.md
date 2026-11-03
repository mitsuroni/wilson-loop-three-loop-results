# Step 3完成：47图全ε连接验收

日期：2026-09-25。[定义桥正文](../../../benchmarks/2ne/theory/generation/remaining-bridges.md)、[全47图生成与I表](../../../benchmarks/2ne/theory/generation/mother-map.md)。本批完成后16图，同时重检首31图的完整ε表达。formal=not_run，未新增Laurent求积分。

## 检查内容

- 原多项式/bounds重新核对47根，所有母族生成项与首批记录相同。
- 47条完整I及后16条完整J与现行修订精确相等。旧F中的ε差分、Q的端点位置、三核crossing、两个通道尺度及κ都显式连接。
- 后16图的命名母积分在ε=2、3及z=2/3、5/4处，与原指定顺序J比较，共32个精确案例。
- 读取旧full-series.json的四个F定义，核对完整求和式及符号单项式差分；不只核对别名。
- 72个双核有限块保留符号ε比较；48个三核有限块覆盖S₃六个排列；32个单字母槽核对反射和交换。12项单元检查包含漏ε、错crossing、错误配对、不支持的字长和表达式语法、有限阶比例冒充全函数关系等控制。
- 旧329系数与92别名只读回代；新文档公式/链接离线渲染。没有重复积分这些系数。

实际退出码、日志和来源/成果哈希见[validation.json](validation.json)，逐图与定义记录见[audit.json](audit.json)，文档检查见[render.json](render.json)。整数ε只用于J及正则母积分，未在Gamma极点处检查I。

## 复算与版本

```powershell
D:/py/python.exe -B -X utf8 workflow/checks/2ne-mothers-complete/run_checks.py
```

`audit.py --write`仅生成本批派生数据[mother-map-r002.json](../../../benchmarks/2ne/theory/generation/mother-map-r002.json)和当前阅读表。旧mother-map.json、两个前阶段检查目录、现行prepared/结果/catalog均保持原字节。变更前四份理论/导航Markdown原件存于`baseline/`，与首批验证报告哈希一致；历史检查按原快照解释，不把旧报告改成47条完成。

本批复用既有母族检查工具及旧F/R块实现，独立路线采用整体多项式展开与有序单项式矩。首31条已有母族计算保持，只重新核对完整表达连接。Solver工作区、Kernel及Lean未修改。

## 边界

任意ε的依据是正文中的积分定义、合法反射/分组及逐块恒等式。有限样点检验实现，不是任意ε插值证明；已有特殊函数和级数的分析求值按原来源复用，不声称本轮独立重算全部积分。Step 4的关系系统、Step 5/6的系数新生成验收、全ε最小性及一般复z仍未完成。

探索读取时曾假设所有结果都有顶层I_compact，实际[2]用SymPy/str的I；最终按现行输入映射取表达并按声明语言解析，未修改旧结果结构。首次未定位到的历史脚本路径仅是读取失败，不计入验证。
