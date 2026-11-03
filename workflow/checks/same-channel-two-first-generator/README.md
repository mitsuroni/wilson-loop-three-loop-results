# [2]六个衍生图的自然全阶生成

日期：2026-09-29。按10/23→12/24→11/25推进，当前六图均有来源明确的全阶生成与任意阶有限词规则。主稿为[[2]族第7节](../../../benchmarks/2ne/families/same-channel-two/three-loop.md)，此前分析版保存为[r001](../../../benchmarks/2ne/families/same-channel-two/three-loop-r001.md)。

完成的解析链：

1. 对10的真实带权二维式展开，得到已解单侧Beta源L和一个混合双端源S；没有沿用无权反对称项积分为零的错误。
2. 对S换元，使z仅进入两条互为倒数的边界；Beta/Cauchy恒等式给出三条上三角源方程。S、E、F在辅助z=−1边界为零，边界由原积分缩放与Gamma倍增公式证明。
3. 从显式矩形3F2源导出一维平方C；其任意阶系数通过shuffle、补字首−1和完整端点极限直接积分。由I10+2I12关系求12/24。
4. 中间插入留下同端源T，其z导数恰为已解HD和初等边界项。T(−1)=0固定原函数；加权纯一维因子等于−F，代回得到11/25。

生成语言是固定Beta/GPL源，尚未压缩成M/N式少指标矩族。没有基底最小性或全阶M/N闭包声明。生成数44/47，其中36/47仍为M/N；[3]20/21/22尚待完成。旧有限阶验收账本56/329保持，formal=not_run。辅助复路径只用于定义解析源和边界，不代表物理复z/i0已验证。

## 实际核验

- `algebra.py` → 42个控制：S的方程和零边界、单侧Beta解、原排序Q01、T到原HD的导数连接、纯Cauchy边界、完整被积函数分解及六分量0/−1矩阵。整数ε控制只比较裸源。
- `seeds.py`、`transform.py`、`assemble.py` → 84个系数/边界控制。先独立生成S至ε⁴及带权Q01，再读旧目标；10/23的14个ε⁻⁶…ε⁰位置完全一致。
- `square.py`、`finish_square.py` → 显式一维平方源ε²…ε⁴的新鲜HyperInt计算和镜像，12/24的14个位置完全一致。原始Maple输入、日志和成功状态存于`runs/`。
- `finish_middle.py` → 20个控制，含中间源的六个伪极点抵消和11/25的14个位置；输入HD来自已经独立验证的单对边生成源，不读取本族目标来拟合。
- 四份数值报告共84个控制：`numeric.json` 39、`full-numeric.json` 15、`exchange-numeric.json` 9、`middle-numeric.json` 21。包含S从零边界沿辅助路径求解、原二维积分、HD驱动导数、镜像和六图完整I；使用三组非整数ε。网格差与ODE比较均为经验检查，不是区间证明。
- `render.mjs`、`verify.py`核对公式渲染、来源和历史快照。旧catalog、prepared、结果及单图derivation均保持原字节。

所有新Maple进程显式加载periodLookups，保持发散检查和遇发散停止。只在完整组合的发散系数抵消后取端点，不使用regHlog删掉发散项。没有进行系数拟合，也没有修改Wilson/Lean。

`finish_middle.py`首次启动缺少Fraction导入，在进入新源组装前报NameError；补齐导入后运行通过，没有接受失败输出。撰写草稿时部分LaTeX反斜杠转义错误在发布前修复，控制字符检查与渲染均单独执行。

## 复验

```powershell
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-first-generator/algebra.py
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-first-generator/seeds.py
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-first-generator/transform.py
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-first-generator/assemble.py
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-first-generator/square.py
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-first-generator/finish_square.py
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-first-generator/finish_middle.py
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-first-generator/numeric.py
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-first-generator/full_numeric.py
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-first-generator/exchange_numeric.py
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-first-generator/middle_numeric.py
node workflow/checks/same-channel-two-first-generator/render.mjs
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-first-generator/verify.py
```

Maple脚本建立新的时间戳目录；再次执行会生成新的检查记录，不能覆盖旧证据来制造通过。`update_docs.py`和`update_navigation.py`是已经执行的一次性发布脚本，不是复验命令。未提交Git。
