# [2]族母积分与三圈端点插入分析

日期：2026-09-29。主稿为[[2]族当前分析](../../../benchmarks/2ne/families/same-channel-two/three-loop.md)。本轮按用户要求转向4e-04及其六个三圈衍生积分，未继续[1,2]或36/39专项。

完成内容：七图原根分类和物理I连接；三种端点差的共同初等二维核；中间插入的反对称投影；4e-04无权反对称消项在加权后的非零余项；FA+FC的一维矩形原函数平方及完整I10+2I12关系。矩形原函数直接复用已有3F2闭式并独立核对。

这些是生成研究的结构推进，不是六图紧凑全阶求值完成。当前38/47生成、36/47 M/N、56/329保持；不声称M/N不可约性、最小基底或额外字母。formal=not_run。

- `algebra.py`：42个控制。重新读取prepared 0.3的七图，核对实际边字/配对与完整物理I；三个整数参数对照原四维有序双核和二维表示，检查一维平方组合与反对称投影；给出加权消项失效的精确2/3控制。没有在Gamma极点上计算完整物理I。
- `numeric.py`：17个控制。三个非整数参数点使用64/112阶映射Gauss求积，比较A/B/C及一维平方，另以35位求积检查两个显式3F2矩形值。网格误差为经验比较，不是区间证明。
- `attempts/`：初版多项式辅助函数复用了积分变量名，导致ε=2的平方关系检查失败；改用独立Dummy后通过。失败尝试保留，未接受失败结果。
- `render.mjs`、`verify.py`：数学渲染、链接及来源字节核对；旧方法、解表、catalog和结果均不修改。

复验：

```powershell
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-descendants/algebra.py
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-descendants/numeric.py
node workflow/checks/same-channel-two-descendants/render.mjs
D:/py/python.exe -B -X utf8 workflow/checks/same-channel-two-descendants/verify.py
```

`update_navigation.py`是一次性发布脚本，不是复验命令。后续从10/23的带权余项和10+2×12的一维平方矩推进，不先把旧Σ或Q系数重命名成新基底。
