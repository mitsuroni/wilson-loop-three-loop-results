# 端点源推广：6e-10/11/12与6e-36

本轮执行4e-04之后的带源推广。正文见[端点源推导](../../../benchmarks/2ne/theory/generation/endpoint-sources.md)，完整非零系数见[新生成解表](../../../benchmarks/2ne/theory/generation/endpoint-source-coefficients.md)。

## 已完成

- 四个原输入的裸组合、crossing参数和物理归一化连接。
- 带源双核的两个二维表示；带上端点项的原函数移位恒等式；保留有序域通量的运动学微分关系；源位置与配对的乘积/反射关系。
- 一个全新Maple会话实际计算六个双核组合系数、三个单核组合系数及τ四阶源矩。未读取旧答案作为CAS输入。
- 由源积分生成τ和E，再精确比较旧定义；四图28位置、9非零全部匹配旧结果。本轮新增6e-11/12的14位置，累计35/329，剩余42图、294位置。
- 以M/N加源生τ、E作有限阶坐标比较；9个最终展示公式全部精确回代为零。未声称全阶闭包或最小性。

## 复验

在仓库根执行：

```powershell
D:/py/python.exe -B -X utf8 workflow/checks/2ne-endpoint-sources/run_checks.py
```

它复用本轮保存的实际CAS输出，核验代码和输入哈希，不重复积分。需要新积分会话时另运行：

```powershell
D:/py/python.exe -B -X utf8 workflow/checks/2ne-endpoint-sources/run_maple.py
```

证据入口：

| 文件 | 内容 |
|---|---|
| [algebra.json](algebra.json) | 四图来源、核微分、原函数移位、四个原域多项式控制及两种配对的字母消元链 |
| [numeric.json](numeric.json) | 六组二维/三维全源对照，四组零源标准函数对照，乘积/反射/源移位检查；三个45位原函数及边界检查 |
| [coefficients.json](coefficients.json) | 全部源矩、28个物理系数、旧结果残差、M/N不相容见证及扩展空间约化 |
| [published.json](published.json) | 9个紧凑展示公式及精确回代状态 |
| [current-run.json](current-run.json) | 当前成功Maple会话位置；每个会话保存脚本、输出、日志及HyperInt/周期表指纹 |
| [render.json](render.json) | Markdown实际解析及KaTeX渲染检查 |
| [validation.json](validation.json) | 本轮命令、保护文件、来源绑定及产物指纹 |

首次系数读取器假设旧JSON有统一外层字段，出现KeyError；修正为兼容两种已存在的结构后通过。未修改旧JSON、积分输出或数学公式来使其通过。原函数边界遗漏、源微分项遗漏及物理crossing奇项遗漏均有负控制。

## 验证边界

数值控制包括双精度变量变换后的Gauss积分，以及45位的一维原函数积分；不是严格误差包围，也没有本轮原四维积分的独立高精度复算。二维/三维比较共享已解析积分的一个端点，另外一个端点用不同方法消去。

带源微分关系引入整数移位源，有限维闭包仍待研究。当前只对有限目标阶数验证M/N/τ/E覆盖；物理复z分支、一般BA/AB权重和三核全阶推广未处理。Wilson MCP未在本会话提供，formal=not_run；Kernel/Lean、旧prepared/results/catalog及历史证据保持。未提交Git。
