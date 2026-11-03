# Π₁ montopar：实现、C/J 与发布检查

日期：2026-09-22。用户要求构建 Π₁ montopar 并输出 C/J，已完成 13/13：2e-Pi1 两图、4e-Pi1 十一图。外部实际 montopar 已回写；不是仅保留暂存实现。

[逐图 C/J](../../../benchmarks/self-energy/pi1/cj-r001.md)、[完整理论](../../../benchmarks/self-energy/pi1/theory-r001.md)、[成果入口](../../../benchmarks/self-energy/pi1/README.md)。schema 1.4，原 s 有序域、完整横向张量、完整颜色和可乘 C 账本；raw 无 Euler/a1/重数。Π₂仍只有已匹配文献来源，生产分支未实现。

## 实现与来源

- `staging/src/montopar/self_energy.py`：严格源解析、内部 δ 连接、颜色迹、几何、显式张量导数、C/J、带 SE 标记的 C4 和发布对象重建验证。复用既有精确多项式工具，未引入运行时 CAS 依赖。
- parser 在 pure 之前识别 SE；未知内部顶点返回 unsupported，Π₂/重复插入明确拒绝。CLI 先筛选后解析，避免其他未实现算符阻断指定 Π₁ 族。
- package_builder 保留原域和算符来源；拒绝分派降级、改坏的 C/J 或断开来源的输出。
- [基线](baseline.json)冻结原代码/测试/输入及 362 个已有输出。`before/` 原字节保留；`source/` 保留完整源 JSON/Markdown。源最外层负号和每图原积分次序都纳入核对。
- [最终外部补丁](upstream-pi1-final.patch)、[回写计划](writeback-plan.json)、[实际回写](writeback.json)：23 文件（8 源码/测试/文档，13 逐图输出，2 汇总包）。418 个其他基线文件保持原字节。未提交或推送 Git。首次回写后的README新增中文段落发现PowerShell管道编码损坏，已单独以UTF-8纠正；[实际补写记录](writeback-doc-fix.json)保留前后哈希，代码和C/J未改，原回写记录未重写。

## 实际检查

| 检查 | 结果与证据 |
|---|---|
| 最终暂存测试 | 370 通过，1 个 slow 外部 Wolfram 集成测试未选；[记录](tests-final.json)、[日志](tests-final.log) |
| 回写后实际外部测试 | 同样 370 通过、1 未选；[记录](tests-external.json)、[日志](tests-external.log) |
| 独立代数和来源连接 | 347 通过；[check_math.py](check_math.py)、[输出](math-checks.json) |
| 二点数值积分 | 8 控制通过；[check_numeric.py](check_numeric.py)、[结果](numeric-checks.json) |
| 旧图生产 CLI 回归 | 改动前后 39 图与 3 汇总包共 42 JSON 逐字节相同；[驱动](check_regression.py)、[记录](regression.json) |
| Pi1 实际上游再生成 | 13 图和 2 汇总包全部与发布结果及外部 workspace 相同；[记录](postwriteback.json) |
| 公式渲染 | 逐图 41、理论 10、筛选表 13，共 64 表达式零错误；[C/J](render-cj.json)、[理论](render-theory.json)、[筛选](render-screening.json) |

347 项涵盖 Fourier Gamma 递推、D 维横向性、完整源相位账本、C₂/C₃ 与未简化泡/FT 的归一化连接，十三图的独立 Gram 几何与轮廓导数、SU(2)/SU(3) 显式生成元迹、C4 尺度变换、原 Markdown 次序、六边对导数矩形和及错误旧邻边因子。13 图的四端点 J 不假设导数和消失。

数值检查实际积分两个二点图的度规/导数控制项，比较 Beta/端点表达与相邻 full/metric。它们是浮点求积及估计误差，不是认证区间，也不冒充四维有序积分求解。

42 个旧 CLI 输出覆盖 3e-3g、4e-4g、4e-3g+3g；其余 pure/mixed 路径由既有测试覆盖，未宣称所有旧类型均新鲜执行了真实颜色 CAS。已有 362 输出的保留另外按字节检查。

## 复验

```powershell
D:/py/python.exe -B -X utf8 workflow/checks/self-energy-pi1/verify.py
D:/py/python.exe -B -X utf8 workflow/checks/self-energy-pi1/verify.py --external
D:/py/python.exe -B -X utf8 workflow/checks/self-energy-pi1/run_tests.py --external
D:/py/python.exe -B -X utf8 workflow/checks/self-energy-pi1/check_math.py
D:/py/python.exe -B -X utf8 workflow/checks/self-energy-pi1/check_numeric.py
```

`stage.py`、`deploy.py` 是一次性操作并拒绝重复覆盖基线/部署记录。`publish.py` 从同一精确 JSON 生成公式，不允许替换已发布 r001 JSON 的内容。`verify.py --external` 使用真实外部源码和上游来源运行 CLI，不写回生产输出，只在本检查目录再生成比较副本。

首次 CLI 调试时，轮廓闭合替换 v4=−v1−v2−v3 未加括号，导致对边来源被误拒；补齐括号后 13 图通过，最终独立相位/几何核验覆盖该分支。没有通过放宽相位检查来接受错误表达。首次 366 测试通过后，新增 4 个来源/汇总包反例，最终为 370。

本轮没有求十三图完整积分或 Laurent 系数，没有运行 Lean/MCP 形式证明、作者 FORM/FORCER 或新增 Feynman 图生成，没有扩大默认 47 条 prepared 根。formal=not_run；一般复运动学延拓仍需后续处理。
