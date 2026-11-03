# 36/39 权重六端点矩：实际核验

2026-10-16。目标是将事先定义的四个独立参数矩求成 Hlog，在同一 M/N 空间内分离旧、新部分，并完整回代物理有限项。[数学推导及结果](../../../benchmarks/2ne/theory/generation/single-opposite-weight-six-20261016/derivation.md)。

结果：21 个旧候选列秩 12，加入四个矩后秩 13。新增商方向可由 E200 代表；E110、E020、E011 的商坐标分别为 1/2、0、1/2。完整 I0 的新增部分为 −E200/16，全部旧部分列在推导中。该结论仅限声明的权重六空间；未证明全阶闭合或二维不可约性。

## 执行与证据

| 步骤 | 实际执行 | 证据 |
|---|---|---|
| 三参数母方程递推至总次数四，边界至次数六 | 160 项精确检查 | [generated.json](generated.json) |
| 七个矩的 z↔1/z 反演与反演两次恢复 | 新 Maple 2025/HyperInt 会话，发散检查开启 | [reflection.json](reflection.json)，其中指向保留的脚本和日志 |
| 完整 M/N 齐次空间精确约化 | 21 列；秩 12→13；所有结果带回代或对偶见证 | [reduction.json](reduction.json) |
| 原裸端点切片、低阶及完整 I | 12 项精确连接，包含两图各 ε⁻²、ε⁻¹、ε⁰ | [reduction.json](reduction.json) |
| 独立稠密精确矩阵与原二维核 | 6 回代，3 完整见证，1 无周期见证，84 列湮灭，5 商关系，14 负控制；14 原核及 crossing 对照 | [verification.json](verification.json) |
| 实际 Wolfbook/MCP | 12 项全通过，无内核消息，操作 cf229691-438a-4d9d-a4ac-c33b72706ad1 | [audit.wl](audit.wl)，[结果](wolfbook-check.json)，[工具返回](wolfbook-result.json) |
| 文档公式和来源审计 | 直接解析正文四个 boxed TeX 公式，再回代精确 Hlog；来源绑定及本地链接 | [artifact-verification.json](artifact-verification.json) |

原核检查使用 z=0.3、0.7、160/280 节点精化，Hlog 用独立的 50 位收敛级数求值；最大绝对差约 2.82e−11。镜像原积分另在 1/z 求值。数值结果不是区间误差证书。

生成脚本不读取物理目标或 M/N 列：先由原 ODE 和 Frobenius 常数写出 generated.json，随后才读取旧结果比较。没有用物理对角线拟合出离对角线矩。符号反演、稠密秩和原核数值各自检查不同环节。

## 可重复执行

在项目根目录，以本机 Python 执行：

~~~powershell
D:/py/python.exe -B -X utf8 workflow/checks/single-opposite-weight-six-20261016/generate.py
D:/py/python.exe -B -X utf8 workflow/checks/single-opposite-weight-six-20261016/reflect.py
D:/py/python.exe -B -X utf8 workflow/checks/single-opposite-weight-six-20261016/reduce.py
D:/py/python.exe -B -X utf8 workflow/checks/single-opposite-weight-six-20261016/verify.py
D:/py/python.exe -B -X utf8 workflow/checks/single-opposite-weight-six-20261016/write_expressions.py
~~~

reflect.py 使用本机 Maple 与 HyperInt 绝对路径，周期表显式加载。每次运行保留独立时间戳目录，不覆盖旧会话日志。

通过 Wolfbook 调用 Get 加载 audit.wl；保存该次真实工具返回为 wolfbook-result.json，再执行：

~~~powershell
D:/py/python.exe -B -X utf8 workflow/checks/single-opposite-weight-six-20261016/verify_artifacts.py
~~~

本轮数学计算无失败或被丢弃的积分块。初次普通沙箱命令因 helper/setup refresh 无法启动，改用已批准的本地执行；一次文档写入调用因模板字符串语法错误未执行，修正分隔符后重发。正文审计先后因解析器不支持 TeX 的单字符分母、多行表达式及报告写出漏导入 Path 而停止；修复审计器后按同一目标重新核验，没有修改公式以绕过检查。工具输出、源哈希和原结果版本分别保留。

Wilson 数学 MCP 当前不可调用，formal=not_run。Wolfbook 和 Maple 是实际 CAS 核验，不构成 Lean 整链接纳。原 C/J、prepared、已发布积分、历史证据、PDF 及 Wilson/Lean 保持；没有 Git 提交。
