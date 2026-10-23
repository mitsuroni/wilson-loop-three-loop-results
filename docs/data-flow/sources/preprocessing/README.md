# 输入准备：原参数 J 与规范 C

2026-10-06新增[原生物理C/J入口](native_physical.py)：`prepare_native_package`要求共同符号契约、可乘C及与圈数一致的Gamma/Euler权重；不转换符号，保留original_C和J，只补一次Euler。17个自能prepared值与现行r002一致，163源Gamma权重通过。此API与冻结自能重放入口分开，默认47个0.3根及命令行不变。[实际检查](../workflow/checks/shared-sign-production/README.md)。

当前版本为 `ordered-color-euler-input/0.3`。此模块在 Kernel 外执行纯 Python 预处理，不调用 Lean、MCP 或证明数据库。

## 默认流程

```text
ordered-source-manifest.json 指定的超立方前积分
  → 检查原参数 J、嵌套边界和变量作用域
  → 完整固定颜色基底重构，提取 a1
  → 验证 Gamma 幂次，加入一次 exp(-n EulerGamma eps)
  → prepared 0.3，作为 Kernel 根输入
```

通常在仓库根目录运行：

```powershell
D:/py/python.exe -m preprocessing
```

输出 [prepared/manifest.json](../benchmarks/2ne/prepared/manifest.json)及 2e/4e/6e 三包，共 47 条。用户接受这些 C/J 为根事实；输入准备质量检查不是数学证明。

## 文件职责

| 文件 | 职责 |
|---|---|
| color.py | 精确颜色代数，在 CA²=1+2CA CF 下完整求坐标后选 a1 |
| ordered.py | 受限算术读取、变量绑定与有序积分 AST；不执行输入代码 |
| prepare.py | 来源哈希、C/Gamma 校验、规范化、可重复发布及 verify_prepared |
| __main__.py | 命令行入口 |
| test_prepare.py、test_ordered.py | 输入、负例、版本完整性回归 |

颜色不是字面搜项或大 N 截断。转换失败不能当作 a1=0。C 的 Gamma 幂必须与圈数一致，已有补偿时拒绝再处理。保存 original_C、完整颜色坐标和 Euler 记录。

J 相对于原参数导出保持不变：传播子因子、上下限和 outer-to-inner 顺序全部保留。s_to_x 与 Jacobian 只作导出追溯，不乘入原参数被积函数；换元、换序或复分支的证明留给下游。

## 更新上游导出

只在框架维护时运行：

```powershell
D:/py/python.exe -m workflow.checks.export_ordered_inputs --montopar-src C:/Users/Lenovo/Documents/learning/wilsonloop/montopar/src
D:/py/python.exe -m preprocessing
```

导出复用本地 montopar 原参数构造函数，核对全部旧超立方表达式，记录代码和数据哈希。外部 montopar 仅是导出依赖，不是日常运行依赖。不得静默改写已发布来源；新来源需审阅版本迁移。当前来源见 [清单](../benchmarks/2ne/ordered-source-manifest.json)。

原始超立方包与 [prepared 0.2 历史快照](../benchmarks/2ne/prepared/history/color-euler-input-0.2/manifest.json)保留。来源包、当前 prepared、历史 prepared 不能相互替代。

## 验证

```powershell
D:/py/python.exe -m unittest preprocessing.test_prepare preprocessing.test_ordered
```

真实结果、耗时和限制统一见 [迁移验收](../workflow/ordered-input-migration.md)。修改准备代码后应重新生成 prepared，因为清单绑定代码哈希。写入先完成全部检查，包文件原子替换，manifest 最后发布；读取必须核对哈希，不能忽略跨文件混合版本。

## 专用自能入口

`python -m preprocessing.self_energy`发布单Pi1、Pi2和双Pi1串联的物理符号2.0修订；`from preprocessing.self_energy import load_current_result`读取当前I，组装系数为+1。见[统一约定](../benchmarks/self-energy/conventions.md)。此专用入口不改变默认47个pure输入、原有prepare.py或Kernel根。旧自能无修订文件名prepared仅作历史重放。
