# 空集族 I 解析结果

目标：13 个 prepared 0.3 空集族积分的完整 epsilon 依赖 I，以及所有极点和有限项 O(epsilon)。
当前表由计算结果生成；解析计算与数值检查状态不代表 Lean 整链完成。

先看 [完整 I 与有限项总表](solutions.md)。

入口：[公共方法](method.md)、[结构化清单](catalog.json)、[推导模板](templates/derivation.md)、[结果模板](templates/result.json)、[结果 schema](templates/result.schema.json)。

| 积分 | n | 状态 | 首项 epsilon 幂次 | 当前结果 |
|---|---:|---|---:|---|
| [2e-01](../../workspaces/2e-01/derivation.md) | 1 | computed_and_checked | -2 | [JSON](../../results/2e-01/r002.json) |
| [4e-01](../../workspaces/4e-01/derivation.md) | 2 | computed_and_checked | -4 | [JSON](../../results/4e-01/r002.json) |
| [4e-03](../../workspaces/4e-03/derivation.md) | 2 | computed_and_checked | -2 | [JSON](../../results/4e-03/r002.json) |
| [6e-01](../../workspaces/6e-01/derivation.md) | 3 | computed_and_checked | -6 | [JSON](../../results/6e-01/r002.json) |
| [6e-02](../../workspaces/6e-02/derivation.md) | 3 | computed_and_checked | -6 | [JSON](../../results/6e-02/r002.json) |
| [6e-03](../../workspaces/6e-03/derivation.md) | 3 | computed_and_checked | -6 | [JSON](../../results/6e-03/r002.json) |
| [6e-07](../../workspaces/6e-07/derivation.md) | 3 | computed_and_checked | -4 | [JSON](../../results/6e-07/r002.json) |
| [6e-08](../../workspaces/6e-08/derivation.md) | 3 | computed_and_checked | -3 | [JSON](../../results/6e-08/r002.json) |
| [6e-09](../../workspaces/6e-09/derivation.md) | 3 | computed_and_checked | -3 | [JSON](../../results/6e-09/r002.json) |
| [6e-17](../../workspaces/6e-17/derivation.md) | 3 | computed_and_checked | -4 | [JSON](../../results/6e-17/r002.json) |
| [6e-18](../../workspaces/6e-18/derivation.md) | 3 | computed_and_checked | -3 | [JSON](../../results/6e-18/r002.json) |
| [6e-19](../../workspaces/6e-19/derivation.md) | 3 | computed_and_checked | -3 | [JSON](../../results/6e-19/r002.json) |
| [6e-35](../../workspaces/6e-35/derivation.md) | 3 | computed_and_checked | -2 | [JSON](../../results/6e-35/r002.json) |

本次运行：[检查汇总](../../../../workflow/checks/empty-family/runs/20260916T083424975899Z/summary.json)。
工作区只有 Markdown；必要脚本和检查输出集中在 workflow/checks/empty-family。原稿和原始参考快照不改写。

同名数学结果由输入逐项重建后确认，仍各自保留 orbit 身份、C、sym 和证据。更高圈及更长混合字不在本批实现范围内。
