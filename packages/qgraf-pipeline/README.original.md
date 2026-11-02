# QGRAF pipeline：当前输出契约

本项目生成／分类 Wilson 四边形的源轨道，当前 1/2/3 圈分别为 2/12/149 个。组合生成器的重放与运行 QGRAF 二进制是不同操作；本次符号迁移重放已有组合生成器。

2026-10-06 起，每条 `orbit_definitions.json` 源声明 `physical_convention.id=wilson-qcd-action-path-fourier/1.0`：Wilson 线 U′=igAU，迹为晚到早，腿标签与动量编号保持；Fourier exp(−ipX)，动量 Wick +i、坐标 Wick −i；完整含 i 的 Π/Γ 包括内部闭圈号，外部 Dyson 系数 +1。Markdown 振幅和 JSON 使用相同字典，双三顶点的压缩源在展示时显式还原五传播子。

下游 Montopar 校验此声明，从原规则直接输出物理 C/J，拒绝未声明的旧振幅。raw 无 Euler/C4；最终 I 的符号来自原生 C，不由求和阶段另补族级负号。数学推导、原生与历史源的关系及使用限制见 [共同符号约定](../montopar/docs/shared-sign-conventions.md)。历史输出应通过显式版本桥读取，不能把新声明附到旧振幅上冒充重新生成。

## 递归审计导出

`qgraf-pipeline 3 --audit-recursion --output-dir outputs/3loop/audit` 单独输出
`recursion-audit.json` 和 `recursion-orbit-definitions.json`，不执行常规生产导出或绘图。
支持 1、2、3 圈的 2e、4e/3e-3g、6e/5e-3g；此模式使用组合生成器，不能与 `--out-file` 合用。

审计保留 2PR 图，新增图使用稳定的 `*-aux-<hash>` 标识，现有生产图号及其定义保持。
轨道定义可供 Montopar 读取；审计 JSON 还保留无标度记录及剔除原因。
它不计算颜色、积分或证明 shuffle 恒等式，不能将 `topology_export_only` 当作校验通过。
不同类型的生成器有不同的边分布预处理，`records` 数不是统一的完整带标签图数。
默认不带审计选项的生产行为保持原样。
