# 高阶系数探索目标的逐项验收

2026-10-16。用户要求用 6e 更高阶 Laurent 系数的新结构辅助分析全阶生成解。本次会话目标同时要求旧阶回归、新物理系数、明确候选空间、原积分来源、跨阶关系与边界，以及一阶保留预测。本页按这些要求验收；最终机器审计为 [r004](verification-r004.json)。

| 要求 | 实际交付与判定 |
|---|---|
| 明确对象、阶数、归一化和基线 | 固定 6e-20、6e-21=22、原 C/I、rho/Gamma/单次 Euler/cosh 与原端点顺序；[基线](runs/r001/baseline.json)及[预测合同](prediction-order4-contract.json)绑定旧字节。完成。 |
| 用新源恢复旧物理零至二阶 | [物理组装](physical-runs/20261016T021005211170Z/assembly-verification.json)实际比较 24 个旧物理/裸源/分子表达式，全部精确差零；遗漏前因子的控制被拒绝。完成。 |
| 新完整 epsilon 三阶、必要时再高一阶 | [完整三阶 Hlog 公式](../../../benchmarks/2ne/theory/generation/pure-three-higher-laurent/coefficients-eps3-hlog.md)及 JSON；[完整四阶 GPL](../../../benchmarks/2ne/theory/generation/pure-three-higher-laurent/coefficients-eps4-gpl.md)及显式词表，共 77,974 行精确回读。泛型源次数与完整物理阶明确区分。完成。 |
| 明确列出 M/N、周期、低源高阶、独立端点插入及乘积 | 指定 Q 系数的 Hlog/周期坐标空间与 crossing；分层秩、全部候选列和精确对偶非成员证书均保存。独立 Euler/Beta/Mellin 仍不足，完整双 J 族覆盖全部零至三阶；[r003](verification-r003.json)复核 399 列独立和 16 个无旧三核种子的完整关系。完成指定空间比较；不声称枚举所有可能函数。 |
| 新结构有原积分解释，而非拟合余项命名 | 从完整 Lambda 原扇区消去 Q_n，识别共同父源；分部积分固定单/双 J 参数曲线。进一步证明物理父源分别为 V_n²/U_n²，完整 I20 只需混合矩，后者由单 J 源驱动。[连续证明](../../../benchmarks/2ne/theory/generation/pure-three-higher-laurent/physical-z-transport.md)。完成。 |
| 全阶生成规则、跨阶关系和明确边界 | 原十矩递推、共享两个投影、41 状态物理 z Fuchsian 系统、全矩阵、固定 z=1 级数边界及正确 epsilon 亚纯延拓。留数次数至多三，给出由前三阶驱动的 Laurent 方程。六种独立双 J 源另有固定 Beta/Gamma 单指标级数。完成。 |
| 一阶未参与构造的预测 | [合同](prediction-order4-contract.json)在完整物理四阶生成前冻结原源规则；原构造脚本和旧系数哈希保持。四个新四阶预测与独立原有序核比较通过，最大观测差 5.79e−11；另有 16 旧阶校准、54 经典 GPL、20 crossing 及四个删宽度项负控制。[实际预测报告](prediction-order4/prediction-verification-r003.json)。完成。 |
| 可复查交付与证据分级 | 连续 Markdown、三/四阶机器数据、自然源系数关系、矩阵、全部脚本与实际 CAS/数值结果；历史版本及失败诊断保留。formal=not_run，无 Git 提交。完成。 |

当前全阶表示是共享端点矩系统加明确一维源。它不是由有限阶拟合定义的“全阶猜想”：矩递推、物理投影与源恒等式都有原积分推导。相反，“所有物理源在全 epsilon 下可由有限个双 J 值直接组合”仍未证明；高阶重建只证明指定阶的有限空间成员关系。数值差也不等同于严格区间证书。

第三检查点曾把“完整独立加权有序 Q 的所有端点连接”列为下一条研究路线。随后发现的直接物理矩投影与 z 系统完成了本目标所需的原积分连接，因而没有为证明当前生成表示而穷举一般 Q。这不是声称 Q 已解决：其完整分类、全阶 BJ-only 公式、最小主积分个数、最小积分维数和 canonical epsilon form 均明确保留为更强的后续问题。用户目标没有要求证明这些最小性结论，执行计划从一开始也明确不预设它们。

保留的限制包括权重十尚未全 z 化为 Hlog/周期规范、HandyG 高权重辅助路线存在 NaN/不一致、以及混合矩流 r001 的辅助奇偶抄写错误；它们各有原始记录，均未用于虚增验收。最终采用独立 GPL 微分求值与原核预测、修正后 r002 的混合矩连接。当前要求已全部有具体交付和验证入口，未以预算或剩余工作难度替代完成判定。
