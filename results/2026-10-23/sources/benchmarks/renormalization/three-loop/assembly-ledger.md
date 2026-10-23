# 三圈逐轨道筛选与组装账本

2026-10-09双TGV系数更正：13/22的ε⁻³至ε⁰已由独立复算修正为r002，[双TGV current](../../double-tgv/4e-3g+3g/current.json)已选入。历史assembly-ledger-r001.json仍绑定两图r001，保留原字节；依赖这些旧系数的部分和/纯CFCA²组装需要刷新，不自动继承旧检查状态。此次未重组完整w3或改动其他颜色结果。[新解析式与组装审计](../../../workflow/checks/double-tgv-independent-review/README.md)。

2026-10-06求和层更新：CFCA TFnf的31图完整求和、普通反项／轮廓减除及RG通过，新[物理current](current.json)选入正式结果。新腿Π₁和Γ₁原源的三槽色序与P相反，物理恢复倍率均为−1，见[独立桥](../../../workflow/checks/w3-ca-nf/physical-bridge.json)。原r001逐轨道账本保留历史；其腿/Γ占位或缺项文字不代表当前状态，原I也未因物理号桥而重写。纯CFCA²仍缺双TGV10图，完整w3未完成。

2026-10-02求解层更新：树级TGV腿Π₁的实际Montopar C/J及12个唯一I均已交付全部极点至有限项，见[现行12轨道](../../self-energy/tgv-pi1/current.json)与[12收口及完整桥](../../self-energy/tgv-pi1/workspaces/3e-3g+Pi1-12/derivation.md)。本页原表保留P0初始快照，其“C/J待生产”不能作为当前求解状态；这12个结果的物理源号、外尺度和C4/web计权接入仍需另行验收，完整w3/RG尚未完成。

2026-10-01腿Π₁生产更新：12条`3e-3g+Pi1`已有实际raw C/J及IBP导出，[整族输入](../../self-energy/tgv-pi1/sources/montopar-r001.json)。下表r001中的“C/J待生产”是原账本历史状态，当前缺口是I求值及prepared/旋转归一化；不手改内容绑定的assembly-ledger-r001.json，也不据此增加已求解计数。两个物理颜色通道的整图系数仍待卷积组装。

2026-10-01后续montopar审计修正：双TGV的26条四外腿记录是有明确解码规则的合法压缩源，已通过当前契约；不能因未显式列出第五条传播子判为错误。原r001告警保留历史语义，由[模块理论审计](../../theory/montopar/double-tgv.md)补充解释。Γ₁四条空占位仍待处理，数学结果和26条未求完计数不变。

实际源、哈希、当前数学包、完整颜色与拓扑比较见[assembly-ledger-r001.json](assembly-ledger-r001.json)。颜色列是现行I的恢复乘子，另乘λ³S3；m、a1与crossing不再乘。非零Lorentz/参数核保留在原C/J中，未知项未作零项删除。

| 轨道 | 边分布 | I或零项 | 完整颜色恢复 | m | 当前动量源 | web或筛选 |
|---|---|---|---|---:|---|---|
| 2e-Pi1+Pi1-01 | [0, 0, 1, 1] | [已解至有限项](../../../benchmarks/self-energy/pi1-chain/results/2e-Pi1+Pi1-01/r002.json) | CF*B1(eps)^2 | 4 | 结构字段已给出 | single_connected_source_color |
| 2e-Pi1+Pi1-02 | [0, 1, 0, 1] | [已解至有限项](../../../benchmarks/self-energy/pi1-chain/results/2e-Pi1+Pi1-02/r002.json) | CF*B1(eps)^2 | 2 | 结构字段已给出 | single_connected_source_color |
| 2e-Pi2-01 | [0, 0, 1, 1] | [已解至有限项](../../../benchmarks/self-energy/pi2/results/2e-Pi2-01/r002.json) | use all three stored operator/color channel coefficients | 4 | 结构字段已给出 | single_connected_source_color |
| 2e-Pi2-02 | [0, 1, 0, 1] | [已解至有限项](../../../benchmarks/self-energy/pi2/results/2e-Pi2-02/r002.json) | use all three stored operator/color channel coefficients | 2 | 结构字段已给出 | single_connected_source_color |
| 3e-3g+Pi1-01 | [0, 0, 1, 2] | C/J待生产 | pending C/J and I normalization | 4 | 结构字段已给出 | pending_source_and_color_normalization |
| 3e-3g+Pi1-02 | [0, 0, 1, 2] | C/J待生产 | pending C/J and I normalization | 4 | 结构字段已给出 | pending_source_and_color_normalization |
| 3e-3g+Pi1-03 | [0, 0, 1, 2] | C/J待生产 | pending C/J and I normalization | 4 | 结构字段已给出 | pending_source_and_color_normalization |
| 3e-3g+Pi1-04 | [0, 0, 2, 1] | C/J待生产 | pending C/J and I normalization | 4 | 结构字段已给出 | pending_source_and_color_normalization |
| 3e-3g+Pi1-05 | [0, 0, 2, 1] | C/J待生产 | pending C/J and I normalization | 4 | 结构字段已给出 | pending_source_and_color_normalization |
| 3e-3g+Pi1-06 | [0, 0, 2, 1] | C/J待生产 | pending C/J and I normalization | 4 | 结构字段已给出 | pending_source_and_color_normalization |
| 3e-3g+Pi1-07 | [0, 1, 0, 2] | C/J待生产 | pending C/J and I normalization | 4 | 结构字段已给出 | pending_source_and_color_normalization |
| 3e-3g+Pi1-08 | [0, 1, 0, 2] | C/J待生产 | pending C/J and I normalization | 4 | 结构字段已给出 | pending_source_and_color_normalization |
| 3e-3g+Pi1-09 | [0, 1, 0, 2] | C/J待生产 | pending C/J and I normalization | 4 | 结构字段已给出 | pending_source_and_color_normalization |
| 3e-3g+Pi1-10 | [0, 1, 1, 1] | C/J待生产 | pending C/J and I normalization | 4 | 结构字段已给出 | pending_source_and_color_normalization |
| 3e-3g+Pi1-11 | [0, 1, 1, 1] | C/J待生产 | pending C/J and I normalization | 4 | 结构字段已给出 | pending_source_and_color_normalization |
| 3e-3g+Pi1-12 | [0, 1, 1, 1] | C/J待生产 | pending C/J and I normalization | 4 | 结构字段已给出 | pending_source_and_color_normalization |
| 3e-Gamma1-01 | [0, 0, 1, 2] | 源占位 | pending tensor/channel C/J; do not assume a common scalar I | 4 | 缺口：central_vertex_placeholder | pending_source_and_color_normalization |
| 3e-Gamma1-02 | [0, 0, 2, 1] | 源占位 | pending tensor/channel C/J; do not assume a common scalar I | 4 | 缺口：central_vertex_placeholder | pending_source_and_color_normalization |
| 3e-Gamma1-03 | [0, 1, 0, 2] | 源占位 | pending tensor/channel C/J; do not assume a common scalar I | 4 | 缺口：central_vertex_placeholder | pending_source_and_color_normalization |
| 3e-Gamma1-04 | [0, 1, 1, 1] | 源占位 | pending tensor/channel C/J; do not assume a common scalar I | 4 | 缺口：central_vertex_placeholder | pending_source_and_color_normalization |
| 4e-3g+3g-01 | [0, 0, 1, 3] | [已解至有限项](../../../workflow/checks/w3-p0/sources/result-dependencies/4e-3g+3g-01/r001.json) | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-02 | [0, 0, 1, 3] | 颜色零 | 0 | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | color零项 |
| 4e-3g+3g-03 | [0, 0, 1, 3] | 未解 | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-04 | [0, 0, 2, 2] | [已解至有限项](../../../benchmarks/double-tgv/4e-3g+3g/results/4e-3g+3g-04/r001.json) | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-05 | [0, 0, 2, 2] | 颜色零 | 0 | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | color零项 |
| 4e-3g+3g-06 | [0, 0, 2, 2] | 部分结果 | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-07 | [0, 0, 3, 1] | [已解至有限项](../../../workflow/checks/w3-p0/sources/result-dependencies/4e-3g+3g-07/r001.json) | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-08 | [0, 0, 3, 1] | 颜色零 | 0 | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | color零项 |
| 4e-3g+3g-09 | [0, 0, 3, 1] | 未解 | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-10 | [0, 1, 0, 3] | [已解至有限项](../../../workflow/checks/w3-p0/sources/result-dependencies/4e-3g+3g-10/r001.json) | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-11 | [0, 1, 0, 3] | 颜色零 | 0 | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | color零项 |
| 4e-3g+3g-12 | [0, 1, 0, 3] | 未解 | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-13 | [0, 1, 1, 2] | [已解至有限项](../../../workflow/checks/w3-p0/sources/result-dependencies/4e-3g+3g-13/r001.json) | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-14 | [0, 1, 1, 2] | 颜色零 | 0 | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | color零项 |
| 4e-3g+3g-15 | [0, 1, 1, 2] | 未解 | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-16 | [0, 1, 2, 1] | 部分结果 | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-17 | [0, 1, 2, 1] | 颜色零 | 0 | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | color零项 |
| 4e-3g+3g-18 | [0, 1, 2, 1] | 未解 | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-19 | [0, 2, 0, 2] | [已解至有限项](../../../workflow/checks/w3-p0/sources/result-dependencies/4e-3g+3g-19/r001.json) | CA^2*CF | 2 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-20 | [0, 2, 0, 2] | 颜色零 | 0 | 2 | 缺口：two_TGV_momentum_definition_missing_or_wrong | color零项 |
| 4e-3g+3g-21 | [0, 2, 0, 2] | 未解 | CA^2*CF | 2 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-22 | [0, 2, 1, 1] | [已解至有限项](../../../workflow/checks/w3-p0/sources/result-dependencies/4e-3g+3g-22/r001.json) | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-23 | [0, 2, 1, 1] | 颜色零 | 0 | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | color零项 |
| 4e-3g+3g-24 | [0, 2, 1, 1] | 未解 | CA^2*CF | 4 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-25 | [1, 1, 1, 1] | 未解 | CA^2*CF | 2 | 缺口：two_TGV_momentum_definition_missing_or_wrong | single_connected_source_color |
| 4e-3g+3g-26 | [1, 1, 1, 1] | 颜色零 | 0 | 1 | 缺口：two_TGV_momentum_definition_missing_or_wrong | color零项 |
| 4e-4g-01 | [0, 0, 1, 3] | Lorentz零 | 0 | 4 | 结构字段已给出 | Lorentz/lightlike零项 |
| 4e-4g-02 | [0, 0, 2, 2] | [已解至有限项](../../../benchmarks/fgv/4e-4g/results/4e-4g-02/r001.json) | CA^2*CF | 4 | 结构字段已给出 | single_connected_source_color |
| 4e-4g-03 | [0, 0, 3, 1] | Lorentz零 | 0 | 4 | 结构字段已给出 | Lorentz/lightlike零项 |
| 4e-4g-04 | [0, 1, 0, 3] | Lorentz零 | 0 | 4 | 结构字段已给出 | Lorentz/lightlike零项 |
| 4e-4g-05 | [0, 1, 1, 2] | [已解至有限项](../../../benchmarks/fgv/4e-4g/results/4e-4g-05/r001.json) | CA^2*CF | 4 | 结构字段已给出 | single_connected_source_color |
| 4e-4g-06 | [0, 1, 2, 1] | [已解至有限项](../../../benchmarks/fgv/4e-4g/results/4e-4g-06/r001.json) | CA^2*CF | 4 | 结构字段已给出 | single_connected_source_color |
| 4e-4g-07 | [0, 2, 0, 2] | [已解至有限项](../../../benchmarks/fgv/4e-4g/results/4e-4g-07/r001.json) | CA^2*CF | 2 | 结构字段已给出 | single_connected_source_color |
| 4e-4g-08 | [0, 2, 1, 1] | [已解至有限项](../../../benchmarks/fgv/4e-4g/results/4e-4g-08/r001.json) | CA^2*CF | 4 | 结构字段已给出 | single_connected_source_color |
| 4e-4g-09 | [1, 1, 1, 1] | [已解至有限项](../../../benchmarks/fgv/4e-4g/results/4e-4g-09/r001.json) | CA^2*CF | 1 | 结构字段已给出 | single_connected_source_color |
| 4e-Pi1-01 | [0, 0, 2, 2] | [已解至有限项](../../../benchmarks/self-energy/pi1/results/4e-Pi1-01/r002.json) | CA*CF*B1(eps) | 4 | 结构字段已给出 | 副本颜色已核对 |
| 4e-Pi1-02 | [0, 0, 2, 2] | [已解至有限项](../../../benchmarks/self-energy/pi1/results/4e-Pi1-02/r002.json) | CA*CF*B1(eps) | 4 | 结构字段已给出 | 副本颜色已核对 |
| 4e-Pi1-03 | [0, 1, 1, 2] | [已解至有限项](../../../benchmarks/self-energy/pi1/results/4e-Pi1-03/r002.json) | CA*CF*B1(eps) | 4 | 结构字段已给出 | 副本颜色已核对 |
| 4e-Pi1-04 | [0, 1, 1, 2] | [已解至有限项](../../../benchmarks/self-energy/pi1/results/4e-Pi1-04/r002.json) | CA*CF*B1(eps) | 4 | 结构字段已给出 | 副本颜色已核对 |
| 4e-Pi1-05 | [0, 1, 2, 1] | [已解至有限项](../../../benchmarks/self-energy/pi1/results/4e-Pi1-05/r002.json) | CA*CF*B1(eps) | 4 | 结构字段已给出 | 副本颜色已核对 |
| 4e-Pi1-06 | [0, 1, 2, 1] | [已解至有限项](../../../benchmarks/self-energy/pi1/results/4e-Pi1-06/r002.json) | CA*CF*B1(eps) | 4 | 结构字段已给出 | 副本颜色已核对 |
| 4e-Pi1-07 | [0, 2, 0, 2] | [已解至有限项](../../../benchmarks/self-energy/pi1/results/4e-Pi1-07/r002.json) | CA*CF*B1(eps) | 2 | 结构字段已给出 | 副本颜色已核对 |
| 4e-Pi1-08 | [0, 2, 0, 2] | [已解至有限项](../../../benchmarks/self-energy/pi1/results/4e-Pi1-08/r002.json) | CA*CF*B1(eps) | 2 | 结构字段已给出 | 副本颜色已核对 |
| 4e-Pi1-09 | [0, 2, 1, 1] | [已解至有限项](../../../benchmarks/self-energy/pi1/results/4e-Pi1-09/r002.json) | CA*CF*B1(eps) | 4 | 结构字段已给出 | 副本颜色已核对 |
| 4e-Pi1-10 | [0, 2, 1, 1] | [已解至有限项](../../../benchmarks/self-energy/pi1/results/4e-Pi1-10/r002.json) | CA*CF*B1(eps) | 4 | 结构字段已给出 | 副本颜色已核对 |
| 4e-Pi1-11 | [1, 1, 1, 1] | [已解至有限项](../../../benchmarks/self-energy/pi1/results/4e-Pi1-11/r002.json) | CA*CF*B1(eps) | 2 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-01 | [0, 0, 2, 3] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-01/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-02 | [0, 0, 2, 3] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-02/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-03 | [0, 0, 2, 3] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-03/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-04 | [0, 0, 2, 3] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-04/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-05 | [0, 0, 3, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-05/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-06 | [0, 0, 3, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-06/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-07 | [0, 0, 3, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-07/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-08 | [0, 0, 3, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-08/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-09 | [0, 1, 1, 3] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-09/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-10 | [0, 1, 1, 3] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-10/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-11 | [0, 1, 1, 3] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-11/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-12 | [0, 1, 1, 3] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-12/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-13 | [0, 1, 2, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-13/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-14 | [0, 1, 2, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-14/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-15 | [0, 1, 2, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-15/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-16 | [0, 1, 2, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-16/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-17 | [0, 1, 2, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-17/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-18 | [0, 1, 3, 1] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-18/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-19 | [0, 1, 3, 1] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-19/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-20 | [0, 1, 3, 1] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-20/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-21 | [0, 1, 3, 1] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-21/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-22 | [0, 2, 0, 3] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-22/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-23 | [0, 2, 0, 3] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-23/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-24 | [0, 2, 0, 3] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-24/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-25 | [0, 2, 0, 3] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-25/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-26 | [0, 2, 1, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-26/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-27 | [0, 2, 1, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-27/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-28 | [0, 2, 1, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-28/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-29 | [0, 2, 1, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-29/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-30 | [0, 2, 1, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-30/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-31 | [0, 2, 2, 1] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-31/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-32 | [0, 2, 2, 1] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-32/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-33 | [0, 2, 2, 1] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-33/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-34 | [0, 2, 2, 1] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-34/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-35 | [0, 2, 2, 1] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-35/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-36 | [0, 3, 1, 1] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-36/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-37 | [0, 3, 1, 1] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-37/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-38 | [0, 3, 1, 1] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-38/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-39 | [0, 3, 1, 1] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-39/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-40 | [1, 1, 1, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-40/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-41 | [1, 1, 1, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-41/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-42 | [1, 1, 1, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-42/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-43 | [1, 1, 1, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-43/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 5e-3g-44 | [1, 1, 1, 2] | [已解至有限项](../../../benchmarks/tgv/5e-3g/results/5e-3g-44/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-01 | [0, 0, 3, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-01/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-02 | [0, 0, 3, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-02/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-03 | [0, 0, 3, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-03/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-04 | [0, 1, 2, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-04/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-05 | [0, 1, 2, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-05/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-06 | [0, 1, 2, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-06/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-07 | [0, 1, 3, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-07/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-08 | [0, 1, 3, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-08/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-09 | [0, 1, 3, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-09/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-10 | [0, 2, 1, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-10/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-11 | [0, 2, 1, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-11/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-12 | [0, 2, 1, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-12/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-13 | [0, 2, 2, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-13/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-14 | [0, 2, 2, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-14/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-15 | [0, 2, 2, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-15/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-16 | [0, 2, 2, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-16/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-17 | [0, 2, 3, 1] | [已解至有限项](../../../benchmarks/2ne/results/6e-17/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-18 | [0, 2, 3, 1] | [已解至有限项](../../../benchmarks/2ne/results/6e-18/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-19 | [0, 2, 3, 1] | [已解至有限项](../../../benchmarks/2ne/results/6e-19/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-20 | [0, 3, 0, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-20/r001.json) | CA^2*CF | 2 | 结构字段已给出 | 副本颜色已核对 |
| 6e-21 | [0, 3, 0, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-21/r001.json) | CA^2*CF | 2 | 结构字段已给出 | 副本颜色已核对 |
| 6e-22 | [0, 3, 0, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-22/r001.json) | CA^2*CF | 2 | 结构字段已给出 | 副本颜色已核对 |
| 6e-23 | [0, 3, 1, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-23/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-24 | [0, 3, 1, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-24/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-25 | [0, 3, 1, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-25/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-26 | [0, 3, 2, 1] | [已解至有限项](../../../benchmarks/2ne/results/6e-26/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-27 | [0, 3, 2, 1] | [已解至有限项](../../../benchmarks/2ne/results/6e-27/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-28 | [0, 3, 2, 1] | [已解至有限项](../../../benchmarks/2ne/results/6e-28/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-29 | [1, 1, 1, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-29/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-30 | [1, 1, 1, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-30/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-31 | [1, 1, 1, 3] | [已解至有限项](../../../benchmarks/2ne/results/6e-31/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-32 | [1, 1, 2, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-32/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-33 | [1, 1, 2, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-33/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-34 | [1, 1, 2, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-34/r001.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-35 | [1, 1, 2, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-35/r002.json) | CA^2*CF | 4 | 结构字段已给出 | 副本颜色已核对 |
| 6e-36 | [1, 2, 1, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-36/r002.json) | CA^2*CF | 2 | 结构字段已给出 | 副本颜色已核对 |
| 6e-37 | [1, 2, 1, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-37/r001.json) | CA^2*CF | 2 | 结构字段已给出 | 副本颜色已核对 |
| 6e-38 | [1, 2, 1, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-38/r001.json) | CA^2*CF | 2 | 结构字段已给出 | 副本颜色已核对 |
| 6e-39 | [1, 2, 1, 2] | [已解至有限项](../../../benchmarks/2ne/results/6e-39/r001.json) | CA^2*CF | 2 | 结构字段已给出 | 副本颜色已核对 |

111已解、12零项、26尚未完整交付。解释与下一步见[P0主稿](p0-assembly.md)。
