# G2 第二里程碑：独立非 cusp 预测与两色原核重建

2026-10-17，WILSON-2.4。完整数学连接见[连续推导](../../../benchmarks/renormalization/theory/matter-independent-20261017/derivation.md)。本目录只增加新证据，不修改既有 I、r001、源包、选择器或第一里程碑。

## 新增结果

1. 给出本阶 RG 的精确盲方向：h3 改 δ、单极点常数改 −δ/3、有限项同时改 −δL+g(ρ)，对偶函数 g 仍保持 RG/crossing。不能由内部一致性推出原积分唯一正确。
2. 实际检索、读取并核对 Falcioni–Gardi–Milloy 的[原论文](https://arxiv.org/abs/1909.00697) §5.2 式 (5.16)。在下述约定桥下，三个物质颜色的 Γ3 预测与项目 2h3 全部精确相同；外部公式来自形状因子/PDF/Drell–Yan 因子化，原文没有声称直接算完三圈闭合 Wilson 积分。
3. 从当前 accepted 原 prepared C/J 重建 `nf2` 与 `Ff` 两色的裸和、反项后各阶、h3 和一般正 z 有限形状，全部与现行结果相同。预测构造不读取发布的单图 I/Laurent 表；程序先写新预测再打开目标比较。研究者此前已知结果，因此不称人类盲测。

## 独立性与方案边界

文献元数据来自 INSPIRE 官方 API，策略为 `open_search`；精确 URL、原 PDF SHA256、印刷页码与公式位置在 [reference.json](reference.json)。原 PDF、完整 API 响应和页面渲染留在 worktree 忽略缓存 `.wilson/matter-independent-literature/`，不纳入交付。按 `search-papers` 与 PDF skill 进行了检索、原文提取及页面目视核对，没有摄入其他知识库。

推导 §2 写明：论文与项目使用不同符号的 Euler/尺度字典，不能把两个 μ 在相同数值下混用。对应同一裸耦合时，μF²=e^(2γE)μP²、λF(μF²)=λP(μP²)，论文自然 G 坐标尺度 μF²/4 等于项目 μx²。共同纯极点耦合 Zλ、无额外有限 Wilson 因子；β0/β1、一圈零非 cusp、两圈两个颜色都是约定控制。论文 μ d/dμ 与项目 μ² 导数之间保留因子二。常数 log shift、有限 f1/f2 Wilson 重整化会改变 Γ 的项也明确写出，没有调整这些量来拟合 h3。

一般四条有限 null 边的应用使用原文式 (5.10) 的有限段局部性，不把本项目一般正 z 解析域称为字面实平行四边形。此处是采用原文因子化关系及明确 MS 字典的外部预测核验，不是该关系的 Lean 证明。

原核重建共享上游 QCD 算符、原 C/J、相同标量化思想和 M 的定义。它提供新的原核积分、展开与实现检查，但不证明图生成完备、算符生产一定正确或一般复运动学连接。`Af` 有限形状没有因此获得独立验证。

## 实际执行与数据

`verify.py` 进行 75 项精确检查，其中 8 项负控制。包括 RG 盲方向、方案变换、原核完整导数项、C 的真实有理/颜色/Gamma 字段、Pi2/chain 几何与原域的逐字段连接、外部 β/非 cusp 预测及两色裸和到有限项。Pi2-Ff 的相消要求保留 Gamma 对数至 ε⁵，包含 ζ5，已经实际保留。

两种精度的直接一维对数差积分在 z=0.3、1、3 计算 H0/H1/H2；没有调用 GPL、MeijerG 或发布的 I。H0 与初等值比较，两个有限形状的 35/55 位提升差小于预设 10⁻²⁵ 阈值。本次最高观测提升差约 3.36×10⁻³⁰；非严格区间证书。

| z | f3,nf2(L=0) | f3,Ff(L=0) |
|---:|---:|---:|
| 0.3 | −0.217534576829153040 | 1.979610310085981091 |
| 1 | −0.223357103229524895 | 1.732799420333618530 |
| 3 | −0.218890447520335831 | 1.937901446536616152 |

实际使用 Python 3.12.10、SymPy 1.14.0、mpmath 1.3.0；最终运行详情见 [verification.json](verification.json)，实际数值、精度和每条残差都保留。数值实现初稿曾把可去点固定截断为 10⁻²⁰，导致 55 位目标出现约 10⁻³⁹ 误差；改为随工作精度变化的 sqrt(mp.eps) 门槛后，在原精度标准下通过，没有放宽阈值。

## 复现

整合后可在新忽略目录独立运行，不覆盖发布 JSON：

```powershell
D:/py/python.exe workflow/checks/matter-independent-20261017/verify.py --root D:/wilson-agent --output-dir D:/wilson-agent/.wilson/matter-independent-review-20261017
```

`--skip-numeric` 只执行精确部分；完整验收使用默认数值检查。主库只读；输出根可以指定。脚本默认从自身目录定位仓库，有任务元数据时读取 `read_only_input_root`。

- [raw-source-predictions.json](raw-source-predictions.json)：源核路线构造的两色预测及独立性说明。
- [verification.json](verification.json)：源哈希、75 项实际检查及直接积分样本。
- [reference.json](reference.json)：精确文献入口、短公式转录、约定位置与原件哈希。
- [verify.py](verify.py)：精确及数值复现入口。
- [literature_search.py](literature_search.py)、[cache_paper.py](cache_paper.py)：可选的原 API/论文缓存工具；不参与数学比较，不下载第三方执行代码。
- [HANDOFF.md](HANDOFF.md)：范围、状态与下一步。

旧 h3 提取脚本亦绑定哈希，原来的 h3 来自单极点，不是事先输入此次文献常数。因此外部非 cusp 比较不是把同一组 I 的两路尺度提取换个名字。

## 未完成范围

`Af` 一般 z 有限形状的独立原核路线仍开放；全 31 图再次积分、本项目图集的独立完备性、一般复域以及形式整链均未做。`formal=not_run`，没有 MCP/Lean 新收据。双 TGV、纯 `CF*CA^2` 及完整三圈完成状态均不在本任务范围内。
