# WILSON-2.4 追加交接 r002

2026-10-17；负责人 `/root/matter_rg` / `@codex-matter-independent`。

worktree `D:/wilson-agent/.wilson/worktrees/WILSON-2.4-matter-independent`，分支 `codex/WILSON-2.4-matter-independent`，基点 `69edac9379bec71feb28ab62e25721ba7e215843`。仅在此前授权两处目录内追加 r002 文件；不修改已交付/集成的十二个 r001 文件、主库、选择器、源包或旧证据。无暂存、提交、推送。

## 修正内容

交叉审阅确认：原文几何日志 Eq5.10 与 r001 所用局部 Q/4 日志的常数连接没有闭合。字面边向量识别下每 cusp 差 log2，必须用 Γpaper=2h−4Klog2。一圈及两圈非零残差已作为实际负控制。不能用一圈 Γ1=0 重新定义未知尺度，也不能以三圈数值相同拟合字典。

因此 r001 外部比较的物理解释降为条件性候选；三个系数算术相等不撤销，但不能称已完成外部物理验证。原 C/J 两色有限形状推导不使用这项外部 Γ，仍有效；Af 有限形状仍开放。

另明确 Rn=1+z^(n eps) 是 crossing 分母，Sn=exp(n eps Lt)Rn。旧式12/13的 R2/R3 应改名𝒞2/𝒞3，为公共 C 前因子，绝不能把它们字面替换成 crossing 分母；新附录给出完整定义和原恒等式。

## 独立复跑

从集成后的主库运行，另存忽略目录：

```powershell
D:/py/python.exe workflow/checks/matter-independent-20261017/verify_conventions_r002.py --root D:/wilson-agent --output-dir D:/wilson-agent/.wilson/matter-conventions-review-r002
```

实际运行 24 项检查通过，其中 7 项负控制。脚本重核原十二文件字节，以及 r001 实际消费的十六个源文件；再检查字面日志偏移、RG所需关系、一/二/三圈非零控制、有限方案变换、两色有限形状和 C2/C3恒等式。`status=passed_with_open_convention_gap` 表示缺口被正确检出，不表示文献物理桥已验证。没有重复运行旧数值求积，也不重跑完整75项或31图。

## 新文件与集成

- 理论目录 `derivation-r002.md`。
- 检查目录 `verify_conventions_r002.py`、`conventions-r002.json`、本 `HANDOFF-r002.md` 和 `files-r002.sha256.json`。

整合窗口应优先链接 r002 附录，并把旧外部对照叙述按本勘误解释；旧记录保持冻结。新 SHA256 清单只登记追加文件，不覆盖原清单。实际校验细节见新 JSON。

## 下一项建议

开独立小任务核查原文参考文献62的几何尺度和 UV有限减除，直接从原定义推导 Eq4.6 与 Eq5.10 之间是否另有 μ/边长常数转换；若无明确来源，继续保留条件性状态。同时可另开 Af 的一般 z 有限形状投影原核验证，不把 nf2/cff 覆盖外推到 Af。

未建立新形式证明或独立上游QCD算符推导，formal=not_run；完整 w3、双 TGV、纯胶子状态保持。
