# 一圈、两圈连续推导标准试行与审核

本轮完成了确定标准、安装 skill、整理14篇逐积分推导及独立审核。当前阅读入口是[积分索引](../../../integrals/README.md)，规范入口是[标准1.0](../../integral-documentation-standard.md)，实际执行指南是[wilson-integral-derivations](../../skills/wilson-integral-derivations/SKILL.md)。用户安装副本位于 `C:/Users/Lenovo/.codex/skills/wilson-integral-derivations`；[安装记录](skill-install.json)绑定与项目规范相同的三个文件。

## 试行结果

一圈2个积分写到 epsilon³，两圈12个写到 epsilon²，共96个Laurent系数位置。每图有 README、derivation、reproduce、validation、精确input/result绑定及code/data/checks入口；按实际用途共享原程序、原语和大输出。当前数据索引逐对象区分原核、已积分子积分、已组装I、函数字典和后置比较表。

[独立逐篇审核](review/README.md)认为14篇正文都已接通，其中9篇明确采用已求引理/计算结果。审核按具体数学式、操作条件、边界常数和系数组装作出结论；完整证据与范围在[逐图audit.json](review/audit.json)。这不是14次新的独立积分计算或正式证明。

两个重点例子可以直接检验标准的实际效果：

- [4e-04](../../../integrals/two-loop/4e-04/derivation.md)：原四维有序积分→两次线性求积→二维Q→K的Euler导数/Beta求和与积分常数→M/N换基→最高要求阶。附录给出端点转换和解析常数。
- [3e-3g-04](../../../integrals/two-loop/3e-3g-04/derivation.md)：完整六面与recoil→D/H/O实际核→Euler减除和共同收敛带→六个完整HLOG jet→历史Laurent组装。当前Python直接采用已组装的归一化I高阶系数，再核表示；它没有重跑D/H/O求积或卷积。

审核发现并修正了M05常数漏项、TGV03递推记号、初始C/J绑定缺失、TGV04代码职责描述及最简单条目的冗余函数材料。审阅者还逐项核对了六张新HLOG表的168个精确坐标，见[实际检查](review/hlog-transcription-check.json)。这些坐标比较不独立验证原HyperInt积分。

## 实际运行与保存的证据

[执行汇总](replay.json)记录了新入口的主运行和真实文件副本重放。14份输出均与主运行及采用result逐字节一致：

| 范围 | 实际通过 | 覆盖边界 |
|---|---|---|
| 一圈 | 80项检查，包含27个数值控制；另记录4个负控制 | 原几何、参数原函数、Gamma/M生成及后置系数比较 |
| 两圈 | 149项检查、2个负控制 | 母函数/Gamma展开、采用I及完整HLOG表示 |
| 两圈单独控制 | 68项 | 原核整指数矩、Pi1混合导数/角点、Beta与M/N/tau及4e04二维系数；不含TGV04新完整五维积分 |

原始主运行输出保存在 [evidence/primary](evidence/primary/)，副本重放的紧凑证据保存在 [evidence/replay](evidence/replay/)。完整临时目录和隔离副本不加入版本控制。首轮新两圈薄适配器缺少同级模块搜索路径，失败日志保留；修复仅作用于新入口，随后两圈重试和全部副本阶段通过。没有为了隐藏失败而改写原运行记录。

从仓库根目录可在新输出目录重放支持的三个阶段：

```powershell
D:/py/python.exe -X utf8 -B workflow/checks/integral-documentation-low-loops-20261103/run_checks.py --output-root replay-output/low-loop-documentation-new
```

各图 reproduce.md 给自身入口、输入输出和原CAS复用范围。历史Maple/Wolfram程序保留原字节及原绝对加载/输出设置；当前命令不会启动它们。

## 对标准的审核结论

[前向试用评估](review/standard-forward-test.md)建议维持1.0，不降低连续推导要求。简单图合并叙述，困难图用明确数学附录；精确共享引理可以采用，但原图的非平凡变换和边界不能以链接代替。正文、实现、复现、数学验证与发布五项分别验收。

本轮的数学结论限于正文给定的正实z及局部epsilon延拓，保留已采用物理根、函数字典和历史计算的原信任范围。没有新增任意复z或z=-1证明、冷CAS全重算、Lean证书或科学选择器切换。

[最终资料检查](artifact-check.json)只负责链接、字节、精确输入绑定与源保护；它不是正文数学审核。4650个受保护原文件、4146个公开文件、主HEAD/暂存及公开7e77ddb保持，实际记录见该检查及[交接](handoff.md)。当前全部是本地整理，未Git提交或同步GitHub。

下一步可按同一标准处理已有三圈6e的39图整批资料。原P2/P5/P6/P7的包装与重放完成仍按历史范围解释，不据此默认其正文已通过本标准；本轮没有自动继续P8。
