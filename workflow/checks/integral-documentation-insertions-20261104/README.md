# 十五个二点／交换插入：连续推导与逐篇审核

本轮按[标准1.0与skill](../../integral-documentation-standard.md)，处理11个4e-Pi1、2个串联Pi1和2个单Pi2。[当前15图入口](../../../integrals/three-loop/README.md#two-point-and-exchange-insertions)和[共享数学](../../../integrals/families/three-loop-insertions/README.md)连接各自完整原核、实际端点约化、已求子积分、全部标量／几何Laurent系数及完整颜色有限项。各图机器结果包含完整七阶颜色级数。

## 实际补齐的推导

每图有README、derivation、reproduce、validation，配套精确input/result绑定和code/data/checks入口。三份共享数学文档和一份明确的函数附录承担必要的共同证明。按数学用途共享代码，保留自身原始域和代入，而非复制十五份假独立求积程序。

- Pi1加交换的六路线：完整横向张量→源几何分子／Fourier参数核→各图端点积分和完整接触面→已求加权H、M及有序T→实际Gamma／crossing与B1颜色卷积。P3、P5、P6的边界不省略；P5展示八个角点和完整rho1 HLOG。
- 串联与Pi2：实际Dyson张量→两端点原核→相邻幂积分或对边四角→完整已求H0→原生C所定符号→全部标量系数和外算符。两个Pi1的B1²与单Pi2的完整A+B/U分别定义。Pi2六个A/B有理式、表观phi极点抵消、U与逆U五阶保守预算均写出，保留Ff通道。

可直接读[4e-Pi1-07](../../../integrals/three-loop/4e-Pi1-07/derivation.md)和[2e-Pi2-02](../../../integrals/three-loop/2e-Pi2-02/derivation.md)检查困难部分。正文给完整必要函数值与最高阶有限求和，没有把待求余积分换名后当作答案。

[逐图科学审核](review/README.md)及[具体证据](review/audit.json)分别评价正文、实现、实际重放、数学验证和发布。审阅者实查21个新打印表达式，包含10个M、rho1、四个chain/Pi2 M式及六个A/B有理式；对应107个非零HLOG坐标，见[实际检查](review/math-transcription-check.json)。这是转录／算术核验，保留原CAS求值与内部Pi2算符的采用范围。

审核修正了一个真正的条件错误：对边横向分子在零角点通常不为零，必须使用D^(-2+2Re eps)与二维径向测度证明严格正Re eps时可积；不能把H矩在零附近的全纯性推广到完整原核。历史Pi1盒装正号也未沿用，当前负号由原生C及connected投影确定。

## 本轮实际运行

[执行汇总](replay.json)记录新主运行：15/15生成成功，8.96秒，无新外部CAS；完整I参考表未被生成读取。105个标量／几何位置及84个已有颜色位置与原表一致。Pi1原表只存标量，故其154个新颜色位置不声称与旧颜色表比较。

全部105标量＋238生成颜色＝343个函数表示，经独立有理shuffle回代通过。另跑原控制：15项源／算符／颜色连接、30项原序与独立cube多项式积分、6个Gamma有限系数双半径Cauchy控制及6个负控制。其余9图未新增原多维finite数值积分，经验误差不称严格区间证书。

[真实复制闭包](evidence/isolation/replay.json)在完整reference目录缺席时生成，保存冻结哈希后再加入比较数据。15图标量和全部颜色数据相同，完整generated.json也与主运行逐字节一致。两次使用相同Python环境和采用函数／算符，不是两次独立解析积分。

数据索引逐图列数学必要M／rho对象，并单列批程序全局加载的字典。P1/P2和相邻两端点图不被误写为依赖全部M。相邻图无需H0，其预算与批程序统一计算H0到二阶、U到五阶的保守实现选择也明确区分。

## 复现与交接

在项目根目录选择新输出目录：

```powershell
D:/py/python.exe -X utf8 -B workflow/checks/integral-documentation-insertions-20261104/run_batch.py --output-root replay-output/insertions-continuous-new
```

生成、冻结、旧表比较、表示回代和范围明确的原控制依次执行；原始输出和日志保存在[evidence](evidence/)。各图薄入口运行相同15项批次，正文与data索引明确说明原Python/Wolfram/Maple资料各自职责。旧绝对路径writer不自动启动。

数学结论限于声明的正实z、初始收敛条带和局部epsilon亚纯延拓。原物理C/J、已求M/T及完整D维Pi2算符保持其原采用身份；无新内部两圈Pi2求积、任意复z／负一边界或formal升级。旧科学结果、四腿工作流、53篇先前条目、公开副本和Git状态的保持见[机械检查](artifact-check.json)；该检查不授予数学证明。

[交接](handoff.md)与[最终文件清单](manifest.json)记录范围。任务WILSON-5.27完成的是连续文档与审核；本轮未提交Git或同步GitHub。文档线下一批为44个5e-3g树三胶子来源；四腿2244的独立研究线保持原计划。
