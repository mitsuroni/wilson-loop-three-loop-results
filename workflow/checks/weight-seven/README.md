# 权重7：低核数比较

本目录承接 [权重≤6统一](../unified-weight-six/README.md)，计算裸双核 q=3、裸三核 q=1，并在明确候选空间中比较。数学说明见 [主稿](../../../benchmarks/2ne/theory/weight-seven.md)。

当前交付：[20260917T052708438941Z](deliveries/20260917T052708438941Z/session-checks.json)。四组新鲜CAS运行成功；8项单元/负例、56项独立数值比较、13组目标表示/不相容见证、5组两R联合比较通过。89项原始输出哈希与316项依赖/快照哈希复验通过。旧47条结果、Kernel/Lean未改写，formal=not_run。

分层秩为16（旧乘积）、22（M/N）、23（六单核）、23（再加裸双核）、28（再加完整两源矩）。Qid的ε³仍可回M/N；两个R的ε¹均可回双核装饰矩，且二者之和已在M/N中。相对M/N，两R合计只新增一维，而不是两个独立方向。这不是全ε可约性或全部单核函数不可约性的证明。

可读输出：[目标精确式](../../../benchmarks/2ne/theory/weight-seven/solutions.md)、[32项两源矩](../../../benchmarks/2ne/theory/weight-seven/source-moments.md)。交付内还保存完整Hlog源矩JSON、28维选定基底、来源快照和实际检查日志。

## 候选范围

- 权重≤6已选基生成的权重7乘积。
- ζ7、M51/M33/M15、N32/N14。
- 既有六个单核母函数的五阶偶系数与 ell 乘四阶奇系数。
- 裸双核三阶及完整两源矩：q=3 的10项偶矩、q=2的6项 ell×奇矩，id/swap 共32项。

最后一项完整性只针对 Q_{π;m,n} 的两源多项式截断，不声称已枚举所有双核端点字。裸项复用全ε恒等式，其余28项逐项新积分。

## 复算入口

先用独立 Maple 会话生成新运行（不会覆盖旧运行）：

~~~powershell
D:/py/python.exe -X utf8 workflow/checks/weight-seven/start_job.py compute_targets
D:/py/python.exe -X utf8 workflow/checks/weight-seven/start_job.py compute_q_fast
D:/py/python.exe -X utf8 workflow/checks/weight-seven/start_job.py compute_single
D:/py/python.exe -X utf8 workflow/checks/weight-seven/start_job.py compute_decorated
~~~

以上是本机隐藏后台子进程；进度和完成状态存入 runs/ 与 active-*.json。分别顺序启动也可以，无需并发；每个进程从 restart 和显式周期表开始，保持发散检查开启。每个任务限时3600秒，Maple路径及HyperInt源码/周期表指纹见各运行provenance.json。

~~~powershell
D:/py/python.exe -X utf8 workflow/checks/weight-seven/progress.py
D:/py/python.exe -X utf8 workflow/checks/weight-seven/numeric_checks.py base
D:/py/python.exe -X utf8 workflow/checks/weight-seven/run_checks.py
~~~

run_checks.py 要求四个CAS任务均成功，执行单元/负例、最终比较、装饰矩数值检查、发布与只读复验；base数值检查单独执行，避免重跑已通过且未变化的目标求积。最终只读复验：

~~~powershell
D:/py/python.exe -X utf8 workflow/checks/weight-seven/verify.py
~~~

新的全轮复算会更新各 current/active 指针；历史原字节快照继续保留，不将旧快照改写成新计算。

## 文件职责

- compute_targets.mpl：原三维双核路线、三核R132/R123的一阶。
- compute_q_fast.mpl：先积两个端点，二维重算双核q=3并核对q=0。
- compute_single.mpl：六单核q=5和六个权重7 M。
- compute_decorated.mpl：完整两源三角形的28个非裸矩。
- hlog7.py：独立四周期ζ2/ζ3/ζ5/ζ7坐标版本；旧 hlog_algebra.py 原字节保留。
- compare.py：分层精确比较、显式转换、不相容见证及两R联合秩。
- numeric_checks.py：Hlog微分方程级数与独立参数求积比较。
- check_primitive.py / test_weight7.py：区间原函数展开、源截断和数据语义负例。
- verify.py：原始输出哈希、不同降维、crossing、已知源矩重构、转换与见证、交付字节检查。
- publish.py：32项源矩表、目标解表、研究主稿和不可变交付快照。
- run_checks.py：实际命令、退出码、耗时和stdout/stderr。

## 实验记录边界

首个目标命令的终端等待在120秒后返回失败，但已启动的Python/Maple进程继续运行；以进程退出、完成标记和status.json判定CAS结果，不把等待超时当成数学失败。随后任务使用显式隐藏后台启动。

第一次装饰矩脚本在积分前遇到Maple受保护变量degree，status.success=false已保留。改名为srcdeg后另开新会话，原失败运行不修改。Python版本适配中的本地启动错误同样不计入成功检查。

数值误差阈值不是严格区间包络。所有“独立”维数都限定于所声明Q-Hlog/形式周期坐标；没有周期超越性、全ε解析最小性或一般高圈闭包证明。
