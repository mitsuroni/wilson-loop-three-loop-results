# [3] 族计算、表示比较与复算

日期：2026-09-17。对应6e-20、6e-21、6e-22；[解表](../../../benchmarks/2ne/families/same-channel-three/solutions.md)、[数学方法](../../../benchmarks/2ne/families/same-channel-three/method.md)。

本轮先运行预处理，47条输入通过且原包未改。新鲜HyperInt给出R123、R132的零阶及crossing分量；原始J经真实配对映射到三标记母核。三个I均无极点，21=22在全ε成立。有限项采用已有M+E，不增加新的独立权重6方向。

## 复算

~~~powershell
D:/py/python.exe -m preprocessing
D:/py/python.exe workflow/checks/same-channel-three/run_maple.py
D:/py/python.exe workflow/checks/same-channel-three/representations.py
D:/py/python.exe workflow/checks/same-channel-three/check_extension.py
D:/py/python.exe workflow/checks/same-channel-three/compact_extension.py
D:/py/python.exe -m unittest discover -s workflow/checks/same-channel-three -p test_*.py
D:/py/python.exe workflow/checks/same-channel-three/check_three.py
D:/py/python.exe workflow/checks/same-channel-three/publish.py
D:/py/python.exe workflow/checks/same-channel-three/verify_results.py
~~~

前三项表示分析输出以独占创建写入当次Maple运行目录；重复执行前应新建Maple运行，而不是覆盖旧报告。只核对已有成果时运行 verify_results.py，无需重算或发布。publish.py 新增修订，不覆盖旧rNNN。

## 实际运行

- 失败启动：runs/20260917T023351452938Z。Maple的D为受保护名称，在任何积分前报错；该运行status.success=false，原日志保留。
- 新鲜HyperInt：runs/20260917T023530229510Z。改用KernelD，新会话加载明确周期表；发散检查全程开启，2个母核、6个raw/plus/minus输出完成，Maple CPU约22.7秒。
- 独立积分检查：runs/20260917T024056933831Z。
- 结果交付：runs/20260917T025201012763Z。
- 指针：[active-run.json](active-run.json)、[current-checks.json](current-checks.json)、[current-delivery.json](current-delivery.json)。
- 范围与清单：[delivery-summary.json](delivery-summary.json)。

## 本轮实际验证范围

- 8个单元/负例：实际配对/颜色、错误负号、跨边排序、重复Euler因子、crossing、单项式矩、无极点和有限项关系。
- 9次原外到内J精确积分，ε=1,2,3、z=2/3；9次终止三重级数与之相等，另核对crossing。只评估J，不在Gamma极点评估I。
- 6个非整数原六维域Sobol案例，每例8×2^17样本，种子、均值、容差保存。
- 4个raw有限母核的完整Hlog/独立降维求积对照，并分别核对其even/odd分量。
- 9个完整ε母核crossing对照。
- 18个全I减有限项的ε减半样本，检查O(ε)缩放。
- 两个有限母核的纯M、M+log、旧τ空间比较，保存6个不相容见证。
- 旧τ+已有E、M+已有E的精确转换及反代；确认R123+R132plus=M01³/36。

所有比较用完整Hlog/ζ坐标；没有用PSLQ构造解析答案。不可表示结论仅针对列明的固定权重有理多项式候选空间。

## 文件职责

reduce_three.py核对原根与C；full_series.py计算已积完的有限多项式块；compute_coefficients.mpl实际求值；check_three.py使用独立的稳定幂差商降维求积与原域数值；三个表示脚本比较函数空间；publish.py生成文档/结果快照；verify_results.py只读重检公式、表示关系、输入和依赖哈希。

已有M、τ、E输出以具体版本和哈希复用并保留副本；没有改写以前的结果或证据。新工作区只有Markdown，复算代码与大JSON留在本目录。

formal=not_run；没有运行Kernel/Lean全链。数值误差是统计或收敛差分估计，不是严格包络。全ε级数按有限块组织收敛，未承诺快速全参数求值；复z/i0与一般正幂仍开放。

本次成果验证器还重检了3个结果、21个系数和687项含快照哈希。[会话检查记录](session-checks.json)汇总实际运行；[导航与47条目录核对](navigation-and-inventory.json)确认50条新文档链接有效及47个当前根均有结果。目录核对不代表重新执行47条积分。旧空集修订按其声明极点阶存储系数（如一圈从−2开始），未为统一格式改写旧结果。
