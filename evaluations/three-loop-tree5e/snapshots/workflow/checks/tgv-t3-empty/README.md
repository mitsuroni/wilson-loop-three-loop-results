# T3[] 六图求解与核验

对象：5e-3g-09/10/26/30/38/39。计算完成：完整ε的已逐项积分级数、42个Laurent位置及18个非零系数，固定z>0、余项O(ε)，formal=not_run。无Kernel、Lean、外部代码改动，未提交Git。

[推导](../../../benchmarks/tgv/families/t3-empty/method.md) · [解表](../../../benchmarks/tgv/families/t3-empty/solutions.md) · [运行](current-run.json) · [冻结交付](current-delivery.json) · [复验](verification.json)

实际检查：

- 6个原始S/P/pure核、完整排序和投影C；3组包含z倒数的精确反射；4负例拒绝漏α交换、漏边序反向、P错号、C翻倍。
- 12个密度多项式矩精确验证；α尺度Jacobian、IBP微分生成式检查。
- 新鲜HyperInt：A/R各4矩，H9/H10/H26各3矩；Maple CPU约16.6秒，发散检查开启。首轮因受保护名称inner失败，新会话改名后全量重算，失败日志保留。
- 42项与现有统一基底精确回代，21个配对系数等式。10/38统一权重；领先三系数在固定Q坐标中秩3，不作为全ε IBP最小性证明。
- 9个65位低维完整I/截断余项样本；12组55位低维J与原七维QMC比较，最大相对差0.0107982061（约1.08%）。8次scramble×2^17点，30bit Sobol加2^-31半格偏移，稳定互补Beta分位，无删点/裁剪；容差为8SE加1e−4相对项，误差为经验估计。
- 15个显式级数有限块与独立四维Gauss–Jacobi求积比较。无穷和的合法性由推导中的支配界给出，未做严格数值尾界。
- 复用旧权重≤6坐标时，329系数/92别名及161哈希复验通过；mixed输入390哈希、二圈TGV28系数及154哈希复验通过。见当前运行reuse.json及日志。

六图I分成09=39、10=38、26=30三组。09只用M；26用M/N；10有限项含既有(E−2τDq4)/96，无新生成元。最高均为二阶极点。

快速离线复验（不启动Maple/数值求积）：

~~~powershell
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t3-empty/verify_results.py
~~~

重算流程：launch.py生成新运行并启动Maple；待maple.log出现T3_EMPTY_COMPLETE，再运行formulas.py、check_math.py、check_numeric.py、full_series.py。长数值任务应后台运行并保存日志。publish.py拒绝覆盖已有r001；更新成果必须使用新修订，不能覆盖冻结文件。freeze.py保存当前依赖字节。

结果JSON内的E是统一基底函数，解析字符串必须提供显式符号表，不能使用SymPy默认Euler常数E。原分类与旧工作区为输入阶段快照，不改写其绑定字节；当前结果由新解表维护。

首次冻结读取二圈旧指针时假定了新格式键delivery，实际是run/provenance.json，因而在写出交付前停止；随后的首次复验也因交付不存在停止。现已按两种真实格式读取，数学结果未改写，失败日志留存。
