# T3[1] 五图求解与核验

本批计算完成：22/23/24/25（13通道）、44（24通道）的完整ε表示及35个Laurent位置、10个非零系数。最高均为简单极点，固定z>0、余项O(ε)，formal=not_run。

[推导](../../../benchmarks/tgv/families/t3-opposite/method.md) · [解表](../../../benchmarks/tgv/families/t3-opposite/solutions.md) · [运行](current-run.json) · [冻结交付](current-delivery.json) · [离线复验](verification.json)

实际完成：

- 5个原核、排序与唯一C/I连接；25↔22、24↔23两组完整反射，4负例拒绝α遗漏、顺序遗漏、分子错号、C翻倍。
- 44的TGV/pure域分离、二圈D→F框架S/P尺度连接与两系数乘积核验。
- T3推前密度低三角/中间矩形的8个多项式矩；3个t/p域分区探针和pure端点原函数。
- 新鲜HyperInt得到H22和K的e¹/e²，4个原结构系数及物理偶组合，发散检查开启，Maple CPU约32秒；完整Hlog及中间结果保留。
- 35个系数精确回代现有基底，14项配对系数恒等式。23/24各系数权重5/6；不新增生成元、不声称全ε IBP最小性。
- 9个余项样本：22用三维Gauss–Jacobi的256/512节点比较，23用二维1024/2048节点比较，节点加密差均小于余项的4%；系数用65位Hlog ODE级数求值。44用65位一维完整表达检查。余项比值符合O(ε)。
- 10组原七维QMC：两点(e,z)=(.65,.4)、(.7,1.3)，各8次scramble×2^17点；最大相对差约0.105%。低维参考另作128/256节点加密；44参考用55位一维求积。
- 12个完整ε级数有限块与独立六维Gauss–Jacobi多项式求积比较。
- 旧统一基底329系数/92别名/161哈希、mixed390哈希、二圈TGV28系数/154哈希、T3[]42系数/174哈希复验通过，见当前运行reuse.json。

数值误差为经验QMC方差或节点加密差，不是严格包围。22–25全ε为已逐项积完的显式收敛级数；44是Gamma/三角函数与既有Meijer G的闭式组合。完整级数未做快速求和或严格数值尾界。

快速离线复验，不启动Maple/数值积分：

~~~powershell
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t3-opposite/verify_results.py
~~~

重算顺序：launch.py创建新运行；等maple.log出现T3_OPPOSITE_COMPLETE后，运行formulas.py、check_math.py、check_numeric.py、full_series.py。check_numeric.py需数分钟，可后台执行并写日志。publish.py拒绝覆盖r001；新成果必须新修订，不能改写已冻结对象。

E是统一基底函数，解析结果字符串必须提供显式符号表。未修改Kernel、Lean、外部montopar或旧输入，未新增形式整链，未提交Git。原分类与旧工作区保持输入阶段快照，由新解表维护当前状态。
