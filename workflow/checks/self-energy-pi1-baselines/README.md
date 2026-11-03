# Π₁输入准备与B_A、B_O求解

2026-09-23。[当前唯一I定义](../../../benchmarks/self-energy/pi1/theory.md)、[B_A](../../../benchmarks/self-energy/pi1/families/b-adjacent/solutions.md)、[B_O](../../../benchmarks/self-energy/pi1/families/b-opposite/solutions.md)。

已完成全部13图的专用颜色/Euler输入准备，并完成B_A/B_O的全εI与14个Laurent位置（7非零）。当前Π₁完成2/13，其余11图尚未求值。

## 归一化准备

`prepare.py`复用`preprocessing.color.project`的完整固定基底重构，分别保存CA和TF nf通道坐标，再提取a₁。二圈a₁=1，三圈a₁=−1/2；原完整Wilson迹及原C保留。整个B₁(D)留在外部颜色前因子中，因此两通道共用一个标量I。

按总圈数检查原Γ(1−ε)幂，确认尚未补偿后只加入一次exp(−n EulerGamma ε)。J、bounds和迭代顺序逐图不变。专用prepared不扩大默认47条Kernel根。[清单](../../../benchmarks/self-energy/pi1/prepared/manifest.json)、[39项准备检查](preparation-checks.json)。

原动量图贡献现记作A=C_source J；唯一I按颜色、耦合、尺度和C4系数提取定义。此前把整个CJ称为目标I的规划已纠正，而不是重算或改变原张量参数化。

## 实际计算与检查

- [compute.py](compute.py)、[30项精确检查](exact-checks.json)：连接完整输入核、二点端点积分、C4与I归一化、BO完整核的导数恒等式、Meijer G母核和全部14系数；精确核对旧Hlog中H₀(0)=M01、2H₀′(0)=M11。
- [numeric.py](numeric.py)、[28项数值/余项控制](numeric-checks.json)：9个BO原完整核的独立二维扇区求积、3个BA端点求积、6个M函数定义积分、4个保留nf=0/4的实际prepared C→I核对，以及6组ε减半的O(ε)余项缩放。
- BO数值使用原完整多项式在两对称三角域上的变换，径向变量消去角点奇性；128/256阶Gauss结果与Meijer G闭式比较。与解析计算中使用的H₀加端点公式是不同求积表示。
- Meijer G、M函数和余项使用55位工作精度；二维Gauss用浮点数和阶数递增比较，不是认证区间或形式证明。

Gamma与Meijer G均为已经求值的标准函数，结果不是未求值积分的改名。有限项只需原有M01、M11及ζ₂、ζ₃；另附经典Li₂/Li₃与完整Hlog坐标。

BO的H₀母核在δ=0附近解析，正实z下D角点的log矩可积；配合Gamma/有理/双曲函数的局部解析性，得到O(ε)余项依据。没有从浮点缩放单独推断余项阶数。

## 证据与迁移

原raw包、生产montopar及旧r001分类未改。`history/`保存被修改的现行文档/导航原字节，[migration.json](migration.json)登记前后语义；新分类为catalog-r002。五对原J反射依旧带必要的z幂，完成同一I提取后五对I均相等。

旧consolidation检查改为重放其原文档快照，并写独立`checks-replay.json`；不把旧CJ对象的检查转授给新I。初次分类构造器加了历史保护，避免覆盖r001；当前分类验证读取r002。

首次BO代数核验时，CAS在未知符号范围下引入了Abs(D)，导致残差未消去。改为在已声明D>0的积分域先按链式法则提出D^(2ε−2)，核对剩余多项式；未放宽恒等式条件。最终30项检查通过。

复验：

```powershell
D:/py/python.exe -B -X utf8 workflow/checks/self-energy-pi1-baselines/prepare.py
D:/py/python.exe -B -X utf8 workflow/checks/self-energy-pi1-baselines/compute.py
D:/py/python.exe -B -X utf8 workflow/checks/self-energy-pi1-baselines/numeric.py
D:/py/python.exe -B -X utf8 workflow/checks/self-energy-pi1-baselines/verify.py
```

`publish.py`是一次性迁移，拒绝重复执行。结果修订和prepared对象由生成器按字节保护，改变定义时应新增版本。未向外部目录写入，未提交Git，未运行新HyperInt或Lean；当前无可调用Wilson MCP，formal=not_run。
