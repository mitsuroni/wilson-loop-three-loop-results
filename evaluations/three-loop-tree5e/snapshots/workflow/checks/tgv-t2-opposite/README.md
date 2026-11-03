# T2[1]：三个原始有序积分的反射求解

对象为5e-3g-16/36/37，完整ε及ε⁻⁶…ε⁰共21位置、11非零，固定z>0、余项O(ε)。[解表](../../../benchmarks/tgv/families/t2-opposite/solutions.md)、[推导](../../../benchmarks/tgv/families/t2-opposite/method.md)。

## 实际检查

- 复用T2[]已检查的通用反射实现，在本轮三个原参数输入上重新执行S/P/pure、端点绑定、原依赖边界、迭代顺序、alpha单纯形、Jacobian及实际C检查。
- 从三个T1候选中逐一匹配，得到唯一对应16→32、36→12、37→11。反射交换13/24通道，z及参考尺度保持，完整J和I均相等。
- 五个负例拒绝：缺alpha置换、缺端点倒序、P错号（源与prepared同步变坏）、C多因子2、错误图配对。
- 21个Laurent系数与原T1表达精确相同，另逐项回代完整Hlog；I36+I37=I16的7条系数关系通过，37系数保持权重6+k。
- 原T1[1]验证器重新执行，21系数及184项源/快照哈希通过。T1的完整ε级数、9个余项及4个尾界案例属于复用证据，本轮没有重算HyperInt或重复这些数值案例。
- 新算6组T2原七维J：(.65,.4)、(.7,1.3)，每图每点8次scramble×2^17点，与已冻结的55位T1完整函数值比较。比较目标是解析J，不是旧QMC均值。最大相对差约0.664%，全部通过预设8倍经验SE加小相对容差；非严格数值包围。
- 使用开立方Sobol半格平移及互补Beta采样，无删点或裁剪。
- 默认47条预处理已执行；TGV专用输入、历史分类、旧结果及已有工作区不改写。

当前运行及交付见[current-run](current-run.json)、[current-delivery](current-delivery.json)。每个新r001保留完整函数定义、原T1修订哈希、逐图反射字典、系数与Hlog坐标。

## 复验

仓库根目录：
~~~powershell
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t2-opposite/verify_results.py
~~~
无需外部montopar、Maple或.wilson；重新核对反射/负例、系数与当前文件/快照，并核对固定的T1交付。

新轮数值复算：
~~~powershell
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t2-opposite/check_math.py
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t2-opposite/check_numeric.py
~~~
前者新建运行目录并固定来源交付；后者约47秒，工具超时不足时使用隐藏后台进程记录日志。publish拒绝覆盖已有r001；重算只复验，需更改结果时另建修订。

没有新增生成元，16只需既有M/周期，36与37有限项的额外方向分别为+(τDq4+10E)/192、−(τDq4+10E)/192。该表示范围继承T1的Q-Hlog坐标分析，不是一般全ε最小主积分定理。

formal=not_run；未更改Kernel/Lean、注册TGV根、证明物理复z/i0延拓或回写外部仓库。本轮未自动提交Git。
