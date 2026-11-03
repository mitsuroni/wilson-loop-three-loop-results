# T1[1]：两通道三个积分的计算交付

完成T1[1]₂₄的5e-3g-11/12和T1[1]₁₃的32，完整ε为已逐项积完的显式收敛级数；ε⁻⁶…ε⁰共21位置、11非零，余项O(ε)，固定z>0。见[解表](../../../benchmarks/tgv/families/t1-opposite/solutions.md)、[原核推导](../../../benchmarks/tgv/families/t1-opposite/method.md)。

## 本次实际运行

- 默认47条预处理运行，TGV仍是专用输入，未扩充Kernel根。
- 原三图的S/P/pure、实际排序与源到prepared连接、alpha换元、Jacobian及唯一C/I规范检查通过。
- 联合xy原函数、Beta分部积分恒等式、6个有序分区/缩放多项式探针通过。
- 完整恒等式J11+J12=J32，且I11+I12=I32；7条系数关系精确通过。
- 新鲜HyperInt：4个有限b矩H0…H3，4个T系数，4个K系数（一次项为0），并作含z^(±ε/2)的物理偶投影。所有发散检查开启，周期表显式加载，运行约3秒Maple CPU时间；不是端到端延迟承诺。
- 四个T系数另从既有2ne的H0/H2完整源矩精确重构。
- 21项与既有统一基底精确回代，无新增生成元。11各系数均匀权重6+k；固定Q-Hlog坐标下，权重6 M/N空间秩12，加入本批后13，额外仅为旧组合τDq4+10E。
- 9个65位降维完整I与Hlog截断余项样本，ε=.03/.015/.0075、z=.4，符合O(ε)。
- 6组原七维J：三图各在(.65,.4)、(.7,1.3)，8次scramble×2^17点，全部通过事先给定的8倍经验SE加小相对容差。最大相对差约1.69%；不是严格包围。沿用上一轮的开立方Sobol半格平移和稳定互补Beta采样，未删点。
- 完整T/K显式级数在两个参数点、n≤32的4项截断与独立一维求积一致于解析尾界内，8→16→32误差递减；界非最紧，数值未作区间算术。
- 复用结果另执行T1[]验证器（49系数/124哈希）、统一权重6验证器（329系数、92转换、161哈希）及mixed输入验证器（44图196项/390哈希）。不把旧HyperInt重算声称为本轮新算。

## 记录和复算

[current-run.json](current-run.json)指向本轮原生Maple文件、Hlog文本、系数及数值报告；[current-delivery.json](current-delivery.json)指向结果和依赖快照。原始输入、旧2ne成果、T1[]/T2[]结果与旧分类文档均保持原字节。当前进度单独更新。

仓库根目录：
~~~powershell
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t1-opposite/verify_results.py
~~~
离线复验原核连接、精确系数/基底和源/快照哈希；不要求外部montopar或.wilson。

重算顺序：
~~~powershell
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t1-opposite/launch.py
# 等待新run的maple.log出现T1_OPPOSITE_COMPLETE且没有Error
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t1-opposite/formulas.py
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t1-opposite/check_math.py
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t1-opposite/check_numeric.py
~~~
launch创建新目录并在后台运行Maple；数值检查约123秒，应给足时间或用隐藏后台进程保存日志。publish只允许创建不存在的r001，不得覆盖已有修订。

首轮Maple误将已为HlogRegInf的列表再次convert而失败，原目录保留；改为按列表/表达式分流后重新启动成功。SymPy第一次未自动消去Beta导数差中的幂，改查除去公共非零密度后的有理恒等式才通过；未将未化简的表达式当作成功。

结果formal=not_run。未更改Kernel/Lean、外部源码或已绑定旧结果，未证明物理复z/i0延拓或全ε主积分最小性，未自动提交Git。

首次离线结果复验把基底符号E误按SymPy默认自然常数解析，因而拒绝。现已使用结果内声明的符号表解析，旧交付目录与失败日志保留；数学结果及其Hlog坐标不变，另冻结新交付。
