# 积分表示框架评审的实际检查

主文：[representation-review.md](../../../benchmarks/2ne/families/same-channel-two/representation-review.md)。本轮不重算积分、不改写数学结果，审视母族定义、冗余、导数闭合与最终展示方案。

[current-review.json](current-review.json)指向实际输出：

- review.json：现有30个结构登记项、8个零项；权重4/5/6自身秩1/3/8；与旧单对边τ空间合并后新增方向0/0/1。
- generator-screen.json：筛选有明确结构指标的新方向；eta^[2],id_{0,1;1}可补足权重6空间，三个Σ可精确表示。选择只是评审候选，未改写旧结果。
- review.py、generator_screen.py：本次执行代码快照。
- provenance.json：实际读取的结构定义、旧基函数与矩数据、既有数值级数实验及代码指纹。
- z导数的链式法则残差为零，给出新增zts/D插入；不由此宣称已建立或不存在有限闭合系统。

首次对一般符号幂直接simplify未能把导数残差归零，未产生成功报告。最终在所用正实域D_i>0上采用明确的exp(logkernel)形式，再用链式求导和有理cancel验证；导数目标式未修改。此处CAS操作不授予形式证明状态。

复现：

~~~powershell
D:/py/python.exe workflow/checks/representation-review/review.py
D:/py/python.exe workflow/checks/representation-review/generator_screen.py
~~~

review.py建立新的时间戳run。两步只读取已有求值结果，输出精确有理Hlog/ζ坐标的秩和等式。坐标秩不是一般超越常数独立性定理。formal: not_run。

本轮另查阅Henn 1304.1806及Panzer 1403.3385的原始论文页面，用于比较微分方程和HyperInt的职责；本地积分关系及秩结论来自本项目定义与实际计算。
