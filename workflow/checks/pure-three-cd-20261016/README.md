# C/D 原排序积分的直接生成与系数检查

日期：2026-10-16。研究对象为 C=R321 与 D=R231=R312 的裸几何积分。完整颜色中它们非零，当前 a1 投影为零；它们可作为反解生产 A/B 的辅助坐标。[数学入口](../../../benchmarks/2ne/theory/generation/pure-three-cd-20261016/derivation.md)。

## 已完成的数学结果

由各自六维有序域，C 降为核乘两个互补矩形，D 降为两个线性幂差乘一个互补矩形。矩形原函数显式为 J_e(x)=−e*x*3F2(1,1,1−e;2,2;x)。四扇区、shuffle 和运动字母 Chen 递归给出任意指定 q 的有限 GPL 词和，字母为 0、1、−z、−1/z，自变量 1。该公式可生成各阶；少数一维母函数的紧凑全 epsilon 闭合仍待研究。

`r001` 完成 epsilon 零至二次；`r002` 完成零至三次。后者直接源词数 Cn2..5 为 128、574、2052、6638，Dn3..6 为 164、688、2412、7498。Maple 单独求值全部新源、卷积已知 rho 因子并写出后，才读取旧 A/B/H/Q 表达式作比较。Python 预先序列化比较文件，比较数据不进入 `direct_words.py` 或新源构造。r002 的 12 个精确回归差全部为零；这是新表示的独立求值与已知阶比较，不宣称预测了此前未知的物理阶。

最终 C 的四阶 Hlog/周期坐标数为 4、24、84、234，D 为 15、57、153、357。完整公式见 [C](../../../benchmarks/2ne/theory/generation/pure-three-cd-20261016/C/coefficients.md)、[D](../../../benchmarks/2ne/theory/generation/pure-three-cd-20261016/D/coefficients.md)及[机器结果](../../../benchmarks/2ne/theory/generation/pure-three-cd-20261016/coefficients.json)。C0=M01³/18，C1 的原 M/N 紧凑表达已按实际 M 定义重新乘积核对。

## 独立检查

- `geometry-verification.json`：epsilon=1,2,3、一般符号 z，原六维有序多项式与二维矩形式逐一一致，错误互补矩形负控制被拒绝。共同归一化因子在这些多项式测试中整体剥离。
- `geometry-wolfbook.json`：实际 Wolfbook/MCP 独立核对 epsilon=2 的两条原域等式及负控制，成功完成，无内核消息。
- `numeric-verification-r002.json`：epsilon=0.23、0.37，z=0.7、1.4，共八个完整非整数原积分比较。C 最大观测差 3.83e−12，D 为 1.54e−12；八个错误矩形负控制拒绝。经验精化，不是区间认证。
- `r002/status.json`、`assembly-verification.json`：实际 HyperInt 日志、源码和输入输出绑定，直接系数权重、旧 Q 回归和原 A/B 关系均通过。
- `r001`／`r002` 的 `mn-comparison.json`、`r001/extension-comparison.json`、两轮 `existing-source-comparison.json`：保留所有有理精确重构与非成员对偶见证。C2 在独立 Euler/Beta/Mellin 空间内，D2、C3、D3 在其外；加入独立双 J 后全部已计算系数可表示。这只涉及明确声明的有限空间。
- [最终重放](verification.json)：读取原报告，重新回代拟合、逐列核对对偶见证、恢复 A/B、检查发布文本和旧 r004 的 35 个绑定文件。

## 复现与历史

已完成的 r001/r002 不覆盖。需要重新积分时使用新目录名，依次运行 `direct_words.py --order 3 --run NEW`、`run_direct.py NEW`；运行前须按本机路径核对 Maple/HyperInt。发布脚本拒绝覆盖当前发布文件。当前结果的无积分只读核验脚本为 `verify.py`；其报告写入一次，重新审计应另存修订版。

首次数值试验中的端点 NaN、随后一个 log1p 端点警告及修正保存在 `numeric-initial-attempt.txt`；第一次完整报告及所用源码另存，接受版本为 r002。没有把失败尝试计为成功，也没有改写它们的历史哈希。

本阶段未修改原 C/J、prepared、旧系数、旧审计或 Wilson/Lean；formal=not_run。全 epsilon 紧凑母函数、一般复 z 分支和最小主积分数不在已完成结论内。
