# 36/39 独立端点生成源：本轮实际核验

2026-10-16。[连续推导与全部公式](../../../benchmarks/2ne/theory/generation/single-opposite-free-endpoints-20261016/derivation.md)。

本轮从原双端点幂矩得到独立 λ/a/b 的非齐次 Gauss 方程、两支边界常数和显式 ₂F₁/₃F₂/₄F₃ 解，再把原物理镜像提前组织为三指标矩。没有修改输入、旧结果、旧证据、PDF 或 Wilson/Lean；未创建新 goal、未提交 Git。formal=not_run。

## 实际通过

| 记录 | 实际检查 |
|---|---|
| [verification.json](verification.json) | 54 精确恒等式；13 错误式控制；12 原二维积分与闭式比较；6 方程残差；6 遗漏边界数值控制；6 一维 M/N 与指标矩比较；3 Taylor 截断缩放 |
| [wolfbook-result.json](wolfbook-result.json) | 原生 Wolfram 14.3 实际运行 11 项通过，包含 4 个错误源控制；无 kernel messages |
| [coefficient-verification.json](coefficient-verification.json) | z=0.4、0.7 的完整 I，ε⁻² 至 ε² 共 10 项对照旧精确 Hlog 向量，最大差 6.95e-11 |
| [euler-verification.json](euler-verification.json) | 非等上参数 Euler 与 Mellin 闭式各 4 项；小 z 两支边界 2 项 |

原二维核与闭式的最大绝对差 5.96e-14。数值结果使用节点精化或高精度比较，不是严格区间证书。高阶系数采用已冻结的 canonical Hlog 向量；多重周期沿用原报告中明确绑定的 zeta_5_3 数值和 Maple 约定，未重算该常数。指标矩仍是候选族，本轮没有精确求出四个权重六矩的 Hlog 字典。

## 重放

在 D:/wilson-agent：

    D:/py/python.exe -B -X utf8 workflow/checks/single-opposite-free-endpoints-20261016/verify.py
    D:/py/python.exe -B -X utf8 workflow/checks/single-opposite-free-endpoints-20261016/check_coefficients.py
    D:/py/python.exe -B -X utf8 workflow/checks/single-opposite-free-endpoints-20261016/check_euler_extension.py

[generators.py](generators.py) 实现全部生成式及直接原核求积。二维网格采用旧 numeric_physical.py 的三次端点映射，单独复制有限网格代码以避免旧模块加载副作用；来源哈希记录在验证文件中。

[algebra.wl](algebra.wl) 经 Wolfbook 的 evaluateExpression 在现有 K1 内执行，Module 局部化全部变量；没有新建、重启或中断其他 kernel。实际 operation ID 为 6c2e6285-10e5-4bc1-b869-39ce42c5d2a2，完整工具返回保存在 wolfbook-result.json。

## 本轮失败与修正

1. 首次 Gauss 规范变换检查停在未通分的有理式；factor 单独没有化为零。改为先 together/cancel 再 factor 后得到严格零，未改变目标恒等式。
2. Euler 数值控制的 B=−0.8 样本存在可积但很尖的端点。45 位直接积分误差为 1.13e-10，未达到原 1e-18 阈值。使用明确换元 t=x^10 后保持同一阈值通过；没有放宽容差或拟合修补。
3. 路径探索中的 PowerShell 通配符和旧目录猜测失败属于只读定位问题；最终脚本使用已找到的实际文件路径。它们没有生成数学成功记录。

## 留下的具体问题

优先为 E200、E110、E020、E011 建立权重六的精确 Hlog/MN 字典，检查对角线之和及离对角线之差。权重九的物理规则已经由固定生成式写在推导式(29)，尚未完成独立 canonical Hlog 检验。全族最小性、跨图全阶闭合和一般复 z 分支没有在本轮证明。
