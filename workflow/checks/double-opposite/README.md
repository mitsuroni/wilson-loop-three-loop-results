# [1,1] 族计算与复现

范围为当前 prepared 0.3 的 4e-06、6e-32、6e-33、6e-34。目标是完整 epsilon 的 I，以及从 epsilon^-6 到有限项的系数；数学依据见 [method.md](../../../benchmarks/2ne/families/double-opposite/method.md)。

```powershell
D:/py/python.exe -m preprocessing
D:/py/python.exe -m unittest discover -s workflow/checks/double-opposite -p test_*.py
D:/py/python.exe workflow/checks/double-opposite/solve_double.py
D:/py/python.exe workflow/checks/double-opposite/publish_tables.py
D:/py/python.exe workflow/checks/double-opposite/verify_results.py
```

- reduce_double.py：读取真实原参数根，识别 D/tilde D、端点次序、邻边权重与 C；输出两个已知核的乘积。
- check_double.py：核对复用源、原顺序整数参数积分、反射后的核和原域数值积分；精确生成 M/tau 表示。
- test_double.py：输入/条件负例、核反射、三圈相等关系、纸面系数及表示回归。
- solve_double.py：逐图计算、检查并新增结果修订；不覆盖旧结果，不修改 Kernel/Lean。
- publish_tables.py、verify_results.py：从 JSON 生成阅读总表，核对公式、依赖和文件。

复用上一批 H0/H1 q=0,1,2 和 M01/M11/M21 的实际 HyperInt 求值及独立 Euler/Cauchy 检查，绑定到具体结果、源文件和哈希；本批不再运行相同的 HyperInt 积分。必要运行文件和代码快照在 runs/，求解工作区仍只有 Markdown。

初次表示导出尝试复用单对边线性 eta 转换，遇到奇部分乘积产生的形式 ell^-2，因检查拒绝而停止。随后改用已证的三项 M→tau 多项式字典并精确反代；没有改动前一族的转换器，也没有将该形式分母保留进答案。原始尝试和已完成的单项保留，完整批次须有 summary.json。

所有完成态结果分别记录 analytic、checks、formal；formal=not_run。Sobol 对照是统计一致性，不是严格数值区间。原 J 的整数 epsilon 检查不评估那些点上的 Gamma 前因子 I。物理复 z/i0 和更高正幂不在本批已完成范围。

当前实际结果与检查汇总：[delivery-summary.json](delivery-summary.json)。
