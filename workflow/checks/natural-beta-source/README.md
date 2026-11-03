# 自然双Beta源子块

正文：[自然生成函数](../../../benchmarks/2ne/theory/generation/natural-beta-source.md)。本轮从实际三核端点差分中的双Beta项直接导出三参数5F4与源方程，没有读取M/N或旧目标系数。

运行：

```powershell
D:/py/python.exe -B -X utf8 workflow/checks/natural-beta-source/check.py
```

[checks.json](checks.json)记录8个符号超几何系数核对、8个方程系数核对、2个终止多项式控制、3个独立Jacobi求积对照，以及遗漏原函数因子的负控制。数值为经验误差检查，非严格区间包围。

该结果仅求完一个真实子块；未证明整个三核已化成5F4，也未建立全族最小性。formal=not_run，未新增Laurent覆盖，未提交Git。
