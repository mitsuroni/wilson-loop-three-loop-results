# 双Π₁串联两图：正式I求解

[CH_A](../../../benchmarks/self-energy/pi1-chain/families/adjacent/solutions.md)和[CH_O](../../../benchmarks/self-energy/pi1-chain/families/opposite/solutions.md)均已计算完成：全ε、14个I位置（9非零）及42个颜色通道位置（27非零）。源/准备输入不变。

`compute.py`实际执行79项精确检查：原核/域、相邻端点、对边完整导数与四角边界、C4、旧Hlog矩、完整G³U至e⁴、14个主I系数及42个颜色卷积系数。Γ余项按固定正实z的局部解析性说明；没有新增M函数或HyperInt运行。

`numeric.py`实际完成48组控制/矩评估：9个完整对边原核扇区求积、3个相邻端点求积、6个M21/M03定义积分、6次实际prepared C与已求J归一化、24组主I/颜色通道余项减半。原核Gauss为双精度非认证积分；其余使用65位mpmath。

首次余项范围1/64…1/512中，02的nf²通道在z=0.5的最后比率为0.46408，未满足距0.5小于0.03；公式与阈值均未改，缩小到1/128…1/2048后全部通过。[首次诊断](first-remainder-range.json)保留，不伪造首次通过。

[历史映射](history-manifest.json)保存准备阶段live Markdown原字节，原catalog-r001和准备记录未改。当前catalog-r002登记求解结果；准备阶段verify.py含“无results”的旧阶段断言，不作为当前交付复验器。当前复验命令为：

```powershell
D:/py/python.exe -B workflow/checks/self-energy-pi1-chain-solve/verify.py
```

无可调用Wilson MCP，未改Kernel/Lean、生产montopar或默认47输入；formal=not_run，一般复运动学延拓未验证。
