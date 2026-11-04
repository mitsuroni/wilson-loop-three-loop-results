# Π₁ P1/P2 求解与核验

本批完成4e-Pi1-01、02、05、06：四个完整I，28个Laurent位置（20非零）。输入及唯一归一化保持不变；当前Pi1累计6/13、42位置（27非零）。

- [P1完整解表](../../../benchmarks/self-energy/pi1/families/p1-adjacent-double/solutions.md)：五阶极点，全阶与z无关。
- [P2完整解表](../../../benchmarks/self-energy/pi1/families/p2-adjacent-middle/solutions.md)：三阶极点，保留cosh尺度比。
- `compute.py`：69项精确检查，实际源核/域/顺序、完整C4组装及28个系数。
- `numeric.py`：46组控制，36个实际四维有序源核求积点（每点8/16阶）、4次prepared C重构（nf=0,4）、6组余项减半。
- `publish.py`：结果文档、catalog-r003及进度入口；旧catalog-r002未覆盖。

运行命令：

```powershell
D:/py/python.exe workflow/checks/self-energy-pi1-adjacent/compute.py
D:/py/python.exe workflow/checks/self-energy-pi1-adjacent/numeric.py
D:/py/python.exe workflow/checks/self-energy-pi1-adjacent/verify.py
```

数值四维求积直接调用prepared integrand：P1用两个三角域的正变量换元；P2对共享三角域r+t>1按r<1/2、t<1/2、中央方形分成三块，端点权重纳入Gauss–Jacobi。不是将闭式重新代入同一路径的比较。浮点求积无认证区间；闭式及余项使用65位mpmath。

完整公式先在正实z、0<Re(eps)<1/2建立，再亚纯延拓至eps=0附近；未验证一般复z的物理延拓。没有可调用的Wilson MCP，未修改Kernel/Lean，formal=not_run。原始source及prepared内容哈希保持；旧报告仅证明原批次。
