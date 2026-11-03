# Π₁ P3计算记录

03/10完整I及14系数位置完成，四阶极点、10非零；累计11/13、77位置46非零。

[解表](../../../benchmarks/self-energy/pi1/families/p3-opposite-se/solutions.md)。`compute.py`实际运行28项精确检查，包含完整核导数恒等式、Beta边界、C4组装、一般M矩公式与旧τ/η的12次Hlog匹配、7个系数的crossing；另外校验源哈希及原域。完整M表达和Hlog向量随结果冻结，旧M文件逐项绑定哈希，未重新运行HyperInt。

`numeric.py`实际运行9组控制：03/10两个带权完整SE核在3个(e,z)点作128/256阶扇区求积，另在3个z上以原M的一维定义积分独立检查有限项O(eps)余项缩放。mpmath为55位；源求积为非认证float64。

完整输入/旧结果不变，当前catalog-r006。一般复z物理延拓未验证，未修改Kernel/Lean，无可调用Wilson MCP，formal=not_run。
