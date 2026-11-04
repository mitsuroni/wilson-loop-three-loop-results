# Π₁ P4计算记录

04/09已得到相同完整I、14系数位置（6非零），实际二阶极点。复用已求值H₀/H₂及τ/η，无新函数。

[解表](../../../benchmarks/self-energy/pi1/families/p4-opposite-bare/solutions.md)。运行 `compute.py` 完成16精确检查及输入/常数Hlog断言；运行 `numeric.py` 完成12数值控制：6次两图降维原核与Meijer G对照、3次04剩余原三角域128/256阶求积、3组独立H矩/余项缩放。mpmath为50位；Gauss为非认证浮点求积。完整Hlog坐标来自源绑定的旧basis-functions.json，未重新运行HyperInt。

起始正实z、0<Re(eps)<1/2，未验证一般复z。无可调用Wilson MCP，formal=not_run。输入、旧结果和旧目录不变，当前catalog-r005。辅助脚本均位于本检查目录。
