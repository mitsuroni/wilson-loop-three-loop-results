# Π₁ P6计算记录

4e-Pi1-11完整I及ε⁻⁶…ε⁰已完成，二阶极点、7位置3非零。复用已有H₀及τ₀₂，没有新函数。当前Π₁累计7/13、49位置30非零。

[解表](../../../benchmarks/self-energy/pi1/families/p6-opposite-distinct/solutions.md) · [推导](../../../benchmarks/self-energy/pi1/workspaces/4e-Pi1-11/derivation.md)

实际运行 `D:/py/python.exe workflow/checks/self-energy-pi1-p6/compute.py`：15项精确检查；输入哈希、sym=2、原顺序和独立边界另作断言。完整Hlog系数直接复用源绑定的basis-functions.json，不重新运行HyperInt。

实际运行 `D:/py/python.exe workflow/checks/self-energy-pi1-p6/numeric.py`：21组控制，9个两核独立二维求积乘积（128/256阶），9个H₀矩数值，3组I与H₀余项缩放。矩的一维积分仅作数值核对，交付系数已保存实际Hlog表达。双精度求积未作认证区间；Meijer G和余项使用55位mpmath。

正实z、0<Re(eps)<1/2成立并向eps=0亚纯延拓；未验证一般复z。formal=not_run，无可调用Wilson MCP，未改Kernel/Lean。prepared及旧目录/结果均保持，新建catalog-r004。
