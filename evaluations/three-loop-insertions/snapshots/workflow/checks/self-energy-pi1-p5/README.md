# Π₁ P5与十三图计算收尾

07/08完整I相等、简单极点；全ε为所有参数已积分的混合指数收敛级数加Meijer G。14位置4非零；Pi1累计13/13、91位置50非零。

[解表](../../../benchmarks/self-energy/pi1/families/p5-opposite-same/solutions.md)。`compute.py`实际完成24项精确检查：源核、反射、完整导数恒等式、八角边界、C4前因子、4个终止级数与原迭代多项式积分、7个Laurent系数的旧Hlog约化。原M及双指标rho_id_0_0_1文件绑定哈希，无新HyperInt运行。

`numeric.py`实际完成9组控制：3组两个源核与降维式对照、3组rho0/rho1与独立有序积分对照、3组有限项余项。55位mpmath计算Meijer G和单个−1字母的Hlog；完整核/有序双核用三维Gauss求积，记录64/128及96/192阶差异。分辨率差不是认证误差，混合无穷级数没有高精度求和或认证尾界。

原source/prepared和旧修订保持。当前catalog-r007；无可调用Wilson MCP，未修改Kernel/Lean，formal=not_run，一般复z物理延拓未验证。
