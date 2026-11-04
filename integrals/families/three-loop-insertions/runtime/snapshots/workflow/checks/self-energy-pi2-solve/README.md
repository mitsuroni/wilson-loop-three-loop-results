# Π₂两图完整I与逐颜色有限项

2026-09-24。用户要求将算符中的pAA/pAf/pFf和Gamma比值也完整展开。本批完成2/2图、6/6颜色系数：42个ε⁻⁶…ε⁰位置，25非零；两图共享几何I保持理论原定义。

入口：[相邻解表](../../../benchmarks/self-energy/pi2/families/adjacent/solutions.md)、[对边解表](../../../benchmarks/self-energy/pi2/families/opposite/solutions.md)、[算符展开](../../../benchmarks/self-energy/pi2/operator-expansion.md)。

## 实际工作

- prepare.py：50项检查，原C/J来源哈希、全颜色固定基底、Ff第二坐标、一次三圈Euler。默认预处理python -m preprocessing另行执行，仍为47根；本批使用专用prepared。
- test_prepare.py：3项unittest，第三项含7种错误输入子例；核验重复Euler、漏Ff、篡改算符/Gamma/域/重数的拒绝行为。
- compute.py：166项精确检查。来源连接、串联J/域/重数匹配、既有Hlog依赖、R到五阶、phi到四阶及42个系数的双路线吻合。
- numeric.py：33组80位控制/矩评估：6个既有M定义矩、3个算符O(eps^5)余项、18组逐图逐颜色O(eps)余项、6组实际prepared C/J完整颜色重构。
- publish.py：从同一精确对象生成显式解表；每图工作区只含derivation.md。没有手抄系数或建立新的数学IR。

共享几何I复用已求串联；新工作是准备Π₂根、核对等价桥和完整算符卷积。算符有表观三阶极点，故R=G10必须保留到ε⁵再合并，保证相邻有限项的ζ₅不丢失。最终表没有未展开p、Gamma或EulerGamma。

## 复算

先执行prepare.py，再compute.py、numeric.py；输入/输出修订若变化，freeze会拒绝覆盖。重复publish会保持历史README，更新展示页。公式采用现有离线VS Code Markdown/KaTeX检查；verify.py连接全部来源、结果、脚本和渲染记录。

实际日志：[preparation-checks.json](preparation-checks.json)、[exact-checks.json](exact-checks.json)、[numeric-checks.json](numeric-checks.json)、[tests.log](tests.log)、[delivery.json](delivery.json)。

原montopar报告绑定的README在[历史快照](history/benchmarks/self-energy/pi2/README.md)保留；原源、理论、C/J以及既有串联结果没有改变。旧报告不覆盖，依赖由当前verify检查。

本次没有新的原核数值求积：Pi2与已验证串联核、域、重数精确相同，复用该证据。新数值比较为非认证高精度支持；不是形式证明。没有可调用Wilson MCP，未改Kernel/Lean，formal=not_run。一般复z、重整化及完整物理图求和不在本批范围。
