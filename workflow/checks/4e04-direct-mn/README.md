# 4e-04直接M/N约化的检查

日期：2026-09-27。对应[未展开积分证明及任意阶公式](../../../benchmarks/2ne/theory/generation/4e04-direct-mn.md)。本轮从矩源的Cauchy方程和解析边界常数建立转换，不使用有限阶拟合。

在仓库根目录执行：

    D:/py/python.exe -B -X utf8 workflow/checks/4e04-direct-mn/run_checks.py

## 实际证据

- [algebra.py](algebra.py) / [algebra.json](algebra.json)：Cauchy分部积分的点态恒等式、完整Q重新组合、两条转换的常数π²/2和0、s→0合并极点的常数账本、七个N源展开控制。没有导入Hlog或候选矩阵。漏掉digamma边界的负控制被检出。
- [numeric.py](numeric.py) / [numeric.json](numeric.json)：45位mpmath直接积分普通可去的Y、Y_c、V、N生成函数，与已有₄F₃式核对。13个参数点覆盖ε=0.17、0.31、0.73、−0.08，以及0.08+0.03i；z为0.4、0.7、1、2.5。两条中间转换及最终Q都检查，最大Q缩放差约1.2×10⁻³⁸。没有数值PV积分、Hlog拟合或严格误差包围。
- [regression.py](regression.py) / [regression.json](regression.json)：先按新任意阶规则生成M/N表达，再读旧矩值与目标做精确回代。Q、Δ和完整I各五阶，共15项零残差；没有rref、PSLQ或候选搜索。
- [render.mjs](render.mjs) / [render.json](render.json)：使用已安装VS Code的Markdown/KaTeX检查正文、R1–R18标签及导航链接。不是编辑器截图检查。
- [validation.json](validation.json)记录实际命令、历史文件保持情况及所有来源/产物绑定；baseline保存本轮更新前的导航MD。

正文附录通过Euler/Gauss连接公式和Mellin核的Gamma展开固定边界常数。它们是解析推导的一部分，数值点仅作辅助排错。没有运行MB.m、HyperInt新积分、Lean或MCP形式整链；formal=not_run。

## 新结论的范围

4e-04的全阶M/N覆盖现在有未展开恒等式及任意阶公式，不再只是五阶经验规律。N权重R=L−log z由两个矩源插入的差直接产生；中间的减除零指标矩会消去。

这不证明任意独立双源f_mn都可约到M/N，也不证明M/N最小独立性、全部6e或3g的闭包。旧结果、catalog及历史报告不改写，原329位置覆盖不重复增加，未提交Git。
