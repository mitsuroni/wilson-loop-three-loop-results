# 4e-04通解主稿的讲解顺序与记号检查

日期：2026-09-27。用户要求参考GF cal中2e-02的推导方式，重写[4e-04主稿](../../../benchmarks/2ne/theory/generation/4e04-direct-mn.md)。本轮只改变讲解组织与等价记号。

## 改写范围

- 正文按原J/C-I约定、一圈M生成函数、四维降维、排序Beta求值、统一矩核、N生成函数、χ边界、完整母函数、展开顺序组织。
- 新记号V_M等于此前裸单核H；V_N等于此前N生成积分除以ε²。最终I公式完全等价，所有标度/Gamma边界保留。
- 参考GF cal文档的结构，没有采用其中手工“丢弃γ_E”、编号混乱或未证明的高j降次声明。当前Z=e^(−γ_E ε)Γ(1−ε)及数学结果保持。
- 旧稿原字节保存在[4e04-direct-mn-r001.md](../../../benchmarks/2ne/theory/generation/4e04-direct-mn-r001.md)。原验证报告绑定该版本；没有修改原报告的哈希或状态来迎合新稿。

## 实际检查

运行：

    D:/py/python.exe -B -X utf8 workflow/checks/4e04-derivation-edit/run_checks.py

[check.py](check.py) / [checks.json](checks.json)核对当前I归一化、V_N的ε²重标度、一圈双曲分解、二维可积拆项、Beta边界合并及前两阶例子；删除边界因子2的负控制被检出。

[render.mjs](render.mjs) / [render.json](render.json)使用已安装VS Code的Markdown/KaTeX检查主稿及导航的数学和链接。原13点数值及15个系数回归只按已绑定报告复用，未重新运行CAS或数值积分。

[validation.json](validation.json)显式登记“旧路径→r001归档”的来源映射并检查其原字节，其他历史输入、结果和证据保持。input.json记录参考稿及旧主稿的来源哈希，baseline保存导航改写前字节。main-body.md仅是本次编排的正文输入，正式入口仍是上述主稿。

未改变M/N覆盖结论、展开阶数、原329位置覆盖或formal状态；未提交Git。
