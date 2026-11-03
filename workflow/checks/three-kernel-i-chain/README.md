# 三核自然生成稿：把辅助函数接回完整物理I

正文：[natural-beta-source.md](../../../benchmarks/2ne/theory/generation/natural-beta-source.md)。

本次按用户反馈改写说明，主线以I₆ₑ₋₂₀为左端，从六维原积分经矩形降维、W/Λ分解、端点幂代入，直到含已求值5F4及明确剩余项的完整等式。另列I₆ₑ₋₂₁=I₆ₑ₋₂₂和两个消元组合。

5F4只对应W的正正区间中的一个双Beta项。它不等于某张完整图，也不足以给出三图的完整全ε通解。本次没有新增积分求值、Laurent覆盖或形式证明。

复现本次检查：

    node workflow/checks/three-kernel-i-chain/render.mjs
    D:/py/python.exe -B -X utf8 workflow/checks/three-kernel-i-chain/check.py

- 核对旧物理归一化、两个组合的消元系数、W四角项的符号、正正区间四项拆分及物理参数代入；带错误组合因子、错误子块符号的负控制。
- 检查正文(I1)–(I6)左端均为完整I₆ₑ₋₂₀；VS Code Markdown/KaTeX离线渲染检查与本地链接检查见[render.json](render.json)，非实时IDE截图。
- [r001](../../../benchmarks/2ne/theory/generation/natural-beta-source-r001.md)逐字节保存原稿。旧Beta报告中原主稿哈希按这个路径验证，其余原文件保持；旧ε²报告的192个绑定保持。
- 原Beta符号/数值求值复用[已有检查](../natural-beta-source/README.md)，没有重新跑CAS积分。形式状态仍为formal=not_run。

[validation.json](validation.json)记录实际结果。update_docs.py为已执行的一次性改写脚本，不属于日常复验入口。未提交Git。
