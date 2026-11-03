# 4e-04：持续追踪Q的主线改写

用户指出，原稿长时间推导辅助对象，难以知道Q到了哪一步。本轮将[主稿](../../../benchmarks/2ne/theory/generation/4e04-direct-mn.md)改成13个完整Q等式，每步注明操作及剩余对象。辅助求值移到附录，并指出返回主线的位置。

- [r003](../../../benchmarks/2ne/theory/generation/4e04-direct-mn-r003.md)保存改写前原字节；r001、r002继续保留。
- 正文新增明确的一维Q式、分部积分后的Q式，以及三个转换分别代入后的Q式。
- 数学结果与归一化保持，未重算积分或提高覆盖数，formal=not_run。

复现：

    D:/py/python.exe -B -X utf8 workflow/checks/4e04-q-chain/run_checks.py

[check.py](check.py)检查13个主线方框左侧均为Q、每步有状态说明，核对原物理结果和辅助恒等式的保留，并逐项检查分部积分、边界合并、三个转换和N生成函数的重标度。没有重新运行旧CAS或数值积分。

[render.mjs](render.mjs)继续检查实际行内数学token、完整KaTeX HTML及普通正文中的漏标公式，提供[HTML预览](preview.html)。这不是实时IDE截图。

[validation.json](validation.json)明确登记旧主稿到r003的来源映射，验证其原字节以及其余历史文件不变。旧检查文件与旧哈希未修改。
