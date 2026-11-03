# 4e-04第四、第五节可读性与行内数学修复

日期：2026-09-27。用户要求拆开两节的推导操作，并修复行内公式。正式入口为[主稿](../../../benchmarks/2ne/theory/generation/4e04-direct-mn.md)。

## 修复内容

- 第四节分为换元、合法拆项、写出排序原函数、求导降层、分部积分、Beta求值和边界回代。
- 第五节按“导数相同→确定常数→换元消去中间量→逐项代回”展开，三条转换单列为BR1、BR2、BR3。
- 全稿62处裸写或局部包裹的数学表达式改成完整数学定界；较长的式子独立成行。
- [r002](../../../benchmarks/2ne/theory/generation/4e04-direct-mn-r002.md)保留改写前原字节，r001继续保留原始研究稿。旧验证报告未修改。

## 检查

在仓库根目录运行：

    D:/py/python.exe -B -X utf8 workflow/checks/4e04-readable-sections/run_checks.py

[check.py](check.py)核对16条原有带标签结果保持完全相同（忽略排版空白），并检查正文新展开的导数、换元/Jacobian和消项步骤。

[render.mjs](render.mjs)使用已安装VS Code的Markdown/KaTeX：

- 分别记录实际解析出的行内与独立公式数，并检查每条行内公式都生成KaTeX HTML。
- 比较已标记行内公式数、实际数学token数与整个HTML的数学节点数。
- 在普通文本token中寻找遗漏的TeX命令、裸下标/上标和数学定界符；以旧r002为负控制，确认检测器会发现之前的问题。
- 检查所有公式、标签和链接，并生成[HTML预览](preview.html)。预览使用本机已安装的KaTeX CSS；这不是实时IDE截图。

[validation.json](validation.json)保存运行记录、旧主稿到r002的来源映射及哈希。sections.md和fix_inline.py是本次文字编辑输入，不是新的求解器。

未重算积分或改变全阶M/N结论，未提升formal状态，未提交Git。
