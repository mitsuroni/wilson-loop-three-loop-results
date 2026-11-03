# Π₁主推导合并

2026-09-23。按用户要求，现行推导统一为[theory.md](../../../benchmarks/self-energy/pi1/theory.md)，只有原动量式定义的I，完整张量直接给出一个C和一个J。主稿不再引入分项I、另一种同名I或正文版本比较。

根README、Pi1 README和设计入口只指向这份主稿；历史修订与审计统一由[历史索引](../../../benchmarks/self-energy/pi1/history.md)进入。数学文件及其证据路径保持原字节，没有删除、移走或重写旧公式。

旧Pi1 README既是导航又被初次发布清单记录了哈希。本次先把其原字节保存到[README快照](history/README-r001.md)，再更新当前导航。原`publication.json`未改；原验证器仅对这一个导航项改为检查原快照，其他数学与数据路径继续原地核验。[导航记录](navigation.json)。这不会把旧收据转授给新的文档或改变旧论文比较结论。

[check.py](check.py)、[checks.json](checks.json)核对主稿的唯一I、完整核与共同C、相邻闭式，以及全部13图的SE/裸线端点和原积分次序。数学内容复用已经核对的推导；本次不重跑生产回归或其余积分。

```powershell
D:/py/python.exe -B -X utf8 workflow/checks/self-energy-pi1-consolidation/check.py
D:/py/python.exe -B -X utf8 workflow/checks/self-energy-pi1/verify.py
```

[公式渲染](render.json)检查现行主稿。旧生产代码、两份JSON输入、既有结果和默认47条prepared均保持；未写入外部目录，formal=not_run。
