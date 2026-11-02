# 3e-3g+Pi1：树级三胶子腿上的单个一圈自能

2026-09-25：新增三圈12个C4轨道，三圈总源图数137→149。四个3e-3g骨架各标记三条腿；同边的两个有序槽位不可合并。镜像只作配对记录，不代表积分相等。

此族与3e-Gamma1（中央三点1PI一圈修正）分开。定义导出四段传播子、一个树TGV、一个Π₁、同一插入腿动量、三个相位、动量δ和外到内有序域。

新族直接采用wilson-self-energy-1pi-insertion/2.0：ΠL=-i(k²g-kk)PL，D0ΠD0连接外系数+1。PL保留内部符号/颜色，不能另加自能外负号。旧2e-Pi1/4e-Pi1/2e-Pi2源仍是legacy，D:/wilson-agent的版本化导入边界已统一它们；本次没有改变旧源及montopar结果。不要把新族再送入旧的翻号适配器。

用组合生成器运行：

```powershell
python -m qgraf_pipeline.cli 3
```

新增PDF位于outputs/3loop/pdfs/3e-3g+Pi1.pdf。所有旧轨道的定义、LaTeX、TikZ及镜像编号保持。

Stage 1已完成；montopar C/J与I尚未实现。当前montopar对新族应拒绝，不能误译为2e/4e-Pi1。完整理论和三圈进度见D:/wilson-agent/benchmarks/self-energy/tgv-pi1/theory.md及../three-loop-status.md。
