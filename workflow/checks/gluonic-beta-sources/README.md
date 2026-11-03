# 纯胶子共同Beta源：实际计算与边界

[共同源推导](../../../benchmarks/theory/gluonic-beta-generators.md) · [逐图清单](../../../benchmarks/theory/gluonic-extra-basis-inventory.md)。本次是积分研究，不改Kernel/Lean。

## 已运行

- `inventory.py`读取47个纯2ne、44个5e-3g、6个非零4e-4g及7个已完成双TGV图的真实系数；30图出现tau/E或其导数。低阶tau别名先使用既有统一转换，不按basis声明误计。
- `algebra.py`：16项logit、乘积变量、谱密度指数、Cauchy根及Gauss方程连接检查。共同收敛域与原物理规范在主稿单独说明，局部代数成功不是形式整链证明。
- `runs/free-beta-w6`：独立Beta指数的权重六源矩；新鲜HyperInt会话。
- `runs/free-beta-w8`：同一源独立重算并扩至权重八；80个单侧矩及镜像。
- `runs/asymmetric-beta-w8`：另补86个含不对称Beta对数的单侧矩及镜像。三个CAS会话均保持发散检测，实际运行日志、脚本哈希及时间单独保存。
- `compare.py`：对称子族的权重六、七、八精确空间；保留最初`comparison.json`。旧VB的对角参数嵌入分别检查13、19、25个已有矩。
- `compare_asymmetric.py`：全部独立指数的完整同权重镜像偶多项式。候选秩14、28、52；加入所测目标后14、28、53。权重八仍缺纯6e-11、12，保留精确非成员证书。这里的“完整”仅指所定义参数矩空间，不指所有可能解析函数。
- `numeric.py`：30项独立控制，含Euler/Gauss、旧VB切片、FGV密度、2ne原二维积分与新谱积分、有序端点原函数、两个关键新矩的直接双重求积及tau/E字典。矩字典最大缩放差约7.2e-12，节点加密最大约3.1e-10；这是经验比较，非严格包围。
- `publish.py`：只从实际系数生成清单及完整字典，`dictionary.json`记录30图的额外部分替换，不擅自删除原M/N、周期或Omega。
- `verify.py`：从CAS输出重建每个矩，回代成员表达、检查非成员泛函、验证低权重独立重算与四胶子Theta tau的完整有理z向量；另记录两个失败目标的共同商空间方向。

## 完整镜像补测

随后补齐G17的偶数L幂，即cosh减1乘两端权重差；其零Beta参数层恰为VN型。`runs/even-beta-w8`新增第四个CAS会话，46个对称、40个不对称单侧矩及镜像；`compare_even.py`给出`comparison-complete.json`。补测后的秩仍为14、28、52，权重八联合秩仍53；故此前缺口并非仅由sinh限制导致。

`numeric_complete.py`另有14项控制，直接检验新矩的镜像符号、节点加密和N型零参数关系；最大缩放差约1.1e-10，细化差约9.8e-9。合计44项数值控制。`validation.json`保留第一阶段863项复验，`verify_complete.py`复用并检查其原数据绑定，再验证所有新偶数矩与最终空间证书。`render-complete.json`绑定最终文档，初稿渲染记录保持。

## 复验

    node workflow/checks/gluonic-beta-sources/render_complete.mjs D:/wilson-agent
    D:/py/python.exe -B -X utf8 workflow/checks/gluonic-beta-sources/verify_complete.py

生成/CAS脚本会写其运行目录，后续修改请另建修订；旧VB与所有物理结果r001保持。三个会话没有关闭发散检查或丢弃发散项。库存解析曾将zeta误作普通符号，已明确为标准函数；局部对数恒等式按正实z假设化简。没有把对称或不对称空间失败改标为通过。

本轮没有新的完整物理图QMC、没有Wilson MCP入口、没有新Lean证明，formal=not_run。结论只包括明示的源恒等式、固定矩、有限阶转换和所测高阶范围；不计新的全epsilon整图完成数或最小基底证明。
