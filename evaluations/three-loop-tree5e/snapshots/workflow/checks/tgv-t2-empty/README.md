# T2[]：精确反射连接与计算核验

本轮完成5e-3g-05/06/07/08/20/21/31的完整ε表达及ε⁻⁶…ε⁰共49位置，40项非零，余项O(ε)，固定z>0。见[解表](../../../benchmarks/tgv/families/t2-empty/solutions.md)、[推导](../../../benchmarks/tgv/families/t2-empty/method.md)。结果formal=not_run。

## 实际路线

边反射(1 2)(3 4)、每条边完整排序反转、alpha1↔alpha3共同给出
05→04、06→03、07→02、08→01、20→19、21→18、31→17。
箭头左侧是T2，右侧是T1。逐图验证原S/P/pure、alpha测度、绝对Jacobian、原有序域与实际C。Gram矩阵及参考尺度不变，所以z不反演，J与I均在完整ε下相等。

本次复用T1已核验的HyperInt矩、Gamma/3F2/显式收敛和及展开余项，没有新跑Maple、改变输入C/J、注册Kernel根或改动外部仓库。默认47条预处理已运行，输出不扩容。

## 已执行检查

- 7组精确反射桥，包含源到prepared的原参数连接和唯一C/I规范。
- 4个负例：缺alpha置换、缺排序反转、P错号、改动边界。
- 49项Laurent传递、49项Hlog回代、49项z反演偶性、21条族内线性系数关系。
- 14组全新T2原七维J检查，在(ε,z)=(0.65,0.4)、(0.7,2.3)，每图每点8次scramble，每次2^17点，对照55位完整解析J。全部满足事先设定的8倍经验标准误差加小相对容差；最大相对差约2.95%，不是高精度或严格包围。
- 旧T1交付7图49系数/124哈希与mixed输入44图196项/390哈希复验通过；旧报告的not_computed仍描述输入阶段，不能覆盖当前解表。
- 当前运行、交付入口见[current-run](current-run.json)、[current-delivery](current-delivery.json)；依赖和快照逐字节核对，不将哈希等同于数学证明。

## 数值失败与修正

保留首次、第二次失败日志和复现点。原Dirichlet采样用1-inverse_beta计算补量时可能舍入到零；本地stable_qmc改为互补参数的inverse_beta恒等式直接算补量，并检查alpha正性与和为1。

第二次失败精确定位到21图一个Sobol坐标等于0的边界采样点（numeric-diagnosis.json），不是解析值矛盾。成功运行对30位Sobol网格统一平移半格2^-31，所有点落在开立方内；未裁剪奇异值、删点或挑选通过的scramble。失败时没有发布结果。

数学检查的首次代码误写5-edge，被域维数比较拒绝；随后改成实际所需的(1 2)(3 4)。第二次用字符串比较等价的指数格式而失败，改为精确符号比较。失败目录保留，最终成功报告重新生成。

## 复验和重算

仓库根目录：
~~~powershell
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t2-empty/verify_results.py
~~~
该命令复算反射桥/负例，核对49项及源/快照哈希；不依赖外部montopar或本机.wilson，不重跑数值积分。

新一轮原核检查：
~~~powershell
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t2-empty/check_math.py
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t2-empty/check_numeric.py
~~~
check_math创建新时间戳目录，check_numeric写到该目录。publish只允许创建尚不存在的r001；复算已有结果使用verify_results，不能覆盖不可变修订。长数值任务应给足执行时间或以隐藏后台进程写日志。

本次没有Lean/MCP正式整链证据；物理复z/i0延拓、一般全ε最小主积分数、其余30个mixed积分均未因此完成。
