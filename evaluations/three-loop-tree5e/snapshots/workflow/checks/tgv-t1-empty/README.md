# T1[] 七图求解与核验

[解表](../../../benchmarks/tgv/families/t1-empty/solutions.md) · [统一推导](../../../benchmarks/tgv/families/t1-empty/method.md) · [结果索引](../../../benchmarks/tgv/families/t1-empty/results-r001.json)。

本次新算5e-3g-01/02/03/04/17/18/19的完整ε表示和ε⁻⁶…ε⁰共49个系数。01–04六阶极点、17/18四阶、19三阶。精确有I02=2I03、I01=I03+I04、I17=I18+I19。系数只含π、ζ3、ζ5、log(z)及乘积；未引入新的M/τ函数。4个裸3e结果保持原样。

## 路线

从唯一prepared C/J及原排序域开始，alpha尺度给1/eps，端点尺度化给有限个Beta因子及比值矩。显式分离ε极点后，剩余F(e,2e)与Z(e)在零点邻域可展开。F是标准3F2；Z具有已逐项积分的正收敛和。18/19全ε仍含这个和，不宣称有限Gamma闭式。

integrate.mpl 在新Maple会话计算F的6项、Z的4项。control.mpl 交换F的积分顺序，并计算可由裸母图Gamma式验证的另一组6项。没有读取旧的T1 mixed系数，没有用PSLQ拟合系数。

## 实际检查

- check_math.py：7图原核/纯线因子/唯一C-I规范、alpha与投影换元；6个反向积分结果、6个Gamma控制；21个线性系数关系和49个crossing系数。
- check_numeric.py：F/Z的6个高精度余项样本，七图I的21个O(eps)样本；七组原始七维J的随机化Sobol（8×2^17），经验误差比较通过，最大相对差约1.05%。
- check_series.py：e=.3时Z的32项显式级数，与独立一维求积比较；差落在解析正尾界内。数值本身非区间算术。
- publish.py：不可变r001，显式49系数与全零字Hlog；定义、参数域、规范和源哈希齐全。
- verify_results.py：离线校验结果、Hlog重建、来源/交付字节，见current-delivery.json。

主要运行见 [current-run.json](current-run.json)。Maple2025、HyperInt及periodLookups路径/哈希保留；数值检查使用mpmath、NumPy、SciPy。首次SymPy核验遇到非整数幂未化简，改用正坐标图因式分解后通过；未关闭HyperInt发散检查。

复验：

    D:/py/python.exe -B -X utf8 workflow/checks/tgv-t1-empty/verify_results.py

重新求积可分别运行 launch.py、control.mpl、check_math.py、check_numeric.py、check_series.py。launch.py创建新run；publish.py拒绝覆盖已有r001。耗时较长的数值任务宜以有日志的后台进程执行，避免前台120秒时限。

原分类README/catalog及原workspace是输入阶段的绑定快照，保留原字节。当前状态另由5e-3g/solutions.md和本次结果索引维护。没有编辑外部montopar或Kernel/Lean，formal=not_run，没有Git提交。
