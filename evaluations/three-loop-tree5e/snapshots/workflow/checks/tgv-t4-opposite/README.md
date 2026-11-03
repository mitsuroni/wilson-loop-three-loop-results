# T4[1] 计算与复验

[五图解表](../../../benchmarks/tgv/families/t4-opposite/solutions.md) · [完整推导](../../../benchmarks/tgv/families/t4-opposite/method.md) · [Ω₆与原基底](../../../benchmarks/tgv/families/t4-opposite/basis-functions.md)

27/28/29/41/42的完整ε表示与ε⁻⁶…ε⁰已计算：35位置12非零。27=29、41=42在全ε成立；28的J满足self-crossing。前三图简单极点，后两图二阶极点。本批完成后5e-3g为44/44、308位置174非零。

## 实际计算与检查

- 五个原C/J、排序、D/F尺度桥、三组反射（含28自身）、五项负例；41/42互补域16个单项矩检查及实际二圈T4的S/P和源结果哈希连接。
- partial.mpl在α和中间边上求积，保留两端点的T4母核；weights.mpl求两个pure矩形的e⁰/e¹。global_ibp.mpl保留cusp/contact/recoil后给出27/28的完整极点及有限项。
- 28的直接点态母核路线与完整边界路线精确一致。27的前一路线出现真实端点发散，转为先保留w=0 cusp；原尝试的V27m1也不复用为正确极点。
- 41/42由J41+J42=JT4·Jpure及反射求解，物理I的1/8系数独立核验；有限项复用二圈T4与H0的实际系数。
- 35项表达精确回代Hlog；权重6相对旧14维空间仅新增Ω₆，增广秩15。有新的有理字坐标见证，不仅是数值拟合。Ω₆系数为27/29的−1/64、28的1/32。
- 12组原七维根采样：其中6组27–29与独立D框架比较；两组41+42因子化、两组42/41反射比较。另有8组全局有序IBP与原核比较。每组8次扰乱、每次2¹⁷个Sobol点。各比较最大相对差分别约0.891%、0.392%、1.204%、0.718%，均通过经验误差准则。
- 两组原六维cusp面与完整Gamma降维对照；三个小ε cusp余项样本符合O(ε)。这只检验cusp，不声称完整I的余项已核验。
- 12个已积分双核级数有限块与独立七维Gauss-Jacobi一致。
- Ω₆在z=0.2、0.4的canonical Hlog与T₆(z)+T₆(1/z)经独立ODE级数及Möbius延拓核对，最大差约2.6×10⁻⁵¹；两个交错周期另由一维积分核对，三个基本Hlog检验通过。工作精度70位、640项，并非严格区间包围。

formal=not_run；没有小ε完整I余项验证、从原积分出发的第二条独立解析路线或物理复z/i0分支证明。Kernel默认47根、Lean和外部montopar均未改。新Ω的数值定义检查不冒充原七维积分有限项的独立求解。

## 全ε与保存格式

27–29为已逐项积分的显式双核收敛级数；41/42为既有T4级数×Meijer G/Gamma/cosh。初始绝对收敛strip为1/2<Re(ε)<1，再亚纯延拓；尚无快速无限求和或严格尾界。

本批结果用计算输出格式wilson-tgv-integral-result/2：Hlog字母可为有理数，并显式存储period_powers顺序[ζ2,ζ3,ζ5,log2,A4,A6]。A4/A6由明确的交错和或实一维积分定义。旧结果修订和Kernel AST保持原字节。

## 复算和离线复验

```powershell
D:/py/python.exe -B -X utf8 workflow/checks/tgv-t4-opposite/verify_results.py
```

只读复验检查冻结依赖、原核连接、35项回代、Ω见证和实际数值记录，不启动Maple/QMC。当前交付由current-delivery.json指定。publish.py只用于首次r001；修改数学结果须新建修订。

launch.py/maple_worker.py每次建立新目录，保存Maple源码、依赖哈希、完整日志和真实退出码；仍需检查错误、完成标记与所需产物，因为Maple内部报错也可能返回进程码0。

check_math.py、full_series.py、check_numeric.py、check_cusp_numeric.py分别保存根/归一化、级数块、原域/IBP/因子化、关键cusp检查。numeric_hlog.py由微分方程构造独立数值级数，不调用HyperInt求值器。wide_hlog.py只扩展本批计算侧坐标，不改变旧全局接口。

## 运行和失败记录

当前使用partial、weights、global_ibp、outer28、omega_basis、basis_reduce、cusp_export各指针所对应的成功运行。global_ibp保留完整边界；outer28是可积的28交叉计算；最初outer.mpl的27点态展开在外层发散检查处终止。

早期global_ibp曾误把乘积RegInf列表交给单字reglim接口；随后改为普通零点极限、标量核反射及完整pure面重算。中间还修正了Maple保留名Zero。对含移动字母的v→0不作直接代入，改为重新求该完整pure面。所有失败运行保留，不授予完成状态。

旧Hlog解析器只接受{0,−1}，无法读取本批真正出现的−1/2、−2及交错周期；因此本批使用独立扩展坐标和显式格式/2。新字母已在两种28边界处理中一致复现，Ω的原函数导数和crossing精确成立。

本轮未提交Git。此前T4[]的八图成果仍保持原字节；本轮接着完成五个T4[1]，当前未提交的T4计算合计13图、91位置44非零。