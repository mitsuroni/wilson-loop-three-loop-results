# 三圈 N=4 总量：闭合有限线的 UV 与坐标尺度连接（r002）

2026-10-19，G5 / WILSON-2.11，Codex，主工作区 master。本文接续[第一阶段](derivation.md)，不修改其打印式直接替换失败的记录。检查脚本和逐项结果在[本阶段报告](../../../../workflow/checks/n4-w3-normalization-20261019/uv-bridge/README.md)。本修订采用研究级文献关系；formal=not_run，绝对有限常数仍未知。

## 先固定要连接的量

取 λ=α_s/π，T=T_F n_f，L=(L_s+L_t)/2，ρ=L_s−L_t=log(s/t)。w_n 是 log W 的系数，F_n 是所有尺度幂一并展开后的严格 ε⁰ 系数。统一
\[
G(\epsilon)=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),\qquad
\mu_x^2={e^{2\gamma_E}\over4}\mu_{\overline{\rm MS}}^2.
\]
基本场按 D=4−2ε 的 DRED 组织：6+2ε 个实标量、两个伴随 Dirac 费米子；内部伴随化不改变外部 Wilson 线的 C_F。采用已有共同裸匹配和耦合约定，不另乘有限 Wilson 重整化。

对于无 running 的 N=4：
\[
w_n=S_n\!\left[-{k_n\over n^2\epsilon^2}-{h_n\over2n\epsilon}\right]
+{k_n\rho^2\over2}+C_n+O(\epsilon),\qquad
S_n=e^{n\epsilon L_s}+e^{n\epsilon L_t},
\]
从而
\[
F_n=-k_n(L^2-\rho^2/4)-h_nL+C_n.
\]
这里的 k、h、C 均可包含颜色；下文 N=4 三圈除以 C_F C_A²。此写法不能把 S_n 外的有限项误当作严格 Laurent 有限项。

## 从裸一圈归一化确定尺度

[Bassetto 等 hep-ph/9303314](https://arxiv.org/pdf/hep-ph/9303314) 的 (3.1) 给出了完整 D 维裸闭合光状平行四边形。他们在轴向 ML 规范独立复现 Feynman 规范结果：
\[
W^{(1)}
=-{g^2C_F\over(2\pi)^{D/2}}\,
 {4\Gamma(D/2-1)\over(4-D)^2}
 \left[(-\varrho^2+i0)^{2-D/2}+(\varrho^2+i0)^{2-D/2}\right],
\quad \varrho^2=\mu_0^2\,n\!\cdot n^*.
\]
此处 \(\varrho\) 是文献的尺度变量，不能与本文的形状 \(\rho\) 混淆。原文矢量满足 \(n\cdot n^*=2LT\)；本项目二次不变量为 \(s=2\,n\cdot n^*\)。作者源的 MS 宏明确是 overline MS。1992 原文的全文本轮未取得；实际读取的是上述 1993 独立复现及 2019 文献，不声称读过未取得的原 PDF。

使用同一个裸耦合
\[
g_B^2\mu_0^{2\epsilon}
=4\pi^2\lambda\left({e^{\gamma_E}\mu_{\overline{\rm MS}}^2\over4\pi}\right)^\epsilon Z_\lambda
\]
后，将原始 π、2π、4π 的 ε 幂统一写成每圈一个 G。得到下列坐标尺度；它们不是三种物理重整化方案：

| 原始表达式 | 提取 G 后的坐标尺度平方 |
|---|---|
| 项目原参数积分的 π 因子 | \(\mu_x^2=e^{2\gamma_E}\mu_{\overline{\rm MS}}^2/4\) |
| (3.1) 已积分闭合轮廓的 2π 因子 | \(\mu_\Box^2=e^{2\gamma_E}\mu_{\overline{\rm MS}}^2/2\) |
| FGM 局部线积分的 4π 因子 | \(\mu_{\rm loc}^2=e^{2\gamma_E}\mu_{\overline{\rm MS}}^2\) |

因而
\[
\mu_\Box^2={\mu_{\rm loc}^2\over2}=2\mu_x^2,\qquad
\mu_\Box^2(n\cdot n^*)=\mu_x^2s.
\]
局部参数积分中的尺度与继承自完整闭合轮廓的几何尺度必须分别解释。FGM 的局部二次式 Q=4uv 也不能在换成边矢量内积时丢掉其 2 的因子。上述字典由裸因子先确定，h₂、h₃ 尚未参与；没有从期望的三圈常数解出一个平移。

一个独立留出是光状平行四边形的 ρ=iπ 分支：对边切向量平行且光状，Feynman 规范下对边一圈收缩为零。余下原始相邻项
\[
-{G(\epsilon)\over\epsilon^2}
\left(e^{\epsilon\ell}+e^{\epsilon(\ell+i\pi)}\right)
\]
在 \(L=\ell+i\pi/2\) 下的有限项为 \(-L^2+\pi^2/12\)，等于本项目 \(F_1/C_F=-L_sL_t+\pi^2/3\) 在同一分支的值。六个 ε>0 的原参数积分数值控制也通过。它不是正 s/t 基点 C₃ 的计算。

若误把 μ_loc 原样复制进闭合几何日志，两条 cusp 日志同时增加 log 2，而 RG 右侧多出 \(-4k_n\log2\)。N=4 两圈、三圈分别得到非零的 \(\pi^2\log2/3\)、\(-11\pi^4\log2/180\)。旧 G2 r002 的字面残差记录保留；本修订提供后续尺度解释，不覆盖其历史证据。

## 有限线的非 cusp 系数

[Falcioni–Gardi–Milloy 1909.00697](https://arxiv.org/pdf/1909.00697) 的 (5.5)、(5.8)–(5.11) 区分有限和无限线：
\[
\mu{d\log W_\Box\over d\mu}
=-2\gamma_{\rm cusp}
[\log\{\mu_\Box^2(n\cdot n^*+i0)\}
+\log\{\mu_\Box^2(-n\cdot n^*+i0)\}]
-\Gamma_\Box ,
\]
\[
\Gamma_\Box=4\Gamma_{\rm co}^{\rm fin}
=4(\Gamma_\cap-\Gamma_\wedge)
=4(\gamma_G-2B_\delta).
\]
不能以四个无限 wedge 的常数代替有限线端点。该文两圈有显式 Wilson 线计算；(5.16) 三圈由因子化与独立 form factor/PDF 结果预测，原文明确并非直接三圈闭合轮廓积分。本阶段按这个文献边界采用它。

在相同不变量日志下，\(\mu d/d\mu=2\partial_L+\beta_\mu\partial_\lambda\)，故 \(h_n=\Gamma_{\Box,n}/2\)。对 N=4 使用[1712.07274](https://arxiv.org/pdf/1712.07274) 所列的 G₀ 与 B，令 \(g_{\rm int}^2=C_A\lambda/4\)：
\[
G_0=-4\zeta_3g_{\rm int}^4+
(32\zeta_5+\tfrac{80}{3}\zeta_2\zeta_3)g_{\rm int}^6,\quad
B=12\zeta_3g_{\rm int}^4-
(80\zeta_5+16\zeta_2\zeta_3)g_{\rm int}^6.
\]
通过三圈的伴随/基本 Casimir 连接给
\[
{h_2^{N=4}\over C_FC_A}=-{7\over2}\zeta_3,\qquad
\boxed{{h_3^{N=4}\over C_FC_A^2}=6\zeta_5+{11\over6}\zeta_2\zeta_3}.
\]
前式独立匹配本地完整两圈结果。这里没有静默改写 1712.07274 的 (2.7) 打印符号；第一阶段 literal substitution 失败的结论仍然有效，也不据此宣布作者勘误。

把新 h₃ 与已知 k₃ 加到[五族反项后差分](../../../../workflow/checks/n4-five-family-20261019/renormalized-difference-r001.json)，四个 QCD 通道 AA、Af、Ff、ff 的 k₃、h₃ 分别与独立 cusp 文献和 FGM (5.16)/2 一致。完整有限尺度导数也用本地 F₁/F₂ 和这些外部 k₃/h₃ 逐色核验。共有 28 条精确检查、6 条原一圈数值控制通过；后者最大归一化误差 3.16e−81。数值控制不覆盖原三圈多重积分。

## 现在总量只余什么

以 \(\zeta_4=\pi^4/90\)，得到
\[
\boxed{\frac{F_3^{N=4}}{C_FC_A^2}
=C_3^{N=4}
-\frac{11\pi^4}{720}(L^2-\rho^2/4)
-\left(6\zeta_5+\frac{11\pi^2\zeta_3}{36}\right)L.}
\]
新[四颜色总量数据](../../../../workflow/checks/n4-w3-normalization-20261019/uv-bridge/qcd-total-up-to-C3-r001.json)给出反项后所有极点和有限函数，AA 的有限常数仍带 C3_N4。它是 research 数据，absolute_F3_complete=false；没有更新共享 current 选择器，也没有完成形式证明或全 QCD AA 有限形状的外部对照。

实际还审计了[2510.07377](https://arxiv.org/pdf/2510.07377) 的作者源与参考文献。(12) 给 correlator/hard/jet/Wilson 因子化，(13) 显示几何日志中的 /16，(14) 再次联系有限线非 cusp 与 threshold 系数；(17)、(21) 的三圈结果仅控制对数项。正文明确留下三圈 regular terms。其 N=4 对照用最高超越权重原理，不能据此填入本地 DRED 的绝对常数。反解因子化需要完整 N=4 correlator 及独立匹配的 N=4 hard/jet 有限函数，当前没有这套有限转换证据。

下一实验被压缩成一个数：
\[
C_3^{N=4}
=\left.[\epsilon^0]\,\frac{(\log W_{N=4})_{\lambda^3}}{C_FC_A^2}\right|_{L=0,\ s/t=1}.
\]
有限常数不受 RG 或形状 Ward 导数约束。四维 HKM 插入函数仍缺最终 D 维积分所需的 evanescent 数据，不能直接积分后称为该常数。

[C3 实验合同 r002](../../../../workflow/checks/n4-w3-normalization-20261019/uv-bridge/C3-extraction-contract-r002.json)区分两种 ε 预算：已有 connected-web/反项路径使用 w₁ 到 ε²、w₂ 到 ε¹；若从未取 log 的 W₃−W₁W₂+W₁³/3 直接组合，则 W₁ 需到 ε⁴、W₂ 需到 ε²，因为 W₂ 的最高极点为 ε⁻⁴。早期 r001 草案只描述反项所需预算，现由 r002 明确补全，保留草案不作当前合同。

同一基点上 \(F_{3,AA}^{QCD}=C_3^{N=4}+\Delta F_{3,AA}\)，现有差分约为 −0.5259594010446349257640。因此，直接计算 N=4 基点或取得经过低圈校准的完整有限转换，才会消去最后这个未知数。G5 继续 active，WILSON-2.11 的有限方案字典及基点验收仍开放。
