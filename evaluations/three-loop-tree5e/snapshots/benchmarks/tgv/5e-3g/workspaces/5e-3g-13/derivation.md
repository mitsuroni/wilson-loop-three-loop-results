# 5e-3g-13

2026-10-08修订：旧r001的单极点与有限项受Laurent组装预算遗漏影响，当前结果为[r002](../../results/5e-3g-13/r002.json)。原核/旧计算记录保留，完整总和递归检验及更正见[连续推导](../../../../recursion-audit/sum-5e.md)。新修订不代表此单图相消方向已独立获证，formal=not_run。

状态：规范原始输入及 IBP 分解已生成并检查；I 的洛朗系数尚未计算。formal=not_run。

输入：[prepared r001](../../prepared/5e-3g-r001.json)。TGV 附着边 [1, 2, 3]，类型 nondegenerate；IBP 共 7 项。

固定三圈颜色坐标为 (−1/4, 1/2, 0)，先完整重构再选择 a₁；Euler 补偿已含一次，禁止重复相乘。α 积分为单纯形，原 s 区域排序保留。

| 项 | 类型 | 相对于共同 C 的系数 | 剩余 s 维数 |
|---|---|---|---|
| boundary_s11_0 | cusp_boundary | (-1)/(4*(1-2*eps)) | 4 |
| boundary_s11_1 | cusp_boundary | (1)/(4*(1-2*eps)) | 4 |
| boundary_s21_0 | cusp_boundary | (1 + z1)/(4*(1-2*eps)) | 4 |
| boundary_s21_s22 | contact_boundary | (-1 - z1)/(4*(1-2*eps)) | 4 |
| boundary_s31_0 | cusp_boundary | (-z1)/(4*(1-2*eps)) | 4 |
| boundary_s31_s32 | contact_boundary | (z1)/(4*(1-2*eps)) | 4 |
| recoil | recoil | 1 | 5 |

下一步：在原有序域上寻找可解析消去的变量，必要时采用已校验的局部超立方映射；合并所有边界及 recoil 后连接唯一 I 定义，再提取 ε⁻⁶…ε⁰。不能由各项的单独发散推断最终极点。
