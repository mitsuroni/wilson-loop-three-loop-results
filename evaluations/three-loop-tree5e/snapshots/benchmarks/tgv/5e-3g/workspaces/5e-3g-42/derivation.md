# 5e-3g-42

状态：规范原始输入及 IBP 分解已生成并检查；I 的洛朗系数尚未计算。formal=not_run。

输入：[prepared r001](../../prepared/5e-3g-r001.json)。TGV 附着边 [1, 2, 3]，类型 nondegenerate；IBP 共 7 项。

固定三圈颜色坐标为 (−1/4, 1/2, 0)，先完整重构再选择 a₁；Euler 补偿已含一次，禁止重复相乘。α 积分为单纯形，原 s 区域排序保留。

| 项 | 类型 | 相对于共同 C 的系数 | 剩余 s 维数 |
|---|---|---|---|
| boundary_s22_1 | cusp_boundary | (-1 - z1)/(4*(1-2*eps)) | 4 |
| boundary_s21_s22 | contact_boundary | (1 + z1)/(4*(1-2*eps)) | 4 |
| boundary_s31_0 | cusp_boundary | (-z1)/(4*(1-2*eps)) | 4 |
| boundary_s31_1 | cusp_boundary | (z1)/(4*(1-2*eps)) | 4 |
| boundary_s11_0 | cusp_boundary | (-1)/(4*(1-2*eps)) | 4 |
| boundary_s11_1 | cusp_boundary | (1)/(4*(1-2*eps)) | 4 |
| recoil | recoil | 1 | 5 |

下一步：在原有序域上寻找可解析消去的变量，必要时采用已校验的局部超立方映射；合并所有边界及 recoil 后连接唯一 I 定义，再提取 ε⁻⁶…ε⁰。不能由各项的单独发散推断最终极点。
