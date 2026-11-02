> **2026-09-19替代说明：本页以下内容为历史审查记录，不能再作为当前实现依据。** 正确推导已统一写入[4e-3g+3g-general-derivation.md](../4e-3g+3g-general-derivation.md)，包含完整D维算子、链式法则、尺度和有序域。当前Python实现及旧C/J尚未修复；历史勾选状态不授予新实现正确性。

---

# 4e-3g+3g 审查：双 TGV 积分 $I = C \times J$

> 目标：独立重推导 + CAS 验证 [4e-3g+3g-general-derivation.md](../4e-3g+3g-general-derivation.md) 的每一节，产出可执行验证脚本和 table-heavy 风格的最终推导文档。

## 任务分解

| 任务 | 范围 | 依赖 |
|------|------|------|
| **小任务**：有效顶点 $V_{2\text{TGV}}$ | §2c Lorentz 分解验证 | [effective-vertices.md](../effective-vertices.md) §3 |
| **大任务**：完整积分翻译 $C \times J$ | §1-§14 全部，含 montopar 代码 | 小任务完成后 |

## 文件夹结构

```
4e-3g+3g-review/
├── README.md                        ← 本文件：审查总索引
├── derivation.md                    ← 最终成品：完整 C×J 推导（table-heavy 风格）
│
├── verify/
│   ├── v01-lorentz-decomposition.m   ← 小任务：独立重算 V_{2TGV} 15 项 Lorentz 展开
│   ├── v02-color-factors.m           ← §3b: 各边分布颜色因子 SUNSimplify
│   ├── v03-operator-table.m          ← §5c/§9c-1: 验证 20 项算子表 ℓ→i∂_X
│   ├── v04-symanzik-param.m          ← §7: 独立推导 B_5, U, 5-点恒等式
│   ├── v05-quadratic-chain-rule.m    ← §8: 验证 F_i^(p), H_ij^(pq), 二次链式法则
│   ├── v06-S-channel-full.m          ← S 道 (3,3,4,4) 完整 C×J 对照
│   ├── v07-T-channel-full.m          ← T 道 (3,4,3,4) 完整 C×J 对照
│   ├── v08-U-channel.m               ← U 道 (3,4,4,3) Jacobi 恒等式
│   └── v09-numerical-crosscheck.m    ← 随机运动学数值代入，全轨道交叉验证
│
└── outputs/                          ← .m 文件的运行输出（.gitignore）
```

## 验证矩阵

| 节 | 内容 | 验证文件 | 方法 | 状态 |
|----|------|---------|------|:---:|
| §1 | 运动学与几何 | — | 与 [basics/02-kinematics.md](../../../basics/02-kinematics.md) 对照 | ⬜ |
| §2c | $V_{2\text{TGV}}$ Lorentz 分解 | `v01-lorentz-decomposition.m` | CAS 独立重算 + 数值交叉验证 | ✅ |
| §3b | 颜色因子 | `v02-color-factors.m` | FeynCalc SUNSimplify | ✅ |
| §4-§5 | $\ell$ 结构分类 + $\ell\to i\partial_X$ | `v03-operator-table.m` | 因子化验证 + 子项计数 | ✅ |
| §6 | Fourier 变换 $z_A,z_B$ 解耦 | `v04-symanzik-param.m` | 独立推导 δ→FT | ✅ |
| §7 | Symanzik 参数化 $B_5,\mathcal{U}$ | `v04-symanzik-param.m` | 独立推导 + 矩阵验证 | ✅ |
| §8 | 二次型链式法则 | `v05-quadratic-chain-rule.m` | 符号求导验证 $F_i^{(p)},H_{ij}^{(pq)}$ | ✅ |
| §9 | $C\times J$ 因子分解 | `v05`-`v06` | Pnum 公式 + 逐轨道对照 | ⬜ |
| §10 | 退化自动处理 | `v06`-`v07` | S 道（双退化）和 T 道（非退化）| ⬜ |
| §11-§13 | $\Delta_{ij}$、Pnum 构造、已知 orbit | `v06`-`v07` | 逐轨道代入对照 | ⬜ |

## 每个 .m 文件的结构约定

```mathematica
(* vXX-xxx.m *)
(* 验证目标：... *)
(* 方法：独立重推导，不引用 general-derivation.md 的中间结果 *)

(* === 1. 设置 === *)
<< FeynCalc`

(* === 2. 定义 === *)
(* 从 effective-vertices.md 或 first principles 出发 *)

(* === 3. 计算 === *)

(* === 4. 对照 === *)
(* 与 general-derivation.md 逐项比对 *)
(* 输出：PASS/FAIL 表格 *)

(* === 5. 数值交叉验证 === *)
(* 随机运动学数值代入 *)
```

## 运行说明

```bash
# 单个文件
wolframscript -f verify/v01-lorentz-decomposition.m

# 全部
for f in verify/v*.m; do wolframscript -f "$f"; done
```

## 与 montopar 代码的对应

| .m 验证文件 | montopar Python 模块 |
|------------|---------------------|
| `v03-operator-table.m` | `double_tgv.py` → `OPERATOR_TABLE` |
| `v04-symanzik-param.m` | `double_tgv.py` → `compute_B5()`, `compute_U()` |
| `v05-quadratic-chain-rule.m` | `double_tgv.py` → `compute_F_i_p()`, `compute_H_ij_pq()` |
| `v06`-`v08` 完整结果 | `double_tgv.py` → golden reference for tests |