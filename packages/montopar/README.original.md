# 2026-10-15 双TGV保留算符出口

`--double-tgv-operator`追加未作用的`J_operator`，保留独立位置母核、D维顶点算符及先微分后限制的次序。`montopar.double_tgv_operator.to_scalar`实际作用导数并恢复原有序标量核；可与`--double-tgv-forms`及`--tgv-ibp`并用，原C与默认输出保持。该模式仍为兼容而生成raw，不宣称省去已有标量计算。见[算符出口与转换](docs/double-tgv-operator.md)。

# 2026-10-09 双TGV紧凑与解析降维初始形式

`--double-tgv-forms`追加可展开的`J_compact`和显式分项`J_reduced`：一般中央参数积成Beta×Gauss，八维到七维；独立双同边有序对输出九项四维Beta边界和。原C/raw与`--tgv-ibp`保持，原始Euler/C4账本不变，无CAS运行依赖。新降维继承正实起始区与共同解析延拓，不宣称任意有限辅助调节下逐点成立。使用、公式及消费合同见[双TGV初始形式](docs/double-tgv-initial-forms.md)。

# 2026-10-06 原生共同符号规则

当前源→C/J 使用 [同一作用量、路径序和 Fourier 字典](docs/shared-sign-conventions.md)，版本 `wilson-qcd-action-path-fourier/1.0`。CLI 要求源声明该契约；所有原生 C 的 `assembly_sign=1`，raw 不含 Euler/C4。FGV、TGV 腿 Π₁、Γ₁的正确符号在 C 生成前由规则得出；旧结果导入的族级符号桥不可用于新输出。下文日期段落是功能开发历史，其旧符号示例以本节为准。

# 2026-10-01 三胶子顶点的Π₁腿修正

2026-10-03: Central one-loop three-gluon vertex `3e-Gamma1` is supported by
schema 1.8: separate gluon/quark C/J, ordered double-simplex raw kernels,
explicit three-sector kernels, and optional complete contour IBP with
`--gamma1-ibp` (alias `--tgv-ibp`). Runtime has no CAS dependency. See
[theory and output contract](docs/gamma1.md). Raw C excludes Euler/C4/crossing;
I is not evaluated. The legacy single-integrand Stage 3 consumer requires an
explicit channel/block adapter.

新增`3e-3g+Pi1`的`tgv_pi1`分派与schema 1.7：12图完整横向张量、原s有序域、两独立α及两项raw C/J，保留两颜色通道；`--tgv-ibp`追加完整端点/接触面和三边recoil，`--package-initial`整族打包。raw阶段无Euler补偿、a1投影或C4重数；I尚未求值。

[完整直接理论及使用](docs/tgv-pi1.md) · [IBP推导](docs/tgv-pi1-ibp.md)。后文“不支持外接TGV”和“IBP尚未完成”等是历史版本描述。既有族的输出协议保持。

# 2026-09-24：双Π₁串联

新增`self_energy_chain`分派与schema 1.5，支持三圈2e-Pi1+Pi1两图的完整C/J、原s域、完整平方颜色和Gamma比U。[理论与运行方式](docs/self-energy-pi1-chain.md)。原单插入schema 1.4输出不变；下文已有能力说明继续适用于各原分支。

# 2026-09-20: double-TGV contour IBP

## 2026-09-22：单个 Π₁ 自能

已支持 `2e-Pi1`（总二圈两图）和 `4e-Pi1`（总三圈十一图）。通过独立 `self_energy` 分派输出 schema 1.4 的完整横向张量 C/J：主表示是原 s 有序域，保留度规/导数分解、完整颜色、一次子图系数和 C4 元数据；raw 不加 Euler、不做 a1 投影、不含重数。未求积分。

```powershell
D:/py/python.exe -m montopar.cli --source ../qgraf-pipeline/outputs/2loop/orbit_definitions.json --type 2e-Pi1 --package-initial
D:/py/python.exe -m montopar.cli --source ../qgraf-pipeline/outputs/3loop/orbit_definitions.json --type 4e-Pi1 --package-initial
```

[推导、C/J 与输出约定](docs/self-energy-pi1.md)。Π₂、重复自能插入及外接 TGV 尚不支持；这些输入明确拒绝。已有 pure/TGV/FGV/mixed/double-TGV 的数学输出保持。


`--tgv-ibp` now also supports `--type "4e-3g+3g"`. It preserves all raw fields and adds complete two-layer ordered IBP C/J terms, including contact-normal corrections and the D-dimensional body remainder. See [output contract](docs/double-tgv-ibp.md) and [complete theory](docs/4e-3g+3g-ibp-theory-r003.md). This exports an equivalent regulated integral representation; evaluation and formal verification remain separate.

# 当前双TGV输出（2026-09-19修复）

4e-3g+3g现为schema 1.3、仅一个规范代表。按修正理论计算完整T1/T2，保留D=4−2ε、链式系数与尺度；原s有序域、五α单纯形、可乘C账本和有向配对显式保存。颜色使用精确四生成元Fierz收缩，无需Mathematica。C4拓扑重数包括相同内部顶点互换；raw阶段不补Euler、不乘重数。

`--type "4e-3g+3g" --package-initial`生成26个workspace输出及汇总，拒绝旧版不完整输入。详见[推导](docs/4e-3g+3g-general-derivation.md)、[逐图C/J](docs/4e-3g+3g-CJ-formulas.md)。本表示未做IBP或积分求值。

# 当前4GV输出（2026-09-18修订）

4e-4g现为schema 1.2，仅一个规范代表；原s有序域、四α单纯形、C可乘账本、分支母式和cube映射均显式保存。零图共用C并标记J=0。新参数 `--fgv-compat-out-dir` 从主结果派生三种兼容输出，修正旧09距离；`--package-initial` 拒绝旧FGV输入。详见[完整推导](docs/4GV-general-derivation.md)、[C/J与命令](docs/4e-4g-CJ-formulas.md)。其他族既有接口保持其版本。

# 当前 TGV / mixed 输出（2026-09-17）

3e-3g、5e-3g 使用单一规范代表及显式 α 单纯形；mixed C 相位已修正，--tgv-ibp 支持完整有序区域。详见 [mixed 有序 IBP](docs/mixed-ordered-ibp.md)。下文若有 A/B 双份输出描述，属于旧版本；当前不生成这两族的 rotated 字段。IBP 完成不代表积分求值完成。

# montopar

montopar 将 qgraf-pipeline 的 orbit 定义直接转换为参数空间中的
$I=C\times J$，结果写入统一的 Wilson-loop workspace。它是 Stage 2 的闭式公式
快速通道，当前支持 pure、TGV、FGV、mixed 与 double-TGV 等已实现 dispatch。

## 安装与基本使用

    cd montopar
    pip install -e .
    montopar --source ../qgraf-pipeline/outputs/2loop/orbit_definitions.json --type 4e

默认输出位置为：

    workspaces/{nloop}loop/{emission_type}/{orbit_id}/montopar.json

处理单个 orbit：

    montopar --source ../qgraf-pipeline/outputs/2loop/orbit_definitions.json --orbit-id 3e-3g-04

## 3e-3g Contour IBP

纯 3e-3g TGV 可以在沿 Wilson 轮廓分部积分后输出积分线性组合：

    montopar --source ../qgraf-pipeline/outputs/2loop/orbit_definitions.json --type 3e-3g --tgv-ibp

该选项不会替换原始表示。输出同时包含：

- `J_raw`：A 组的原始三维积分；
- `J_raw_rotated`：C4+1 后 B 组的原始积分；
- `J_ibp`：A 组的 IBP 后积分和；
- `J_ibp_rotated`：B 组重新计算得到的 IBP 后积分和。

求导前星型核使用公共前因子

$$
\boxed{
C_\star
=\frac{C_A C_F(4\pi\alpha_s)^2\Gamma(1-2\varepsilon)}
{256\pi^{4-2\varepsilon}}
(-2v_{14}+i0)^{2\varepsilon}
=\frac{C_{\rm raw}}{4(1-2\varepsilon)}.
}
$$

边界项引用 $C_\star$，非退化 recoil 项继续引用 $C_{\rm raw}$。

### 已实现的分解

| orbit | 类型 | A 组项数 | B 组项数 | IBP 结构 |
|---|---|---:|---:|---|
| 3e-3g-01 | 两腿同边 | 3 | 3 | 2 cusp + 1 contact |
| 3e-3g-02 | 两腿同边 | 3 | 3 | 2 cusp + 1 contact |
| 3e-3g-03 | 两腿同边 | 3 | 3 | 2 cusp + 1 contact |
| 3e-3g-04 | 三腿互异 | 7 | 7 | 6 cusp + 1 recoil |

退化情形的 IBP 必须在原始有序三角域上完成。同边参数的对角面是 contact
边界，不能按普通超立方体面丢弃。非退化情形在三个独立的 $[0,1]$ 变量上产生六个
有向边界面，余下部分保留为具有较短分子的三维 recoil 积分。

### `J_ibp` 数据结构

每个 IBP 表示是可独立处理的 term 列表，而不是单个大 integrand：

    {
      "J_ibp": {
        "schema_version": "1.0",
        "method": "contour_ibp",
        "status": "complete",
        "case": "nondegenerate",
        "prefactors": {
          "C_raw": "...",
          "C_star": "...",
          "relation": "C_star = C_raw/(4*(1-2*eps))"
        },
        "operator": {
          "derivative_vector": ["z1", "-1 - z1", "1"]
        },
        "terms": [
          {
            "term_id": "boundary_s21_0",
            "kind": "cusp_boundary",
            "coefficient": "-z1",
            "C_ref": "C_star",
            "combined_C": "...",
            "fixed_parameters": {"s21": "0"},
            "identified_parameters": {},
            "integrand": "...",
            "variables": ["x1", "x2"],
            "bounds": {"x1": [0, 1], "x2": [0, 1]},
            "denom_power": "1 - 2*eps"
          }
        ],
        "certificates": {
          "operator_identity": true,
          "boundary_orientation": true,
          "pnum_identity": true
        }
      }
    }

若域结构无法可靠识别或代数证书失败，`status` 为
`unsupported`，并保留 `reason`；原始 `J_raw` 不受影响。

## 下游边界

当前 j-integral/dimreduce 默认读取单个 `J_raw.integrand`。montopar 已为每个
`J_ibp.terms[*]` 提供完整系数、积分核、变量与边界，但 Stage 3 尚未自动逐项
初始化和汇总这些 term。不要在没有显式多 term bridge 的情况下把 `J_ibp` 冒充
`J_raw`。

## 测试

    cd montopar
    python -m pytest -k "not slow" -q

Contour IBP 的主要测试位于 `tests/test_tgv_ibp.py`，并由 CLI、emitter 与
package-builder 测试覆盖序列化和透传。

## 推导与设计

- `docs/3e-3g-contour-ibp-sumCJ.md`：四个 orbit 的 A/B $\sum_iC_iJ_i$ 汇总；
- `docs/3e-3g-02-contour-ibp-derivation.md`：退化 02 的手动推导与闭式结果；
- `docs/3e-3g-04-contour-ibp-derivation.md`：非退化 04 的 A/B 分解；
- `docs/3e-3g-CJ-formulas.md`：3e-3g 原始表示与 IBP 总览；
- `../docs/superpowers/specs/2026-09-07-montopar-3e-3g-contour-ibp-design.md`：
  schema、算法和测试契约。

## 关键约定

- 参考标度必须写为 $(\pm v_{ij}+i0)^{n\varepsilon}$，`i0` 位于幂的括号内；
- A/B 组使用同一算法，B 组必须从 C4+1 旋转后的 orbit 重新计算；
- 边界约束先作用于原始 $s$ 变量，再建立该 term 自己的 $s\to x$ 映射；
- `J_raw` 是兼容接口和等价性证书，不因启用 IBP 而删除。

Double-TGV source validation was tightened on 2026-09-19: see [the closed source contract](docs/double-tgv-source-contract.md). Extra topology and duplicate endpoints are rejected; schema 1.3 C/J for the current 26 diagrams is unchanged.

## Bare two-loop 1PI self energy

`2e-Pi2` is supported at total three loops (schema 1.6, dispatch `self_energy_pi2`). Both original ordered C/J outputs retain AA/Af/Ff, including the CF^2 TF nf second coordinate. Full-D coefficients and the common-I convention are in [theory](docs/self-energy-pi2.md). No Euler factor or color projection is applied in raw output.

```powershell
montopar --source ../qgraf-pipeline/outputs/3loop/orbit_definitions.json --type 2e-Pi2 --package-initial
```
