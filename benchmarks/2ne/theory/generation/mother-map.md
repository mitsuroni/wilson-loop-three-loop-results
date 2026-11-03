# Step 3：逐图母族生成表

由 `workflow/checks/2ne-mothers-complete/audit.py --write` 从原多项式、bounds及通道配对重建。
定义和归一化见 [母族理论](mothers.md)。上横线表示同一母积分取值于1/z；所有全边常数用R2、R12表示。
本表覆盖全部47图的J生成与现行完整ε的I连接。[Q/R/F定义桥](remaining-bridges.md)给出后16图的逐式推导；没有用有限阶关系代替全ε关系。

$$J_G=z^{b\epsilon}\epsilon^{-P}B_G(\epsilon,z),\qquad b=a+r/2,\quad P=2(n-r).$$

| 图 | 族 | b | P | κ | B_G |
|---|---|---:|---:|---:|---|
| 2e-01 | [] | 1 | 2 | -1 | $1$ |
| 2e-02 | [1] | 1/2 | 0 | -1/2 | $\mathcal T^{[1]}_{1;\varnothing;\varnothing}$ |
| 4e-01 | [] | 2 | 4 | -1/4 | $\frac{1}{4}$ |
| 4e-02 | [1] | 3/2 | 2 | -1/4 | $\mathcal T^{[1]}_{1;\varnothing;\varnothing} - \mathcal T^{[1]}_{1;\varnothing;A}$ |
| 4e-03 | [] | 1 | 4 | -1/4 | $1 - R_{2}$ |
| 4e-04 | [2] | 1 | 0 | -1/8 | $\mathcal T^{[2]}_{1,2;\varnothing,\varnothing;\varnothing,\varnothing}$ |
| 4e-05 | [1] | 1/2 | 2 | -1/4 | $\mathcal T^{[1]}_{1;B;\varnothing}$ |
| 4e-06 | [1,1] | 1 | 0 | -1/16 | $\overline{\mathcal T^{[1]}_{1;\varnothing;\varnothing}} \mathcal T^{[1]}_{1;\varnothing;\varnothing}$ |
| 6e-01 | [] | 3 | 6 | -1/8 | $\frac{1}{36}$ |
| 6e-02 | [] | 3 | 6 | -1/16 | $\frac{1}{36}$ |
| 6e-03 | [] | 3 | 6 | -1/16 | $\frac{1}{36}$ |
| 6e-04 | [1] | 5/2 | 4 | -1/8 | $\frac{\mathcal T^{[1]}_{1;\varnothing;\varnothing}}{4} - \frac{\mathcal T^{[1]}_{1;\varnothing;A}}{2} + \frac{\mathcal T^{[1]}_{1;\varnothing;AA}}{2}$ |
| 6e-05 | [1] | 5/2 | 4 | -1/16 | $\frac{\mathcal T^{[1]}_{1;\varnothing;\varnothing}}{4} - \frac{\mathcal T^{[1]}_{1;\varnothing;A}}{2} + \frac{\mathcal T^{[1]}_{1;\varnothing;AA}}{2}$ |
| 6e-06 | [1] | 5/2 | 4 | -1/16 | $\frac{\mathcal T^{[1]}_{1;\varnothing;A}}{2} - \mathcal T^{[1]}_{1;\varnothing;AA}$ |
| 6e-07 | [] | 2 | 6 | -1/16 | $- \frac{R_{12}}{2} + \frac{R_{2}}{2}$ |
| 6e-08 | [] | 2 | 6 | -1/8 | $\frac{R_{12}}{4} - \frac{R_{2}}{2} + \frac{1}{4}$ |
| 6e-09 | [] | 2 | 6 | -1/16 | $\frac{R_{12}}{4} - \frac{R_{2}}{2} + \frac{1}{4}$ |
| 6e-10 | [2] | 2 | 2 | -1/8 | $\mathcal T^{[2]}_{1,2;\varnothing,\varnothing;\varnothing,\varnothing} - \mathcal T^{[2]}_{1,2;\varnothing,\varnothing;\varnothing,A}$ |
| 6e-11 | [2] | 2 | 2 | -1/16 | $\mathcal T^{[2]}_{1,2;\varnothing,\varnothing;\varnothing,A} - \mathcal T^{[2]}_{1,2;\varnothing,\varnothing;A,\varnothing}$ |
| 6e-12 | [2] | 2 | 2 | -1/16 | $\mathcal T^{[2]}_{2,1;\varnothing,\varnothing;\varnothing,\varnothing} - \mathcal T^{[2]}_{2,1;\varnothing,\varnothing;\varnothing,A}$ |
| 6e-13 | [1] | 3/2 | 4 | -1/16 | $R_{2} \mathcal T^{[1]}_{1;B;\varnothing} - R_{2} \mathcal T^{[1]}_{1;B;A}$ |
| 6e-14 | [1] | 3/2 | 4 | -1/8 | $\mathcal T^{[1]}_{1;B;\varnothing} \left(1 - R_{2}\right) + \mathcal T^{[1]}_{1;B;A} \left(R_{2} - 1\right)$ |
| 6e-15 | [1] | 3/2 | 4 | -1/16 | $\mathcal T^{[1]}_{1;B;A} \left(1 - R_{2}\right)$ |
| 6e-16 | [1] | 3/2 | 4 | -1/16 | $\mathcal T^{[1]}_{1;\varnothing;\varnothing} \left(1 - R_{2}\right) + \mathcal T^{[1]}_{1;\varnothing;A} \left(R_{2} - 1\right) + \mathcal T^{[1]}_{1;B;\varnothing} \left(R_{2} - 1\right) + \mathcal T^{[1]}_{1;B;A} \left(1 - R_{2}\right)$ |
| 6e-17 | [] | 1 | 6 | -1/16 | $- \frac{R_{12}}{2} + \frac{R_{2}}{2}$ |
| 6e-18 | [] | 1 | 6 | -1/8 | $\frac{R_{12}}{4} - \frac{R_{2}}{2} + \frac{1}{4}$ |
| 6e-19 | [] | 1 | 6 | -1/16 | $\frac{R_{12}}{4} - \frac{R_{2}}{2} + \frac{1}{4}$ |
| 6e-20 | [3] | 3/2 | 0 | -1/16 | $\mathcal T^{[3]}_{1,2,3;\varnothing,\varnothing,\varnothing;\varnothing,\varnothing,\varnothing}$ |
| 6e-21 | [3] | 3/2 | 0 | -1/32 | $\mathcal T^{[3]}_{1,3,2;\varnothing,\varnothing,\varnothing;\varnothing,\varnothing,\varnothing}$ |
| 6e-22 | [3] | 3/2 | 0 | -1/32 | $\mathcal T^{[3]}_{2,1,3;\varnothing,\varnothing,\varnothing;\varnothing,\varnothing,\varnothing}$ |
| 6e-23 | [2] | 1 | 2 | -1/8 | $\mathcal T^{[2]}_{1,2;B,\varnothing;\varnothing,\varnothing}$ |
| 6e-24 | [2] | 1 | 2 | -1/16 | $\mathcal T^{[2]}_{2,1;B,\varnothing;\varnothing,\varnothing}$ |
| 6e-25 | [2] | 1 | 2 | -1/16 | $\mathcal T^{[2]}_{1,2;\varnothing,B;\varnothing,\varnothing} - \mathcal T^{[2]}_{1,2;B,\varnothing;\varnothing,\varnothing}$ |
| 6e-26 | [1] | 1/2 | 4 | -1/8 | $\frac{\mathcal T^{[1]}_{1;BB;\varnothing}}{2}$ |
| 6e-27 | [1] | 1/2 | 4 | -1/16 | $\frac{\mathcal T^{[1]}_{1;B;\varnothing}}{2} - \mathcal T^{[1]}_{1;BB;\varnothing}$ |
| 6e-28 | [1] | 1/2 | 4 | -1/16 | $\frac{\mathcal T^{[1]}_{1;BB;\varnothing}}{2}$ |
| 6e-29 | [1] | 3/2 | 4 | -1/8 | $\mathcal T^{[1]}_{1;\varnothing;B} - \mathcal T^{[1]}_{1;\varnothing;BA} - \mathcal T^{[1]}_{1;\varnothing;AB}$ |
| 6e-30 | [1] | 3/2 | 4 | -1/16 | $\mathcal T^{[1]}_{1;\varnothing;BA}$ |
| 6e-31 | [1] | 3/2 | 4 | -1/16 | $\mathcal T^{[1]}_{1;\varnothing;\varnothing} \left(1 - R_{2}\right) - \mathcal T^{[1]}_{1;\varnothing;B} + \mathcal T^{[1]}_{1;\varnothing;AB}$ |
| 6e-32 | [1,1] | 2 | 2 | -1/16 | $\overline{\mathcal T^{[1]}_{1;\varnothing;\varnothing}} \mathcal T^{[1]}_{1;\varnothing;\varnothing} - \overline{\mathcal T^{[1]}_{1;\varnothing;\varnothing}} \mathcal T^{[1]}_{1;\varnothing;A} - \mathcal T^{[1]}_{1;\varnothing;\varnothing} \overline{\mathcal T^{[1]}_{1;\varnothing;B}} + \mathcal T^{[1]}_{1;\varnothing;A} \overline{\mathcal T^{[1]}_{1;\varnothing;B}}$ |
| 6e-33 | [1,1] | 2 | 2 | -1/8 | $\overline{\mathcal T^{[1]}_{1;\varnothing;B}} \mathcal T^{[1]}_{1;\varnothing;\varnothing} - \overline{\mathcal T^{[1]}_{1;\varnothing;B}} \mathcal T^{[1]}_{1;\varnothing;A}$ |
| 6e-34 | [1,1] | 2 | 2 | -1/16 | $\overline{\mathcal T^{[1]}_{1;\varnothing;B}} \mathcal T^{[1]}_{1;\varnothing;A}$ |
| 6e-35 | [] | 1 | 6 | -1/16 | $R_{2}^{2} - 2 R_{2} + 1$ |
| 6e-36 | [1] | 5/2 | 4 | -1/32 | $\mathcal T^{[1]}_{1;\varnothing;\varnothing} - \mathcal T^{[1]}_{1;\varnothing;A} - \mathcal T^{[1]}_{1;A;\varnothing} + \mathcal T^{[1]}_{1;A;A}$ |
| 6e-37 | [1,2] | 3/2 | 0 | -1/16 | $\overline{\mathcal T^{[1]}_{1;\varnothing;\varnothing}} \mathcal T^{[2]}_{1,2;\varnothing,\varnothing;\varnothing,\varnothing}$ |
| 6e-38 | [1,2] | 3/2 | 0 | -1/32 | $\overline{\mathcal T^{[1]}_{1;\varnothing;\varnothing}} \mathcal T^{[2]}_{2,1;\varnothing,\varnothing;\varnothing,\varnothing}$ |
| 6e-39 | [1] | 1/2 | 4 | -1/32 | $\mathcal T^{[1]}_{1;B;B}$ |

## 已闭合的47条完整I连接

G在下列自动生成公式中表示补偿因子Z；R2、R12为理论正文中的Gamma比值。H字典完全固定；Qid01/Qid10表示较大/较小端点的A插入，R123等记录配对。横线/后缀bar均表示z→1/z。完整定义见[连接推导](remaining-bridges.md)。
每条公式从本页B_G经统一crossing投影产生，与旧修订的紧凑全ε式精确相等；这里没有重新积分Laurent系数。

### 2e-01

$$I=- \frac{G}{eps^{2}}$$

[现行结果](../../../../benchmarks/2ne/results/2e-01/r002.json)；全ε符号残差为0。

### 2e-02

$$I=- \frac{G H_{0}}{2 \cosh{\left(\frac{\ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/2e-02/r002.json)；全ε符号残差为0。

### 4e-01

$$I=- \frac{G^{2}}{16 eps^{4}}$$

[现行结果](../../../../benchmarks/2ne/results/4e-01/r002.json)；全ε符号残差为0。

### 4e-02

$$I=- \frac{G^{2} \left(\frac{H_{0} e^{\frac{\ell eps}{2}}}{2} + \frac{H_{0} e^{- \frac{\ell eps}{2}}}{2} - \frac{H_{1} e^{\frac{\ell eps}{2}}}{2} - \frac{\bar{H1} e^{- \frac{\ell eps}{2}}}{2}\right)}{4 eps^{2} \cosh{\left(\ell eps \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/4e-02/r002.json)；全ε符号残差为0。

### 4e-03

$$I=- \frac{G^{2} \left(1 - R_{2}\right)}{4 eps^{4} \cosh{\left(\ell eps \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/4e-03/r002.json)；全ε符号残差为0。

### 4e-04

$$I=- \frac{G^{2} Qid}{8 \cosh{\left(\ell eps \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/4e-04/r002.json)；全ε符号残差为0。

### 4e-05

$$I=- \frac{G^{2} \left(\frac{H_{0} e^{\frac{\ell eps}{2}}}{2} + \frac{H_{0} e^{- \frac{\ell eps}{2}}}{2} - \frac{H_{1} e^{\frac{\ell eps}{2}}}{2} - \frac{\bar{H1} e^{- \frac{\ell eps}{2}}}{2}\right)}{4 eps^{2} \cosh{\left(\ell eps \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/4e-05/r001.json)；全ε符号残差为0。

### 4e-06

$$I=- \frac{G^{2} H_{0}^{2}}{16 \cosh{\left(\ell eps \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/4e-06/r002.json)；全ε符号残差为0。

### 6e-01

$$I=- \frac{G^{3}}{288 eps^{6}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-01/r002.json)；全ε符号残差为0。

### 6e-02

$$I=- \frac{G^{3}}{576 eps^{6}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-02/r002.json)；全ε符号残差为0。

### 6e-03

$$I=- \frac{G^{3}}{576 eps^{6}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-03/r002.json)；全ε符号残差为0。

### 6e-04

$$I=- \frac{G^{3} \left(\frac{H_{0} e^{\ell eps}}{8} + \frac{H_{0} e^{- \ell eps}}{8} - \frac{H_{1} e^{\ell eps}}{4} - \frac{\bar{H1} e^{- \ell eps}}{4} + \frac{H_{2} e^{\ell eps}}{8} + \frac{\bar{H2} e^{- \ell eps}}{8}\right)}{8 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-04/r002.json)；全ε符号残差为0。

### 6e-05

$$I=- \frac{G^{3} \left(\frac{H_{0} e^{\ell eps}}{8} + \frac{H_{0} e^{- \ell eps}}{8} - \frac{H_{1} e^{\ell eps}}{4} - \frac{\bar{H1} e^{- \ell eps}}{4} + \frac{H_{2} e^{\ell eps}}{8} + \frac{\bar{H2} e^{- \ell eps}}{8}\right)}{16 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-05/r001.json)；全ε符号残差为0。

### 6e-06

$$I=- \frac{G^{3} \left(\frac{H_{1} e^{\ell eps}}{4} + \frac{\bar{H1} e^{- \ell eps}}{4} - \frac{H_{2} e^{\ell eps}}{4} - \frac{\bar{H2} e^{- \ell eps}}{4}\right)}{16 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-06/r001.json)；全ε符号残差为0。

### 6e-07

$$I=- \frac{G^{3} \left(- \frac{R_{12}}{2} + \frac{R_{2}}{2}\right) \cosh{\left(\frac{\ell eps}{2} \right)}}{16 eps^{6} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-07/r002.json)；全ε符号残差为0。

### 6e-08

$$I=- \frac{G^{3} \left(\frac{R_{12}}{4} - \frac{R_{2}}{2} + \frac{1}{4}\right) \cosh{\left(\frac{\ell eps}{2} \right)}}{8 eps^{6} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-08/r002.json)；全ε符号残差为0。

### 6e-09

$$I=- \frac{G^{3} \left(\frac{R_{12}}{4} - \frac{R_{2}}{2} + \frac{1}{4}\right) \cosh{\left(\frac{\ell eps}{2} \right)}}{16 eps^{6} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-09/r002.json)；全ε符号残差为0。

### 6e-10

$$I=- \frac{G^{3} \left(\left(Qid - Qid_{01}\right) e^{\frac{\ell eps}{2}} + \left(Qid - \bar{Qid01}\right) e^{- \frac{\ell eps}{2}}\right)}{16 eps^{2} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-10/r002.json)；全ε符号残差为0。

### 6e-11

$$I=- \frac{G^{3} \left(\left(Qid_{01} - Qid_{10}\right) e^{\frac{\ell eps}{2}} + \left(\bar{Qid01} - \bar{Qid10}\right) e^{- \frac{\ell eps}{2}}\right)}{32 eps^{2} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-11/r002.json)；全ε符号残差为0。

### 6e-12

$$I=- \frac{G^{3} \left(\left(Qswap - Qswap_{01}\right) e^{\frac{\ell eps}{2}} + \left(Qswap - \bar{Qswap01}\right) e^{- \frac{\ell eps}{2}}\right)}{32 eps^{2} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-12/r002.json)；全ε符号残差为0。

### 6e-13

$$I=- \frac{G^{3} \left(H_{0} R_{2} - H_{1} R_{2} - \bar{H1} R_{2} + HX R_{2}\right)}{16 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-13/r001.json)；全ε符号残差为0。

### 6e-14

$$I=- \frac{G^{3} \left(- H_{0} R_{2} + H_{0} + H_{1} R_{2} - H_{1} + \bar{H1} R_{2} - \bar{H1} - HX R_{2} + HX\right)}{8 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-14/r001.json)；全ε符号残差为0。

### 6e-15

$$I=- \frac{G^{3} \left(- \frac{H_{1} R_{2}}{2} + \frac{H_{1}}{2} - \frac{\bar{H1} R_{2}}{2} + \frac{\bar{H1}}{2} + HX R_{2} - HX\right)}{16 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-15/r001.json)；全ε符号残差为0。

### 6e-16

$$I=- \frac{G^{3} \left(- \frac{H_{1} R_{2}}{2} + \frac{H_{1}}{2} - \frac{\bar{H1} R_{2}}{2} + \frac{\bar{H1}}{2} + HX R_{2} - HX\right)}{16 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-16/r001.json)；全ε符号残差为0。

### 6e-17

$$I=- \frac{G^{3} \left(- \frac{R_{12}}{2} + \frac{R_{2}}{2}\right) \cosh{\left(\frac{\ell eps}{2} \right)}}{16 eps^{6} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-17/r002.json)；全ε符号残差为0。

### 6e-18

$$I=- \frac{G^{3} \left(\frac{R_{12}}{4} - \frac{R_{2}}{2} + \frac{1}{4}\right) \cosh{\left(\frac{\ell eps}{2} \right)}}{8 eps^{6} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-18/r002.json)；全ε符号残差为0。

### 6e-19

$$I=- \frac{G^{3} \left(\frac{R_{12}}{4} - \frac{R_{2}}{2} + \frac{1}{4}\right) \cosh{\left(\frac{\ell eps}{2} \right)}}{16 eps^{6} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-19/r002.json)；全ε符号残差为0。

### 6e-20

$$I=- \frac{G^{3} R_{123}}{16 \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-20/r001.json)；全ε符号残差为0。

### 6e-21

$$I=- \frac{G^{3} \left(R_{132} + R_{213}\right)}{64 \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-21/r001.json)；全ε符号残差为0。

### 6e-22

$$I=- \frac{G^{3} \left(R_{132} + R_{213}\right)}{64 \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-22/r001.json)；全ε符号残差为0。

### 6e-23

$$I=- \frac{G^{3} \left(\left(Qid - Qid_{01}\right) e^{\frac{\ell eps}{2}} + \left(Qid - \bar{Qid01}\right) e^{- \frac{\ell eps}{2}}\right)}{16 eps^{2} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-23/r002.json)；全ε符号残差为0。

### 6e-24

$$I=- \frac{G^{3} \left(\left(Qswap - Qswap_{01}\right) e^{\frac{\ell eps}{2}} + \left(Qswap - \bar{Qswap01}\right) e^{- \frac{\ell eps}{2}}\right)}{32 eps^{2} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-24/r002.json)；全ε符号残差为0。

### 6e-25

$$I=- \frac{G^{3} \left(\left(Qid_{01} - Qid_{10}\right) e^{\frac{\ell eps}{2}} + \left(\bar{Qid01} - \bar{Qid10}\right) e^{- \frac{\ell eps}{2}}\right)}{32 eps^{2} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-25/r002.json)；全ε符号残差为0。

### 6e-26

$$I=- \frac{G^{3} \left(\frac{H_{0} e^{\ell eps}}{8} + \frac{H_{0} e^{- \ell eps}}{8} - \frac{H_{1} e^{\ell eps}}{4} - \frac{\bar{H1} e^{- \ell eps}}{4} + \frac{H_{2} e^{\ell eps}}{8} + \frac{\bar{H2} e^{- \ell eps}}{8}\right)}{8 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-26/r001.json)；全ε符号残差为0。

### 6e-27

$$I=- \frac{G^{3} \left(\frac{H_{1} e^{\ell eps}}{4} + \frac{\bar{H1} e^{- \ell eps}}{4} - \frac{H_{2} e^{\ell eps}}{4} - \frac{\bar{H2} e^{- \ell eps}}{4}\right)}{16 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-27/r001.json)；全ε符号残差为0。

### 6e-28

$$I=- \frac{G^{3} \left(\frac{H_{0} e^{\ell eps}}{8} + \frac{H_{0} e^{- \ell eps}}{8} - \frac{H_{1} e^{\ell eps}}{4} - \frac{\bar{H1} e^{- \ell eps}}{4} + \frac{H_{2} e^{\ell eps}}{8} + \frac{\bar{H2} e^{- \ell eps}}{8}\right)}{16 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-28/r001.json)；全ε符号残差为0。

### 6e-29

$$I=- \frac{G^{3} \left(H_{0} - H_{1} - \bar{H1} + HW\right)}{8 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-29/r001.json)；全ε符号残差为0。

### 6e-30

$$I=- \frac{G^{3} \left(- \frac{H_{0} R_{2}}{2} + \frac{H_{1}}{2} + \frac{\bar{H1}}{2} - \frac{HW}{2}\right)}{16 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-30/r001.json)；全ε符号残差为0。

### 6e-31

$$I=- \frac{G^{3} \left(- \frac{H_{0} R_{2}}{2} + \frac{H_{1}}{2} + \frac{\bar{H1}}{2} - \frac{HW}{2}\right)}{16 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-31/r001.json)；全ε符号残差为0。

### 6e-32

$$I=- \frac{G^{3} \left(\left(H_{0} H_{1} - H_{1}^{2}\right) e^{\frac{\ell eps}{2}} + \left(H_{0} \bar{H1} - \bar{H1}^{2}\right) e^{- \frac{\ell eps}{2}}\right)}{32 eps^{2} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-32/r001.json)；全ε符号残差为0。

### 6e-33

$$I=- \frac{G^{3} \left(\left(H_{0}^{2} - 2 H_{0} H_{1} + H_{1}^{2}\right) e^{\frac{\ell eps}{2}} + \left(H_{0}^{2} - 2 H_{0} \bar{H1} + \bar{H1}^{2}\right) e^{- \frac{\ell eps}{2}}\right)}{16 eps^{2} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-33/r001.json)；全ε符号残差为0。

### 6e-34

$$I=- \frac{G^{3} \left(\left(H_{0} H_{1} - H_{1}^{2}\right) e^{\frac{\ell eps}{2}} + \left(H_{0} \bar{H1} - \bar{H1}^{2}\right) e^{- \frac{\ell eps}{2}}\right)}{32 eps^{2} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-34/r001.json)；全ε符号残差为0。

### 6e-35

$$I=- \frac{G^{3} \left(R_{2}^{2} - 2 R_{2} + 1\right) \cosh{\left(\frac{\ell eps}{2} \right)}}{16 eps^{6} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-35/r002.json)；全ε符号残差为0。

### 6e-36

$$I=- \frac{G^{3} \left(\frac{H_{0} e^{\ell eps}}{2} + \frac{H_{0} e^{- \ell eps}}{2} - H_{1} e^{\ell eps} - \bar{H1} e^{- \ell eps} + \frac{HD e^{\ell eps}}{2} + \frac{\bar{HD} e^{- \ell eps}}{2}\right)}{32 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-36/r002.json)；全ε符号残差为0。

### 6e-37

$$I=- \frac{G^{3} H_{0} Qid}{16 \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-37/r001.json)；全ε符号残差为0。

### 6e-38

$$I=- \frac{G^{3} H_{0} Qswap}{32 \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-38/r001.json)；全ε符号残差为0。

### 6e-39

$$I=- \frac{G^{3} \left(\frac{H_{0} e^{\ell eps}}{2} + \frac{H_{0} e^{- \ell eps}}{2} - H_{1} e^{\ell eps} - \bar{H1} e^{- \ell eps} + \frac{HD e^{\ell eps}}{2} + \frac{\bar{HD} e^{- \ell eps}}{2}\right)}{32 eps^{4} \cosh{\left(\frac{3 \ell eps}{2} \right)}}$$

[现行结果](../../../../benchmarks/2ne/results/6e-39/r001.json)；全ε符号残差为0。

完整机器记录见 [mother-map-r002.json](mother-map-r002.json)，实际检查见 [复算说明](../../../../workflow/checks/2ne-mothers-complete/README.md)。
