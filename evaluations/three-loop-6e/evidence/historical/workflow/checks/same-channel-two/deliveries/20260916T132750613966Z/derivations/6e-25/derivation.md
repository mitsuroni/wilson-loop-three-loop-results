# 6e-25：[2] 族解析推导

状态：全 ε 显式收敛级数及 ε⁻⁶…ε⁰ 系数已计算并核验；formal: not_run。当前根为 prepared 0.3；z>0，实 0<ε<1 的完整 I 及 ε=0 邻域展开。

## 根与原顺序

原外到内顺序：s31, s23, s42, s22, s41, s21。

~~~text
((-(1+z1))^2*1) * ((1-s21)*(1-s41) + z1*s21*s41)^(-1+eps) * ((1-s22)*s31*1)^(-1+eps) * ((1-s23)*(1-s42) + z1*s23*s42)^(-1+eps)
~~~

原边界：
~~~json
{
  "s21": [
    "0",
    "s22"
  ],
  "s22": [
    "0",
    "s23"
  ],
  "s23": [
    "0",
    "1"
  ],
  "s31": [
    "0",
    "1"
  ],
  "s41": [
    "0",
    "s42"
  ],
  "s42": [
    "0",
    "1"
  ]
}
~~~

## 从根降维

| 原边低到高 | 端点字 |
|---|---|
| s21, s22, s23 | OBO |
| s31 | A |
| s41, s42 | OO |

配对：id；原邻边积分后的分子权重：(1-v)^eps-(1-w)^eps。应用同时反射后归入 id_middle。

a1=1/4，sym=4，κ=-1/16。C 的 Euler 补偿逐图核对，只用一次。

标准 t 低高：s41, s42；s 低高：s21, s23。可积性、换序、核反射、切片与完整求和定义见[族方法](../../families/same-channel-two/method.md)。

## 完整 J 与 I

$$J=\frac{z^{\varepsilon} F_{id middle}{\left(\varepsilon,\frac{1}{z} \right)}}{\varepsilon}.$$

$$I=- \frac{\left(z^{\frac{\varepsilon}{2}} F_{id middle}{\left(\varepsilon,z \right)} + z^{- \frac{\varepsilon}{2}} F_{id middle}{\left(\varepsilon,\frac{1}{z} \right)}\right) e^{- 3 \gamma \varepsilon} \Gamma^{3}\left(1 - \varepsilon\right)}{32 \varepsilon \cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)}}.$$

F_id_middle 是[显式收敛级数](../../families/same-channel-two/full-series.json)，每个有限块已积完；实现见 workflow/checks/same-channel-two/full_series.py。它不是未求值积分的别名，也不宣称有限特殊函数闭式。

## Laurent 系数

| 幂次 | 系数 |
|---:|---|
| -6 | $0$ |
| -5 | $0$ |
| -4 | $0$ |
| -3 | $0$ |
| -2 | $0$ |
| -1 | $\frac{M{\left(0,1,z \right)} M{\left(1,1,z \right)}}{64} + \frac{\pi^{2} M{\left(1,1,z \right)}}{96} + \frac{M{\left(1,3,z \right)}}{96} - \frac{N{\left(1,2,z \right)}}{64}$ |
| 0 | $SigmaB$ |

余项 O(ε)。Hlog 全式和 Σ 的定义见[系数函数](../../families/same-channel-two/basis-functions.md)。三圈有限项保留 crossing 奇部分乘 ln z；G 与 sech 的一次项为零。

## 检查与输出

结果修订：[ r002 ](../../results/6e-25/r002.json)。

原顺序整数 J、非整数原域 Sobol、独立对数核积分、余项缩放、精确基底转换及收敛级数检查见[检查目录](../../../../workflow/checks/same-channel-two/README.md)。没有执行 Kernel/Lean 整链或物理复分支延拓。
