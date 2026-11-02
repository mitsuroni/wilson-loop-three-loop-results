# 4e-3g+3g：修正后的原参数 C/J（schema 1.3）

2026-09-19。26个规范轨道的raw C/J已重新生成：17个非零、9个颜色零。仅维护一个规范代表；下列公式尚未求积分，尚未做颜色a1投影、Euler补偿或crossing组装。

令 $\epsilon=e$、$A_{\rm ref}=-2v_{14}+i0$。共同标量因子为

$$C_* = C_A^2C_F\left(\frac{\alpha_s}{\pi}\right)^3
\frac{\pi^{3e}\Gamma(2-3e)}{64}A_{\rm ref}^{3e}.$$

下表的C列为 $C/C_*$；零图的有效J置零，但完整几何核、域和零原因仍保留在JSON。旋转重数另列，**没有乘入C或J**。

| 图号 | 顶点有向腿序及位置 | C/C* | C4重数 | 状态 |
|---|---|---:|---:|---|
| 01 | `(s31,s41; s42,s43)` | -1 | 4 | 非零核 |
| 02 | `(s31,s42; s41,s43)` | 0 | 4 | 颜色零 |
| 03 | `(s31,s43; s41,s42)` | 1 | 4 | 非零核 |
| 04 | `(s31,s32; s41,s42)` | -1 | 4 | 非零核 |
| 05 | `(s31,s41; s32,s42)` | 0 | 4 | 颜色零 |
| 06 | `(s31,s42; s32,s41)` | 1 | 4 | 非零核 |
| 07 | `(s31,s32; s33,s41)` | -1 | 4 | 非零核 |
| 08 | `(s31,s33; s32,s41)` | 0 | 4 | 颜色零 |
| 09 | `(s31,s41; s32,s33)` | 1 | 4 | 非零核 |
| 10 | `(s21,s41; s42,s43)` | -1 | 4 | 非零核 |
| 11 | `(s21,s42; s41,s43)` | 0 | 4 | 颜色零 |
| 12 | `(s21,s43; s41,s42)` | 1 | 4 | 非零核 |
| 13 | `(s21,s31; s41,s42)` | -1 | 4 | 非零核 |
| 14 | `(s21,s41; s31,s42)` | 0 | 4 | 颜色零 |
| 15 | `(s21,s42; s31,s41)` | 1 | 4 | 非零核 |
| 16 | `(s21,s31; s32,s41)` | -1 | 4 | 非零核 |
| 17 | `(s21,s32; s31,s41)` | 0 | 4 | 颜色零 |
| 18 | `(s21,s41; s31,s32)` | 1 | 4 | 非零核 |
| 19 | `(s21,s22; s41,s42)` | -1 | 2 | 非零核 |
| 20 | `(s21,s41; s22,s42)` | 0 | 2 | 颜色零 |
| 21 | `(s21,s42; s22,s41)` | 1 | 2 | 非零核 |
| 22 | `(s21,s22; s31,s41)` | -1 | 4 | 非零核 |
| 23 | `(s21,s31; s22,s41)` | 0 | 4 | 颜色零 |
| 24 | `(s21,s41; s22,s31)` | 1 | 4 | 非零核 |
| 25 | `(s11,s21; s31,s41)` | -1 | 2 | 非零核 |
| 26 | `(s11,s31; s21,s41)` | 0 | 1 | 颜色零 |

每图原s上下限与外到内次序均显式保存在输出。设 $d\mu_\alpha=d^5\alpha\,\delta(1-\sum_{i=1}^5\alpha_i)$，全部α非负，

$$\mathcal U=(\alpha_1+\alpha_2+\alpha_5)(\alpha_3+\alpha_4+\alpha_5)-\alpha_5^2.$$

17个非零图共同采用

$$\boxed{J=\int_{\mathcal D}d^4s\int d\mu_\alpha\,
\frac{\prod_i\alpha_i^{-e}}{\mathcal U^{1+2e}}
\left[\frac{P_1}{\widehat B_F^{3-3e}}+\frac{P_2}{\widehat B_F^{2-3e}}\right],
\quad P_1=-\frac{2-3e}{8}\widehat T_1,\quad P_2=-\frac14\widehat T_2.}$$

$\widehat B_F$继承共同Feynman处方，在正实起始区为 $\widehat B+i0\mathcal U$；一般复运动学延拓尚未数值验证。

完整的 $P_1,P_2,\widehat T_1,\widehat T_2$ 多项式已在JSON逐图展开。为便于阅读，下文列原域、B和较短的T2结构；T1使用统一的向量公式，所有量由腿位置确定：

$$\xi_i=\sum_{r<e_i}E_r+s_iE_{e_i},\qquad
E_r\cdot_*E_t=\frac{v_r\cdot v_t}{v_{14}},\qquad t_i=E_{e_i},$$

$$g_i=-2\alpha_i\sum_{j\ne i}\alpha_jc_{ij}(\xi_i-\xi_j),$$

$$V_L=t_{12}(g_2-g_1)+t_1(2t_2\cdot_*g_1+t_2\cdot_*g_2)
-t_2(t_1\cdot_*g_1+2t_1\cdot_*g_2),$$

$$V_R=t_{34}(g_4-g_3)+t_3(2t_4\cdot_*g_3+t_4\cdot_*g_4)
-t_4(t_3\cdot_*g_3+2t_3\cdot_*g_4),\qquad
\widehat T_1=-V_L\cdot_*V_R.$$

这里 $t_{ij}=t_i\cdot_*t_j$，$c_{12}=\alpha_3+\alpha_4+\alpha_5$、$c_{34}=\alpha_1+\alpha_2+\alpha_5$、其余跨顶点 $c_{ij}=\alpha_5$。第二个分子由下文K±给出：

$$\widehat T_2=-2\alpha_5\big[(\alpha_1\alpha_3+\alpha_2\alpha_4)K_+
+(\alpha_1\alpha_4+\alpha_2\alpha_3)K_-\big].$$

K±已保留 $D=4-2e$，不是四维近似。非零图19、21、25的C4重数为2；其他非零图为4。零图20、26的重数分别为2、1。

## 4e-3g+3g-01

外到内：`s31 → s43 → s42 → s41`。原界：`s31: 0…1`，`s41: 0…s42`，`s42: 0…s43`，`s43: 0…1`。

$$\widehat B=- \alpha_{1} z_{1} \left(s_{31} - 1\right) \left(\alpha_{2} \alpha_{3} s_{41} + \alpha_{2} \alpha_{4} s_{41} + \alpha_{2} \alpha_{5} s_{41} + \alpha_{3} \alpha_{5} s_{42} + \alpha_{4} \alpha_{5} s_{43}\right).$$

$$K_+=0,\qquad K_-=0.$$

## 4e-3g+3g-02

外到内：`s31 → s43 → s42 → s41`。原界：`s31: 0…1`，`s42: 0…s43`，`s41: 0…s42`，`s43: 0…1`。

$C=0$，颜色筛选后的有效 $J=0$。未求积的几何多项式另保留在JSON。

## 4e-3g+3g-03

外到内：`s31 → s43 → s42 → s41`。原界：`s31: 0…1`，`s43: 0…1`，`s41: 0…s42`，`s42: 0…s43`。

$$\widehat B=- \alpha_{1} z_{1} \left(s_{31} - 1\right) \left(\alpha_{2} \alpha_{3} s_{43} + \alpha_{2} \alpha_{4} s_{43} + \alpha_{2} \alpha_{5} s_{43} + \alpha_{3} \alpha_{5} s_{41} + \alpha_{4} \alpha_{5} s_{42}\right).$$

$$K_+=0,\qquad K_-=0.$$

## 4e-3g+3g-04

外到内：`s32 → s31 → s42 → s41`。原界：`s31: 0…s32`，`s32: 0…1`，`s41: 0…s42`，`s42: 0…1`。

$$\widehat B=- \alpha_{5} z_{1} \left(\alpha_{3} s_{41} + \alpha_{4} s_{42}\right) \left(\alpha_{1} s_{31} - \alpha_{1} + \alpha_{2} s_{32} - \alpha_{2}\right).$$

$$K_+=z_{1}^{2},\qquad K_-=- z_{1}^{2}.$$

## 4e-3g+3g-05

外到内：`s32 → s31 → s42 → s41`。原界：`s31: 0…s32`，`s41: 0…s42`，`s32: 0…1`，`s42: 0…1`。

$C=0$，颜色筛选后的有效 $J=0$。未求积的几何多项式另保留在JSON。

## 4e-3g+3g-06

外到内：`s32 → s31 → s42 → s41`。原界：`s31: 0…s32`，`s42: 0…1`，`s32: 0…1`，`s41: 0…s42`。

$$\widehat B=- z_{1} \left(\alpha_{1} \alpha_{2} \alpha_{3} s_{31} s_{42} - \alpha_{1} \alpha_{2} \alpha_{3} s_{42} + \alpha_{1} \alpha_{2} \alpha_{4} s_{31} s_{42} - \alpha_{1} \alpha_{2} \alpha_{4} s_{42} + \alpha_{1} \alpha_{2} \alpha_{5} s_{31} s_{42} - \alpha_{1} \alpha_{2} \alpha_{5} s_{42} + \alpha_{1} \alpha_{3} \alpha_{4} s_{32} s_{41} - \alpha_{1} \alpha_{3} \alpha_{4} s_{41} + \alpha_{1} \alpha_{4} \alpha_{5} s_{31} s_{41} - \alpha_{1} \alpha_{4} \alpha_{5} s_{41} + \alpha_{2} \alpha_{3} \alpha_{4} s_{32} s_{41} - \alpha_{2} \alpha_{3} \alpha_{4} s_{41} + \alpha_{2} \alpha_{3} \alpha_{5} s_{32} s_{42} - \alpha_{2} \alpha_{3} \alpha_{5} s_{42} + \alpha_{3} \alpha_{4} \alpha_{5} s_{32} s_{41} - \alpha_{3} \alpha_{4} \alpha_{5} s_{41}\right).$$

$$K_+=- 2 z_{1}^{2} \left(eps + 1\right),\qquad K_-=z_{1}^{2} \left(2 eps - 7\right).$$

## 4e-3g+3g-07

外到内：`s33 → s32 → s31 → s41`。原界：`s31: 0…s32`，`s32: 0…s33`，`s33: 0…1`，`s41: 0…1`。

$$\widehat B=- \alpha_{4} s_{41} z_{1} \left(\alpha_{1} \alpha_{3} s_{33} - \alpha_{1} \alpha_{3} + \alpha_{1} \alpha_{5} s_{31} - \alpha_{1} \alpha_{5} + \alpha_{2} \alpha_{3} s_{33} - \alpha_{2} \alpha_{3} + \alpha_{2} \alpha_{5} s_{32} - \alpha_{2} \alpha_{5} + \alpha_{3} \alpha_{5} s_{33} - \alpha_{3} \alpha_{5}\right).$$

$$K_+=0,\qquad K_-=0.$$

## 4e-3g+3g-08

外到内：`s33 → s32 → s31 → s41`。原界：`s31: 0…s32`，`s33: 0…1`，`s32: 0…s33`，`s41: 0…1`。

$C=0$，颜色筛选后的有效 $J=0$。未求积的几何多项式另保留在JSON。

## 4e-3g+3g-09

外到内：`s33 → s32 → s31 → s41`。原界：`s31: 0…s32`，`s41: 0…1`，`s32: 0…s33`，`s33: 0…1`。

$$\widehat B=- \alpha_{2} s_{41} z_{1} \left(\alpha_{1} \alpha_{3} s_{31} - \alpha_{1} \alpha_{3} + \alpha_{1} \alpha_{4} s_{31} - \alpha_{1} \alpha_{4} + \alpha_{1} \alpha_{5} s_{31} - \alpha_{1} \alpha_{5} + \alpha_{3} \alpha_{5} s_{32} - \alpha_{3} \alpha_{5} + \alpha_{4} \alpha_{5} s_{33} - \alpha_{4} \alpha_{5}\right).$$

$$K_+=0,\qquad K_-=0.$$

## 4e-3g+3g-10

外到内：`s21 → s43 → s42 → s41`。原界：`s21: 0…1`，`s41: 0…s42`，`s42: 0…s43`，`s43: 0…1`。

$$\widehat B=\alpha_{1} \left(\alpha_{2} \alpha_{3} s_{21} s_{41} z_{1} + \alpha_{2} \alpha_{3} s_{21} s_{41} - \alpha_{2} \alpha_{3} s_{21} - \alpha_{2} \alpha_{3} s_{41} + \alpha_{2} \alpha_{3} + \alpha_{2} \alpha_{4} s_{21} s_{41} z_{1} + \alpha_{2} \alpha_{4} s_{21} s_{41} - \alpha_{2} \alpha_{4} s_{21} - \alpha_{2} \alpha_{4} s_{41} + \alpha_{2} \alpha_{4} + \alpha_{2} \alpha_{5} s_{21} s_{41} z_{1} + \alpha_{2} \alpha_{5} s_{21} s_{41} - \alpha_{2} \alpha_{5} s_{21} - \alpha_{2} \alpha_{5} s_{41} + \alpha_{2} \alpha_{5} + \alpha_{3} \alpha_{5} s_{21} s_{42} z_{1} + \alpha_{3} \alpha_{5} s_{21} s_{42} - \alpha_{3} \alpha_{5} s_{21} - \alpha_{3} \alpha_{5} s_{42} + \alpha_{3} \alpha_{5} + \alpha_{4} \alpha_{5} s_{21} s_{43} z_{1} + \alpha_{4} \alpha_{5} s_{21} s_{43} - \alpha_{4} \alpha_{5} s_{21} - \alpha_{4} \alpha_{5} s_{43} + \alpha_{4} \alpha_{5}\right).$$

$$K_+=0,\qquad K_-=0.$$

## 4e-3g+3g-11

外到内：`s21 → s43 → s42 → s41`。原界：`s21: 0…1`，`s42: 0…s43`，`s41: 0…s42`，`s43: 0…1`。

$C=0$，颜色筛选后的有效 $J=0$。未求积的几何多项式另保留在JSON。

## 4e-3g+3g-12

外到内：`s21 → s43 → s42 → s41`。原界：`s21: 0…1`，`s43: 0…1`，`s41: 0…s42`，`s42: 0…s43`。

$$\widehat B=\alpha_{1} \left(\alpha_{2} \alpha_{3} s_{21} s_{43} z_{1} + \alpha_{2} \alpha_{3} s_{21} s_{43} - \alpha_{2} \alpha_{3} s_{21} - \alpha_{2} \alpha_{3} s_{43} + \alpha_{2} \alpha_{3} + \alpha_{2} \alpha_{4} s_{21} s_{43} z_{1} + \alpha_{2} \alpha_{4} s_{21} s_{43} - \alpha_{2} \alpha_{4} s_{21} - \alpha_{2} \alpha_{4} s_{43} + \alpha_{2} \alpha_{4} + \alpha_{2} \alpha_{5} s_{21} s_{43} z_{1} + \alpha_{2} \alpha_{5} s_{21} s_{43} - \alpha_{2} \alpha_{5} s_{21} - \alpha_{2} \alpha_{5} s_{43} + \alpha_{2} \alpha_{5} + \alpha_{3} \alpha_{5} s_{21} s_{41} z_{1} + \alpha_{3} \alpha_{5} s_{21} s_{41} - \alpha_{3} \alpha_{5} s_{21} - \alpha_{3} \alpha_{5} s_{41} + \alpha_{3} \alpha_{5} + \alpha_{4} \alpha_{5} s_{21} s_{42} z_{1} + \alpha_{4} \alpha_{5} s_{21} s_{42} - \alpha_{4} \alpha_{5} s_{21} - \alpha_{4} \alpha_{5} s_{42} + \alpha_{4} \alpha_{5}\right).$$

$$K_+=0,\qquad K_-=0.$$

## 4e-3g+3g-13

外到内：`s21 → s31 → s42 → s41`。原界：`s21: 0…1`，`s31: 0…1`，`s41: 0…s42`，`s42: 0…1`。

$$\widehat B=- \alpha_{1} \alpha_{2} \alpha_{3} s_{21} s_{31} + \alpha_{1} \alpha_{2} \alpha_{3} s_{31} - \alpha_{1} \alpha_{2} \alpha_{4} s_{21} s_{31} + \alpha_{1} \alpha_{2} \alpha_{4} s_{31} - \alpha_{1} \alpha_{2} \alpha_{5} s_{21} s_{31} + \alpha_{1} \alpha_{2} \alpha_{5} s_{31} + \alpha_{1} \alpha_{3} \alpha_{5} s_{21} s_{41} z_{1} + \alpha_{1} \alpha_{3} \alpha_{5} s_{21} s_{41} - \alpha_{1} \alpha_{3} \alpha_{5} s_{21} - \alpha_{1} \alpha_{3} \alpha_{5} s_{41} + \alpha_{1} \alpha_{3} \alpha_{5} + \alpha_{1} \alpha_{4} \alpha_{5} s_{21} s_{42} z_{1} + \alpha_{1} \alpha_{4} \alpha_{5} s_{21} s_{42} - \alpha_{1} \alpha_{4} \alpha_{5} s_{21} - \alpha_{1} \alpha_{4} \alpha_{5} s_{42} + \alpha_{1} \alpha_{4} \alpha_{5} - \alpha_{2} \alpha_{3} \alpha_{5} s_{31} s_{41} z_{1} + \alpha_{2} \alpha_{3} \alpha_{5} s_{41} z_{1} - \alpha_{2} \alpha_{4} \alpha_{5} s_{31} s_{42} z_{1} + \alpha_{2} \alpha_{4} \alpha_{5} s_{42} z_{1}.$$

$$K_+=- z_{1} \left(z_{1} + 1\right),\qquad K_-=z_{1} \left(z_{1} + 1\right).$$

## 4e-3g+3g-14

外到内：`s21 → s31 → s42 → s41`。原界：`s21: 0…1`，`s41: 0…s42`，`s31: 0…1`，`s42: 0…1`。

$C=0$，颜色筛选后的有效 $J=0$。未求积的几何多项式另保留在JSON。

## 4e-3g+3g-15

外到内：`s21 → s31 → s42 → s41`。原界：`s21: 0…1`，`s42: 0…1`，`s31: 0…1`，`s41: 0…s42`。

$$\widehat B=\alpha_{1} \alpha_{2} \alpha_{3} s_{21} s_{42} z_{1} + \alpha_{1} \alpha_{2} \alpha_{3} s_{21} s_{42} - \alpha_{1} \alpha_{2} \alpha_{3} s_{21} - \alpha_{1} \alpha_{2} \alpha_{3} s_{42} + \alpha_{1} \alpha_{2} \alpha_{3} + \alpha_{1} \alpha_{2} \alpha_{4} s_{21} s_{42} z_{1} + \alpha_{1} \alpha_{2} \alpha_{4} s_{21} s_{42} - \alpha_{1} \alpha_{2} \alpha_{4} s_{21} - \alpha_{1} \alpha_{2} \alpha_{4} s_{42} + \alpha_{1} \alpha_{2} \alpha_{4} + \alpha_{1} \alpha_{2} \alpha_{5} s_{21} s_{42} z_{1} + \alpha_{1} \alpha_{2} \alpha_{5} s_{21} s_{42} - \alpha_{1} \alpha_{2} \alpha_{5} s_{21} - \alpha_{1} \alpha_{2} \alpha_{5} s_{42} + \alpha_{1} \alpha_{2} \alpha_{5} - \alpha_{1} \alpha_{3} \alpha_{4} s_{31} s_{41} z_{1} + \alpha_{1} \alpha_{3} \alpha_{4} s_{41} z_{1} - \alpha_{1} \alpha_{3} \alpha_{5} s_{21} s_{31} + \alpha_{1} \alpha_{3} \alpha_{5} s_{31} + \alpha_{1} \alpha_{4} \alpha_{5} s_{21} s_{41} z_{1} + \alpha_{1} \alpha_{4} \alpha_{5} s_{21} s_{41} - \alpha_{1} \alpha_{4} \alpha_{5} s_{21} - \alpha_{1} \alpha_{4} \alpha_{5} s_{41} + \alpha_{1} \alpha_{4} \alpha_{5} - \alpha_{2} \alpha_{3} \alpha_{4} s_{31} s_{41} z_{1} + \alpha_{2} \alpha_{3} \alpha_{4} s_{41} z_{1} - \alpha_{2} \alpha_{3} \alpha_{5} s_{31} s_{42} z_{1} + \alpha_{2} \alpha_{3} \alpha_{5} s_{42} z_{1} - \alpha_{3} \alpha_{4} \alpha_{5} s_{31} s_{41} z_{1} + \alpha_{3} \alpha_{4} \alpha_{5} s_{41} z_{1}.$$

$$K_+=2 z_{1} \left(eps + 1\right) \left(z_{1} + 1\right),\qquad K_-=- z_{1} \left(2 eps - 7\right) \left(z_{1} + 1\right).$$

## 4e-3g+3g-16

外到内：`s21 → s32 → s31 → s41`。原界：`s21: 0…1`，`s31: 0…s32`，`s32: 0…1`，`s41: 0…1`。

$$\widehat B=- \alpha_{1} \alpha_{2} \alpha_{3} s_{21} s_{31} + \alpha_{1} \alpha_{2} \alpha_{3} s_{31} - \alpha_{1} \alpha_{2} \alpha_{4} s_{21} s_{31} + \alpha_{1} \alpha_{2} \alpha_{4} s_{31} - \alpha_{1} \alpha_{2} \alpha_{5} s_{21} s_{31} + \alpha_{1} \alpha_{2} \alpha_{5} s_{31} - \alpha_{1} \alpha_{3} \alpha_{4} s_{32} s_{41} z_{1} + \alpha_{1} \alpha_{3} \alpha_{4} s_{41} z_{1} - \alpha_{1} \alpha_{3} \alpha_{5} s_{21} s_{32} + \alpha_{1} \alpha_{3} \alpha_{5} s_{32} + \alpha_{1} \alpha_{4} \alpha_{5} s_{21} s_{41} z_{1} + \alpha_{1} \alpha_{4} \alpha_{5} s_{21} s_{41} - \alpha_{1} \alpha_{4} \alpha_{5} s_{21} - \alpha_{1} \alpha_{4} \alpha_{5} s_{41} + \alpha_{1} \alpha_{4} \alpha_{5} - \alpha_{2} \alpha_{3} \alpha_{4} s_{32} s_{41} z_{1} + \alpha_{2} \alpha_{3} \alpha_{4} s_{41} z_{1} - \alpha_{2} \alpha_{4} \alpha_{5} s_{31} s_{41} z_{1} + \alpha_{2} \alpha_{4} \alpha_{5} s_{41} z_{1} - \alpha_{3} \alpha_{4} \alpha_{5} s_{32} s_{41} z_{1} + \alpha_{3} \alpha_{4} \alpha_{5} s_{41} z_{1}.$$

$$K_+=- z_{1} \left(2 eps - 7\right),\qquad K_-=2 z_{1} \left(eps + 1\right).$$

## 4e-3g+3g-17

外到内：`s21 → s32 → s31 → s41`。原界：`s21: 0…1`，`s32: 0…1`，`s31: 0…s32`，`s41: 0…1`。

$C=0$，颜色筛选后的有效 $J=0$。未求积的几何多项式另保留在JSON。

## 4e-3g+3g-18

外到内：`s21 → s32 → s31 → s41`。原界：`s21: 0…1`，`s41: 0…1`，`s31: 0…s32`，`s32: 0…1`。

$$\widehat B=\alpha_{1} \alpha_{2} \alpha_{3} s_{21} s_{41} z_{1} + \alpha_{1} \alpha_{2} \alpha_{3} s_{21} s_{41} - \alpha_{1} \alpha_{2} \alpha_{3} s_{21} - \alpha_{1} \alpha_{2} \alpha_{3} s_{41} + \alpha_{1} \alpha_{2} \alpha_{3} + \alpha_{1} \alpha_{2} \alpha_{4} s_{21} s_{41} z_{1} + \alpha_{1} \alpha_{2} \alpha_{4} s_{21} s_{41} - \alpha_{1} \alpha_{2} \alpha_{4} s_{21} - \alpha_{1} \alpha_{2} \alpha_{4} s_{41} + \alpha_{1} \alpha_{2} \alpha_{4} + \alpha_{1} \alpha_{2} \alpha_{5} s_{21} s_{41} z_{1} + \alpha_{1} \alpha_{2} \alpha_{5} s_{21} s_{41} - \alpha_{1} \alpha_{2} \alpha_{5} s_{21} - \alpha_{1} \alpha_{2} \alpha_{5} s_{41} + \alpha_{1} \alpha_{2} \alpha_{5} - \alpha_{1} \alpha_{3} \alpha_{5} s_{21} s_{31} + \alpha_{1} \alpha_{3} \alpha_{5} s_{31} - \alpha_{1} \alpha_{4} \alpha_{5} s_{21} s_{32} + \alpha_{1} \alpha_{4} \alpha_{5} s_{32} - \alpha_{2} \alpha_{3} \alpha_{5} s_{31} s_{41} z_{1} + \alpha_{2} \alpha_{3} \alpha_{5} s_{41} z_{1} - \alpha_{2} \alpha_{4} \alpha_{5} s_{32} s_{41} z_{1} + \alpha_{2} \alpha_{4} \alpha_{5} s_{41} z_{1}.$$

$$K_+=z_{1},\qquad K_-=- z_{1}.$$

## 4e-3g+3g-19

外到内：`s22 → s21 → s42 → s41`。原界：`s21: 0…s22`，`s22: 0…1`，`s41: 0…s42`，`s42: 0…1`。

$$\widehat B=\alpha_{5} \left(\alpha_{1} \alpha_{3} s_{21} s_{41} z_{1} + \alpha_{1} \alpha_{3} s_{21} s_{41} - \alpha_{1} \alpha_{3} s_{21} - \alpha_{1} \alpha_{3} s_{41} + \alpha_{1} \alpha_{3} + \alpha_{1} \alpha_{4} s_{21} s_{42} z_{1} + \alpha_{1} \alpha_{4} s_{21} s_{42} - \alpha_{1} \alpha_{4} s_{21} - \alpha_{1} \alpha_{4} s_{42} + \alpha_{1} \alpha_{4} + \alpha_{2} \alpha_{3} s_{22} s_{41} z_{1} + \alpha_{2} \alpha_{3} s_{22} s_{41} - \alpha_{2} \alpha_{3} s_{22} - \alpha_{2} \alpha_{3} s_{41} + \alpha_{2} \alpha_{3} + \alpha_{2} \alpha_{4} s_{22} s_{42} z_{1} + \alpha_{2} \alpha_{4} s_{22} s_{42} - \alpha_{2} \alpha_{4} s_{22} - \alpha_{2} \alpha_{4} s_{42} + \alpha_{2} \alpha_{4}\right).$$

$$K_+=\left(z_{1} + 1\right)^{2},\qquad K_-=- \left(z_{1} + 1\right)^{2}.$$

## 4e-3g+3g-20

外到内：`s22 → s21 → s42 → s41`。原界：`s21: 0…s22`，`s41: 0…s42`，`s22: 0…1`，`s42: 0…1`。

$C=0$，颜色筛选后的有效 $J=0$。未求积的几何多项式另保留在JSON。

## 4e-3g+3g-21

外到内：`s22 → s21 → s42 → s41`。原界：`s21: 0…s22`，`s42: 0…1`，`s22: 0…1`，`s41: 0…s42`。

$$\widehat B=\alpha_{1} \alpha_{2} \alpha_{3} s_{21} s_{42} z_{1} + \alpha_{1} \alpha_{2} \alpha_{3} s_{21} s_{42} - \alpha_{1} \alpha_{2} \alpha_{3} s_{21} - \alpha_{1} \alpha_{2} \alpha_{3} s_{42} + \alpha_{1} \alpha_{2} \alpha_{3} + \alpha_{1} \alpha_{2} \alpha_{4} s_{21} s_{42} z_{1} + \alpha_{1} \alpha_{2} \alpha_{4} s_{21} s_{42} - \alpha_{1} \alpha_{2} \alpha_{4} s_{21} - \alpha_{1} \alpha_{2} \alpha_{4} s_{42} + \alpha_{1} \alpha_{2} \alpha_{4} + \alpha_{1} \alpha_{2} \alpha_{5} s_{21} s_{42} z_{1} + \alpha_{1} \alpha_{2} \alpha_{5} s_{21} s_{42} - \alpha_{1} \alpha_{2} \alpha_{5} s_{21} - \alpha_{1} \alpha_{2} \alpha_{5} s_{42} + \alpha_{1} \alpha_{2} \alpha_{5} + \alpha_{1} \alpha_{3} \alpha_{4} s_{22} s_{41} z_{1} + \alpha_{1} \alpha_{3} \alpha_{4} s_{22} s_{41} - \alpha_{1} \alpha_{3} \alpha_{4} s_{22} - \alpha_{1} \alpha_{3} \alpha_{4} s_{41} + \alpha_{1} \alpha_{3} \alpha_{4} + \alpha_{1} \alpha_{4} \alpha_{5} s_{21} s_{41} z_{1} + \alpha_{1} \alpha_{4} \alpha_{5} s_{21} s_{41} - \alpha_{1} \alpha_{4} \alpha_{5} s_{21} - \alpha_{1} \alpha_{4} \alpha_{5} s_{41} + \alpha_{1} \alpha_{4} \alpha_{5} + \alpha_{2} \alpha_{3} \alpha_{4} s_{22} s_{41} z_{1} + \alpha_{2} \alpha_{3} \alpha_{4} s_{22} s_{41} - \alpha_{2} \alpha_{3} \alpha_{4} s_{22} - \alpha_{2} \alpha_{3} \alpha_{4} s_{41} + \alpha_{2} \alpha_{3} \alpha_{4} + \alpha_{2} \alpha_{3} \alpha_{5} s_{22} s_{42} z_{1} + \alpha_{2} \alpha_{3} \alpha_{5} s_{22} s_{42} - \alpha_{2} \alpha_{3} \alpha_{5} s_{22} - \alpha_{2} \alpha_{3} \alpha_{5} s_{42} + \alpha_{2} \alpha_{3} \alpha_{5} + \alpha_{3} \alpha_{4} \alpha_{5} s_{22} s_{41} z_{1} + \alpha_{3} \alpha_{4} \alpha_{5} s_{22} s_{41} - \alpha_{3} \alpha_{4} \alpha_{5} s_{22} - \alpha_{3} \alpha_{4} \alpha_{5} s_{41} + \alpha_{3} \alpha_{4} \alpha_{5}.$$

$$K_+=- 2 \left(eps + 1\right) \left(z_{1} + 1\right)^{2},\qquad K_-=\left(2 eps - 7\right) \left(z_{1} + 1\right)^{2}.$$

## 4e-3g+3g-22

外到内：`s22 → s21 → s31 → s41`。原界：`s21: 0…s22`，`s22: 0…1`，`s31: 0…1`，`s41: 0…1`。

$$\widehat B=- \alpha_{1} \alpha_{3} \alpha_{4} s_{31} s_{41} z_{1} + \alpha_{1} \alpha_{3} \alpha_{4} s_{41} z_{1} - \alpha_{1} \alpha_{3} \alpha_{5} s_{21} s_{31} + \alpha_{1} \alpha_{3} \alpha_{5} s_{31} + \alpha_{1} \alpha_{4} \alpha_{5} s_{21} s_{41} z_{1} + \alpha_{1} \alpha_{4} \alpha_{5} s_{21} s_{41} - \alpha_{1} \alpha_{4} \alpha_{5} s_{21} - \alpha_{1} \alpha_{4} \alpha_{5} s_{41} + \alpha_{1} \alpha_{4} \alpha_{5} - \alpha_{2} \alpha_{3} \alpha_{4} s_{31} s_{41} z_{1} + \alpha_{2} \alpha_{3} \alpha_{4} s_{41} z_{1} - \alpha_{2} \alpha_{3} \alpha_{5} s_{22} s_{31} + \alpha_{2} \alpha_{3} \alpha_{5} s_{31} + \alpha_{2} \alpha_{4} \alpha_{5} s_{22} s_{41} z_{1} + \alpha_{2} \alpha_{4} \alpha_{5} s_{22} s_{41} - \alpha_{2} \alpha_{4} \alpha_{5} s_{22} - \alpha_{2} \alpha_{4} \alpha_{5} s_{41} + \alpha_{2} \alpha_{4} \alpha_{5} - \alpha_{3} \alpha_{4} \alpha_{5} s_{31} s_{41} z_{1} + \alpha_{3} \alpha_{4} \alpha_{5} s_{41} z_{1}.$$

$$K_+=- z_{1} - 1,\qquad K_-=z_{1} + 1.$$

## 4e-3g+3g-23

外到内：`s22 → s21 → s31 → s41`。原界：`s21: 0…s22`，`s31: 0…1`，`s22: 0…1`，`s41: 0…1`。

$C=0$，颜色筛选后的有效 $J=0$。未求积的几何多项式另保留在JSON。

## 4e-3g+3g-24

外到内：`s22 → s21 → s31 → s41`。原界：`s21: 0…s22`，`s41: 0…1`，`s22: 0…1`，`s31: 0…1`。

$$\widehat B=\alpha_{1} \alpha_{2} \alpha_{3} s_{21} s_{41} z_{1} + \alpha_{1} \alpha_{2} \alpha_{3} s_{21} s_{41} - \alpha_{1} \alpha_{2} \alpha_{3} s_{21} - \alpha_{1} \alpha_{2} \alpha_{3} s_{41} + \alpha_{1} \alpha_{2} \alpha_{3} + \alpha_{1} \alpha_{2} \alpha_{4} s_{21} s_{41} z_{1} + \alpha_{1} \alpha_{2} \alpha_{4} s_{21} s_{41} - \alpha_{1} \alpha_{2} \alpha_{4} s_{21} - \alpha_{1} \alpha_{2} \alpha_{4} s_{41} + \alpha_{1} \alpha_{2} \alpha_{4} + \alpha_{1} \alpha_{2} \alpha_{5} s_{21} s_{41} z_{1} + \alpha_{1} \alpha_{2} \alpha_{5} s_{21} s_{41} - \alpha_{1} \alpha_{2} \alpha_{5} s_{21} - \alpha_{1} \alpha_{2} \alpha_{5} s_{41} + \alpha_{1} \alpha_{2} \alpha_{5} - \alpha_{1} \alpha_{3} \alpha_{4} s_{22} s_{31} + \alpha_{1} \alpha_{3} \alpha_{4} s_{31} - \alpha_{1} \alpha_{4} \alpha_{5} s_{21} s_{31} + \alpha_{1} \alpha_{4} \alpha_{5} s_{31} - \alpha_{2} \alpha_{3} \alpha_{4} s_{22} s_{31} + \alpha_{2} \alpha_{3} \alpha_{4} s_{31} + \alpha_{2} \alpha_{3} \alpha_{5} s_{22} s_{41} z_{1} + \alpha_{2} \alpha_{3} \alpha_{5} s_{22} s_{41} - \alpha_{2} \alpha_{3} \alpha_{5} s_{22} - \alpha_{2} \alpha_{3} \alpha_{5} s_{41} + \alpha_{2} \alpha_{3} \alpha_{5} - \alpha_{2} \alpha_{4} \alpha_{5} s_{31} s_{41} z_{1} + \alpha_{2} \alpha_{4} \alpha_{5} s_{41} z_{1} - \alpha_{3} \alpha_{4} \alpha_{5} s_{22} s_{31} + \alpha_{3} \alpha_{4} \alpha_{5} s_{31}.$$

$$K_+=2 \left(eps + 1\right) \left(z_{1} + 1\right),\qquad K_-=- \left(2 eps - 7\right) \left(z_{1} + 1\right).$$

## 4e-3g+3g-25

外到内：`s11 → s21 → s31 → s41`。原界：`s11: 0…1`，`s21: 0…1`，`s31: 0…1`，`s41: 0…1`。

$$\widehat B=- \alpha_{1} \alpha_{2} \alpha_{3} s_{11} s_{21} z_{1} + \alpha_{1} \alpha_{2} \alpha_{3} s_{21} z_{1} - \alpha_{1} \alpha_{2} \alpha_{4} s_{11} s_{21} z_{1} + \alpha_{1} \alpha_{2} \alpha_{4} s_{21} z_{1} - \alpha_{1} \alpha_{2} \alpha_{5} s_{11} s_{21} z_{1} + \alpha_{1} \alpha_{2} \alpha_{5} s_{21} z_{1} - \alpha_{1} \alpha_{3} \alpha_{4} s_{31} s_{41} z_{1} + \alpha_{1} \alpha_{3} \alpha_{4} s_{41} z_{1} + \alpha_{1} \alpha_{3} \alpha_{5} s_{11} s_{31} z_{1} + \alpha_{1} \alpha_{3} \alpha_{5} s_{11} s_{31} - \alpha_{1} \alpha_{3} \alpha_{5} s_{11} z_{1} - \alpha_{1} \alpha_{3} \alpha_{5} s_{31} z_{1} + \alpha_{1} \alpha_{3} \alpha_{5} z_{1} - \alpha_{1} \alpha_{4} \alpha_{5} s_{11} s_{41} + \alpha_{1} \alpha_{4} \alpha_{5} s_{11} - \alpha_{2} \alpha_{3} \alpha_{4} s_{31} s_{41} z_{1} + \alpha_{2} \alpha_{3} \alpha_{4} s_{41} z_{1} - \alpha_{2} \alpha_{3} \alpha_{5} s_{21} s_{31} + \alpha_{2} \alpha_{3} \alpha_{5} s_{31} + \alpha_{2} \alpha_{4} \alpha_{5} s_{21} s_{41} z_{1} + \alpha_{2} \alpha_{4} \alpha_{5} s_{21} s_{41} - \alpha_{2} \alpha_{4} \alpha_{5} s_{21} - \alpha_{2} \alpha_{4} \alpha_{5} s_{41} + \alpha_{2} \alpha_{4} \alpha_{5} - \alpha_{3} \alpha_{4} \alpha_{5} s_{31} s_{41} z_{1} + \alpha_{3} \alpha_{4} \alpha_{5} s_{41} z_{1}.$$

$$K_+=- 2 eps z_{1}^{2} + 7 z_{1}^{2} + 10 z_{1} + 1,\qquad K_-=2 eps z_{1}^{2} + 2 z_{1}^{2} + 8 z_{1} - 1.$$

## 4e-3g+3g-26

外到内：`s11 → s21 → s31 → s41`。原界：`s11: 0…1`，`s31: 0…1`，`s21: 0…1`，`s41: 0…1`。

$C=0$，颜色筛选后的有效 $J=0$。未求积的几何多项式另保留在JSON。

## 使用与状态

从montopar目录可用以下命令重新生成实际输出和汇总包：

~~~powershell
D:/py/python.exe -m montopar.cli --source ../qgraf-pipeline/outputs/3loop/orbit_definitions.json --type "4e-3g+3g" --package-initial
~~~

标准安装后也可直接使用`montopar`命令。原参数是`J_raw_s`，可选超立方是`J_raw`；α单纯形、变量顺序、颜色有向配对、旋转重数和分支处方全部自包含。C组件严格相乘等于combined，不附加I^6或v14²假账本。

本版为原始参数表示，不是IBP约化式或积分解析解；formal=not_run。正确性核验包含全部26图与独立顶点/二阶算子的对照。未修改旧的审计快照和积分结果。
