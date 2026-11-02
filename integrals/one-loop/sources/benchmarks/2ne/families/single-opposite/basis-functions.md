# [1] 族实际系数函数

定义和表示选择见 [basis.md](basis.md)。以下全由本批重新计算的 HyperInt 输出生成；tau_(0/1/2,q) 保留旧阶乘规范，其他 tau/eta 的精确含义见该文 §6。
Hlog 为标准 Goncharov 迭代积分，字母为 0、−1；这里给出实际公式，不只留下系数提取符号。

## tau0q0

来源：H0，q=0，通道=plus，自然权重=2。

$$
\tau_{0,0}=- \frac{\pi^{2}}{2} - \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}
$$

## tau0q1

来源：H0，q=1，通道=plus，自然权重=3。

$$
\tau_{0,1}=\zeta\left(3\right) + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{2} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{4} + \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)} - \frac{\operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{2}
$$

## tau0q2

来源：H0，q=2，通道=plus，自然权重=4。

$$
\tau_{0,2}=- \frac{\pi^{4}}{30} - \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)} + \frac{\zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{2} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)}}{2} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{4} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{4} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{8} - \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)} + \frac{\operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)}}{2} + \frac{\operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)}}{2} - \frac{\operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{2}
$$

## tau0q3

来源：H0，q=3，通道=plus，自然权重=5。

$$
\tau_{0,3}=\zeta\left(5\right) + \frac{\pi^{4} \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{30} - \frac{\pi^{4} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{60} + \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)} - \frac{\zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{2} - \frac{\zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{2} + \frac{\zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{4} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1\right) \right)}}{2} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0\right) \right)}}{4} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1\right) \right)}}{4} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)}}{8} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1\right) \right)}}{4} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)}}{8} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)}}{8} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{16} + \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0, \  0\right) \right)} - \frac{\operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0\right) \right)}}{2} - \frac{\operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0\right) \right)}}{2} + \frac{\operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0\right) \right)}}{2} - \frac{\operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0\right) \right)}}{2} + \frac{\operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0\right) \right)}}{4} + \frac{\operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0\right) \right)}}{4} - \frac{\operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0\right) \right)}}{4}
$$

## tau0q4

来源：H0，q=4，通道=plus，自然权重=6。

$$
\tau_{0,4}=- \frac{\pi^{6}}{315} - \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)} + \frac{\zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{2} - \frac{\pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)}}{30} + \frac{\pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{60} + \frac{\pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{60} - \frac{\pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{120} - \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1\right) \right)} + \frac{\zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0\right) \right)}}{2} + \frac{\zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1\right) \right)}}{2} - \frac{\zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)}}{4} + \frac{\zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1\right) \right)}}{2} - \frac{\zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)}}{4} - \frac{\zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)}}{4} + \frac{\zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{8} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  -1\right) \right)}}{2} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0\right) \right)}}{4} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  -1\right) \right)}}{4} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)}}{8} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  -1\right) \right)}}{4} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0\right) \right)}}{8} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1\right) \right)}}{8} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)}}{16} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  -1\right) \right)}}{4} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0\right) \right)}}{8} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1\right) \right)}}{8} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)}}{16} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1\right) \right)}}{8} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0\right) \right)}}{16} + \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1\right) \right)}}{16} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{32} - \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  -1, \  0, \  0\right) \right)} + \frac{\operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0, \  0, \  0\right) \right)}}{2} + \frac{\operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  -1, \  0, \  0\right) \right)}}{2} - \frac{\operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0, \  0\right) \right)}}{2} + \frac{\operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  -1, \  0, \  0\right) \right)}}{2} - \frac{\operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0, \  0\right) \right)}}{4} - \frac{\operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1, \  0, \  0\right) \right)}}{4} + \frac{\operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0, \  0\right) \right)}}{4} + \frac{\operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  -1, \  0, \  0\right) \right)}}{2} - \frac{\operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0, \  0\right) \right)}}{4} - \frac{\operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1, \  0, \  0\right) \right)}}{4} + \frac{\operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0, \  0\right) \right)}}{4} - \frac{\operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1, \  0, \  0\right) \right)}}{4} + \frac{\operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0, \  0\right) \right)}}{8} + \frac{\operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1, \  0, \  0\right) \right)}}{8} - \frac{3 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0, \  0\right) \right)}}{16}
$$

## tau1q0

来源：H1，q=0，通道=plus，自然权重=2。

$$
\tau_{1,0}=- \frac{\pi^{2}}{2} - \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}
$$

## tau1q1

来源：H1，q=1，通道=plus，自然权重=3。

$$
\tau_{1,1}=2 \zeta\left(3\right) + \pi^{2} \operatorname{Hlog}{\left(z,\left( -1,\right) \right)} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{2} + 2 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)} - \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}
$$

## tau1q2

来源：H1，q=2，通道=plus，自然权重=4。

$$
\tau_{1,2}=- \frac{2 \pi^{4}}{15} - 4 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)} + 2 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)} - 2 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)} + \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)} + \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)} - \frac{11 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{24} - 4 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)} + 2 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)} + 2 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)} - \frac{3 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{2}
$$

## tau1q3

来源：H1，q=3，通道=plus，自然权重=5。

$$
\tau_{1,3}=8 \zeta\left(5\right) + \frac{4 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{15} - \frac{113 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{720} + 8 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)} - 4 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)} - 4 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)} + \frac{\zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{2} + 4 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1\right) \right)} - 2 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0\right) \right)} - 2 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1\right) \right)} + \frac{11 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)}}{12} - 2 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1\right) \right)} + \frac{13 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)}}{12} + \frac{5 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)}}{4} - \frac{5 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{8} + 8 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0, \  0\right) \right)} - 4 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0\right) \right)} - 4 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0\right) \right)} + 3 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0\right) \right)} - 4 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0\right) \right)} + \frac{5 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0\right) \right)}}{2} + \frac{5 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0\right) \right)}}{2} - 2 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0\right) \right)}
$$

## tau1q4

来源：H1，q=4，通道=plus，自然权重=6。

$$
\tau_{1,4}=- \frac{16 \pi^{6}}{315} - 16 \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)} - \frac{\pi^{2} \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{6} + 7 \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)} - \frac{8 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)}}{15} + \frac{113 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{360} + \frac{113 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{360} - \frac{17 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{90} - 16 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1\right) \right)} + 8 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0\right) \right)} + 8 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1\right) \right)} - \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)} + 8 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1\right) \right)} - 3 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)} - 5 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)} + \frac{\zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{2} - 8 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  -1\right) \right)} + 4 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0\right) \right)} + 4 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  -1\right) \right)} - \frac{11 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)}}{6} + 4 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  -1\right) \right)} - \frac{13 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0\right) \right)}}{6} - \frac{5 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1\right) \right)}}{2} + \frac{5 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)}}{4} + 4 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  -1\right) \right)} - \frac{13 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0\right) \right)}}{6} - \frac{5 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1\right) \right)}}{2} + \frac{5 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)}}{4} - \frac{5 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1\right) \right)}}{2} + \frac{17 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0\right) \right)}}{12} + \frac{7 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1\right) \right)}}{4} - \frac{83 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{96} - 16 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  -1, \  0, \  0\right) \right)} + 8 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0, \  0, \  0\right) \right)} + 8 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  -1, \  0, \  0\right) \right)} - 6 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0, \  0\right) \right)} + 8 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  -1, \  0, \  0\right) \right)} - 5 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0, \  0\right) \right)} - 5 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1, \  0, \  0\right) \right)} + 4 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0, \  0\right) \right)} + 8 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  -1, \  0, \  0\right) \right)} - 5 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0, \  0\right) \right)} - 5 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1, \  0, \  0\right) \right)} + 4 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0, \  0\right) \right)} - 5 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1, \  0, \  0\right) \right)} + \frac{7 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0, \  0\right) \right)}}{2} + \frac{7 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1, \  0, \  0\right) \right)}}{2} - \frac{47 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0, \  0\right) \right)}}{16}
$$

## eta1q0

来源：H1，q=0，通道=ell_times_minus，自然权重=3。

$$
\eta_{1,0}=0
$$

## eta1q1

来源：H1，q=1，通道=ell_times_minus，自然权重=4。

$$
\eta_{1,1}=\frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{6} + 2 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}
$$

## eta1q2

来源：H1，q=2，通道=ell_times_minus，自然权重=5。

$$
\eta_{1,2}=\frac{17 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{360} + 2 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)} - \frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)}}{3} - \frac{2 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)}}{3} - \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)} + \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)} - 4 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0\right) \right)} - 4 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0\right) \right)} - 2 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0\right) \right)} + 5 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0\right) \right)}
$$

## eta1q3

来源：H1，q=3，通道=ell_times_minus，自然权重=6。

$$
\eta_{1,3}=2 \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)} + \frac{\pi^{2} \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{3} - \frac{17 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{180} - \frac{17 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{180} + \frac{8 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{45} - 4 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)} + 4 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)} + \frac{2 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)}}{3} + \frac{4 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0\right) \right)}}{3} + 2 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1\right) \right)} - 2 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)} + \frac{4 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0\right) \right)}}{3} + 2 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1\right) \right)} - 2 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)} + 2 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1\right) \right)} - \frac{7 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0\right) \right)}}{3} - 3 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1\right) \right)} + \frac{25 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{12} + 8 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0, \  0\right) \right)} + 8 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0, \  0\right) \right)} + 4 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1, \  0, \  0\right) \right)} - 10 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0, \  0\right) \right)} + 8 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0, \  0\right) \right)} + 4 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1, \  0, \  0\right) \right)} - 10 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0, \  0\right) \right)} + 4 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1, \  0, \  0\right) \right)} - 10 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0, \  0\right) \right)} - 6 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1, \  0, \  0\right) \right)} + \frac{21 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0, \  0\right) \right)}}{2}
$$

## tau2q0

来源：H2，q=0，通道=plus，自然权重=2。

$$
\tau_{2,0}=- \frac{\pi^{2}}{2} - \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}
$$

## tau2q1

来源：H2，q=1，通道=plus，自然权重=3。

$$
\tau_{2,1}=3 \zeta\left(3\right) + \frac{3 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{2} - \frac{3 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{4} + 3 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)} - \frac{3 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{2}
$$

## tau2q2

来源：H2，q=2，通道=plus，自然权重=4。

$$
\tau_{2,2}=- \frac{41 \pi^{4}}{120} - 9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)} + \frac{9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{2} - \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)}}{2} + \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{4} + \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{4} - \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{8} - 9 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)} + \frac{9 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)}}{2} + \frac{9 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)}}{2} - \frac{7 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{2}
$$

## tau2q3

来源：H2，q=3，通道=plus，自然权重=5。

$$
\tau_{2,3}=\frac{\pi^{2} \zeta\left(3\right)}{2} + 30 \zeta\left(5\right) + \frac{41 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{40} - \frac{157 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{240} + 27 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)} - \frac{27 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{2} - \frac{27 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{2} + \frac{3 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{4} + \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1\right) \right)}}{2} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0\right) \right)}}{4} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1\right) \right)}}{4} + \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)}}{8} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1\right) \right)}}{4} + \frac{31 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)}}{8} + \frac{39 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)}}{8} - \frac{43 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{16} + 27 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0, \  0\right) \right)} - \frac{27 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0\right) \right)}}{2} - \frac{27 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0\right) \right)}}{2} + \frac{21 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0\right) \right)}}{2} - \frac{27 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0\right) \right)}}{2} + \frac{39 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0\right) \right)}}{4} + \frac{39 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0\right) \right)}}{4} - \frac{33 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0\right) \right)}}{4}
$$

## tau2q4

来源：H2，q=4，通道=plus，自然权重=6。

$$
\tau_{2,4}=- \frac{283 \pi^{6}}{840} - \frac{9 \zeta^{2}\left(3\right)}{2} - 90 \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)} - \frac{3 \pi^{2} \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{2} - \frac{3 \pi^{2} \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{4} + 36 \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)} - \frac{123 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)}}{40} + \frac{157 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{80} + \frac{157 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{80} - \frac{649 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{480} - 81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1\right) \right)} + \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0\right) \right)}}{2} + \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1\right) \right)}}{2} - \frac{9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)}}{4} + \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1\right) \right)}}{2} - \frac{45 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)}}{4} - \frac{117 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)}}{4} + \frac{9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{8} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  -1\right) \right)}}{2} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0\right) \right)}}{4} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  -1\right) \right)}}{4} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)}}{8} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  -1\right) \right)}}{4} - \frac{93 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0\right) \right)}}{8} - \frac{117 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1\right) \right)}}{8} + \frac{129 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)}}{16} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  -1\right) \right)}}{4} - \frac{93 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0\right) \right)}}{8} - \frac{117 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1\right) \right)}}{8} + \frac{129 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)}}{16} - \frac{117 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1\right) \right)}}{8} + \frac{141 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0\right) \right)}}{16} + \frac{189 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1\right) \right)}}{16} - \frac{209 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{32} - 81 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  -1, \  0, \  0\right) \right)} + \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0, \  0, \  0\right) \right)}}{2} + \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  -1, \  0, \  0\right) \right)}}{2} - \frac{63 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0, \  0\right) \right)}}{2} + \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  -1, \  0, \  0\right) \right)}}{2} - \frac{117 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0, \  0\right) \right)}}{4} - \frac{117 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1, \  0, \  0\right) \right)}}{4} + \frac{99 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0, \  0\right) \right)}}{4} + \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  -1, \  0, \  0\right) \right)}}{2} - \frac{117 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0, \  0\right) \right)}}{4} - \frac{117 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1, \  0, \  0\right) \right)}}{4} + \frac{99 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0, \  0\right) \right)}}{4} - \frac{117 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1, \  0, \  0\right) \right)}}{4} + \frac{189 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0, \  0\right) \right)}}{8} + \frac{189 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1, \  0, \  0\right) \right)}}{8} - \frac{327 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0, \  0\right) \right)}}{16}
$$

## eta2q0

来源：H2，q=0，通道=ell_times_minus，自然权重=3。

$$
\eta_{2,0}=0
$$

## eta2q1

来源：H2，q=1，通道=ell_times_minus，自然权重=4。

$$
\eta_{2,1}=\frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{3} + 4 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}
$$

## eta2q2

来源：H2，q=2，通道=ell_times_minus，自然权重=5。

$$
\eta_{2,2}=\frac{17 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{120} + 6 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)} - \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)} - 2 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)} - 3 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)} + 3 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)} - 12 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0\right) \right)} - 12 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0\right) \right)} - 6 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0\right) \right)} + 15 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0\right) \right)}
$$

## eta2q3

来源：H2，q=3，通道=ell_times_minus，自然权重=6。

$$
\eta_{2,3}=9 \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)} + \frac{3 \pi^{2} \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{2} - \frac{17 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{40} - \frac{17 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{40} + \frac{311 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{360} - 18 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)} + 18 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)} + 3 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)} + 6 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0\right) \right)} + 9 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1\right) \right)} - 9 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)} + 6 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0\right) \right)} + 9 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1\right) \right)} - 9 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)} + 9 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1\right) \right)} - \frac{21 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0\right) \right)}}{2} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1\right) \right)}}{2} + \frac{61 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{6} + 36 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0, \  0\right) \right)} + 36 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0, \  0\right) \right)} + 18 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1, \  0, \  0\right) \right)} - 45 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0, \  0\right) \right)} + 36 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0, \  0\right) \right)} + 18 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1, \  0, \  0\right) \right)} - 45 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0, \  0\right) \right)} + 18 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1, \  0, \  0\right) \right)} - 45 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0, \  0\right) \right)} - 27 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1, \  0, \  0\right) \right)} + 48 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0, \  0\right) \right)}
$$

## tauWq0

来源：HW，q=0，通道=plus，自然权重=2。

$$
\tau^{+}_{W,0}=- \frac{\pi^{2}}{2} - \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}
$$

## tauWq1

来源：HW，q=1，通道=plus，自然权重=3。

$$
\tau^{+}_{W,1}=3 \zeta\left(3\right) + \frac{3 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{2} - \frac{3 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{4} + 3 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)} - \frac{3 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{2}
$$

## tauWq2

来源：HW，q=2，通道=plus，自然权重=4。

$$
\tau^{+}_{W,2}=- \frac{13 \pi^{4}}{60} - 9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)} + \frac{9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{2} - \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)}}{2} + \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{4} + \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{4} - \frac{19 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{24} - 9 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)} + \frac{9 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)}}{2} + \frac{9 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)}}{2} - \frac{5 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{2}
$$

## tauWq3

来源：HW，q=3，通道=plus，自然权重=5。

$$
\tau^{+}_{W,3}=- \pi^{2} \zeta\left(3\right) + 21 \zeta\left(5\right) + \frac{13 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{20} - \frac{13 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{40} + 27 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)} - \frac{27 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{2} - \frac{27 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{2} + \frac{3 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{4} + \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1\right) \right)}}{2} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0\right) \right)}}{4} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1\right) \right)}}{4} + \frac{19 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)}}{8} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1\right) \right)}}{4} + \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)}}{8} + \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)}}{8} - \frac{19 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{16} + 27 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0, \  0\right) \right)} - \frac{27 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0\right) \right)}}{2} - \frac{27 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0\right) \right)}}{2} + \frac{15 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0\right) \right)}}{2} - \frac{27 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0\right) \right)}}{2} + \frac{27 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0\right) \right)}}{4} + \frac{27 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0\right) \right)}}{4} - \frac{15 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0\right) \right)}}{4}
$$

## tauWq4

来源：HW，q=4，通道=plus，自然权重=6。

$$
\tau^{+}_{W,4}=- \frac{41 \pi^{6}}{420} + 9 \zeta^{2}\left(3\right) - 63 \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)} + 3 \pi^{2} \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)} - \frac{3 \pi^{2} \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{2} + \frac{63 \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{2} - \frac{39 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)}}{20} + \frac{39 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{40} + \frac{39 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{40} - \frac{263 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{720} - 81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1\right) \right)} + \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0\right) \right)}}{2} + \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1\right) \right)}}{2} - \frac{9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)}}{4} + \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1\right) \right)}}{2} - \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)}}{4} - \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)}}{4} + \frac{9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{8} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  -1\right) \right)}}{2} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0\right) \right)}}{4} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  -1\right) \right)}}{4} - \frac{57 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)}}{8} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  -1\right) \right)}}{4} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0\right) \right)}}{8} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1\right) \right)}}{8} + \frac{57 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)}}{16} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  -1\right) \right)}}{4} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0\right) \right)}}{8} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1\right) \right)}}{8} + \frac{57 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)}}{16} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1\right) \right)}}{8} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0\right) \right)}}{16} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1\right) \right)}}{16} - \frac{163 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{96} - 81 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  -1, \  0, \  0\right) \right)} + \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0, \  0, \  0\right) \right)}}{2} + \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  -1, \  0, \  0\right) \right)}}{2} - \frac{45 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0, \  0\right) \right)}}{2} + \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  -1, \  0, \  0\right) \right)}}{2} - \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0, \  0\right) \right)}}{4} - \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1, \  0, \  0\right) \right)}}{4} + \frac{45 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0, \  0\right) \right)}}{4} + \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  -1, \  0, \  0\right) \right)}}{2} - \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0, \  0\right) \right)}}{4} - \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1, \  0, \  0\right) \right)}}{4} + \frac{45 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0, \  0\right) \right)}}{4} - \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1, \  0, \  0\right) \right)}}{4} + \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0, \  0\right) \right)}}{8} + \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1, \  0, \  0\right) \right)}}{8} - \frac{91 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0, \  0\right) \right)}}{16}
$$

## tauDq0

来源：HD，q=0，通道=plus，自然权重=2。

$$
\tau^{+}_{D,0}=- \frac{\pi^{2}}{2} - \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}
$$

## tauDq1

来源：HD，q=1，通道=plus，自然权重=3。

$$
\tau^{+}_{D,1}=3 \zeta\left(3\right) + \frac{3 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{2} - \frac{3 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{4} + 3 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)} - \frac{3 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{2}
$$

## tauDq2

来源：HD，q=2，通道=plus，自然权重=4。

$$
\tau^{+}_{D,2}=- \frac{19 \pi^{4}}{80} - 9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)} + \frac{9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{2} - \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)}}{2} + \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{4} + \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{4} - \frac{7 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{8} - 9 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)} + \frac{9 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)}}{2} + \frac{9 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)}}{2} - 3 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}
$$

## tauDq3

来源：HD，q=3，通道=plus，自然权重=5。

$$
\tau^{+}_{D,3}=\frac{7 \pi^{2} \zeta\left(3\right)}{12} + \frac{29 \zeta\left(5\right)}{2} + \frac{57 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{80} - \frac{89 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{288} + 27 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)} - \frac{27 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{2} - \frac{27 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{2} + \frac{9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{4} + \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1\right) \right)}}{2} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0\right) \right)}}{4} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1\right) \right)}}{4} + \frac{21 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)}}{8} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1\right) \right)}}{4} + \frac{89 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)}}{24} + \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)}}{8} - \frac{71 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{48} + 27 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0, \  0\right) \right)} - \frac{27 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0\right) \right)}}{2} - \frac{27 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0\right) \right)}}{2} + 9 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0\right) \right)} - \frac{27 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0\right) \right)}}{2} + \frac{35 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0\right) \right)}}{4} + \frac{27 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0\right) \right)}}{4} - \frac{11 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0\right) \right)}}{2}
$$

## tauDq4

来源：HD，q=4，通道=plus，自然权重=6。

$$
\tau^{+}_{D,4}=- \frac{811 \pi^{6}}{7560} - \frac{5 \zeta^{2}\left(3\right)}{4} - \frac{87 \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{2} - \frac{7 \pi^{2} \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{4} - \frac{5 \pi^{2} \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{8} + \frac{67 \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{4} - \frac{171 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)}}{80} + \frac{89 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{96} + \frac{581 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{480} - \frac{85 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{192} - 81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1\right) \right)} + \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0\right) \right)}}{2} + \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1\right) \right)}}{2} - \frac{27 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)}}{4} + \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1\right) \right)}}{2} - \frac{73 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)}}{4} - \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)}}{4} + \frac{19 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{8} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  -1\right) \right)}}{2} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0\right) \right)}}{4} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  -1\right) \right)}}{4} - \frac{63 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)}}{8} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  -1\right) \right)}}{4} - \frac{89 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0\right) \right)}}{8} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1\right) \right)}}{8} + \frac{71 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)}}{16} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  -1\right) \right)}}{4} - \frac{89 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0\right) \right)}}{8} - \frac{97 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1\right) \right)}}{8} + \frac{87 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)}}{16} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1\right) \right)}}{8} + \frac{97 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0\right) \right)}}{16} + \frac{97 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1\right) \right)}}{16} - \frac{77 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{32} - 81 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  -1, \  0, \  0\right) \right)} + \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0, \  0, \  0\right) \right)}}{2} + \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  -1, \  0, \  0\right) \right)}}{2} - 27 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0, \  0\right) \right)} + \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  -1, \  0, \  0\right) \right)}}{2} - \frac{105 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0, \  0\right) \right)}}{4} - \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1, \  0, \  0\right) \right)}}{4} + \frac{33 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0, \  0\right) \right)}}{2} + \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  -1, \  0, \  0\right) \right)}}{2} - \frac{105 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0, \  0\right) \right)}}{4} - \frac{97 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1, \  0, \  0\right) \right)}}{4} + \frac{37 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0, \  0\right) \right)}}{2} - \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1, \  0, \  0\right) \right)}}{4} + \frac{129 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0, \  0\right) \right)}}{8} + \frac{97 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1, \  0, \  0\right) \right)}}{8} - \frac{161 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0, \  0\right) \right)}}{16}
$$

## etaDq0

来源：HD，q=0，通道=ell_times_minus，自然权重=3。

$$
\eta_{D,0}=0
$$

## etaDq1

来源：HD，q=1，通道=ell_times_minus，自然权重=4。

$$
\eta_{D,1}=\frac{\pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{3} + 4 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}
$$

## etaDq2

来源：HD，q=2，通道=ell_times_minus，自然权重=5。

$$
\eta_{D,2}=\frac{17 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{240} + 2 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)} - \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)} - \frac{3 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)}}{2} - 2 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)} + \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{4} - 12 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0\right) \right)} - 9 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0\right) \right)} - 4 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0\right) \right)} + \frac{25 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0\right) \right)}}{2}
$$

## etaDq3

来源：HD，q=3，通道=ell_times_minus，自然权重=6。

$$
\eta_{D,3}=\frac{9 \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{2} + \frac{3 \pi^{2} \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{4} - \frac{17 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{80} - \frac{17 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{80} + \frac{293 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{720} - 6 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)} + 3 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)} + 12 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)} - \frac{9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{2} + 3 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)} + \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0\right) \right)}}{2} + 6 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1\right) \right)} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)}}{4} + \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0\right) \right)}}{2} + 6 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1\right) \right)} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)}}{4} + 6 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1\right) \right)} - \frac{15 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0\right) \right)}}{2} - 9 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1\right) \right)} + \frac{16 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{3} + 36 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0, \  0\right) \right)} + 27 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0, \  0\right) \right)} + 12 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1, \  0, \  0\right) \right)} - \frac{75 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0, \  0\right) \right)}}{2} + 27 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0, \  0\right) \right)} + 12 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1, \  0, \  0\right) \right)} - \frac{75 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0, \  0\right) \right)}}{2} + 12 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1, \  0, \  0\right) \right)} - 33 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0, \  0\right) \right)} - 18 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1, \  0, \  0\right) \right)} + \frac{63 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0, \  0\right) \right)}}{2}
$$

## tauXq0

来源：HX，q=0，通道=plus，自然权重=2。

$$
\tau^{+}_{X,0}=- \frac{\pi^{2}}{2} - \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}
$$

## tauXq1

来源：HX，q=1，通道=plus，自然权重=3。

$$
\tau^{+}_{X,1}=3 \zeta\left(3\right) + \frac{3 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{2} - \frac{3 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{4} + 3 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)} - \frac{3 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{2}
$$

## tauXq2

来源：HX，q=2，通道=plus，自然权重=4。

$$
\tau^{+}_{X,2}=- \frac{3 \pi^{4}}{10} - 9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)} + \frac{9 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{2} - \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)}}{2} + \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{4} + \frac{9 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{4} - \frac{23 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{24} - 9 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)} + \frac{9 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)}}{2} + \frac{9 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)}}{2} - \frac{5 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{2}
$$

## tauXq3

来源：HX，q=3，通道=plus，自然权重=5。

$$
\tau^{+}_{X,3}=\frac{3 \pi^{2} \zeta\left(3\right)}{2} + 21 \zeta\left(5\right) + \frac{9 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{10} - \frac{9 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{20} + 27 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)} - \frac{27 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{2} - \frac{27 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{2} + \frac{19 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{4} + \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1\right) \right)}}{2} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0\right) \right)}}{4} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1\right) \right)}}{4} + \frac{23 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)}}{8} - \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1\right) \right)}}{4} + \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)}}{8} + \frac{27 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)}}{8} - \frac{23 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{16} + 27 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0, \  0\right) \right)} - \frac{27 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0\right) \right)}}{2} - \frac{27 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0\right) \right)}}{2} + \frac{15 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0\right) \right)}}{2} - \frac{27 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0\right) \right)}}{2} + \frac{27 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0\right) \right)}}{4} + \frac{27 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0\right) \right)}}{4} - \frac{15 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0\right) \right)}}{4}
$$

## tauXq4

来源：HX，q=4，通道=plus，自然权重=6。

$$
\tau^{+}_{X,4}=- \frac{353 \pi^{6}}{1680} - 3 \zeta^{2}\left(3\right) - 63 \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)} - \frac{9 \pi^{2} \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1,\right) \right)}}{2} + \frac{9 \pi^{2} \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{4} + \frac{63 \zeta\left(5\right) \operatorname{Hlog}{\left(z,\left( 0,\right) \right)}}{2} - \frac{27 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  -1\right) \right)}}{10} + \frac{27 \pi^{4} \operatorname{Hlog}{\left(z,\left( -1, \  0\right) \right)}}{20} + \frac{27 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  -1\right) \right)}}{20} - \frac{13 \pi^{4} \operatorname{Hlog}{\left(z,\left( 0, \  0\right) \right)}}{20} - 81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1\right) \right)} + \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0\right) \right)}}{2} + \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1\right) \right)}}{2} - \frac{57 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0\right) \right)}}{4} + \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1\right) \right)}}{2} - \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0\right) \right)}}{4} - \frac{81 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1\right) \right)}}{4} + \frac{57 \zeta\left(3\right) \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0\right) \right)}}{8} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  -1\right) \right)}}{2} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0\right) \right)}}{4} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  -1\right) \right)}}{4} - \frac{69 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0\right) \right)}}{8} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  -1\right) \right)}}{4} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0\right) \right)}}{8} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1\right) \right)}}{8} + \frac{69 \pi^{2} \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0\right) \right)}}{16} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  -1\right) \right)}}{4} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0\right) \right)}}{8} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1\right) \right)}}{8} + \frac{69 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0\right) \right)}}{16} - \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1\right) \right)}}{8} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0\right) \right)}}{16} + \frac{81 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1\right) \right)}}{16} - \frac{203 \pi^{2} \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0\right) \right)}}{96} - 81 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  -1, \  0, \  0\right) \right)} + \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  -1, \  0, \  0, \  0\right) \right)}}{2} + \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  -1, \  0, \  0\right) \right)}}{2} - \frac{45 \operatorname{Hlog}{\left(z,\left( -1, \  -1, \  0, \  0, \  0, \  0\right) \right)}}{2} + \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  -1, \  0, \  0\right) \right)}}{2} - \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  -1, \  0, \  0, \  0\right) \right)}}{4} - \frac{81 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  -1, \  0, \  0\right) \right)}}{4} + \frac{45 \operatorname{Hlog}{\left(z,\left( -1, \  0, \  0, \  0, \  0, \  0\right) \right)}}{4} + \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  -1, \  0, \  0\right) \right)}}{2} - \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  -1, \  0, \  0, \  0\right) \right)}}{4} - \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  -1, \  0, \  0\right) \right)}}{4} + \frac{45 \operatorname{Hlog}{\left(z,\left( 0, \  -1, \  0, \  0, \  0, \  0\right) \right)}}{4} - \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  -1, \  0, \  0\right) \right)}}{4} + \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  -1, \  0, \  0, \  0\right) \right)}}{8} + \frac{81 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  -1, \  0, \  0\right) \right)}}{8} - \frac{91 \operatorname{Hlog}{\left(z,\left( 0, \  0, \  0, \  0, \  0, \  0\right) \right)}}{16}
$$

