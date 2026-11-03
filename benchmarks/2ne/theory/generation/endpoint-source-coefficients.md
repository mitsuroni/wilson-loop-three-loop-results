# 带端点源四图：本轮重建的完整系数

本表由[源积分推导](endpoint-sources.md)的新鲜计算生成。沿用唯一C/I归一化，所有旧结果保持不变。下列系数之外，既定 ε⁻⁶…ε⁰ 范围内的系数全部为零。

M/N与τ、E的定义见[统一权重六基底](../unified-weight-six.md)。本表的τ、E均由本轮源矩重新积分后与旧定义核对。坐标约化在积分生成之后执行，不作为源的定义或全阶闭包证明。

## 6e-10

$$
\begin{aligned}
I_{6e-10}^{(-1)}= {}&- \frac{M_{0,1} \zeta\left(3\right)}{8} - \frac{M_{1,3}}{96}.
\end{aligned}
$$

$$
\begin{aligned}
I_{6e-10}^{(0)}= {}&- \frac{E}{12} + \frac{263 M_{0,1}^{3}}{34560} + \frac{\tau^{+}_{D,4}}{24}\\
&{}+ \frac{47 \pi^{2} M_{0,1}^{2}}{34560} + \frac{M_{0,1} M_{2,1}}{24} + \frac{59 \pi^{4} M_{0,1}}{8640}\\
&{}- \frac{5 M_{1,1}^{2}}{192} - \frac{7 M_{1,1} \zeta\left(3\right)}{16} + \frac{\pi^{2} M_{2,1}}{96}\\
&{}- \frac{13 M_{2,3}}{1536} - \frac{9 M_{4,1}}{1024} - \frac{N_{2,2}}{64}.
\end{aligned}
$$

## 6e-11

$$
\begin{aligned}
I_{6e-11}^{(-1)}= {}&\frac{M_{0,1} M_{1,1}}{64} + \frac{\pi^{2} M_{1,1}}{96} + \frac{M_{1,3}}{96}\\
&{}- \frac{N_{1,2}}{64}.
\end{aligned}
$$

$$
\begin{aligned}
I_{6e-11}^{(0)}= {}&\frac{7 E}{96} - \frac{601 M_{0,1}^{3}}{138240} + \frac{\tau^{+}_{D,4}}{96}\\
&{}- \frac{139 \pi^{2} M_{0,1}^{2}}{138240} - \frac{13 M_{0,1} M_{2,1}}{768} - \frac{\pi^{4} M_{0,1}}{864}\\
&{}+ \frac{11 M_{1,1}^{2}}{384} + \frac{M_{1,1} \zeta\left(3\right)}{16} + \frac{\pi^{2} M_{2,1}}{192}\\
&{}+ \frac{47 M_{2,3}}{6144} - \frac{9 M_{4,1}}{4096} - \frac{N_{2,2}}{128}.
\end{aligned}
$$

## 6e-12

$$
\begin{aligned}
I_{6e-12}^{(-1)}= {}&- \frac{M_{0,1} M_{1,1}}{64} + \frac{M_{0,1} \zeta\left(3\right)}{16} - \frac{\pi^{2} M_{1,1}}{48}\\
&{}- \frac{M_{1,3}}{64} + \frac{N_{1,2}}{32}.
\end{aligned}
$$

$$
\begin{aligned}
I_{6e-12}^{(0)}= {}&\frac{E}{96} + \frac{77 M_{0,1}^{3}}{138240} - \frac{5 \tau^{+}_{D,4}}{96}\\
&{}+ \frac{623 \pi^{2} M_{0,1}^{2}}{138240} - \frac{19 M_{0,1} M_{2,1}}{768} + \frac{\pi^{4} M_{0,1}}{17280}\\
&{}+ \frac{5 M_{1,1}^{2}}{384} + \frac{M_{1,1} \zeta\left(3\right)}{32} - \frac{5 \pi^{2} M_{2,1}}{384}\\
&{}- \frac{67 M_{2,3}}{6144} + \frac{45 M_{4,1}}{4096} + \frac{5 N_{2,2}}{256}.
\end{aligned}
$$

## 6e-36

$$
\begin{aligned}
I_{6e-36}^{(-2)}= {}&\frac{\pi^{2} M_{0,1}}{384} - \frac{M_{2,1}}{128}.
\end{aligned}
$$

$$
\begin{aligned}
I_{6e-36}^{(-1)}= {}&- \frac{M_{0,1} M_{1,1}}{64} + \frac{M_{0,1} \zeta\left(3\right)}{8} - \frac{7 \pi^{2} M_{1,1}}{384}\\
&{}+ \frac{M_{1,3}}{768} - \frac{M_{3,1}}{128} + \frac{3 N_{1,2}}{128}.
\end{aligned}
$$

$$
\begin{aligned}
I_{6e-36}^{(0)}= {}&- \frac{43 M_{0,1}^{3}}{46080} + \frac{91 \pi^{2} M_{0,1}^{2}}{15360} - \frac{\tau^{+}_{D,4}}{32}\\
&{}- \frac{17 M_{0,1} M_{2,1}}{512} + \frac{401 \pi^{4} M_{0,1}}{92160} - \frac{19 \pi^{2} M_{2,1}}{1024}\\
&{}- \frac{29 M_{2,3}}{6144} + \frac{31 M_{4,1}}{12288} + \frac{5 N_{2,2}}{512}.
\end{aligned}
$$

完整Hlog、源矩、精确残差及不相容见证见[机器记录](../../../../workflow/checks/2ne-endpoint-sources/coefficients.json)。formal=not_run。
