# P5：07/08完整I及有限项（现行r002）

现行数学约定为[自能2.0](../../../conventions.md)，输入和结果见[当前清单](../../../current.json)。本页旧过程的原字节保存在[迁移前快照](../../../../../workflow/checks/self-energy-signs/history/benchmarks/self-energy/pi1/families/p5-opposite-same/solutions.md)。本项目在导入边界完成规范化，外部montopar/qgraf源码及输出仍保持原版。


以下I已采用物理1PI插入约定，组装系数为+1。G、U和各母函数保持原定义；J及有序域不变。系数表覆盖ε⁻⁶至ε⁰，未显示的位置为零，余项O(ε)。

## 4e-Pi1-07

[机器结果](../../results/4e-Pi1-07/r002.json)；[现行prepared](../../prepared/4e-Pi1-r002.json)。

$$
I=- \frac{\left(Tmixed_{sum} + \frac{H_{03} - 2 H2m \sinh{\left(\varepsilon \log{\left(z \right)} \right)} + \left(H_{01} - 2 H2p\right) \cosh{\left(\varepsilon \log{\left(z \right)} \right)}}{2 \varepsilon \left(\varepsilon + 1\right)}\right) e^{- 3 \gamma_E \varepsilon} \Gamma^{3}\left(1 - \varepsilon\right)}{128 \varepsilon \left(1 - 2 \varepsilon\right) \cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)}}
$$

上式缩写及显式已积分级数沿用以下绑定：

```text
Tmixed_sum = T(2eps,eps,z)+T(eps,2eps,z)
H03 = H0(3eps,z)
H01 = H0(eps,z)
H2p = (H2(eps,z)+H2(eps,1/z))/2
H2m = (H2(eps,z)-H2(eps,1/z))/2
T(a,b,z) = ((1+z)/sqrt(z))^(a+b)*Sum[(1-a)_k*(1-b)_l/(k!*l!)*B(k,l,z),{k,0,infinity},{l,0,infinity}]
B(k,l,z) = Sum[p_k(i,j;z)*p_l(m,n;z)/((i+1)*(i+m+2)*(j+1)*(j+n+2)),{i,j=0..k},{m,n=0..l}]
p_k(i,j;z) = Sum[k!*(-1)^h*z^(k-i-j+h)/((k-i-j+h)!*(i-h)!*(j-h)!*h!*(1+z)^(k-h)),{h,max(0,i+j-k),min(i,j)}]
order = each finite B block first, then the double sum; no arbitrary rearrangement of signed monomials
domain = z>0, initially 0<Re(eps)<1/2 with a=2eps,b=eps; locally normally convergent near eps=0
```

| ε幂次 | I 的系数 |
|---:|---|
| -1 | $- \frac{M_{01}^{2}}{384} + \frac{M_{11}}{256}$ |
| 0 | $- \frac{L^{2} M_{01}}{512} - \frac{M_{01}^{2}}{192} + \frac{M_{03}}{384} + \frac{M_{11}}{256} + \frac{M_{21}}{256} - \frac{3 \rho_{1}}{128}$ |

## 4e-Pi1-08

[机器结果](../../results/4e-Pi1-08/r002.json)；[现行prepared](../../prepared/4e-Pi1-r002.json)。

$$
I=- \frac{\left(Tmixed_{sum} + \frac{H_{03} - 2 H2m \sinh{\left(\varepsilon \log{\left(z \right)} \right)} + \left(H_{01} - 2 H2p\right) \cosh{\left(\varepsilon \log{\left(z \right)} \right)}}{2 \varepsilon \left(\varepsilon + 1\right)}\right) e^{- 3 \gamma_E \varepsilon} \Gamma^{3}\left(1 - \varepsilon\right)}{128 \varepsilon \left(1 - 2 \varepsilon\right) \cosh{\left(\frac{3 \varepsilon \log{\left(z \right)}}{2} \right)}}
$$

上式缩写及显式已积分级数沿用以下绑定：

```text
Tmixed_sum = T(2eps,eps,z)+T(eps,2eps,z)
H03 = H0(3eps,z)
H01 = H0(eps,z)
H2p = (H2(eps,z)+H2(eps,1/z))/2
H2m = (H2(eps,z)-H2(eps,1/z))/2
T(a,b,z) = ((1+z)/sqrt(z))^(a+b)*Sum[(1-a)_k*(1-b)_l/(k!*l!)*B(k,l,z),{k,0,infinity},{l,0,infinity}]
B(k,l,z) = Sum[p_k(i,j;z)*p_l(m,n;z)/((i+1)*(i+m+2)*(j+1)*(j+n+2)),{i,j=0..k},{m,n=0..l}]
p_k(i,j;z) = Sum[k!*(-1)^h*z^(k-i-j+h)/((k-i-j+h)!*(i-h)!*(j-h)!*h!*(1+z)^(k-h)),{h,max(0,i+j-k),min(i,j)}]
order = each finite B block first, then the double sum; no arbitrary rearrangement of signed monomials
domain = z>0, initially 0<Re(eps)<1/2 with a=2eps,b=eps; locally normally convergent near eps=0
```

| ε幂次 | I 的系数 |
|---:|---|
| -1 | $- \frac{M_{01}^{2}}{384} + \frac{M_{11}}{256}$ |
| 0 | $- \frac{L^{2} M_{01}}{512} - \frac{M_{01}^{2}}{192} + \frac{M_{03}}{384} + \frac{M_{11}}{256} + \frac{M_{21}}{256} - \frac{3 \rho_{1}}{128}$ |

由旧定义到现行定义的线性变换、颜色和系数连接已核对；未重新积分，不继承新形式证明。formal=not_run。
