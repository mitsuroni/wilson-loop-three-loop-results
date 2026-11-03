# [2] 七图：以统一母族系数表达 I

此表与[既有Σ结果](solutions.md)数学相同，使用从通用有序 J 导出的母族系数。推导见[unified-structure.md](unified-structure.md)，所有基函数的已求值 Hlog 见[结构基函数](structural-basis-functions.md)。z>0；formal: not_run。

定义 ρ=τ^{[2],π}_{m,n;q}=[ε^q]Q^+，χ=η^{[2],π}_{m,n;q}=ln z·[ε^q]Q^-。括号内的所有函数具有相同 z。对于η，q只计母核的ε阶，外乘ln z使权重再加1。

全体 ε⁻⁶…ε⁻² 系数为零；4e-04 的 ε⁻¹ 也为零。下式余项均为 O(ε)，而非完整 ε 函数。

## 4e-04

$$I=- \frac{\tau^{[2],\mathrm{id}}_{0,0;0}}{8}+O(\varepsilon).
$$

## 6e-10 = 6e-23

$$I=- \frac{\eta^{[2],\mathrm{id}}_{0,0;1}}{16} + \frac{\eta^{[2],\mathrm{id}}_{0,1;1}}{16} - \frac{\tau^{[2],\mathrm{id}}_{0,0;2}}{8} + \frac{\tau^{[2],\mathrm{id}}_{0,1;2}}{8} + \frac{- \frac{\tau^{[2],\mathrm{id}}_{0,0;1}}{8} + \frac{\tau^{[2],\mathrm{id}}_{0,1;1}}{8}}{\varepsilon}+O(\varepsilon).
$$

## 6e-11 = 6e-25

$$I=- \frac{\eta^{[2],\mathrm{id}}_{0,1;1}}{32} + \frac{\eta^{[2],\mathrm{id}}_{1,0;1}}{32} - \frac{\tau^{[2],\mathrm{id}}_{0,1;2}}{16} + \frac{\tau^{[2],\mathrm{id}}_{1,0;2}}{16} + \frac{- \frac{\tau^{[2],\mathrm{id}}_{0,1;1}}{16} + \frac{\tau^{[2],\mathrm{id}}_{1,0;1}}{16}}{\varepsilon}+O(\varepsilon).
$$

## 6e-12 = 6e-24

$$I=- \frac{\eta^{[2],\mathrm{swap}}_{0,0;1}}{32} + \frac{\eta^{[2],\mathrm{swap}}_{0,1;1}}{32} - \frac{\tau^{[2],\mathrm{swap}}_{0,0;2}}{16} + \frac{\tau^{[2],\mathrm{swap}}_{0,1;2}}{16} + \frac{- \frac{\tau^{[2],\mathrm{swap}}_{0,0;1}}{16} + \frac{\tau^{[2],\mathrm{swap}}_{0,1;1}}{16}}{\varepsilon}+O(\varepsilon).
$$

## Σ 的位置

ΣA/B/C 就是上述三个三圈式的有限项。旧名字可以继续用于读表，但不作为新的基础积分定义。新系数均从同一个Q族产生，可以按m,n,q推广；低阶、基底覆盖和独立生成元数应分开判断。
