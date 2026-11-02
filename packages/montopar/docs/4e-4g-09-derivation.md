# 4e-4g-09：修正后的几何与C/J

2026-09-18。取 $(a,b,c,d)=(1,2,3,4)$，$s_{11}=x_1,s_{21}=x_2,s_{31}=x_3,s_{41}=x_4$，Jacobian为1。通用推导见[4GV修订稿](4GV-general-derivation.md)。

首尾边属于 $4\to1$ 相邻：$V_4=-v_4$，因此

$$X_1-X_4=x_1v_1+(1-x_4)v_4,\qquad d_{14}=x_1(1-x_4).$$

旧文档和旧三份导出中的 $(1-x_1)x_4$ 错误，不能由“无向距离”解释，因为两式不是整体点差反号。

$$\begin{aligned}
S={}&\alpha_1\alpha_2z(1-x_1)x_2
+\alpha_1\alpha_3[z(1-x_1)(1-x_3)+x_1x_3]
+\alpha_1\alpha_4x_1(1-x_4)\\
&+\alpha_2\alpha_3(1-x_2)x_3
+\alpha_2\alpha_4[(1-x_2)(1-x_4)+zx_2x_4]
+\alpha_3\alpha_4z(1-x_3)x_4.
\end{aligned}$$

颜色通道为 $C_A^2C_F(-1/4,0,1/4)$，合并的无量纲分子为

$$\ell_{09}(z)=-(1+4z+z^2)/4.$$

$$C_{\rm raw}=\frac{C_A^2C_F}{16}(\alpha_s/\pi)^3\Gamma(2-3\epsilon)
\pi^{3\epsilon}(-2v_{14}+i0)^{3\epsilon},$$

$$J_{09}=\int_{[0,1]^4}d^4x\int_{\alpha_i\ge0}d^4\alpha\,
\delta(1-\sum\alpha)\ell_{09}(z)(\prod\alpha_i)^{-\epsilon}S^{-2+3\epsilon}.$$

原s表示与cube表示均保存单纯形。物理分支按通用推导从 $z>0,v_{14}<0$ 区域延拓，不能单独给S任意选择分支。raw未补Euler，也未乘orbit重数；四点星图的C4几何重数为1，仍不包含在这里的C内。

新 `final.json`、`state.json`、`export.dimreduce.json` 必须从本轮主输出重新生成。逐图几何、归一化和数据契约检查不能由本图单独代替。
