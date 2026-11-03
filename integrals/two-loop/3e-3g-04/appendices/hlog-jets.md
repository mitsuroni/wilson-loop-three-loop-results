# Complete evaluated HLOG jets used at positive epsilon orders

Every function in this appendix is on real $z>0$. Let $G_{(a_1,\ldots,a_n)}=G(a_1,\ldots,a_n;z)$ with the first-letter kernel $dt/(t-a_1)$, and $G_{\varnothing}=1$. Every row below is an explicit coefficient in the finite formula
\[
X(z)=\sum_{\mathbf w\in\mathcal W_X}c_{X,\mathbf w}(z)G_{\mathbf w}(z).
\]
The table lists the entire finite set $\mathcal W_X$, including the empty word when its period coefficient is nonzero; all omitted words have zero coefficient. This is a complete analytic expression, not a link to a final result or a definition by an unknown integral. The coefficients are exact rational functions of $z$ and zeta periods. The native original outputs are the six files `D2.txt`, `D3.txt`, `H1.txt`, `H2.txt`, `O1.txt`, `O2.txt` in the immutable [positive-epsilon run](../../../../workflow/checks/tgv04-high-order/runs/positive-epsilon/). These tables transcribe those outputs into first-letter HLOG coordinates.

The kernels, domains, actual variable orders and necessary source orders are specified in the main text. The adopted computational evaluation lemma is: integrating those six displayed kernels over those domains in the divergence-checked regularization, with real positive $z$, gives exactly the finite coordinate tables below. The historic HyperInt session, endpoint regularization, fibration and coordinate checks are reused; the current Python entry adopts the already assembled normalized I coefficients and does not independently recompute this lemma or live-convolve the subintegral jets.

## $D_2$

| Word $\mathbf w$ | Coefficient $c_{D2,\mathbf w}$ |
|---|---|
| $\varnothing$ | $\frac{2 \zeta_{2} \zeta_{3}}{z + 1} + \frac{2 \zeta_{5}}{z + 1}$ |
| $(0)$ | $\frac{7 \zeta_{2}^{2}}{2 \left(z + 1\right)}$ |
| $(0,-1,0)$ | $\frac{4 \zeta_{2}}{z + 1}$ |
| $(-1,0,0,0,0)$ | $\frac{10}{z + 1}$ |
| $(-1,-1,0,0,0)$ | $- \frac{4}{z + 1}$ |
| $(0,-1,0,0,0)$ | $\frac{4}{z + 1}$ |
| $(0,0)$ | $\frac{6 \zeta_{3}}{z + 1}$ |
| $(-1,0)$ | $- \frac{4 \zeta_{3}}{z + 1}$ |
| $(0,0,0)$ | $- \frac{6 \zeta_{2}}{z + 1}$ |
| $(0,0,0,0,0)$ | $- \frac{17}{z + 1}$ |
| $(-1)$ | $- \frac{17 \zeta_{2}^{2}}{5 \left(z + 1\right)}$ |
| $(-1,-1,0)$ | $- \frac{4 \zeta_{2}}{z + 1}$ |
| $(-1,0,0)$ | $\frac{6 \zeta_{2}}{z + 1}$ |

## $D_3$

| Word $\mathbf w$ | Coefficient $c_{D3,\mathbf w}$ |
|---|---|
| $\varnothing$ | $\frac{92 \zeta_{2}^{3}}{35 z + 35} - \frac{2 \zeta_{3}^{2}}{z + 1}$ |
| $(-1,-1,0,0,0,0)$ | $- \frac{20}{z + 1}$ |
| $(-1,0,0,0,0,0)$ | $\frac{34}{z + 1}$ |
| $(0,-1,0,0,0,0)$ | $\frac{20}{z + 1}$ |
| $(-1,-1,-1,0,0,0)$ | $\frac{8}{z + 1}$ |
| $(0,-1,-1,0,0,0)$ | $- \frac{8}{z + 1}$ |
| $(-1,0,-1,0,0,0)$ | $- \frac{8}{z + 1}$ |
| $(0,0,-1,0,0,0)$ | $\frac{8}{z + 1}$ |
| $(0,0,0,0,0,0)$ | $- \frac{49}{z + 1}$ |
| $(-1,0,-1,0)$ | $- \frac{8 \zeta_{2}}{z + 1}$ |
| $(0,-1,-1,0)$ | $- \frac{8 \zeta_{2}}{z + 1}$ |
| $(-1)$ | $- \frac{4 \left(\zeta_{2} \zeta_{3} + \zeta_{5}\right)}{z + 1}$ |
| $(0,-1,0,0)$ | $\frac{12 \zeta_{2}}{z + 1}$ |
| $(-1,-1)$ | $\frac{34 \zeta_{2}^{2}}{5 \left(z + 1\right)}$ |
| $(0,0,-1,0)$ | $\frac{8 \zeta_{2}}{z + 1}$ |
| $(0,0,0)$ | $\frac{16 \zeta_{3}}{z + 1}$ |
| $(-1,0)$ | $- \frac{7 \zeta_{2}^{2}}{z + 1}$ |
| $(0,0,0,0)$ | $- \frac{10 \zeta_{2}}{z + 1}$ |
| $(-1,-1,0)$ | $\frac{8 \zeta_{3}}{z + 1}$ |
| $(-1,0,0)$ | $- \frac{12 \zeta_{3}}{z + 1}$ |
| $(-1,-1,0,0)$ | $- \frac{12 \zeta_{2}}{z + 1}$ |
| $(0,-1)$ | $- \frac{34 \zeta_{2}^{2}}{5 \left(z + 1\right)}$ |
| $(0,-1,0)$ | $- \frac{8 \zeta_{3}}{z + 1}$ |
| $(0,0)$ | $\frac{71 \zeta_{2}^{2}}{10 \left(z + 1\right)}$ |
| $(0)$ | $\frac{2 \left(\zeta_{2} \zeta_{3} + 5 \zeta_{5}\right)}{z + 1}$ |
| $(-1,0,0,0)$ | $\frac{12 \zeta_{2}}{z + 1}$ |
| $(-1,-1,-1,0)$ | $\frac{8 \zeta_{2}}{z + 1}$ |

## $H_1$

| Word $\mathbf w$ | Coefficient $c_{H1,\mathbf w}$ |
|---|---|
| $(-1,0,-1,0,0)$ | $\frac{2}{z}$ |
| $(-1,0,-1)$ | $\frac{6 \zeta_{2}}{z}$ |
| $(0,-1)$ | $- \frac{4 \zeta_{3}}{z}$ |
| $(-1,-1,0,0,0)$ | $\frac{6}{z}$ |
| $(-1)$ | $\frac{41 \zeta_{2}^{2}}{5 z}$ |
| $(-1,-1,0)$ | $\frac{4 \zeta_{2}}{z}$ |
| $(-1,-1,-1,0,0)$ | $\frac{2}{z}$ |
| $(0,-1,-1)$ | $- \frac{12 \zeta_{2}}{z}$ |
| $(-1,-1)$ | $\frac{2 \zeta_{3}}{z}$ |
| $(-1,-1,-1)$ | $\frac{6 \zeta_{2}}{z}$ |
| $(0,-1,0)$ | $\frac{4 \zeta_{2}}{z}$ |
| $(-1,0,0)$ | $- \frac{2 \zeta_{2}}{z}$ |
| $(0,-1,-1,0,0)$ | $- \frac{4}{z}$ |

## $H_2$

| Word $\mathbf w$ | Coefficient $c_{H2,\mathbf w}$ |
|---|---|
| $(-1,0,0,0)$ | $- \frac{2 \zeta_{2}}{z}$ |
| $(-1,-1,-1)$ | $\frac{2 \zeta_{3}}{z}$ |
| $(0,-1)$ | $- \frac{82 \zeta_{2}^{2}}{5 z}$ |
| $(-1,0,-1,-1)$ | $- \frac{18 \zeta_{2}}{z}$ |
| $(-1,-1)$ | $\frac{56 \zeta_{2}^{2}}{5 z}$ |
| $(-1)$ | $- \frac{2 \left(5 \zeta_{2} \zeta_{3} - \zeta_{5}\right)}{z}$ |
| $(-1,-1,-1,0)$ | $\frac{4 \zeta_{2}}{z}$ |
| $(-1,-1,0,-1)$ | $\frac{6 \zeta_{2}}{z}$ |
| $(-1,-1,-1,-1)$ | $\frac{6 \zeta_{2}}{z}$ |
| $(-1,0,-1,0)$ | $\frac{12 \zeta_{2}}{z}$ |
| $(0,0,-1,0)$ | $- \frac{8 \zeta_{2}}{z}$ |
| $(-1,-1,0,0)$ | $\frac{6 \zeta_{2}}{z}$ |
| $(0,0,-1)$ | $\frac{8 \zeta_{3}}{z}$ |
| $(0,-1,0,-1)$ | $- \frac{12 \zeta_{2}}{z}$ |
| $(-1,0,-1)$ | $- \frac{6 \zeta_{3}}{z}$ |
| $(-1,0,0,-1)$ | $\frac{6 \zeta_{2}}{z}$ |
| $(0,-1,-1,-1)$ | $- \frac{12 \zeta_{2}}{z}$ |
| $(0,-1,0,0)$ | $\frac{4 \zeta_{2}}{z}$ |
| $(-1,0)$ | $\frac{33 \zeta_{2}^{2}}{5 z}$ |
| $(0,-1,-1,0)$ | $- \frac{8 \zeta_{2}}{z}$ |
| $(0,-1,-1)$ | $- \frac{4 \zeta_{3}}{z}$ |
| $(0,0,-1,-1)$ | $\frac{24 \zeta_{2}}{z}$ |
| $(0,0,-1,-1,0,0)$ | $\frac{8}{z}$ |
| $(0,-1,-1,-1,0,0)$ | $- \frac{4}{z}$ |
| $(-1,0,0,-1,0,0)$ | $\frac{2}{z}$ |
| $(-1,0,-1,0,0,0)$ | $\frac{6}{z}$ |
| $(-1,-1,0,-1,0,0)$ | $\frac{2}{z}$ |
| $(-1,-1,0,0,0,0)$ | $\frac{14}{z}$ |
| $(0,-1,-1,0,0,0)$ | $- \frac{12}{z}$ |
| $(0,-1,0,-1,0,0)$ | $- \frac{4}{z}$ |
| $(-1,-1,-1,0,0,0)$ | $\frac{6}{z}$ |
| $(-1,-1,-1,-1,0,0)$ | $\frac{2}{z}$ |
| $(-1,0,-1,-1,0,0)$ | $- \frac{6}{z}$ |

## $O_1$

| Word $\mathbf w$ | Coefficient $c_{O1,\mathbf w}$ |
|---|---|
| $\varnothing$ | $\frac{44 \zeta_{2}^{2}}{5} - \frac{25 \zeta_{2} \zeta_{3}}{2} - 2 \zeta_{3} - \frac{13 \zeta_{5}}{2}$ |
| $(-1,-1,0,0,0)$ | $10$ |
| $(-1,-1,0,0)$ | $6$ |
| $(-1,-1,-1,0,0)$ | $-4$ |
| $(0,0,0,0,0)$ | $7$ |
| $(0,-1,0,0)$ | $-4$ |
| $(0,-1,0,0,0)$ | $-9$ |
| $(-1,0,0,0,0)$ | $- \frac{17}{2}$ |
| $(-1,0,-1,0,0)$ | $5$ |
| $(0,0,0,0)$ | $2$ |
| $(-1,0,0,0)$ | $-5$ |
| $(0,0,-1,0,0)$ | $-6$ |
| $(0,-1,-1,0,0)$ | $5$ |
| $(0,0,-1)$ | $- 18 \zeta_{2}$ |
| $(0,0)$ | $6 \left(\zeta_{2} - \zeta_{3}\right)$ |
| $(0,-1,-1)$ | $15 \zeta_{2}$ |
| $(0,0,0)$ | $11 \zeta_{2}$ |
| $(-1,-1)$ | $2 \left(9 \zeta_{2} - 2 \zeta_{3}\right)$ |
| $(0,-1,0)$ | $- 14 \zeta_{2}$ |
| $(-1,-1,-1)$ | $- 12 \zeta_{2}$ |
| $(-1,0,0)$ | $- 2 \left(5 \zeta_{2} + 1\right)$ |
| $(-1,-1,0)$ | $14 \zeta_{2}$ |
| $(-1,0,-1)$ | $15 \zeta_{2}$ |
| $(-1,0)$ | $- 11 \zeta_{2} + 4 \zeta_{3}$ |
| $(0,-1)$ | $- 12 \zeta_{2} + 5 \zeta_{3}$ |
| $(0)$ | $\frac{34 \zeta_{2}^{2} + 10 \zeta_{2} - 15 \zeta_{3}}{5}$ |
| $(-1)$ | $\frac{7 \zeta_{2}^{2} - 120 \zeta_{2} + 120 \zeta_{3}}{20}$ |

## $O_2$

| Word $\mathbf w$ | Coefficient $c_{O2,\mathbf w}$ |
|---|---|
| $\varnothing$ | $- \frac{314 \zeta_{2}^{3}}{105} + \frac{88 \zeta_{2}^{2}}{5} - 25 \zeta_{2} \zeta_{3} + \frac{7 \zeta_{3}^{2}}{2} - 4 \zeta_{3} - 13 \zeta_{5}$ |
| $(-1,-1,0,-1,0,0)$ | $-9$ |
| $(0,0,-1,-1,0,0)$ | $16$ |
| $(-1,-1,-1,0,0)$ | $-8$ |
| $(-1,-1,-1,0,0,0)$ | $-15$ |
| $(0,0,0,0,0)$ | $14$ |
| $(-1,0,0,-1)$ | $27 \zeta_{2}$ |
| $(0,0,0)$ | $2 \left(11 \zeta_{2} - 10 \zeta_{3}\right)$ |
| $(0,-1,-1,-1)$ | $- 27 \zeta_{2}$ |
| $(-1,-1,0)$ | $7 \left(4 \zeta_{2} - \zeta_{3}\right)$ |
| $(-1,-1,-1)$ | $- 2 \left(12 \zeta_{2} - 5 \zeta_{3}\right)$ |
| $(0,0,-1,-1)$ | $48 \zeta_{2}$ |
| $(-1,-1)$ | $\frac{61 \zeta_{2}^{2} + 180 \zeta_{2} - 40 \zeta_{3}}{5}$ |
| $(0,-1,0,-1,0,0)$ | $8$ |
| $(0,0,-1,0,0,0)$ | $-26$ |
| $(-1,0,0,-1,0,0)$ | $9$ |
| $(0,-1,0,0,0,0)$ | $-27$ |
| $(-1,0,-1,-1,0,0)$ | $-9$ |
| $(-1,-1,0,0,0)$ | $20$ |
| $(-1,0,-1,0,0)$ | $10$ |
| $(0,0,0,-1,0,0)$ | $-16$ |
| $(0,0,-1,0,0)$ | $-12$ |
| $(0,-1,0,0,0)$ | $-18$ |
| $(0,0,0,0,0,0)$ | $31$ |
| $(-1,-1,-1,-1,0,0)$ | $10$ |
| $(-1,0,0,0,0)$ | $-17$ |
| $(-1,-1,0,0,0,0)$ | $28$ |
| $(0,-1,-1,0,0)$ | $10$ |
| $(-1,0,-1,0,0,0)$ | $17$ |
| $(-1,0,0,0,0,0)$ | $- \frac{49}{2}$ |
| $(0,-1,-1,-1,0,0)$ | $-9$ |
| $(0,-1,-1,0,0,0)$ | $16$ |
| $(0,-1,0,-1)$ | $24 \zeta_{2}$ |
| $(0)$ | $\frac{136 \zeta_{2}^{2} - 275 \zeta_{2} \zeta_{3} + 40 \zeta_{2} - 60 \zeta_{3} - 155 \zeta_{5}}{10}$ |
| $(0,-1,0)$ | $- 28 \zeta_{2} + 9 \zeta_{3}$ |
| $(-1,0,-1)$ | $3 \left(10 \zeta_{2} - 3 \zeta_{3}\right)$ |
| $(-1)$ | $\frac{7 \zeta_{2}^{2} + 20 \zeta_{2} \zeta_{3} - 120 \zeta_{2} + 120 \zeta_{3} + 100 \zeta_{5}}{10}$ |
| $(-1,-1,-1,0)$ | $- 25 \zeta_{2}$ |
| $(0,-1,-1,0)$ | $25 \zeta_{2}$ |
| $(-1,0,-1,0)$ | $26 \zeta_{2}$ |
| $(-1,-1,-1,-1)$ | $30 \zeta_{2}$ |
| $(0,0,-1,0)$ | $- 42 \zeta_{2}$ |
| $(0,-1,-1)$ | $3 \left(10 \zeta_{2} - 3 \zeta_{3}\right)$ |
| $(0,0,-1)$ | $- 4 \left(9 \zeta_{2} - 4 \zeta_{3}\right)$ |
| $(-1,0)$ | $- \frac{99 \zeta_{2}^{2} + 440 \zeta_{2} - 160 \zeta_{3}}{20}$ |
| $(0,0,0,-1)$ | $- 48 \zeta_{2}$ |
| $(-1,-1,0,-1)$ | $- 27 \zeta_{2}$ |
| $(-1,0,0)$ | $- 2 \left(10 \zeta_{2} - 5 \zeta_{3} + 2\right)$ |
| $(-1,-1,0,0)$ | $2 \left(13 \zeta_{2} + 6\right)$ |
| $(0,0,0,0)$ | $4 \left(8 \zeta_{2} + 1\right)$ |
| $(0,-1)$ | $- \frac{2 \left(8 \zeta_{2}^{2} + 60 \zeta_{2} - 25 \zeta_{3}\right)}{5}$ |
| $(-1,0,-1,-1)$ | $- 27 \zeta_{2}$ |
| $(0,-1,0,0)$ | $- 23 \zeta_{2} - 8$ |
| $(-1,0,0,0)$ | $- 2 \left(9 \zeta_{2} + 5\right)$ |
| $(0,0)$ | $\frac{19 \zeta_{2}^{2} + 24 \zeta_{2} - 24 \zeta_{3}}{2}$ |

## The finite recoil jet

The remaining lower jet, required in the same assembly, is
\[
\begin{aligned}
O_0={}&-G_{-1,0,0}+(-3\zeta_2+3\zeta_3)G_{-1}+22\zeta_2^2/5-\zeta_3
+(\zeta_2-3\zeta_3/2)G_0-2G_{0,-1,0,0}+3\zeta_2G_{0,0}\\
&-6\zeta_2G_{0,-1}-11\zeta_2G_{-1,0}/2+G_{0,0,0,0}+3G_{-1,-1,0,0}
+9\zeta_2G_{-1,-1}-5G_{-1,0,0,0}/2.
\end{aligned}
\]
Together with $D_0,D_1,H_0,O_{-1}$ in the main text, the six coordinate tables supply every jet needed through $[e^2]I$. The explicit period-only rows fix endpoint constants; derivative or symbol comparisons alone would not determine them. Equations (8)-(9) in the main text now assemble a fully specified HLOG expression before the optional $M,N,\tau$ coordinate compression.
