# Complete physical epsilon-four coefficients

The [explicit coefficient file](coefficients-eps4-gpl.json) contains the full weight-ten coefficients of `6e-20` and `6e-21=6e-22`. It uses exact rational coefficients and finite products of explicitly listed GPL words. All Gamma/Euler, rho and cosh normalization factors are included. This representation has not been reduced to a weight-ten Hlog/period basis.

For `z>0`, every GPL has argument 1 and letters `0,1,Z=-z,IZ=-1/z`, with kernel `dt/(t-a)` on the real branch. The JSON `GPL_words` array is an explicit word table. A row

```text
[c,[i1,...,ik],p,q,r,s]
```

means

\[
c\,\prod_{j=1}^kG(\mathrm{GPL\_words}[i_j];1)
\log^p z\,\log^q\!\frac{z}{(1+z)^2}\,\zeta_2^r\zeta_3^s.
\]

Indices are zero based. Summing the 38,894 rows for I20 or 39,080 rows for I21 gives its complete epsilon-four coefficient. No row denotes an unevaluated source integral or a fitted remainder. This is a computational representation, not a minimal function basis.

To display the normalization compactly, let `r_j,q` be the degree `q+3` coefficient of the physical source numerator after the factors `epsilon^-3`, `P` and `rho^(3 epsilon/2)` are removed. The exact source definitions and signs are in the [continuous derivation](derivation.md). Write

\[
a=\frac32\log\rho,\qquad
p_2=\frac32\zeta_2-\frac98\log^2z,\qquad
p_4=\frac{57}{40}\zeta_2^2-\frac{27}{16}\zeta_2\log^2z+\frac{135}{128}\log^4z.
\]

The fixed convolution used to generate the file is

\[
I_j^{(4)}=r_{j,4}+a r_{j,3}
+\left(\frac{a^2}{2}+p_2\right)r_{j,2}
+\left(\frac{a^3}{6}+ap_2+\zeta_3\right)r_{j,1}
+\left(\frac{a^4}{24}+\frac{a^2p_2}{2}+a\zeta_3+p_4\right)r_{j,0}.
\]

The source rule and this normalization were frozen before generating the new source degree seven. The output exactly preserves all previously exported source coefficients. At `z=0.4,0.7`, the new predictions agree with independent original ordered-kernel quadratures with maximum observed discrepancy `5.79e-11`; reciprocal checks also pass. These are numerical tests with empirical refinement, not certified bounds or an exact all-z identity check. The [actual prediction audit](../../../../../workflow/checks/pure-three-higher-laurent-20261016/prediction-order4/prediction-verification-r003.json) records the method, precision controls and limitations. Formal proof status is `not_run`.
