# Validation of the two-loop result

The calculation connects twelve source diagrams to ordered parameter integrals, evaluated Laurent series, and a renormalized Wilson-loop exponent. The checks of these objects have different scopes.

| Check | Conclusion |
| --- | --- |
| Source and parameterization | The diagram labels, ordered domains, normalization factors and color/orbit assignments connect the source integrals to the normalized I values |
| Integral expansions | The catalogue contains all twelve scalar series through epsilon squared, together with the two self-energy color channels |
| Bare assembly | The fourth-order pole cancels in the exponent sum; the finite M_(2,1) and M_(0,1) squared structures cancel after the stated identities |
| Reference counterterm | The two color channels and all three surviving pole orders agree with the prediction from the known beta function and cusp/non-cusp anomalous dimensions |
| Finite RG equation | The scale dependence of F_2 agrees with the reference anomalous dimensions when the running of the coupling is included |
| Complete finite comparison | The full finite function, including its scale-independent shape dependence and constants, agrees with the published expression after matching the observable, coupling, logarithms and branches |

## Pole cancellation and its independent test

The bare sums S_A and S_f contain poles. The combination that becomes finite is w_2-hat + delta w_2 - Z_(4,2), where delta w_2 = -beta_0 w_1-hat/(4 epsilon). Extracting Z_(4,2) from those same poles and subtracting it defines the renormalized result. The independent check uses a counterterm predicted from reference anomalous dimensions:

\[
Z_{4,2}^{\rm ref}
=\frac{3C_F\beta_0}{8\varepsilon^3}
+\frac{C_F\beta_0L/4-K_2^{\rm ref}/2}{\varepsilon^2}
-\frac{K_2^{\rm ref}L+\gamma_2^{\rm ref}/4}{\varepsilon}.
\]

Here K_F = Gamma_cusp^F = C_F lambda + K_2 lambda squared + ... is the lightlike fundamental cusp anomalous dimension, and gamma_4 = gamma_2 lambda squared + ... is the non-cusp term in the report's convention. Every predicted pole matches the calculated sum.

## What the finite RG equation fixes

At fixed geometry, dL/d log(mu squared) = 1 and z is constant. Coupling running gives

\[
\partial_L F_2-\frac{\beta_0}{4}F_1=-2K_2L-\frac12\gamma_2.
\]

This relation tests the finite scale logarithms. Adding a function h(z) independent of L leaves it unchanged at this order. The complete finite comparison therefore provides additional information: it tests the M_(1,1) contribution and the finite constants as well as the logarithms.

The reference is the renormalized QCD Wilson loop supplied in the authors' `WRen.m` ancillary file for Chen et al., archived at [Zenodo 17200473](https://doi.org/10.5281/zenodo.17200473). Write the file as W_R = 1 + a_anc W_1 + a_anc^2 W_2 + ..., where a_anc = alpha_s/(4 pi) = lambda/4. Then F_2^ref = (W_2 - W_1^2/2)/16 after the variable and branch conversion. The factor 1/16 comes from the coupling conversion; the quadratic subtraction comes from taking log W_R. The one-loop term fixes log(s_a^2 mu_a^2) = L + ell/2, log(t_a^2 mu_a^2) = L - ell/2 and s_a^2/t_a^2 = z. The analytic difference from the calculated F_2 is zero on the matched positive-real-z branch.

## Individual integrals and higher epsilon orders

Agreement of the physical sum constrains combinations of diagrams. Compensating errors in separate diagrams could leave that sum unchanged. The individual integral source and expansion checks therefore remain necessary. The finite physical comparison also does not certify every diagram's epsilon and epsilon-squared coefficients. Those coefficients are supported by the exact formulas, saved analytic results and high-order supplements used to construct the catalogue.

The archived report verification contains 427 catalogue/assembly checks, 66 momentum-to-parameter checks, and 52 checks of the renormalization chain with five negative controls. The external finite comparison and numerical controls are retained from the original physical analysis, with 36 analytic checks and 50 numerical controls. These are separate collections with overlapping purposes, rather than a count of independent proofs. Their source bindings are recorded in `provenance.json`.

The supplied `check-results.wl` checks the exported data and the physical assembly in a fresh evaluation. The source integrations and full literature comparison are retained evidence. The scope is fixed z > 0 with the stated scale prescription; general complex continuation and a complete machine-checked formal proof are outside these results.

## References

The cusp coefficient is compared in the lambda = alpha_s/pi convention with A. Grozin, J. M. Henn, G. P. Korchemsky and P. Marquard, *The three-loop cusp anomalous dimension in QCD and its supersymmetric extensions*, JHEP 01 (2016) 140, [arXiv:1510.07803](https://arxiv.org/abs/1510.07803), Eq. (6.12).

The non-cusp and complete finite comparisons use H. Chen, P. F. Monni, Z. Pang, G. Vita and H. X. Zhu, *Correlation Function/Wilson Loop Duality in Gauge Theory from EFT*, [arXiv:2510.07377](https://arxiv.org/abs/2510.07377), and the authors' `ADs.m` and `WRen.m` [ancillary files](https://doi.org/10.5281/zenodo.17200473).


The diagram organization and the separate supersymmetric comparison also use
J. M. Drummond, J. Henn, G. P. Korchemsky and E. Sokatchev,
*On planar gluon amplitudes/Wilson loops duality*, Nucl. Phys. B 795 (2008)
52-68, [arXiv:0709.2368](https://arxiv.org/abs/0709.2368), Sections 2.2-2.3.
That comparison requires the N=4 matter insertion and dimensional-reduction
conventions. The four-point finite constant is taken with the subsequent
correction in C. Anastasiou et al., *Two-Loop Polygon Wilson Loops in N=4 SYM*,
JHEP 05 (2009) 115, [arXiv:0902.2245](https://arxiv.org/abs/0902.2245),
Eq. (3.15) and the text after Eq. (3.16). This supersymmetric comparison and
the QCD F_2^ref from Chen et al. have different theory and normalization inputs.
