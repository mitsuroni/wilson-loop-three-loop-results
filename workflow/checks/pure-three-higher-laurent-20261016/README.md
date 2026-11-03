# Higher Laurent coefficients as constraints on the complete generator

2026-10-16. The user explicitly established this exploration goal. Its deliverables and requirement-by-requirement acceptance are complete. [Scope and milestones](../../pure-three-higher-laurent-plan.md), [continuous mathematical record](../../../benchmarks/2ne/theory/generation/pure-three-higher-laurent/derivation.md), [acceptance audit](goal-acceptance-audit.md).

Complete new physical orders three and four of 6e-20 and 6e-21=6e-22 have been generated. The frozen generator predicts reserved order four. Original-integral identities reduce both physical coupled parts to the existing ten quadratic moments, with fixed one-dimensional J insertions; their physical-z transport is now explicit. The natural weighted two-J family reconstructs all physical/bare orders zero through three without using any old three-kernel coefficient as a generator. BJ-only all-epsilon formulas and the full ordered-Q endpoint classification remain stronger open questions. [Final audit](verification-r004.json), [natural-source stage](verification-r003.json), [reserved-order stage](verification-r002.json), [initial complete-order-three stage](verification.json).

## Initial actual run

[start_exploration.py](start_exploration.py) ran successfully and created a fresh [r001](runs/r001/verification.json) directory in approximately 10.2 seconds. Its source constructors were imported from the previous research directory without running their report-writing main functions. The old inputs, scripts and reports retained their hashes.

- [baseline.json](runs/r001/baseline.json) checks 1456 stored physical/bare Hlog-period coordinates with exact rational coefficients, correct weights, no duplicate coordinates, and I21=I22. It binds the existing epsilon-zero/one/two JSON by SHA256. This is an audit of stored results, not a source-to-physical regression or a new proof that negative Laurent powers vanish.
- [W-parent-words.json](runs/r001/W-parent-words.json) explicitly stores all 56 generic parent-source coefficients through total source order six, including auxiliary low-order terms that cancel in the physical projection. Every rational coefficient and GPL word can be read directly from the file. This is not an assembled W or physical I table.
- [Lambda-sector-words.json](runs/r001/Lambda-sector-words.json) explicitly stores complete generic sector coefficients at source orders three through six. The term counts are 86,342,1136,3405, with weights six through nine. Exact serialization roundtrips, weights and the A=1 endpoint-log cancellation checks passed at each order. The two A=1 sectors share the same coalescence pattern; this is not two independent numerical integrations.
- [new-source-numeric.json](runs/r001/new-source-numeric.json) checks the new sixth source coefficient at (A,B,C,d)=(-1/3,-1/4,-1/5,-1/6). Its 3405-term GPL expression agrees with direct quadrature of the separately expanded original complete kernel, with absolute difference 7.05e-56 at 60 working digits. GPL sums use 140 terms; the independent quadrature uses 30 nodes per dimension and a 100-term evaluation of J's coefficient. This is an observed discrepancy, not an interval or tail-error certificate.

Both source files specify their letter maps and the exact product convention. The Lambda boundary-J degree is an additional explicitly evaluated GPL factor, not an undefined period. The export does not substitute the three physical sectors, add the mirror, or convolve the rho/Gamma/cosh factors. Therefore no new complete physical epsilon-cubed result or new function-space rank is claimed by this run. `formal=not_run`.

## Reproduction

```powershell
D:/py/python.exe -B -X utf8 workflow/checks/pure-three-higher-laurent-20261016/start_exploration.py --run r003 --max-source-order 7
```

The default maximum source order is six. Each run uses a new directory and refuses to overwrite an old run. The actual [r002](runs/r002/verification.json) extends to source degree seven after the prediction rule was frozen; the command above requests a fresh reproduction directory. Inputs and outputs are hash-bound in each run.

## Complete physical order three

The fresh [physical run](physical-runs/20261016T021005211170Z/status.json) reads the explicit source exports and uses [physical_sources.py](physical_sources.py) to assemble every oriented sector and corner term. Endpoint logarithms are retained until their exact coefficients cancel. W uses the coupled parent source and all three corner terms; Lambda uses its complete sector. No old target coefficient enters the calculation.

[run_physical.py](run_physical.py) starts a fresh Maple/HyperInt session with the period table explicitly loaded, divergence checks enabled and no blockwise regHlog. It integrates all physical numerators through source degree six, adds the mirror and rho factor, constructs H independently from its J primitive, and restores the full Gamma/Euler/cosh prefactor. The [actual log](physical-runs/20261016T021005211170Z/maple.log) reports Maple CPU time about 498 seconds; this is not a cold-start wall-time benchmark.

[assemble_data.py](assemble_data.py) and the [exact assembly report](physical-runs/20261016T021005211170Z/assembly-verification.json) verify 24 exact regressions: six independent old physical expressions, twelve old bare/source expressions, and six source numerators. Four missing-normalization negative controls are rejected. The complete new weight-nine physical results contain 356 and 357 Hlog/period coordinates. W3 and Lambda3 contain 357 and 230. These are coordinate counts, not master-integral counts.

Published artifacts are [the full formulas](../../../benchmarks/2ne/theory/generation/pure-three-higher-laurent/coefficients-eps3-hlog.md) and [machine-readable coefficients](../../../benchmarks/2ne/theory/generation/pure-three-higher-laurent/coefficients-through-eps3.json). The latter is a byte copy of the newly calculated [coefficient file](physical-runs/20261016T021005211170Z/coefficients.json), not a copy of the old baseline.

[numeric_physical.py](numeric_physical.py) compares the new order-three coefficients with the original ordered three/four-dimensional endpoint integrals at z=0.4,0.7. All eight comparisons and four missing-interval-variance negative controls pass. The [actual report](physical-runs/20261016T021005211170Z/numeric-verification.json) records maximum observed discrepancy 2.19e-13. The original-coordinate calculation and 55-digit Hlog ODE series are independent of the new word-endpoint integration. Quadrature errors are empirical.

## Declared weight-nine spaces

[run_m9.py](run_m9.py) computes all eight weight-nine M moments from their defining natural source in a fresh CAS session. [The run](candidate-runs/20261016T021400575141Z/status.json) and [independent defining-integral checks](candidate-runs/20261016T021400575141Z/numeric-verification.json) pass; the maximum observed discrepancy is 3.50e-42.

[compare_weight9.py](compare_weight9.py) declares the exact finite spaces and saves all columns before elimination. The accepted [r002 comparison](physical-runs/20261016T021005211170Z/r002-comparison-weight9.json) includes an explicit normalization-equivalence gate and exact dual witnesses. Earlier exploratory comparison data and its executed source snapshot are retained separately.

| Space | Rank | Rank after both physical targets |
|---|---:|---:|
| M/N, registered periods and products | 59 | 61 |
| Old weight-six extras and products added | 63 | 65 |
| Known normalization logarithms added | 156 | 158 |
| Old W/Lambda orders zero through two added | 161 | 163 |

## Existing Gauss/Mellin natural sources

[gauss_family.py](gauss_family.py) expands the existing P source through the exact Gauss ODE and Gamma ratio, then generates Mellin moments by prepending zero letters. The family is the one defined in the single-opposite source note: its Gauss function has equal upper parameters a,a. It is not the fully independent a,b family.

[The fresh family run](gauss-runs/20261016T023349031430Z/status.json) exports 36 P and 84 Mellin-P moments through weight nine and calculates both crossing components. [Independent numeric checks](gauss-runs/20261016T023349031430Z/numeric-verification.json) compare eight high-weight defining Euler integrals and four Pochhammer/Mellin sums with the recurrence output; all pass. The maximum Euler discrepancy is 4.36e-54. The original lower identities P10/P20/P11 are also checked exactly.

[compare_gauss9.py](compare_gauss9.py) closes the declared source algebra under all homogeneous products through weight nine, keeping both crossing parities at intermediate weights. The [full comparison](physical-runs/20261016T021005211170Z/comparison-gauss-weight9.json) gives even rank 373, increasing to 375 after the two new physical coefficients. W3/Lambda3 give the same quotient dimension. Exact column data and nonmembership witnesses are preserved.

[period_filtration.py](period_filtration.py) tests where the obstruction starts. At maximal Hlog word length nine, after discarding all positive-weight period decorations, the rank is already 225 to 227. [The complete filtration](physical-runs/20261016T021005211170Z/period-filtration.json) retains two new directions at every level. Therefore the deficit in this declared family is present in the ordinary maximal symbol; it is not solely a missing period or known normalization correction.

These are finite Q-linear source/product spaces. The tests do not cover arbitrary rational-z coefficients or every independently weighted lower-kernel source. The independent a,b,c Euler extension is tested separately below.

## Three independent Euler indices

[free_euler_family.py](free_euler_family.py) derives all coefficients from the general Gauss ODE for `2F1(a,b;1+c;-z)`, with the Gamma ratio retained. It generates 120 independent Euler moments and 210 moments with an independent Mellin index through weight nine. The [fresh run](free-euler-runs/20261016T030846505040Z/status.json) computes both crossing components in a fresh Maple session. Its 120 old diagonal-source regressions, 12 [independent defining-integral/Pochhammer checks](free-euler-runs/20261016T030846505040Z/numeric-verification.json), and 240 crossing-component regressions pass.

[compare_free_euler9.py](compare_free_euler9.py) closes all homogeneous products of the declared generators, preserving both parities at intermediate weights. [The exact comparison](physical-runs/20261016T021005211170Z/comparison-free-euler-weight9.json) gives rank 373 to 375, with two independent full-function obstructions. Previously constructed maximal-symbol duals also annihilate the enlarged candidate columns. Thus releasing the Gauss upper-parameter binding does not explain the missing directions in this registered finite space. Arbitrary independently weighted two-kernel sources are not thereby excluded.

## Reserved order and additional checks

[prediction-order4-contract.json](prediction-order4-contract.json) froze the source rules, physical mappings, source-degree budget, prefactor and original-interval-moment formulas before generating physical order four. It fixes z=0.4,0.7 as numerical validation points. The source was kept fixed when executing this prediction.

[normalization-moments.wl](normalization-moments.wl) was actually evaluated through the connected Wolfbook MCP kernel. Six exact normalization/interval identities and one omitted-moment negative control pass. [Actual tool result](normalization-moments-wolfbook.json), operation `540008a4-0075-4413-b039-233e99a47852`. All symbols are locally scoped. These are CAS checks, not Lean proof admission.

The four simple HandyG stdin checks passed. The subsequent high-weight attempt exposed the installed `geval` parser's 20-character line limit. A [workspace Fortran adapter](handyg_batch.f90) uses the unchanged installed library with explicit numeric inputs. Of 100 requested samples, only 26 returned before a prolonged 27th-word calculation was stopped: five were NaN, three finite values failed the predeclared comparison tolerance, and 18 agreed. This attempt is retained as a [diagnostic](prediction-order4/handyg-high-weight-diagnostic.json), not counted as a complete physical prediction check.

[predict_order4.py](predict_order4.py) supplies longer exact exports to the frozen physical mapper in memory, preserving its code and all old records. It includes the common endpoint cancellation, all sectors and corners, the mirror, H cubed, rho and the full Gamma/Euler/cosh convolution. [The full finite-GPL output](prediction-order4/finite-gpl-coefficients.json) contains physical orders zero through four. [Published order-four data](../../../benchmarks/2ne/theory/generation/pure-three-higher-laurent/coefficients-eps4-gpl.json) contains 38,894/39,080 exact rows for I20/I21; its explicit GPL word table avoids repeated long words. [All 77,974 rows roundtrip exactly](published-eps4.json).

[numeric_order4_original.py](numeric_order4_original.py) independently inserts the fourth logarithmic moment into the original ordered kernels. Its linear primitive retains `center^4/24+center^2*width^2/48+width^4/1920`, followed by three/four-dimensional tensor quadrature. Four missing-fourth-moment controls are rejected. [The raw original-kernel result](prediction-order4/original-quadrature.json) reads no source GPLs.

[evaluate_prediction_ode.py](evaluate_prediction_ode.py) evaluates all GPLs simultaneously through their shared suffix ODE on `x=logistic(s)`. The [precision audit](precision_prediction_ode.py) controls small-x absolute errors that otherwise propagate through long zero prefixes. The [accepted prediction report](prediction-order4/prediction-verification-r003.json) records four new physical comparisons with maximum observed discrepancy 5.79e-11, 16 old-order calibrations, and 54 classical-polylogarithm controls. Earlier less precise reports remain as separate records. [Reciprocal checks](prediction-order4/crossing-verification.json) cover 20 order-zero-through-four cases. No target-fitting parameter was introduced. Errors are empirical; no interval certificate, exact all-z weight-ten Hlog reduction, or Lean proof is claimed.

## Subsequent research questions

The physical-slice connection is now established below through a shared moment system. A stronger finite BJ-only all-epsilon replacement and an exhaustive freely weighted ordered-Q endpoint classification remain open. The finite reconstructions do not prove either statement. These are further research questions beyond the completed high-coefficient-guided exploration; no minimality conclusion is asserted.

## Independent double-Beta extension

[beta_convolution_family.py](beta_convolution_family.py) derives the independent `3F2(a,b,c;1+d,1+e;-z)` jets from its third-order ODE and origin boundary. These are zero-parameter jets of this specific family. Upper S3 and lower S2 symmetry reduce 462 raw jets to 97; one independent Mellin index gives 120 further canonical jets. [The fresh run](beta-convolution-runs/20261016T034627077157Z/status.json) and 16 [independent Pochhammer/original-Beta controls](beta-convolution-runs/20261016T034627077157Z/numeric-verification.json) pass.

The full product algebra remains [373 to 375](beta-convolution-runs/20261016T034627077157Z/comparison.json), or [372 to 374](beta-convolution-runs/20261016T034627077157Z/primitive-comparison.json) after removing all old three-kernel and E/tau generators. This negative result is for the declared parameter-origin family and products, not for all hypergeometric functions, integer parameter shifts, Horn sources or one-dimensional integral kernels.

## Complete Lambda reduction

[The derivation](../../../benchmarks/2ne/theory/generation/pure-three-higher-laurent/lambda-one-dimensional-insertions.md) eliminates the old Q endpoint moment, identifies the only bilinear endpoint term with minus the existing W Horn parent, and reconstructs every remaining term as an explicit one-dimensional Y/Z/J primitive. [Actual Wolfbook/MCP algebra](lambda-moment-wolfbook.json), [three generic polynomial original-integral controls](lambda-compression-polynomial.json), and [two noninteger full-sector quadratures](lambda-compression-numeric.json) pass. Wrong-parent-sign and missing-zero-mode controls fail as expected.

[check_lambda_compression.py](check_lambda_compression.py) starts a fresh HyperInt session and checks all 21 physical-sector differences at source degrees one through seven. [All are exactly zero](lambda-compression-fibration-N7/status.json); this includes the degree needed for physical epsilon four. The log reports Maple CPU time about 2184 seconds, not a wall-time benchmark. The earlier same-argument GPL word comparison had nonzero unsimplified differences; these vanished after full fibration. It is retained as an insufficient representation test, not reported as a failed integral identity.

## Natural weighted bilinear-J family

[bilinear_j_family.py](bilinear_j_family.py) generates raw independent coefficients of `int (1-X_s*x)^a J_b(X_t*x) J_c(X_u*x) dx/x`, with X0=1, X1=-z. The six argument patterns follow directly from the Lambda reduction; the coefficients are exact three-word shuffles with an outer zero. [The fresh run](bilinear-j-runs/20261016T040428164992Z/status.json) evaluates all 158 moments through weight nine and both crossing parts. Sixteen [original one-dimensional controls](bilinear-j-runs/20261016T040428164992Z/numeric-verification.json) pass.

[bilinear_j_single_series.py](bilinear_j_single_series.py) gives all six full parameter functions as explicit single-index sums. Their coefficients use finite Beta/Gamma moments, first-order rational recurrences and fixed boundary values; no coefficient integral is left unspecified. [Actual tests](bilinear-j-single-series-verification.json) include rational convolution identities, six noninteger source integrals and two z=2 Cayley continuations. The unit-period sum can be slow near zero exponents; numerical observations are not certified remainder bounds.

Adding the full weighted family gives [399 to 399](bilinear-j-runs/20261016T040428164992Z/comparison.json). The [unweighted J-product subfamily](bilinear-j-runs/20261016T040428164992Z/subset-plain/comparison.json) remains 381 to 383, so the independent elementary endpoint exponent matters in this comparison. The first new source directions occur in odd weight seven, [95 to 97 with exact dual certificates](bilinear-j-runs/20261016T040428164992Z/new-weight7-directions.json).

The stronger [unseeded comparison](bilinear-j-runs/20261016T040428164992Z/unseeded/comparison.json) removes every known three-kernel coefficient and old E/tau generator before building candidates. Its even ranks at weights six through nine are 49,99,199,399; all 16 physical/bare coefficients through order three reconstruct exactly. [publish_bilinear_reconstruction.py](publish_bilinear_reconstruction.py) expands every basis alias into actual source moments and [replays their products exactly](published-natural-sources.json). [Published natural-source relations and conventions](../../../benchmarks/2ne/theory/generation/pure-three-higher-laurent/natural-source-reconstruction.md).

The finite relations are long coordinate expressions, not a compact all-epsilon physical formula or minimal master basis. The current reliable whole-epsilon representation still uses the common Horn parent plus the derived one-dimensional insertions. The strict [direct parser](hlog9_direct_parser.py) is independently checked against the unchanged SymPy parser and rejects unsupported letters, periods, nonlinear words and syntax.

## Shared physical generator and z transport

[shared_physical_generator.py](shared_physical_generator.py) proves the physical Lambda parent sum is the V_n square of the existing complete moment system, and its reciprocal is the U_n square. It derives complete normalized I20/I21 projections and fixed one-dimensional bilinear-J curves for the remainder. [Exact results](shared-physical-generator/exact-verification.json) include both n parities, the two physical formulas, a wrong-mixed-sign negative and all 21 remainder coefficient comparisons through source degree seven. [Explicit remainder coefficients](shared-physical-generator/compact-one-dimensional-R-numerators.json) are auxiliary source numerators, not a new physical table.

[physical_z_system.py](physical_z_system.py) derives the z flow of the original complete endpoint moments. Four singular leading rows of the x equation have exactly zero theta_x coefficient in that flow, so closure does not invert a singular matrix. Scaling uv, v² and vZ by z at each jet level yields [the full 41-state Fuchsian connection](physical-z-transport/connection.json), with 73/82 sparse entries, only 0 and -1 finite poles, and epsilon degree at most three. These are realization states, not a master count or canonical epsilon form.

[Exact controls](physical-z-transport/exact-verification.json) cover generic symmetric-square identities, row alignment and five all-z integer cases, with forcing/weight-derivative negatives. [Actual Wolfbook/MCP](physical-z-transport/flatness-wolfbook.json) verifies all-parameter discrete/z compatibility. [Eight noninteger original-kernel comparisons](physical-z-transport/numeric-verification.json) at e=8.37,12.37 and z=0.7,1.4 have maximum scaled observed difference 5.89e-13. [Residue spectra](physical-z-transport/spectral-verification.json) have affine epsilon eigenvalues with integer offsets; this alone is not a proof of an epsilon form or MPL completeness.

[numeric_shared_physical_generator.py](numeric_shared_physical_generator.py) uses only the three F0 projections at z=1. At e=0.23,0.37, four complete physical comparisons with independent original ordered three/four-dimensional kernels pass, with maximum observed difference 1.02e-12; ten primitive-J controls and two reversed-mixed-sign negatives also pass. [Actual results](shared-physical-generator/numeric-verification.json). This does not evaluate the higher auxiliary x jets termwise near e=0.

[The mixed projection flow](mixed_projection_flow.py) further expresses theta_z M entirely in one-J Y sources. [Accepted r002](shared-physical-generator/mixed-flow-verification-r002.json) checks connected parity substitutions and four derivatives of original parent integrals. The r001 auxiliary V parity transcription and its correction are [retained explicitly](shared-physical-generator/mixed-flow-initial-check.txt); r001 is not used as the connected symbolic proof.

[The continuous proof](../../../benchmarks/2ne/theory/generation/pure-three-higher-laurent/physical-z-transport.md) specifies the z=1 boundary sums, initial Re(e)>4 endpoint operation, meromorphic continuation, noncommuting x/epsilon limits, and the Laurent recursion driven by the previous three orders. All epsilon identities arise from the original integrals; finite CAS and empirical numerical controls are distinguished from Lean proof admission. `formal=not_run`.
