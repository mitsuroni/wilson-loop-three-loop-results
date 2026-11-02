# Current one-to-three-loop Wilson-loop result

The current calculation provides complete one/two-loop results and all four
three-loop color channels through epsilon0, with direct native-source assembly
and RG validation. Read the exact [one-loop](../../benchmarks/current/results/1loop.json),
[two-loop](../../benchmarks/current/results/2loop.json) and
[three-loop](../../benchmarks/current/results/3loop.json) files, or run:

```text
python -m workflow.current_calculation result 3
python -m workflow.current_calculation replay --output-dir replay-output/current-new
```

With lambda=alpha_s_MSbar/pi and T=TF nf,

$$\log W=\sum_{n\ge1}\lambda^n w_n,\qquad
w_3=C_F\left[C_A^2w_{AA}+C_ATw_{Af}+C_FTw_{Ff}+T^2w_{ff}\right].$$

Each channel contains every power from epsilon^-6 through epsilon0. The files
separate U, bare B=S3 U, coupling CT, post-counterterm P=B+CT, pole logZ, finite
F, cusp k3 and noncusp h3. Keep L=(Ls+Lt)/2 distinct from rho=log(z). The finite
functions belong to log W_R, not the ordinary Wilson-loop coefficient.

The three-loop counterterm is

$$CT_3=-\frac{2b_0}{\epsilon}B_2+
\left(\frac{b_0^2}{\epsilon^2}-\frac{b_1}{2\epsilon}\right)B_1,$$
$$b_0=\frac{11C_A-4T}{12},\qquad
b_1=\frac{34C_A^2-20C_AT-12C_FT}{48}.$$

B1 and B2 keep their own S1/S2. The full finite relation is

$$\partial_LF_3-2b_0F_2-b_1F_1=-2k_3L-h_3,$$

with the corresponding color extraction. Compact function-family coefficients
and independently constructed native common-z HLOG expressions are both stored.
All published full-color coefficients reduce to M/N, zeta and physical scale L;
the complete channel needs no independent tau/E or explicit rational-z coefficient.
This does not assert the same closure for every individual source or positive order.

Af is directly reconstructed from31 current native orbit integrals. AA reads83
tree orbits,31 insertion orbits and35 four-leg sources grouped into nine joint
objects. Extra5e letters/periods are retained before their exact cancellation.
AA's MN coefficients are independently obtained by rational basis conversion;
old totals are holdouts, not construction inputs.

The [source catalogue](../../benchmarks/current/catalog.json),
[source/joint ledger](source-ledger.json), [matter ledger](matter-ledger.json),
[single convention](../../workflow/current-conventions.md) and
[actual validation](../../validation/current-2026-11-02/README.md) explain the
connections. Existing I2344, new admitted I1234 and the other original integral
derivations/dictionaries retain their prior research acceptance and reuse scopes.
This release is not a fresh integration of every original kernel or a formal proof.
RG cannot fix an arbitrary L-independent finite function.

The [23 October archive](../2026-10-23/README.md) and
[low-loop pilot](../2026-11-02-low-loops/README.md) preserve their original dated
scope. Use the current entrypoint for new reading and assembly; detailed English
three-loop evaluation articles remain the next publication layer.
