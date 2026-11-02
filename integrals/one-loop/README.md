# One loop: from the two exchange graphs to evaluated integrals

A null quadrilateral gives the simplest place to see every normalization
boundary before attempting two or three loops. The adjacent exchange has
two endpoint poles, while the opposite exchange is finite at epsilon=0.
Both start with the same Fourier propagator and the same scalar prefactor;
their different behavior comes from the geometry of the integration domain.

| Orbit | Geometry / C4 multiplicity | Continuous derivation | Exact normalized I | Fresh checks |
| --- | --- | --- | --- | --- |
| [2e-01](2e-01/README.md) | Adjacent edges3/4; m=4 | [Derivation](2e-01/derivation.md) | [JSON](2e-01/result.json) | [Record](2e-01/checks.json) |
| [2e-02](2e-02/README.md) | Opposite edges2/4; m=2 | [Derivation](2e-02/derivation.md) | [JSON](2e-02/result.json) | [Record](2e-02/checks.json) |

Each result contains the **normalized physical I**, including its C4 orbit,
crossing combination and exactly one Euler compensation. CF, the coupling
lambda=alpha_s/pi and the common scale factor S1 are removed. In particular,
do not multiply either I by m again. At this order W-1 and log W agree.
The full coefficient is restored as CF lambda S1(I_adjacent+I_opposite), with
S1 retaining its own one-loop scale convention.

The new execution starts with the frozen graph, native raw C/J and original
ordered/prepared inputs. It checks their phase/Gram geometry and exact C
connection, evaluates the original parameter primitives, and generates the
Gamma/M-moment convolution through epsilon3. Predictions are saved before
the old coefficient tables are opened. The published fields are ordinary
Laurent coefficients without factorials, including explicit pole zeros.

The extra positive orders have a purpose: one-loop epsilon2 enters the
three-loop finite counterterm, and one-loop epsilon3 is available for the
two-loop post-counterterm positive-order budget. These orders were actually
generated and checked, not appended as a requested but uncomputed budget.
An all-epsilon Gamma answer and a standard all-epsilon Meijer G answer are
also provided. The documented M generating formula determines any fixed
higher order, but this release's explicit checked table stops at epsilon3.

The [shared reproducer](reproduce.py) needs only Python, SymPy and mpmath;
[versions](requirements.txt) are pinned. Run from the repository root:

```text
python integrals/one-loop/reproduce.py --output-dir one-loop-replay --order 3 --dps 40
```

The output directory must be new and outside the input/result trees. It
contains `2e-01/result.json`, `2e-02/result.json`, a generated-before-holdout
snapshot, numerical readouts and the complete check record. It never invokes
the archived source scripts or writes into the packet/data/packages trees.
See the individual [adjacent](2e-01/reproduce.md) and
[opposite](2e-02/reproduce.md) reproduction notes.

The actual successful execution is preserved in
[validation/checks.json](validation/checks.json) and
[validation/numeric.json](validation/numeric.json).80 checks, including four
negative controls, and27 numerical controls passed. The largest observed
scaled numerical discrepancy was about8.88e-27 at40 working digits. This is
an empirical precision control, not a certified interval. Source contracts,
the Fourier propagator rule and standard Euler/Mellin identities remain
explicit mathematical inputs; no integral Lean proof is claimed.

The [historical source manifest](source-manifest.json) binds unchanged
research derivations, original result objects, old check records and relevant
calculation scripts. They establish provenance and comparison targets;
they were not run in place or used as the source of the new coefficient
table. General physical complex-z continuation remains outside the claimed
domain. `formal=not_run`.
