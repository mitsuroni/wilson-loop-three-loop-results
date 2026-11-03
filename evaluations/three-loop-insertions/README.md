# P5: fifteen two-point and exchange-insertion integrals

This packet contains eleven4e-Pi1, two serial-Pi1 and two single-Pi2 source entries, with their own English articles, original graph/raw/native C/J bindings, complete scalar Laurent coefficients and full external colour convolution. The six Pi1 bodies and two chain/Pi2 geometries share code and continuous mathematical explanations.

[The catalogue](catalogue.json) lists all15 articles. Start with [the common interpretation](methods/README.md), [Pi1 plus exchange](methods/pi1.md), [chains and Pi2](methods/chains-pi2.md) and [independent controls](methods/independent-controls.md). Each own entry gives the original ordered kernel/bounds, normalization, reduction and order budget, every scalar coefficient epsilon^-6 through0, each coloured finite coefficient, and machine-readable complete coloured series.

## What is computed and reused

The active scripts read original sources and existing evaluated M/HLOG function dictionaries, then repeat the root reduction, Gamma/rational expansion, crossing/orbit assembly and complete operator convolution. They do not read the copied final result tables during generation. Original internal-bubble rules, evaluated moments and the adopted full-D Pi2 literature operator are reused explicitly; this is not a second independent analytic calculation of every original loop graph.

The adapter is limited to these frozen source roots. Its six Pi1 route assignments reuse the original family derivations and verify the bound source objects; they are not presented as a new general topology classifier or automatic integration method for arbitrary diagrams.

The current native C determines signs before integration. Euler compensation is consumed once. The unprojected full trace remains in the input; connected I and external colour factors have the exact scope stated in the family method. The Pi1 branch's absent Ff/ff connected channels are not statements that its unprojected CF-squared trace vanishes. A serial Pi1 chain and one Pi2 insertion have distinct operators even when sharing the endpoint geometry; the current geometric sign relation follows from their actual C factors.

The result has105 complete scalar/geometric positions and238 explicit nonzero-channel positions, including structural epsilon zeros within those channel tables. Every one of these343 expressions has both complete0,-1HLOG/zeta coordinates and a named-function expression. M01/M11/M21/M03 abbreviations are identified with the original M moments; the distinct rho_id_0_0_0/1 functions retain their complete evaluated dictionary, and rho means log(z).

## Run the declared batch

Use CPython3.12 and the versions in [requirements-lock.txt](requirements-lock.txt). The actual run used CPython3.12.10, SymPy1.14.0 and mpmath1.3.0 on Windows; no fresh package installation or second operating system was tested. From the packet root, choose a NEW output directory:

```text
python -B scripts/verify_manifest.py
python -B scripts/evaluate_batch.py --output runs/new-batch
python -B scripts/verify_batch.py --results runs/new-batch/generated.json --out runs/new-batch/verification.json
python -B scripts/check_representations.py --results runs/new-batch/generated.json --out runs/new-batch/representations.json
python -B scripts/independent_controls.py --results runs/new-batch/generated.json --out runs/new-batch/independent-controls.json
```

No installed Maple, HyperInt, private C-drive source tree or proof database is required. Snapshot files preserve original bytes. Historical writer scripts and Markdown in snapshots are archival evidence; their old local paths/statuses do not become current commands. The enclosing original workspace's staging/preservation scripts are local-audit tools; scripts in this packet are the portable entrypoints.

The source inventory in [inputs.json](inputs.json) separates native/root sources from final reference bindings. Existing final references are in references/ and are opened only by checks after generation. Preflight hashes original historical material, so this is not an answer-free filesystem or human-blind discovery claim. The generator's mathematical read closure excludes final result JSON and the old Pi2 operator-expansion table; Pi2 multipliers are regenerated from its original A/B source.

## Actual evidence and its limits

[Generation](evidence/generation/run.json) completed15/15 in approximately7.49seconds, with no external CAS and no access denials. The [candidate freeze](evidence/generation/generation-freeze.json) preceded final reference export. The [exact comparison](evidence/generation/verification.json) passes105 scalar positions and84 explicitly stored historical coloured positions. The Pi1 references store scalar I; its154 coloured positions are newly reconstructed with the full source operator rather than falsely labelled old coloured-table entries.

[Independent rational-shuffle back-substitution](evidence/generation/representations.json) passes all343 scalar/colour function expressions. [Source/operator controls](evidence/generation/independent-controls-r002.json) independently check all15 tensor numerators, source operator/colour connections and238 coloured convolutions,30 paired original-order/cube polynomial integrations, six adjacent Gamma coefficient Cauchy controls at two radii, and six negative examples. The largest observed Cauchy residual is approximately2.4e-71; this is an empirical comparison, not a certified numerical enclosure.

The other nine roots have no newly run original-multidimensional finite numerical control in this batch. Positive-integer-epsilon source identities alone do not determine Laurent finite terms. No new complex i0 continuation, direct z=-1 calculation or integral-level formal proof is claimed. formal=not_run throughout.

[The actual copied-closure replay](evidence/isolation/report.json) starts with118 declared code/source/dictionary files and no final reference directory. After generation, final references/source-reading snapshots and checkers are added. All15 generated results are byte-identical to the primary candidate, and105 scalar plus343 representation checks pass again in that copy. The run uses the same installed Python environment. Python audit hooks are application-level checks, not an OS sandbox; this finite in-process batch adds no separately enforced memory/wall bound.

The packet changes no adopted scientific coefficient, selector, normalization convention or public repository. Its status at construction is local curation/replay complete, with commit/publication controlled by a later explicit user request.
