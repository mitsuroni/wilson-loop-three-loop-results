# Three-loop 6e: 39 complete integral entries

This packet connects the original labelled graphs and ordered C/J roots to all Laurent coefficients of their physical I through epsilon^0. There are 39 individual English articles and six shared mathematical routes. The completed 6e-01 pilot is reused; this batch supplies the other 38 entries and repeats its inexpensive Gamma calculation as a consistency check.

The result is a reproducible curation and reassembly of existing research mathematics. Ten adjacent-chain integrals are evaluated and expanded afresh. The remaining routes reuse actual archived CAS integrals after exact request/script/execution/raw-output/parser checks, then repeat the reduction, order budgeting and physical assembly. No new CAS integration or formal proof is claimed.

## Follow the data

1. [The catalogue](catalogue.json) selects an individual article. Each entry has its own graph/source bindings, ordered J with dependent bounds and original order, normalization, subintegral budget, seven coefficient positions and checks.
2. [The method index](methods/README.md) explains the six routes, their shared subintegrals and the precise scope of aliases. Start with [common normalization](methods/common-normalization.md) for the C/J-to-I connection and function conventions.
3. [The source manifest](inputs/source-manifest.json) binds the graph, original ordered source and native full-colour cube presentation separately. Sources below inputs/source are unchanged whole-package bytes; readable per-orbit bindings point to actual fields. The accepted projected ordered roots remain in [the engine input](engine/benchmarks/2ne/prepared/6e.json).
4. [The frozen generated candidate](evidence/generation/generated.json) was produced without reading final tables. Its [run](evidence/generation/run.json), [actual primitive consumption](evidence/generation/primitive-replay.json) and [pre-comparison freeze](evidence/generation/generation-freeze.json) preserve the execution record.
5. [Exact verification](evidence/generation/verification.json) covers 273 Laurent positions, 273 function substitutions, every orbit's source bridge, full-colour reconstruction and normalization. [Independent controls](evidence/generation/independent-controls-r002.json) and [their scope](methods/independent-controls.md) are recorded separately.

| Route | Individual entries | Active computation |
| --- | ---: | --- |
| Adjacent chains | 10 | New chain/Gamma expansion |
| Single opposite | 15 | Archived centered moments, new original-root assembly |
| Same-channel two opposite | 6 | Archived ordered gap jets, new assembly |
| Same-channel three opposite | 3 | Archived original zero-order kernels, new assembly |
| Double opposite | 3 | Archived gap jets, new assembly |
| Mixed opposite | 2 | Archived original zero-order kernels, new assembly |

All result positions epsilon^-6 through epsilon^0 are explicit, including established zeros. The complete HLOG answers use 0,-1 and zeta periods; the compact representations use the supplied M/N and, when needed, tau/E definitions. An exact independent rational shuffle verifies each compact formula. Naming a function is not used as a substitute for solving its coefficients.

## Reproduce the batch

Use CPython 3.12 with the versions in [requirements-lock.txt](requirements-lock.txt). The actual validation used CPython 3.12.10 on Windows and the already installed declared libraries. No online installation or second operating system was tested. The 27 copied project modules are unchanged; the four original manifest/input files are retained because the original reader validates all 47 roots in that manifest.

From this directory, choose a new output directory:

```powershell
python -B scripts/verify_manifest.py
python -B scripts/evaluate_batch.py --output runs/replay-new
python -B scripts/verify_batch.py --results runs/replay-new/generated.json --out runs/replay-new/verification.json
python -B scripts/independent_controls.py --input engine/benchmarks/2ne/prepared/6e.json --results runs/replay-new/generated.json --out runs/replay-new/independent-controls.json
```

The default batch needs no installed Maple or HyperInt. Backend paths and hashes in the raw archive are historical provenance used to regenerate and match the original scripts; these paths are not executed or opened. The packet includes all 388 original primitive/source files, about 11.2 MB. Admission covers 71 scripts; this 39-orbit execution consumes 69 units and 174 distinct raw channel files. Two admitted low-loop units are unused. These numbers are not counts of independent master integrals.

The evaluator rejects reads outside its declared code/data/library/output closure and forbids external processes and network connections. These Python audit checks are not an operating-system sandbox. The in-process finite Gamma expansion does not add an enforced memory/time bound. Historical raw scripts and old orchestration writers remain evidence, rather than current production commands.

[The actual isolated replay](evidence/isolation/report.json) initially contained 425 generation files and no final references or results. After generation and freeze, reference/source data were added for verification. All 39 coefficient payloads were identical and all 273 positions passed again. This is a copied-closure replay in the same installed Python environment; absolute provenance locations and whole result JSON hashes can differ after relocation.

## What the checks establish

All 39 source colours, orbit factors and normalizations are checked individually. Fourteen paired polynomial source identities at epsilon=1,2 compare the original iterated integral with a separately constructed cube/Jacobian monomial integration for seven representatives spanning all six routes. These identities do not determine the Laurent finite terms of divergent integrals.

The holomorphic sources 20/21/22/37/38 additionally have new original six-dimensional finite-integral controls at z=1/2: eight randomized Halton shifts of 16384 samples each, compared with a separately implemented HLOG series evaluator. They pass the declared empirical tolerance; observed relative residuals are approximately 0.028% to 0.532%. Neither empirical errors nor series convergence constitute rigorous error enclosures. The existing 6e-01 pilot has its own independent controls outside this packet.

The equality 37=38 applies only to c0. Direct exact source integration gives (I37/Gcal^3-I38/Gcal^3) at epsilon=2,z=1/2 equal to 1161/266240, excluding a full-epsilon identity. Other shared full-I statements retain their original derivations and are not proved merely by matching seven coefficients.

Candidate generation precedes reference export in the recorded controller sequence. Two checker scripts were completed after that freeze, while the generator/provider and input/primitive bytes stayed unchanged; the final manifest binds their final versions. The freeze is preserved as historical evidence, rather than silently updated to a later checker. This is not a trusted timestamp or a claim of human blindness.

Positive real z and local meromorphic continuation at epsilon=0 are the mathematical scope. A new direct z=-1 boundary calculation, arbitrary complex i0 continuation, a second independent analytic integration of all 39 diagrams and an integral-level Lean proof are outside this packet. Formal state remains not_run. No scientific selector or upstream source is changed by the packet, and it has not been committed or uploaded.
