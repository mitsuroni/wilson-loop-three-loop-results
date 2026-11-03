# Three-loop tree-three-gluon evaluation packet

This local P6 packet curates all44 `5e-3g` sources in eight actual evaluation families. Each entry connects its original momentum components and native/ordered C/J to all seven physical Laurent functions epsilon^-6 through epsilon^0. Shared derivations and programs avoid repeating the same integration route44 times.

Start with the [common momentum/C/J/I method](methods/README.md), the [31 T1/T2/T3 routes](methods/t123.md), the [13 T4 routes](methods/t4.md), [independent source/finite controls](methods/independent-controls.md), and [function back-substitution](methods/representations.md). The [per-source index](integrals/three-loop/README.md) links every own article and machine-readable result. [catalogue.json](catalogue.json), [inputs.json](inputs.json) and [manifest.json](manifest.json) bind the portable data.

Existing HyperInt primitive evaluations are explicitly reused. The present run recomputes exact assembly; it does not claim44 independent CAS integrations. In particular27/28 reuse original fully integrated global-IBP output, whereas the corrected T4-empty results are reassembled from raw Q jets with fresh inversion and the correct absolute epsilon budget. Individual27/28/29 extra letters and periods remain present through the finite result.

With the packages in requirements-lock.txt installed, choose a new output directory:

```text
python -B scripts/evaluate_batch.py --output runs/new-batch
python -B scripts/verify_batch.py --results runs/new-batch/generated.json --out runs/new-batch/verification.json
python -B scripts/check_representations.py --candidate runs/new-batch/generated.json --references references/all-results.json --output runs/new-batch/representations.json
python -B scripts/independent_controls.py --results runs/new-batch/generated.json --output runs/new-batch/independent-controls.json
python -B scripts/verify_manifest.py
```

The final checker reads copied reference objects after generation. Original scripts/logs/raw outputs and evaluated dictionaries are present under snapshots; a cold Maple/HyperInt run requires the original external CAS/packages and is outside this Python replay. Active paths are packet-relative. Source/profile/old packets and public GitHub bytes remain unchanged. Positive-real-z scope and formal=not_run are explicit; see evidence/generation and evidence/isolation for actual execution records. This packet has not been committed or uploaded.
