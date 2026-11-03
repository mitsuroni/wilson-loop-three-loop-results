# Reproduce 5e-3g-35

From the packet directory with the locked Python packages installed, choose a new output directory:

```text
python -B scripts/evaluate_batch.py --output runs/new-batch
python -B scripts/verify_batch.py --results runs/new-batch/generated.json --out runs/new-batch/verification.json
python -B scripts/check_representations.py --candidate runs/new-batch/generated.json --references references/all-results.json --output runs/new-batch/representations.json
python -B scripts/independent_controls.py --results runs/new-batch/generated.json --output runs/new-batch/independent-controls.json
python -B scripts/verify_manifest.py
```

Read the `5e-3g-35` member. Generation reuses archived primitive evaluations; the final checker then reads the copied reference objects. The original Maple scripts/output/logs are retained in snapshots, but a cold Maple/HyperInt rerun and package installation are outside this Python replay.
