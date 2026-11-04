# Reproduce 4e-Pi1-03

From the packet directory, use a new output directory:

```text
python -B scripts/evaluate_batch.py --output runs/new-batch
python -B scripts/verify_batch.py --results runs/new-batch/generated.json --out runs/new-batch/verification.json
python -B scripts/check_representations.py --results runs/new-batch/generated.json --out runs/new-batch/representations.json
python -B scripts/independent_controls.py --results runs/new-batch/generated.json --out runs/new-batch/independent-controls.json
```

Read the `4e-Pi1-03` member. The batch shares six Pi1 bodies and two chain/Pi2 geometries while retaining each original source and distinct operator. Install the versions in requirements-lock.txt. It requires no CAS or private source path; the original local paths in archived provenance are historical observations. The final checker reads the copied final tables only after generation.
