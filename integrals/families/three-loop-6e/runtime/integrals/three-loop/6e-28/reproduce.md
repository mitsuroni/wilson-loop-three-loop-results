# Reproduce 6e-28

From the package root, run the shared 39-orbit batch into a NEW directory:

```powershell
python -B scripts/evaluate_batch.py --output runs/replay-new
python -B scripts/verify_batch.py --results runs/replay-new/generated.json --out runs/replay-new/verification.json
python -B scripts/independent_controls.py --results runs/replay-new/generated.json --input engine/benchmarks/2ne/prepared/6e.json --out runs/replay-new/independent-controls.json
```

Read the `6e-28` member of generated.json. All seven requested positions are known; positions below the declared lower_bound are established zero. One batch shares primitive requests across graphs. The default command replays archived original CAS text and freshly assembles I. It does not require Maple and does not constitute another cold CAS integration. Install the versions declared by requirements-lock.txt. Preserve original primitive scripts/logs and their manifest hashes; an altered or missing primitive is rejected.
