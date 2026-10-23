# Reproduce the checks

Open index.html directly for the offline report. Keep the assets and sources directories beside it. No web connection is needed to read or typeset the formulas.

For the mathematical checks, use Python 3.12 with SymPy 1.14.0 and mpmath 1.3.0. From the extracted package directory run:

```text
python reproduce/reproduce.py --output-dir ../matter-replay
```

The command verifies the source hashes, then writes new outputs outside the frozen report. It replays the existing 452 checks, rebuilds the compact tables and diagram ledger, compares the external UV data and evaluates numerical samples. It does not launch Maple or Mathematica. The mixed r-to-z HyperInt assembly witness is reused with its original source bindings.

The fresh check output copied from the October 17 checker retains that checker's historical `date` field and conservative capability flags. The execution is recorded here as 23 October 2026; the newer UV and finite-source comparisons are separate evidence in uv-checks.json. A passing replay does not create a new formal proof or independently integrate the Af finite shape.

`coefficients.json` uses exact SymPy expression strings, including `M_0_1` and `zeta(3)`. The conversion dictionary fixes all symbols and source GPL maps. Source snapshots are byte-preserving copies; package/source hashes do not confer mathematical proof status.

Data: [coefficients](coefficients.json), [diagram ledger](diagram-ledger.json), [conversion dictionary](conversion-dictionary.json), [source manifest](source-manifest.json), [compact checks](coefficient-checks.json), [UV checks](uv-checks.json), [numerical samples](numeric-samples.json), [historical checker replay](rg-replay/checks.json).
