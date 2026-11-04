# Reproduce 4e-Pi1-11

From the repository root, use CPython3.12 with the [declared requirements](../../families/three-loop-insertions/runtime/requirements.txt), SymPy1.14.0 and mpmath1.3.0. The configured local interpreter is D:/py/python.exe. Choose a NEW output directory:

```powershell
D:/py/python.exe -X utf8 -B integrals/three-loop/4e-Pi1-11/code/reproduce.py --output replay-output/4e-Pi1-11-new
```

The own adapter runs the whole fifteen-source batch. Select this member in generated.json; it also emits other members. It preserves the byte-bound historical sources and writes into the chosen new directory.

The original-root source identities are checked by [pi1.py](../../families/three-loop-insertions/runtime/scripts/pi1.py), which expands the six precisely adopted endpoint/face reduction formulas, Gamma/rational factors, solved M/T coefficients and full B1 colour convolution. It does not derive all six reductions afresh from arbitrary topology and does not perform a new integration of the internal gluon/ghost/quark bubble or the adopted M/rho1 functions. The mathematical derivation exposes those source-to-formula bridges explicitly.

Original Wolfram/native work supplies the accepted C/J; original Maple/HyperInt evaluated the stated moments; Python performs the supported finite-series reconstruction. Historical writers keep old absolute paths and are not launched by this command. [data/index.json](data/index.json) specifies actual sources, integrated subobjects, assembly code and post-generation comparisons separately.

The actual new batch passed105 scalar/geometric positions,84 stored historical colour positions and343 representation checks covering all238 generated colour positions. These are totals for15 entries, not independent tests of this one original integral. See the [actual run](../../../workflow/checks/integral-documentation-insertions-20261104/evidence/primary/generation/run.json), [verification](../../../workflow/checks/integral-documentation-insertions-20261104/evidence/primary/verification.json) and [representation checks](../../../workflow/checks/integral-documentation-insertions-20261104/evidence/primary/representations.json). Generation reads no whole-I reference table and opens no new CAS process. The follow-up checker files are separate from generation.

The [batch report](../../../workflow/checks/integral-documentation-insertions-20261104/replay.json) records exact commands, copied-closure checks and limitations. Use a fresh directory on rerun; preserve failed output for diagnosis. No new cold-CAS checkpoint/resume mechanism is implied by this finite Python assembly. Scientific continuity is audited independently of these execution receipts.
