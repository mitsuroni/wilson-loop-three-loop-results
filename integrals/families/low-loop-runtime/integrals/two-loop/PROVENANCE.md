# Which calculation is rebuilt, and which integration is reused?

The twelve `result.json` files are in one native physical convention. Their
current individual source records were checked coefficient by coefficient
through the finite term; their positive orders were checked against the
preserved higher-order report and complete available HLOG data.

| Object | Construction in this publication | Earlier mathematical work adopted |
|---|---|---|
| Six pure-exchange coefficients | New Gamma/V/VN coefficient recurrence and physical prefactor convolution, all seven powers | Ordered endpoint integrations, evaluated V/VN and Q identities, fixed M/N HLOG dictionary |
| TGV01/02 | New series from their full Gamma/Beta expression | Same-edge path/simplex integrations and the actual reflection map |
| TGV03 | New contact-preserving single-kernel/Beta/Gamma expansion | Ordered-domain IBP, contact-face identity and fixed moment source |
| TGV04 | Exact current-source and complete HLOG-to-MN/tau comparisons | Original six-face plus subtracted-recoil HyperInt integrations, including D2/H1/O1 and D3/H2/O2 |
| Pi1 adjacent/opposite | New expansion of native complete-I formulas and scalar/color sign checks | Full transverse insertion/Fourier derivation, elementary/corner and opposite-kernel evaluations |

`source-manifest.json` binds the71 ordinary source/dictionary/reference files
needed by `reproduce.py`. `CAS-archive-manifest.json` separately binds21
additional original Maple save files and preparation scripts needed to inspect
or rerun the archived TGV04 sessions. Both manifests preserve original bytes;
their scopes differ. The archive files are not secretly executed by the Python
replay.

The baseline `canonical_source` in the machine catalogue is the individual
accepted Laurent result through epsilon0. The epsilon¹/epsilon² values are
newly expanded from the evaluated mother formula for eleven diagrams, with
the previous report serving as a comparison. TGV04's positive orders instead
consume its explicitly declared historical integrated compact/HLOG results.
The exact reference pointers are recorded in `representation-sources.json`.

The new checks are deliberately broader than copying the table: original
ordered pure kernels at positive integer epsilon, full Pi1 mixed derivatives
and corners, the surviving original TGV01 Beta integral, original real M/N/tau
definitions, and 4e04's original two-variable Q0/Q1/Q2 coefficient kernels are
actually evaluated. They support the retained derivation chain. They do not
claim that all twelve original noninteger-epsilon multidimensional integrals
were independently reintegrated in this release.

The native Pi1 r002 I already carries the physical negative sign. Its assembly
multiplier is +1 with the full epsilon-dependent B1 weight. The historical
r001-to-r002 multiplier -1 is provenance metadata, not another operation to
apply to the released coefficients.

Original report builders and CAS files under `sources/` retain old hard-coded
roots and writers. They are archival scientific sources. Use the packet's
`reproduce.py` and `check_integrals.py` for the portable checks. Do not execute
an old writer against an immutable snapshot. No Maple/Mathematica run, new
whole-integral Lean proof, or new global complex-domain continuation is claimed.
See `VALIDATION.json` for the actual fresh checks and scope.
