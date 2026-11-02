# Validation of the complete low-loop example

The validation layers answer different questions. The actual records below keep
mathematical evaluation, representation conversion, assembly, numerical controls
and publication portability distinct.

| Layer | Actual evidence | Scope |
| --- | --- | --- |
| One-loop original geometry, parameter primitives, Gamma/moment expansion | [80 checks including 4 negative controls](../../../integrals/one-loop/validation/checks.json), [27 numerical controls](../../../integrals/one-loop/validation/numeric.json) | Fresh coefficient generation through epsilon3; standard Fourier/Euler/Mellin inputs stated in the articles |
| Two-loop mother/Gamma series and full HLOG coordinates | [149 checks and 2 negative controls](../../../integrals/two-loop/verification.json) | 12 complete tables, 84 positions; historical source integrations and dictionaries reused explicitly |
| Two-loop original/reduced-kernel and moment controls | [68 controls](../../../integrals/two-loop/integral-checks.json) | Listed ordered polynomial, Pi1, Beta, Q-kernel and M/N/tau controls; empirical mesh/quadrature, not new full raw TGV04 integration |
| Source and scalar sign migration | [Independent version bridge](version-bridge-checks.json) | Native Pi1 r002 versus legacy r001 at C/J/domain/operator and individual coefficient levels |
| Full-color/scales/counterterms/RG and author finite result | [68 checks including 5 negative controls](../checks.json) | Exact color/pole/finite comparisons, separate cusp/noncusp and RG-blind finite direction |
| Actual integral engines before downstream assembly | [Fresh regeneration run](fresh-regeneration.json) | Executes both declared child routes; not a new fourteen-integral CAS derivation |
| Finite function reading | [12 new two-precision samples](../numeric-samples.json) | Original M11 integral versus polylogs, external finite comparison and crossing at 40/60 digits |
| Independent isolated reproduction | [Read-only review and runs](isolated-reproduction.json) | Copied declared closure, isolated imports, blocked network and private source reads, new output paths |

The [source manifest](../source-manifest.json) binds unmodified historical totals
and self-authored comparison programs. They are holdouts/provenance, not sources
for constructing the new individual tables. [Reference provenance](../reference/provenance.json)
attributes the external derived fixture and identifies optional original author
files. Those downloaded files are not distributed here.

The complete [publication manifest](../../../LOW-LOOP-PUBLICATION.json) binds the
immutable low-loop files. Mutable reader navigation is checked separately and
versioned by Git. The original upstream full-tree manifest remains tied to its
original commit; older mathematical evidence is not silently rebound to new files.

All listed checks ran in their stated scope. A successful run or checksum is not
a proof of the original QFT calculation. In particular the standard graph rules,
Abelian shuffle identity and explicitly identified historical integrations remain
research inputs. General complex-z continuation and a Lean integral chain are
outside this packet: `formal=not_run`.
