# Gamma1 source contract (2026-10-03)

The four 3e-Gamma1 orbits now export a complete effective, bare one-loop 1PI three-gluon vertex source. Each source contains three external propagators and phases, l1+l2+l3 conservation, ordered contour bounds, closed Lorentz/colour indices, total g0^6 order, C4 metadata and a versioned full-D operator definition.

The operator comes from Davydychev, Osland and Tarasov, hep-ph/9605348v3, DOI 10.1103/PhysRevD.54.4087, with the 1999 erratum DOI 10.1103/PhysRevD.59.109901. Use (2.4), (3.4)-(3.9), (3.17)-(3.22), D=4-2eps, Feynman gauge xi=0 and Tr_Dirac(1)=4. The common g0^2 eta/(4pi)^(D/2), CA or TF nf, and scalar master definitions are explicitly separated in gamma1_operator.py; the full definition is embedded in the JSON sources and hash-bound.

The project convention is V1=i Gamma_DOT, matched to the tree rule g0 f times the kinematic tensor. This is an explicit convention dictionary; independent complete one-loop phase validation remains open. The source contains the one-loop increment, excludes external Pi1 insertions and counterterms, and must not be multiplied by another internal loop measure.

central_vertex.validate_definition checks the canonical source contract, including full operator content and its reference. Changes to momenta, phases, indices, bounds, ordering, operator coefficients or normalization are rejected.

The dedicated three-point route bypasses the legacy shared-momentum template. Other higher-point/higher-loop effective vertices are outside this change. Enumeration and mirror metadata remain unchanged; the four source orbits are not merged.

This completes Stage 1 in the stated operator convention. Montopar has no Gamma1 C/J implementation yet and explicitly rejects these sources instead of routing them to pure/tree TGV. I remains unsolved.
