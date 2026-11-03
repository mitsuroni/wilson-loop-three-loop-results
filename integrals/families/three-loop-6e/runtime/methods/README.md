# The six evaluation routes for the 39 three-loop exchange integrals

The common problem is to evaluate each accepted ordered parameter integral `J`, restore its own colour and symmetry normalization, and obtain every physical `I` coefficient from epsilon power -6 through 0. The diagrams share six kinds of mathematics, so they share derivations and programs while retaining individual source and result records.

Start with [the normalization and interpretation](common-normalization.md), then follow the route selected by the actual propagators and bounds.

| Route | 6e orbit numbers | Count | Subintegrals needed through the finite term |
| --- | --- | ---: | --- |
| [Adjacent chains](empty.md) | 01,02,03,07,08,09,17,18,19,35 | 10 | Beta/Gamma words and scalar Taylor coefficients through order 6 |
| [One opposite propagator](single-opposite.md) | 04,05,06,13,14,15,16,26,27,28,29,30,31,36,39 | 15 | Centered two-variable moments through Taylor order 4 |
| [Two opposite propagators in one channel](same-channel-two.md) | 10,11,12,23,24,25 | 6 | Ordered gap kernel through Taylor order 2; its order-0 cancellation is retained |
| [Three opposite propagators in one channel](same-channel-three.md) | 20,21,22 | 3 | Original ordinary integral at epsilon=0 |
| [One opposite propagator in each channel, with an adjacent pair](double-opposite.md) | 32,33,34 | 3 | Ordered gap kernel through Taylor order 2, or the equivalent H0/H1 product convolution |
| [One and two opposite propagators](mixed-opposite.md) | 37,38 | 2 | Original ordinary integral at epsilon=0, equivalently H0 times Qid/Qswap at zero |

The portable replay recomputes chain integrals and all physical assembly. For routes that require HyperInt, it reuses the archived primitive execution only after checking the mathematical request, generated script, successful execution record, raw output and parser. This is a replay of existing integration evidence. It is not another independent CAS integration, and agreement with the old final table is a separate check performed after candidate generation.

The following equalities have different meanings. They do not remove the individual entries from the batch.

| Relation | Established scope | Reason |
| --- | --- | --- |
| 10=23, 11=25, 12=24 | Complete physical I | Their original J functions are related by crossing, with the same checked normalization |
| 21=22 | Complete physical I | Crossing exchanges R132 and R213; their original J functions generally differ |
| 32=34 | Original J and complete physical I | The two separately integrated channel factors exchange; both kappa factors are -1/16 |
| 37=38 | Finite coefficient only | Different Q kernels and different kappa factors coincide only at the displayed order |

The precise selected coefficients, source bindings, subintegral records and checks belong to the per-orbit data. A family identity or a successful replay does not by itself supply a Lean proof of an original integral. This batch has `formal=not_run`.
