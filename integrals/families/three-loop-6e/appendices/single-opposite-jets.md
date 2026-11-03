# Complete evaluated single-opposite HLOG jets through order4

For each H_i in the [shared derivation](../single-opposite.md), write
H_i(epsilon,z)=sum_{q=0}^4 h_iq(z)epsilon^q+O(epsilon⁵).
These are ordinary coefficients: the1/q! belongs in the integrated Taylor
kernel(S5), and no further factorial multiplies these answers. In this
appendix G(w;z) is the HLOG with ordered word w in letters0,-1, G(empty;z)=1.
The real z>0 branch and exact zeta periods are fixed throughout.

Every displayed function is a complete already-integrated coefficient, not
an implicit definition by coefficient extraction. The native HyperInt
raw.txt was parsed into exact arithmetic and typeset; the archived file is
linked for the immutable original syntax, request, compute.mpl and log.
A link alone is not the evaluation: the complete expression follows it.

For reflected jets use h_iq(1/z) on this same real branch. The runtime also
retains native reflected.txt with its exact positive-z fibration. H0,HW,HX
are crossing even by the actual variable substitutions. H1,H2,HD need both
channels in physical pairing. HX's native profile(0,1;1,0) equals the defined
(1,0;0,1) by t<->s, with no scale or factorial change.

## H0: parameters (0, 0, 0, 0)

Order0; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-27c948aee63e43ff915d1d184389af99/order-0/raw.txt).

$$h_{0,0}(z)=- G(0,0;z) - \frac{\pi^{2}}{2}.$$

Order1; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-0272749e13a64b9583e0491f9f6cd666/order-1/raw.txt).

$$h_{0,1}(z)=\frac{\pi^{2} G(-1;z)}{2} - \frac{G(0,0,0;z)}{2} - \frac{\pi^{2} G(0;z)}{4} + G(-1,0,0;z) + \zeta\left(3\right).$$

Order2; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-0272749e13a64b9583e0491f9f6cd666/order-2/raw.txt).

$$h_{0,2}(z)=- \frac{G(0,0,0,0;z)}{2} - \frac{\pi^{2} G(-1,-1;z)}{2} + \frac{\pi^{2} G(0,-1;z)}{4} + \frac{G(0;z) \zeta\left(3\right)}{2} + \frac{\pi^{2} G(-1,0;z)}{4} - \frac{\pi^{2} G(0,0;z)}{8} - G(-1,-1,0,0;z) + \frac{G(-1,0,0,0;z)}{2} + \frac{G(0,-1,0,0;z)}{2} - G(-1;z) \zeta\left(3\right) - \frac{\pi^{4}}{30}.$$

Order3; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-037e74ca8da848638b3464ab8cefdd71/order-3/raw.txt).

$$h_{0,3}(z)=- \frac{G(-1,-1,0,0,0;z)}{2} - \frac{\pi^{2} G(0,-1,-1;z)}{4} + \frac{G(0,0,-1,0,0;z)}{4} + \frac{\pi^{2} G(-1,0,0;z)}{8} - \frac{G(0,-1,-1,0,0;z)}{2} - \frac{\pi^{4} G(0;z)}{60} - \frac{G(-1,0,-1,0,0;z)}{2} - \frac{\pi^{2} G(0,0,0;z)}{16} + \frac{G(0,-1,0,0,0;z)}{4} + \frac{G(0,0;z) \zeta\left(3\right)}{4} - \frac{\pi^{2} G(-1,-1,0;z)}{4} + \frac{\pi^{4} G(-1;z)}{30} + \frac{\pi^{2} G(-1,-1,-1;z)}{2} + G(-1,-1,-1,0,0;z) + \frac{\pi^{2} G(0,0,-1;z)}{8} - \frac{\pi^{2} G(-1,0,-1;z)}{4} - \frac{G(0,0,0,0,0;z)}{4} - \frac{G(0,-1;z) \zeta\left(3\right)}{2} + \frac{G(-1,0,0,0,0;z)}{2} + \frac{\pi^{2} G(0,-1,0;z)}{8} + G(-1,-1;z) \zeta\left(3\right) - \frac{G(-1,0;z) \zeta\left(3\right)}{2} + \zeta\left(5\right).$$

Order4; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-037e74ca8da848638b3464ab8cefdd71/order-4/raw.txt).

$$h_{0,4}(z)=\frac{G(0,0,0,-1,0,0;z)}{8} + \frac{G(0,-1,-1,-1,0,0;z)}{2} - \frac{G(0,-1,0,-1,0,0;z)}{4} + \frac{G(0,0,-1,0,0,0;z)}{8} - G(-1,-1,-1,-1,0,0;z) - \frac{G(-1,0,-1,0,0,0;z)}{4} + \frac{G(-1,0,0,0,0,0;z)}{4} - \frac{G(-1,0,0,-1,0,0;z)}{4} - \frac{\pi^{2} G(-1,0,-1,0;z)}{8} - \frac{\pi^{4} G(-1,-1;z)}{30} + \frac{\pi^{2} G(0,-1,0,0;z)}{16} - G(-1;z) \zeta\left(5\right) - \frac{G(-1,-1,0,0,0,0;z)}{2} + \frac{\pi^{2} G(-1,-1,-1,0;z)}{4} + \frac{G(-1,-1,0;z) \zeta\left(3\right)}{2} - \frac{\pi^{2} G(0,0,-1,-1;z)}{8} - \frac{\pi^{2} G(-1,-1,0,0;z)}{8} - \frac{\pi^{2} G(-1,0,0,-1;z)}{8} - \frac{G(0,0,-1;z) \zeta\left(3\right)}{4} + \frac{G(-1,0,-1;z) \zeta\left(3\right)}{2} + \frac{\pi^{2} G(-1,-1,0,-1;z)}{4} - \frac{G(0,-1,0;z) \zeta\left(3\right)}{4} + \frac{\pi^{2} G(-1,0,-1,-1;z)}{4} + \frac{G(-1,-1,-1,0,0,0;z)}{2} + \frac{G(0,-1,-1;z) \zeta\left(3\right)}{2} - \frac{\pi^{4} G(0,0;z)}{120} + \frac{\pi^{2} G(0,0,-1,0;z)}{16} + \frac{\pi^{2} G(0,-1,-1,-1;z)}{4} + \frac{\pi^{2} G(0,0,0,-1;z)}{16} + \frac{\pi^{4} G(-1,0;z)}{60} - \frac{\pi^{2} G(0,-1,-1,0;z)}{8} + \frac{\pi^{2} G(-1,0,0,0;z)}{16} - \frac{G(-1,0,0;z) \zeta\left(3\right)}{4} + \frac{G(0,0,0;z) \zeta\left(3\right)}{8} - \frac{G(0,0,-1,-1,0,0;z)}{4} - G(-1,-1,-1;z) \zeta\left(3\right) - \frac{\pi^{2} G(-1,-1,-1,-1;z)}{2} + \frac{G(0;z) \zeta\left(5\right)}{2} - \frac{\pi^{2} G(0,-1,0,-1;z)}{8} - \frac{\pi^{2} G(0,0,0,0;z)}{32} + \frac{\pi^{4} G(0,-1;z)}{60} + \frac{G(-1,-1,0,-1,0,0;z)}{2} + \frac{G(0,-1,0,0,0,0;z)}{4} + \frac{G(-1,0,-1,-1,0,0;z)}{2} - \frac{3 G(0,0,0,0,0,0;z)}{16} - \frac{G(0,-1,-1,0,0,0;z)}{4} - \frac{\pi^{6}}{315}.$$

## H1: parameters (0, 0, 1, 0)

Order0; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-cb202a7b22a4409db3fe08c778f6ee64/order-0/raw.txt).

$$h_{1,0}(z)=- G(0,0;z) - \frac{\pi^{2}}{2}.$$

Order1; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-cb202a7b22a4409db3fe08c778f6ee64/order-1/raw.txt).

$$h_{1,1}(z)=\pi^{2} G(-1;z) - \frac{G(0,0,0;z)}{2} + 2 G(-1,0,0;z) - \frac{5 \pi^{2} G(0;z)}{12} + 2 \zeta\left(3\right).$$

Order2; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-cb202a7b22a4409db3fe08c778f6ee64/order-2/raw.txt).

$$h_{1,2}(z)=\frac{\pi^{2} G(0,-1;z)}{2} + \frac{5 \pi^{2} G(-1,0;z)}{6} - \frac{\pi^{2} G(0,0;z)}{8} - 4 G(-1;z) \zeta\left(3\right) + G(-1,0,0,0;z) - 2 \pi^{2} G(-1,-1;z) + 3 G(0;z) \zeta\left(3\right) - 4 G(-1,-1,0,0;z) - \frac{G(0,0,0,0;z)}{2} + G(0,-1,0,0;z) - \frac{31 \pi^{4}}{360}.$$

Order3; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-8a0ff292fb25440bbcbf462a5c927cbb/order-3/raw.txt).

$$h_{1,3}(z)=\frac{5 \pi^{2} G(0,-1,0;z)}{12} + \frac{G(0,0;z) \zeta\left(3\right)}{2} - \frac{5 \pi^{2} G(0,0,0;z)}{48} - \pi^{2} G(-1,0,-1;z) - 2 G(-1,-1,0,0,0;z) + 4 \pi^{2} G(-1,-1,-1;z) - \frac{5 \pi^{2} G(-1,-1,0;z)}{3} + 8 G(-1,-1,-1,0,0;z) - 2 G(0,-1;z) \zeta\left(3\right) + \frac{31 \pi^{4} G(-1;z)}{180} - \pi^{2} G(0,-1,-1;z) - 2 G(0,-1,-1,0,0;z) + \frac{G(0,-1,0,0,0;z)}{2} - \frac{49 \pi^{4} G(0;z)}{720} - \frac{G(0,0,0,0,0;z)}{4} + \frac{G(0,0,-1,0,0;z)}{2} - 2 G(-1,0,-1,0,0;z) + G(-1,0,0,0,0;z) + 8 G(-1,-1;z) \zeta\left(3\right) - 6 G(-1,0;z) \zeta\left(3\right) + \frac{\pi^{2} G(0,0,-1;z)}{4} + \frac{\pi^{2} G(-1,0,0;z)}{4} + \frac{\pi^{2} \zeta\left(3\right)}{3} + 10 \zeta\left(5\right).$$

Order4; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-8a0ff292fb25440bbcbf462a5c927cbb/order-4/raw.txt).

$$h_{1,4}(z)=\frac{5 \pi^{2} G(0,0,-1,0;z)}{24} + 12 G(-1,-1,0;z) \zeta\left(3\right) + \frac{\pi^{2} G(0,-1,0,0;z)}{8} - \frac{\pi^{2} G(-1,0,0,-1;z)}{2} + 4 G(0,-1,-1;z) \zeta\left(3\right) + \frac{49 \pi^{4} G(-1,0;z)}{360} + \frac{10 \pi^{2} G(-1,-1,-1,0;z)}{3} + \frac{5 \pi^{2} G(-1,0,0,0;z)}{24} - 16 G(-1,-1,-1;z) \zeta\left(3\right) + \frac{\pi^{2} G(0,0,0,-1;z)}{8} + 2 \pi^{2} G(-1,-1,0,-1;z) - \frac{\pi^{2} G(0,-1,0,-1;z)}{2} - 8 \pi^{2} G(-1,-1,-1,-1;z) - \frac{\pi^{2} G(0;z) \zeta\left(3\right)}{6} + 11 G(0;z) \zeta\left(5\right) - \frac{\pi^{2} G(0,0,0,0;z)}{32} - 3 G(0,-1,0;z) \zeta\left(3\right) + \frac{31 \pi^{4} G(0,-1;z)}{360} - G(0,0,-1;z) \zeta\left(3\right) - \frac{5 \pi^{2} G(-1,0,-1,0;z)}{6} - \frac{31 \pi^{4} G(0,0;z)}{1440} - \frac{\pi^{2} G(-1,-1,0,0;z)}{2} - G(-1,0,0;z) \zeta\left(3\right) - \frac{5 \pi^{2} G(0,-1,-1,0;z)}{6} + 2 \pi^{2} G(0,-1,-1,-1;z) + 4 G(-1,-1,0,-1,0,0;z) - G(0,-1,0,-1,0,0;z) - \frac{3 G(0,0,0,0,0,0;z)}{16} - 16 G(-1,-1,-1,-1,0,0;z) + 4 G(0,-1,-1,-1,0,0;z) + \frac{G(0,-1,0,0,0,0;z)}{2} - G(0,-1,-1,0,0,0;z) + 4 G(-1,0,-1,-1,0,0;z) + \frac{G(0,0,0,-1,0,0;z)}{4} - G(-1,0,0,-1,0,0;z) + \frac{3 G(0,0,0;z) \zeta\left(3\right)}{4} - G(-1,0,-1,0,0,0;z) - 2 G(-1,-1,0,0,0,0;z) + \frac{G(0,0,-1,0,0,0;z)}{4} + \frac{G(-1,0,0,0,0,0;z)}{2} + 4 G(-1,-1,-1,0,0,0;z) - G(0,0,-1,-1,0,0;z) - \frac{\pi^{2} G(0,0,-1,-1;z)}{2} - \frac{31 \pi^{4} G(-1,-1;z)}{90} + 2 \pi^{2} G(-1,0,-1,-1;z) + 4 G(-1,0,-1;z) \zeta\left(3\right) - 20 G(-1;z) \zeta\left(5\right) - \frac{2 \pi^{2} G(-1;z) \zeta\left(3\right)}{3} - \frac{31 \pi^{6}}{1008} - 2 \zeta^{2}\left(3\right).$$

## H2: parameters (0, 0, 2, 0)

Order0; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-6c2197cf1ffe43e2acbe37fa539159c1/order-0/raw.txt).

$$h_{2,0}(z)=- G(0,0;z) - \frac{\pi^{2}}{2}.$$

Order1; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-6c2197cf1ffe43e2acbe37fa539159c1/order-1/raw.txt).

$$h_{2,1}(z)=- \frac{7 \pi^{2} G(0;z)}{12} - \frac{G(0,0,0;z)}{2} + 3 G(-1,0,0;z) + \frac{3 \pi^{2} G(-1;z)}{2} + 3 \zeta\left(3\right).$$

Order2; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-6c2197cf1ffe43e2acbe37fa539159c1/order-2/raw.txt).

$$h_{2,2}(z)=\frac{3 G(0,-1,0,0;z)}{2} + \frac{3 G(-1,0,0,0;z)}{2} + \frac{3 \pi^{2} G(0,-1;z)}{4} - \frac{9 \pi^{2} G(-1,-1;z)}{2} - 9 G(-1,-1,0,0;z) - \frac{G(0,0,0,0;z)}{2} - 9 G(-1;z) \zeta\left(3\right) + \frac{15 G(0;z) \zeta\left(3\right)}{2} - \frac{\pi^{2} G(0,0;z)}{8} + \frac{7 \pi^{2} G(-1,0;z)}{4} - \frac{\pi^{4}}{5}.$$

Order3; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-6c2197cf1ffe43e2acbe37fa539159c1/order-3/raw.txt).

$$h_{2,3}(z)=\frac{3 G(0,-1,0,0,0;z)}{4} + \frac{3 \pi^{2} G(-1,0,0;z)}{8} - \frac{2 \pi^{4} G(0;z)}{9} - \frac{21 \pi^{2} G(-1,-1,0;z)}{4} + \frac{3 \pi^{4} G(-1;z)}{5} + \frac{3 G(0,0;z) \zeta\left(3\right)}{4} - \frac{9 G(-1,0,-1,0,0;z)}{2} + \frac{3 G(-1,0,0,0,0;z)}{2} + \frac{7 \pi^{2} G(0,-1,0;z)}{8} - \frac{45 G(-1,0;z) \zeta\left(3\right)}{2} - \frac{9 G(-1,-1,0,0,0;z)}{2} - \frac{9 \pi^{2} G(0,-1,-1;z)}{4} - \frac{7 \pi^{2} G(0,0,0;z)}{48} - \frac{9 \pi^{2} G(-1,0,-1;z)}{4} + \frac{27 \pi^{2} G(-1,-1,-1;z)}{2} + 27 G(-1,-1,-1,0,0;z) - \frac{G(0,0,0,0,0;z)}{4} - \frac{9 G(0,-1,-1,0,0;z)}{2} + 27 G(-1,-1;z) \zeta\left(3\right) + \frac{3 \pi^{2} G(0,0,-1;z)}{8} + \frac{3 G(0,0,-1,0,0;z)}{4} - \frac{9 G(0,-1;z) \zeta\left(3\right)}{2} + 2 \pi^{2} \zeta\left(3\right) + 39 \zeta\left(5\right).$$

Order4; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-6c2197cf1ffe43e2acbe37fa539159c1/order-4/raw.txt).

$$h_{2,4}(z)=- \frac{9 \pi^{2} G(0,-1,0,-1;z)}{8} - \frac{21 \pi^{2} G(-1,0,-1,0;z)}{8} - 81 G(-1,-1,-1;z) \zeta\left(3\right) - \frac{\pi^{2} G(0,0,0,0;z)}{32} + \frac{27 \pi^{2} G(-1,0,-1,-1;z)}{4} + \frac{3 \pi^{4} G(0,-1;z)}{10} + \frac{27 G(0,-1,-1;z) \zeta\left(3\right)}{2} - \frac{9 \pi^{2} G(0,0,-1,-1;z)}{8} - \frac{9 \pi^{2} G(-1,0,0,-1;z)}{8} + \frac{7 \pi^{2} G(-1,0,0,0;z)}{16} - \frac{\pi^{4} G(0,0;z)}{20} - \frac{9 \pi^{2} G(-1,-1,0,0;z)}{8} + \frac{27 G(-1,0,-1;z) \zeta\left(3\right)}{2} - \frac{9 G(0,0,-1;z) \zeta\left(3\right)}{4} - \pi^{2} G(0;z) \zeta\left(3\right) + \frac{123 G(0;z) \zeta\left(5\right)}{2} + \frac{2 \pi^{4} G(-1,0;z)}{3} - \frac{9 G(-1,0,0;z) \zeta\left(3\right)}{4} - \frac{9 \pi^{4} G(-1,-1;z)}{5} + \frac{27 \pi^{2} G(0,-1,-1,-1;z)}{4} + \frac{7 \pi^{2} G(0,0,-1,0;z)}{16} + \frac{3 \pi^{2} G(0,-1,0,0;z)}{16} - \frac{45 G(0,-1,0;z) \zeta\left(3\right)}{4} - \frac{81 \pi^{2} G(-1,-1,-1,-1;z)}{2} + \frac{15 G(0,0,0;z) \zeta\left(3\right)}{8} - \frac{3 G(0,0,0,0,0,0;z)}{16} - \frac{9 G(-1,0,-1,0,0,0;z)}{4} - 81 G(-1,-1,-1,-1,0,0;z) - \frac{9 G(-1,-1,0,0,0,0;z)}{2} + \frac{3 G(0,-1,0,0,0,0;z)}{4} + \frac{27 G(0,-1,-1,-1,0,0;z)}{2} - \frac{9 G(0,-1,0,-1,0,0;z)}{4} + \frac{27 G(-1,-1,-1,0,0,0;z)}{2} - \frac{9 G(0,-1,-1,0,0,0;z)}{4} - \frac{9 G(0,0,-1,-1,0,0;z)}{4} - \frac{21 \pi^{2} G(0,-1,-1,0;z)}{8} + \frac{3 G(0,0,0,-1,0,0;z)}{8} - \frac{9 G(-1,0,0,-1,0,0;z)}{4} + \frac{3 G(0,0,-1,0,0,0;z)}{8} + \frac{3 G(-1,0,0,0,0,0;z)}{4} + \frac{27 G(-1,-1,0,-1,0,0;z)}{2} + \frac{27 G(-1,0,-1,-1,0,0;z)}{2} - 117 G(-1;z) \zeta\left(5\right) - 6 \pi^{2} G(-1;z) \zeta\left(3\right) + \frac{27 \pi^{2} G(-1,-1,0,-1;z)}{4} + \frac{63 \pi^{2} G(-1,-1,-1,0;z)}{4} + \frac{3 \pi^{2} G(0,0,0,-1;z)}{16} + \frac{135 G(-1,-1,0;z) \zeta\left(3\right)}{2} - \frac{41 \pi^{6}}{210} - 18 \zeta^{2}\left(3\right).$$

## HW: parameters (0, 0, 1, 1)

Order0; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-da088b2b40c84b2ea248f4486ed1593d/order-0/raw.txt).

$$h_{W,0}(z)=- G(0,0;z) - \frac{\pi^{2}}{2}.$$

Order1; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-da088b2b40c84b2ea248f4486ed1593d/order-1/raw.txt).

$$h_{W,1}(z)=\frac{3 \pi^{2} G(-1;z)}{2} + 3 G(-1,0,0;z) - \frac{3 \pi^{2} G(0;z)}{4} - \frac{3 G(0,0,0;z)}{2} + 3 \zeta\left(3\right).$$

Order2; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-da088b2b40c84b2ea248f4486ed1593d/order-2/raw.txt).

$$h_{W,2}(z)=- \frac{5 G(0,0,0,0;z)}{2} + \frac{9 \pi^{2} G(-1,0;z)}{4} - 9 G(-1;z) \zeta\left(3\right) - \frac{9 \pi^{2} G(-1,-1;z)}{2} + \frac{9 G(0;z) \zeta\left(3\right)}{2} - 9 G(-1,-1,0,0;z) - \frac{19 \pi^{2} G(0,0;z)}{24} + \frac{9 G(0,-1,0,0;z)}{2} + \frac{9 G(-1,0,0,0;z)}{2} + \frac{9 \pi^{2} G(0,-1;z)}{4} - \frac{13 \pi^{4}}{60}.$$

Order3; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-da088b2b40c84b2ea248f4486ed1593d/order-3/raw.txt).

$$h_{W,3}(z)=- \frac{27 G(0,-1;z) \zeta\left(3\right)}{2} - \frac{27 G(-1,-1,0,0,0;z)}{2} + \frac{27 \pi^{2} G(0,0,-1;z)}{8} + 27 G(-1,-1,-1,0,0;z) + \frac{3 G(0,0;z) \zeta\left(3\right)}{4} - \frac{15 G(0,0,0,0,0;z)}{4} - \frac{27 \pi^{2} G(-1,-1,0;z)}{4} - \frac{27 \pi^{2} G(0,-1,-1;z)}{4} + \frac{15 G(-1,0,0,0,0;z)}{2} - \frac{27 \pi^{2} G(-1,0,-1;z)}{4} + 27 G(-1,-1;z) \zeta\left(3\right) - \frac{27 G(-1,0,-1,0,0;z)}{2} + \frac{27 \pi^{2} G(0,-1,0;z)}{8} + \frac{27 \pi^{2} G(-1,-1,-1;z)}{2} + \frac{19 \pi^{2} G(-1,0,0;z)}{8} - \frac{27 G(-1,0;z) \zeta\left(3\right)}{2} + \frac{27 G(0,-1,0,0,0;z)}{4} + \frac{27 G(0,0,-1,0,0;z)}{4} + \frac{13 \pi^{4} G(-1;z)}{20} - \frac{13 \pi^{4} G(0;z)}{40} - \frac{19 \pi^{2} G(0,0,0;z)}{16} - \frac{27 G(0,-1,-1,0,0;z)}{2} - \pi^{2} \zeta\left(3\right) + 21 \zeta\left(5\right).$$

Order4; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-da088b2b40c84b2ea248f4486ed1593d/order-4/raw.txt).

$$h_{W,4}(z)=- \frac{81 \pi^{2} G(-1,0,0,-1;z)}{8} + \frac{81 \pi^{2} G(-1,-1,-1,0;z)}{4} - \frac{81 \pi^{2} G(0,-1,0,-1;z)}{8} - \frac{163 \pi^{2} G(0,0,0,0;z)}{96} + \frac{81 \pi^{2} G(-1,0,-1,-1;z)}{4} - 63 G(-1;z) \zeta\left(5\right) + 3 \pi^{2} G(-1;z) \zeta\left(3\right) - \frac{81 G(0,-1,0;z) \zeta\left(3\right)}{4} - \frac{3 \pi^{2} G(0;z) \zeta\left(3\right)}{2} + \frac{63 G(0;z) \zeta\left(5\right)}{2} + \frac{81 \pi^{2} G(0,0,-1,0;z)}{16} - \frac{9 G(-1,0,0;z) \zeta\left(3\right)}{4} - \frac{81 \pi^{2} G(0,-1,-1,0;z)}{8} + \frac{39 \pi^{4} G(-1,0;z)}{40} - \frac{81 \pi^{2} G(-1,-1,-1,-1;z)}{2} - \frac{81 \pi^{2} G(0,0,-1,-1;z)}{8} - \frac{81 \pi^{2} G(-1,0,-1,0;z)}{8} - \frac{39 \pi^{4} G(-1,-1;z)}{20} - \frac{57 \pi^{2} G(-1,-1,0,0;z)}{8} + \frac{81 G(-1,-1,0;z) \zeta\left(3\right)}{2} - 81 G(-1,-1,-1;z) \zeta\left(3\right) + \frac{39 \pi^{4} G(0,-1;z)}{40} + \frac{9 G(0,0,0;z) \zeta\left(3\right)}{8} + \frac{81 \pi^{2} G(0,-1,-1,-1;z)}{4} + \frac{57 \pi^{2} G(-1,0,0,0;z)}{16} + \frac{81 G(0,-1,-1;z) \zeta\left(3\right)}{2} + \frac{45 G(0,-1,0,0,0,0;z)}{4} - 81 G(-1,-1,-1,-1,0,0;z) - \frac{81 G(0,0,-1,-1,0,0;z)}{4} - \frac{81 G(0,-1,-1,0,0,0;z)}{4} + \frac{45 G(-1,0,0,0,0,0;z)}{4} - \frac{81 G(0,-1,0,-1,0,0;z)}{4} - \frac{81 G(-1,0,0,-1,0,0;z)}{4} - \frac{81 G(-1,0,-1,0,0,0;z)}{4} + \frac{81 G(-1,-1,0,-1,0,0;z)}{2} + \frac{81 G(-1,0,-1,-1,0,0;z)}{2} - \frac{81 G(0,0,-1;z) \zeta\left(3\right)}{4} + \frac{81 G(0,0,0,-1,0,0;z)}{8} + \frac{81 G(0,-1,-1,-1,0,0;z)}{2} - \frac{45 G(-1,-1,0,0,0,0;z)}{2} - \frac{91 G(0,0,0,0,0,0;z)}{16} + \frac{81 G(-1,-1,-1,0,0,0;z)}{2} + \frac{81 G(0,0,-1,0,0,0;z)}{8} + \frac{57 \pi^{2} G(0,-1,0,0;z)}{16} - \frac{263 \pi^{4} G(0,0;z)}{720} + \frac{81 \pi^{2} G(0,0,0,-1;z)}{16} + \frac{81 G(-1,0,-1;z) \zeta\left(3\right)}{2} + \frac{81 \pi^{2} G(-1,-1,0,-1;z)}{4} - \frac{41 \pi^{6}}{420} + 9 \zeta^{2}\left(3\right).$$

## HD: parameters (1, 0, 1, 0)

Order0; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-a59fe28fcef34593922acada910e65cd/order-0/raw.txt).

$$h_{D,0}(z)=- G(0,0;z) - \frac{\pi^{2}}{2}.$$

Order1; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-a59fe28fcef34593922acada910e65cd/order-1/raw.txt).

$$h_{D,1}(z)=3 G(-1,0,0;z) + \frac{3 \pi^{2} G(-1;z)}{2} - \frac{7 \pi^{2} G(0;z)}{12} - \frac{G(0,0,0;z)}{2} + 3 \zeta\left(3\right).$$

Order2; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-a59fe28fcef34593922acada910e65cd/order-2/raw.txt).

$$h_{D,2}(z)=\frac{3 G(-1,0,0,0;z)}{2} - \frac{\pi^{2} G(0,0;z)}{8} - 9 G(-1;z) \zeta\left(3\right) - 9 G(-1,-1,0,0;z) - \frac{G(0,0,0,0;z)}{2} + \frac{5 \pi^{2} G(0,-1;z)}{4} + \frac{7 \pi^{2} G(-1,0;z)}{4} - \frac{9 \pi^{2} G(-1,-1;z)}{2} + \frac{11 G(0;z) \zeta\left(3\right)}{2} + \frac{5 G(0,-1,0,0;z)}{2} - \frac{\pi^{4}}{6}.$$

Order3; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-a59fe28fcef34593922acada910e65cd/order-3/raw.txt).

$$h_{D,3}(z)=- \frac{33 G(-1,0;z) \zeta\left(3\right)}{2} + \frac{35 \pi^{2} G(0,-1,0;z)}{24} + \frac{3 \pi^{2} G(-1,0,0;z)}{8} - \frac{15 G(-1,0,-1,0,0;z)}{2} + \frac{5 G(0,-1,0,0,0;z)}{4} - \frac{15 G(0,-1;z) \zeta\left(3\right)}{2} - \frac{7 \pi^{2} G(0,0,0;z)}{48} + \frac{3 G(0,0,-1,0,0;z)}{4} + \frac{27 \pi^{2} G(-1,-1,-1;z)}{2} - \frac{15 \pi^{2} G(0,-1,-1;z)}{4} - \frac{G(0,0,0,0,0;z)}{4} - \frac{15 G(0,-1,-1,0,0;z)}{2} - \frac{15 \pi^{2} G(-1,0,-1;z)}{4} + \frac{\pi^{4} G(-1;z)}{2} + 27 G(-1,-1;z) \zeta\left(3\right) - \frac{19 \pi^{4} G(0;z)}{180} + 27 G(-1,-1,-1,0,0;z) - \frac{9 G(-1,-1,0,0,0;z)}{2} + \frac{3 G(-1,0,0,0,0;z)}{2} + \frac{3 G(0,0;z) \zeta\left(3\right)}{4} + \frac{3 \pi^{2} G(0,0,-1;z)}{8} - \frac{21 \pi^{2} G(-1,-1,0;z)}{4} + \frac{4 \pi^{2} \zeta\left(3\right)}{3} + 19 \zeta\left(5\right).$$

Order4; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-a59fe28fcef34593922acada910e65cd/order-4/raw.txt).

$$h_{D,4}(z)=\frac{63 \pi^{2} G(-1,-1,-1,0;z)}{4} - \frac{55 G(0,-1,0;z) \zeta\left(3\right)}{4} - \frac{35 \pi^{2} G(-1,0,-1,0;z)}{8} + \frac{45 \pi^{2} G(0,-1,-1,-1;z)}{4} + \frac{5 \pi^{2} G(0,-1,0,0;z)}{16} + \frac{5 \pi^{2} G(0,0,0,-1;z)}{16} + \frac{45 \pi^{2} G(-1,-1,0,-1;z)}{4} - \frac{3 \pi^{4} G(-1,-1;z)}{2} + \frac{45 \pi^{2} G(-1,0,-1,-1;z)}{4} - \frac{9 \pi^{2} G(-1,0,0,-1;z)}{8} - \frac{2 \pi^{2} G(0;z) \zeta\left(3\right)}{3} + \frac{43 G(0;z) \zeta\left(5\right)}{2} - \frac{25 \pi^{2} G(0,-1,0,-1;z)}{8} - \frac{35 \pi^{2} G(0,-1,-1,0;z)}{8} - \frac{9 G(0,0,-1;z) \zeta\left(3\right)}{4} - \frac{9 \pi^{2} G(-1,-1,0,0;z)}{8} + \frac{45 G(-1,0,-1;z) \zeta\left(3\right)}{2} - \frac{9 \pi^{2} G(0,0,-1,-1;z)}{8} + \frac{7 \pi^{2} G(0,0,-1,0;z)}{16} - \frac{\pi^{2} G(0,0,0,0;z)}{32} + \frac{19 \pi^{4} G(-1,0;z)}{60} - \frac{9 G(-1,0,0;z) \zeta\left(3\right)}{4} - \frac{81 \pi^{2} G(-1,-1,-1,-1;z)}{2} + \frac{45 G(0,-1,-1;z) \zeta\left(3\right)}{2} + \frac{7 \pi^{2} G(-1,0,0,0;z)}{16} + \frac{45 G(-1,0,-1,-1,0,0;z)}{2} - \frac{15 G(-1,0,-1,0,0,0;z)}{4} + \frac{3 G(0,0,-1,0,0,0;z)}{8} + \frac{27 G(-1,-1,-1,0,0,0;z)}{2} - \frac{9 G(-1,-1,0,0,0,0;z)}{2} - \frac{3 G(0,0,0,0,0,0;z)}{16} + \frac{5 G(0,-1,0,0,0,0;z)}{4} + \frac{45 G(0,-1,-1,-1,0,0;z)}{2} - 81 G(-1,-1,-1,-1,0,0;z) - \frac{15 G(0,-1,-1,0,0,0;z)}{4} + \frac{11 G(0,0,0;z) \zeta\left(3\right)}{8} - \frac{9 G(-1,0,0,-1,0,0;z)}{4} + \frac{5 G(0,0,0,-1,0,0;z)}{8} + \frac{3 G(-1,0,0,0,0,0;z)}{4} - \frac{25 G(0,-1,0,-1,0,0;z)}{4} + \frac{45 G(-1,-1,0,-1,0,0;z)}{2} - \frac{9 G(0,0,-1,-1,0,0;z)}{4} + \frac{99 G(-1,-1,0;z) \zeta\left(3\right)}{2} - 81 G(-1,-1,-1;z) \zeta\left(3\right) + \frac{5 \pi^{4} G(0,-1;z)}{12} - 57 G(-1;z) \zeta\left(5\right) - 4 \pi^{2} G(-1;z) \zeta\left(3\right) - \frac{\pi^{4} G(0,0;z)}{24} - \frac{\pi^{6}}{15} - 8 \zeta^{2}\left(3\right).$$

## HX: parameters (0, 1, 1, 0)

Order0; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-1bfaa1b10aca494cb620ea93305ed781/order-0/raw.txt).

$$h_{X,0}(z)=- G(0,0;z) - \frac{\pi^{2}}{2}.$$

Order1; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-1bfaa1b10aca494cb620ea93305ed781/order-1/raw.txt).

$$h_{X,1}(z)=3 G(-1,0,0;z) - \frac{3 G(0,0,0;z)}{2} + \frac{3 \pi^{2} G(-1;z)}{2} - \frac{3 \pi^{2} G(0;z)}{4} + 3 \zeta\left(3\right).$$

Order2; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-1bfaa1b10aca494cb620ea93305ed781/order-2/raw.txt).

$$h_{X,2}(z)=\frac{9 G(0,-1,0,0;z)}{2} - \frac{9 \pi^{2} G(-1,-1;z)}{2} + \frac{9 G(-1,0,0,0;z)}{2} - 9 G(-1,-1,0,0;z) + \frac{9 \pi^{2} G(-1,0;z)}{4} + \frac{9 \pi^{2} G(0,-1;z)}{4} - \frac{23 \pi^{2} G(0,0;z)}{24} - \frac{5 G(0,0,0,0;z)}{2} - 9 G(-1;z) \zeta\left(3\right) + \frac{9 G(0;z) \zeta\left(3\right)}{2} - \frac{3 \pi^{4}}{10}.$$

Order3; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-1bfaa1b10aca494cb620ea93305ed781/order-3/raw.txt).

$$h_{X,3}(z)=\frac{23 \pi^{2} G(-1,0,0;z)}{8} - \frac{27 G(-1,0,-1,0,0;z)}{2} + \frac{27 G(0,-1,0,0,0;z)}{4} + \frac{19 G(0,0;z) \zeta\left(3\right)}{4} + 27 G(-1,-1;z) \zeta\left(3\right) + \frac{27 G(0,0,-1,0,0;z)}{4} - \frac{27 \pi^{2} G(0,-1,-1;z)}{4} - \frac{27 G(0,-1;z) \zeta\left(3\right)}{2} + \frac{27 \pi^{2} G(0,0,-1;z)}{8} - \frac{9 \pi^{4} G(0;z)}{20} - \frac{27 G(-1,0;z) \zeta\left(3\right)}{2} - \frac{27 G(-1,-1,0,0,0;z)}{2} + \frac{9 \pi^{4} G(-1;z)}{10} - \frac{23 \pi^{2} G(0,0,0;z)}{16} - \frac{15 G(0,0,0,0,0;z)}{4} + 27 G(-1,-1,-1,0,0;z) + \frac{15 G(-1,0,0,0,0;z)}{2} - \frac{27 \pi^{2} G(-1,0,-1;z)}{4} - \frac{27 \pi^{2} G(-1,-1,0;z)}{4} + \frac{27 \pi^{2} G(-1,-1,-1;z)}{2} - \frac{27 G(0,-1,-1,0,0;z)}{2} + \frac{27 \pi^{2} G(0,-1,0;z)}{8} + \frac{3 \pi^{2} \zeta\left(3\right)}{2} + 21 \zeta\left(5\right).$$

Order4; [native evaluated raw coefficient](../runtime/primitives/objects/cas-runs/single-opposite-1bfaa1b10aca494cb620ea93305ed781/order-4/raw.txt).

$$h_{X,4}(z)=- \frac{57 G(-1,0,0;z) \zeta\left(3\right)}{4} - \frac{69 \pi^{2} G(-1,-1,0,0;z)}{8} + \frac{81 \pi^{2} G(0,0,0,-1;z)}{16} - \frac{13 \pi^{4} G(0,0;z)}{20} - 63 G(-1;z) \zeta\left(5\right) - \frac{9 \pi^{2} G(-1;z) \zeta\left(3\right)}{2} - \frac{81 \pi^{2} G(0,-1,0,-1;z)}{8} - \frac{81 \pi^{2} G(-1,0,-1,0;z)}{8} - \frac{203 \pi^{2} G(0,0,0,0;z)}{96} + \frac{81 \pi^{2} G(-1,-1,0,-1;z)}{4} + \frac{81 G(-1,-1,0;z) \zeta\left(3\right)}{2} - \frac{81 \pi^{2} G(-1,-1,-1,-1;z)}{2} - \frac{27 \pi^{4} G(-1,-1;z)}{10} + \frac{81 G(-1,0,-1;z) \zeta\left(3\right)}{2} + \frac{81 \pi^{2} G(-1,0,-1,-1;z)}{4} - \frac{81 \pi^{2} G(-1,0,0,-1;z)}{8} + \frac{69 \pi^{2} G(-1,0,0,0;z)}{16} + \frac{81 G(0,-1,-1;z) \zeta\left(3\right)}{2} + \frac{81 \pi^{2} G(-1,-1,-1,0;z)}{4} + \frac{27 \pi^{4} G(-1,0;z)}{20} + \frac{81 \pi^{2} G(0,-1,-1,-1;z)}{4} + \frac{27 \pi^{4} G(0,-1;z)}{20} + \frac{9 \pi^{2} G(0;z) \zeta\left(3\right)}{4} + \frac{63 G(0;z) \zeta\left(5\right)}{2} - \frac{81 G(0,-1,0;z) \zeta\left(3\right)}{4} - 81 G(-1,-1,-1;z) \zeta\left(3\right) - 81 G(-1,-1,-1,-1,0,0;z) - \frac{81 G(0,-1,0,-1,0,0;z)}{4} - \frac{81 G(-1,0,-1,0,0,0;z)}{4} - \frac{81 G(0,0,-1,-1,0,0;z)}{4} + \frac{81 G(-1,-1,-1,0,0,0;z)}{2} + \frac{81 G(-1,-1,0,-1,0,0;z)}{2} + \frac{45 G(-1,0,0,0,0,0;z)}{4} + \frac{45 G(0,-1,0,0,0,0;z)}{4} - \frac{45 G(-1,-1,0,0,0,0;z)}{2} - \frac{81 G(-1,0,0,-1,0,0;z)}{4} + \frac{57 G(0,0,0;z) \zeta\left(3\right)}{8} + \frac{81 G(0,0,-1,0,0,0;z)}{8} + \frac{81 G(0,0,0,-1,0,0;z)}{8} + \frac{81 G(0,-1,-1,-1,0,0;z)}{2} - \frac{91 G(0,0,0,0,0,0;z)}{16} + \frac{81 G(-1,0,-1,-1,0,0;z)}{2} - \frac{81 G(0,-1,-1,0,0,0;z)}{4} + \frac{69 \pi^{2} G(0,-1,0,0;z)}{16} - \frac{81 \pi^{2} G(0,0,-1,-1;z)}{8} - \frac{81 \pi^{2} G(0,-1,-1,0;z)}{8} + \frac{81 \pi^{2} G(0,0,-1,0;z)}{16} - \frac{81 G(0,0,-1;z) \zeta\left(3\right)}{4} - \frac{353 \pi^{6}}{1680} - 3 \zeta^{2}\left(3\right).$$

