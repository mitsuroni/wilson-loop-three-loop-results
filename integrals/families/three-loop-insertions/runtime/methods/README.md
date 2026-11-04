# Fifteen two-point and exchange-insertion integrals

P5 covers eleven Pi1-plus-exchange sources, two serial-Pi1 sources and two single-Pi2 sources. Each orbit retains its own labelled momentum/source object, ordered kernel, normalization and seven Laurent positions. Shared mathematical methods reduce repetition without dropping separate source weights.

Read [Pi1 with a bare exchange](pi1.md) for the six P1–P6 bodies and [serial Pi1 and single Pi2](chains-pi2.md) for the adjacent/opposite geometries and their distinct full external operators. The [independent-control matrix](independent-controls.md) states what is checked beyond reassembly.

The physical convention is wilson-qcd-action-path-fourier/1.0. Native C/J already contains the physical sign and exactly one Euler compensation exp(-3 gamma_E epsilon); the original full colour is retained. The connected scalar I uses the documented fixed-basis a1 projection and crossing/orbit extraction. That scalar, its external colour polynomial, and the original unprojected source are different objects. In particular a Pi2 geometry reused from a chain is not a replacement for the Pi2 operator.

Let rho=log(z), z>0, and Gcal=exp(-gamma_E epsilon) Gamma(1-epsilon). For a documented R factor, the scalar construction has the form

\[
 I(z,\epsilon)=\frac{m a_1}{2}R(\epsilon)
 \frac{J(z,\epsilon)+z^{3\epsilon}J(1/z,\epsilon)}{1+z^{3\epsilon}}.
\]

The family articles specify what is included in R and what remains external. A common equation alone does not identify the family normalization. Pi2 retains its geometric I and each coefficient of its full operator separately; the output contains their actual epsilon convolution before any sum of colour channels.

The target for every scalar or nonzero coloured contribution is the complete series from epsilon^-6 through epsilon^0 with O(epsilon), including established zeros. If I starts at epsilon^-p and a regular external multiplier is phi=sum(phi_j epsilon^j), its finite contribution is sum(phi_j c_-j,j=0,...,p). The whole required jet is retained, rather than substituting phi(0).

All HLOG coordinates use dt/(t-a), base point0 and positive endpointz, with letters0,-1 and zeta periods. Zero-letter words use the tangential log convention. M and N are the original unnormalised moments, N_kj=M_k,j+1-rho M_kj; no factorial redefinition is introduced. Other coefficient functions retain the exact definitions in their bound dictionary and complete HLOG expansions.

The original source rules, original evaluated moment/HLOG dictionaries and the full-D Pi2 literature operator are explicitly reused. The batch reexecutes root reductions, Gamma jets, physical assembly and operator convolutions, with exact comparisons performed only after candidate freeze. There is no new independent integration of all internal self-energy graphs, arbitrary complex i0 continuation or formal proof; formal=not_run.
