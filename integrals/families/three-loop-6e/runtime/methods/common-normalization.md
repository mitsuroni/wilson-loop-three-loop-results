# From an ordered J to the physical I

Each three-loop input supplies its own numerator, three propagator polynomials, six parameter bounds and an explicit outer-to-inner integration order. These are accepted input facts. The evaluator first checks those objects and derives the endpoint structure; a diagram label alone never selects an answer.

Write epsilon for the dimensional regulator and let z=z1>0. Throughout these articles logarithms of positive arguments are real. The source may contain an adjacent kernel q(1-u)v, with q in {1,z}, or one of the opposite kernels

\[
D_z(t,s)=(1-t)(1-s)+zts,\qquad
\widetilde D_z(t,s)=z(1-t)(1-s)+ts=zD_{1/z}(t,s).
\]

The positive interior permits the adjacent factorization

\[
q\,[q(1-u)v]^{\epsilon-1}
=q^\epsilon(1-u)^{\epsilon-1}v^{\epsilon-1}.
\]

It does not permit arbitrary changes of the root's integration order. For the present exchange graphs every endpoint belongs to exactly one propagator. At fixed z>0, an opposite kernel vanishes only at two corners, where it is comparable to the sum of the two corner distances. An adjacent endpoint has an integrable power for Re epsilon>0. Enlarging the ordered domain to the full endpoint cube therefore provides an absolute-integrability bound for small positive Re epsilon. This establishes the specific Fubini regrouping used below. The original order remains part of each input record and of the independent root checks.

## Colour, symmetry and the unique Gamma factor

Let n=3 and let a1 be the first coordinate after reconstructing the full ordered colour basis under CA^2=1+2 CA CF. This is a fixed-basis projection. Let sym be the orbit multiplicity and put

\[
\mathcal G(\epsilon)=e^{-\gamma_E\epsilon}\Gamma(1-\epsilon),\qquad
\kappa=\frac{\mathrm{sym}\,a_1}{2(-2)^n}.
\]

The accepted prepared C already contains exactly one factor exp(-n gamma_E epsilon), with the validated power n of Gamma(1-epsilon). The mathematical expression \(\mathcal G^n\) below denotes that factor as consumed from C; it is not an instruction to add a second Euler compensation. The common scale and coupling factors are extracted according to the current C/I convention. The resulting definition is

\[
I(\epsilon,z)=\kappa\mathcal G(\epsilon)^n
\frac{J(\epsilon,z)+z^{n\epsilon}J(\epsilon,1/z)}
{1+z^{n\epsilon}}.
\tag{1}
\]

The two J terms are evaluated before being combined. A source J need not be crossing invariant. The crossing and normalization in (1) apply to the physical I, which is why some simplifications are legal only after assembly. No additional W versus log W conversion or RG subtraction is performed inside this per-integral packet.

Equation (1) is first connected to the original integral for z>0 and small positive Re epsilon, then continued meromorphically to epsilon=0. At fixed positive z its denominator is nonzero in a sufficiently small disk around zero. A complex physical i0 continuation, including a route to z=-1, requires a separately specified branch connection and is not established by these positive-real-z articles.

## Expansion and function conventions

The target is

\[
I(\epsilon,z)=\sum_{p=-6}^{0}c_p(z)\epsilon^p+O(\epsilon).
\]

A coefficient below a proved lower pole bound is zero, rather than missing. A cancellation within the computed jet can improve the actual pole order further. The six methods state their required Taylor order before multiplying explicit poles; taking only the apparent finite part of a subintegral can otherwise omit contributions to c0.

All complete HLOG results use

\[
G(a,w;z)=\int_0^z\frac{G(w;t)}{t-a}\,dt,
\quad G(\varnothing;z)=1,
\quad G(\underbrace{0,\ldots,0}_{m};z)=\frac{\log^m z}{m!}.
\]

The letters are 0 and -1. The all-zero words use the displayed tangential-basepoint convention; the positive endpoint is ordinary. Products are expanded by shuffle when comparing word coordinates. The empty word represents constants, including products of zeta values. Thus a Gamma-only answer is also a complete HLOG answer, with empty words and, when present, zero-letter logarithmic words.

For a compact representation define

\[
r=\frac{z}{1+z},\quad W(\beta)=\log[\beta(1-\beta)],\quad
R(\beta)=\log\frac{1-\beta}{\beta},\quad \Lambda(\beta)=\rho+R(\beta),\quad \rho=\log z,
\]

\[
M_{k,j}=\int_0^1\frac{W^k\Lambda^j}{\beta-r}\,d\beta,\qquad
N_{k,j}=\int_0^1\frac{W^kR\Lambda^j}{\beta-r}\,d\beta,\qquad j\ge1.
\]

At beta=r, Lambda=0, so the quotient is read by continuity; it does not require a principal value. Useful exact relations are

\[
N_{k,j}=M_{k,j+1}-\log z\,M_{k,j},\qquad
M_{0,1}=-\frac{\log^2z+\pi^2}{2},\qquad M_{0,3}=-M_{0,1}^2.
\]

The tau/eta functions are coefficient functions of explicitly defined centered kernels in [the one-opposite method](single-opposite.md); E is defined in [the three-opposite method](same-channel-three.md). Their complete saved HLOG values are what makes them usable as a result representation. Introducing a name for an unevaluated integral would not establish a new solved coefficient. A compact-family conversion is checked in complete HLOG/period coordinates and is a different operation from evaluating the original integral.
