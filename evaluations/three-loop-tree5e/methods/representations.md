# Exact back-substitution of the compact function views

The package regenerates the Laurent coefficients in explicit iterated-integral
coordinates before consulting the selected old answers.  A compact expression
in the project's \(M,N,\tau,E,\Omega_6\) names is useful only if it expands
back into that coefficient.  [The independent checker](../scripts/check_representations.py) performs this
separate check after candidate freezing.  It does not import the producer,
its shuffle library, SymPy, or the producer's HyperInt parser.

The common exact coordinate space has basis elements

\[
G(a_1,\ldots,a_n;z)\,
\zeta_2^{p_2}\zeta_3^{p_3}\zeta_5^{p_5}
\log(2)^{p_L}A_4^{p_4}A_6^{p_6},\qquad
a_i\in\{0,-1,-\tfrac12,-2\}.
\]

The six period exponents are stored in the displayed order.  Three-zeta
records are embedded by appending three zeros.  Original T1/T2 rows with
`Hlog_word`, `pi_power`, `zeta3_power` and `zeta5_power` are embedded with
\(\pi^{2m}=6^m\zeta_2^m\); their saved rational coefficient already includes
the factorial converting \(\log^n z\) to \(n!G(0^n;z)\).
All coefficients use Python's exact `Fraction` arithmetic.

The checker reads coefficient expressions through a restricted Python AST.
It admits rational arithmetic, nonnegative integer powers, explicit Hlog
words at endpoint \(z\), \(\log z\), \(\log2\), \(\zeta_n\), and named
generators with explicit admitted dictionaries.  Even periods are reduced
by \(\zeta_4=2\zeta_2^2/5\) and
\(\zeta_6=8\zeta_2^3/35\).  Unsupported letters, unbound names, and division
by a function are rejected rather than guessed or numerically fitted.

Each named generator is replaced by its saved Hlog vector.  Products use a
separate recursion counting all order-preserving interleavings, including
repeated letters.  For example
\(G(0;z)^2=2G(0,0;z)\); failing to retain this multiplicity would change
finite coefficients.  A dictionary with duplicate coordinate rows is
summed exactly before comparison.  Equality is exact equality in the
declared word/period coordinates, without PSLQ or numerical tolerances.

The M dictionary is independently read from the 15 raw `M_*_*.txt` outputs
in the snapshotted single-opposite run
`workflow/checks/single-opposite/runs/20260916T090538777022Z`.  If a selected
reference also supplies any of these M generators in `basis_definitions`,
the two dictionaries must agree exactly.  Other explicit reference
`basis_definitions`, including N, tau, E and Omega6, are admitted reused
evaluated dictionaries.  Their reuse is listed for each orbit in the report.
This check does not claim to reintegrate the defining functions or prove
their independence.  In particular, an Omega6 dictionary can carry the wider
letters and the fixed periods \(A_4,A_6\) without being silently replaced by
a \(\{0,-1\}\) answer.

For every stored Laurent position the named expression is expanded and
compared with the frozen regenerated coefficient.  When the reference also
stores `hlog_coefficients`, those coordinates are compared separately with
the same candidate.  Thus a matching named expression is not accepted merely
because it matches its own copied coordinate table.  The comparison report
retains both outcomes, the generators used, raw M source hashes and the
number of checked positions.

The standalone self-check exercises repeated-letter and mixed-letter
shuffles, even-zeta conversion, rational wide letters, a changed rational
coefficient, and rejection of an unknown generator, unexpected letter and
function denominator.  These controls check the arithmetic/parser boundary;
they are not a proof of the original integrals.  The post-freeze run checked
all 44 reference records: all 308 named expressions and all 308 separately
stored Hlog coordinate lists agree exactly with the frozen regenerated
candidates.  All eight arithmetic/parser controls passed.  These observed
comparisons do not follow merely from the self-checks.  The mathematical
conclusions remain computational and the formal status remains `not_run`.
