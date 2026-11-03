# Current analytic results

**[Read the four-color result article](four-color-results.md).** It displays the full three-loop Laurent coefficients and finite analytic functions in one place, with definitions of $M,N$, the color decomposition and the scale variable $L$.

| Object | Exact data |
|---|---|
| One-loop series, including positive orders | [1loop.json](../../benchmarks/current/results/1loop.json) |
| Two-loop series, including positive orders | [2loop.json](../../benchmarks/current/results/2loop.json) |
| Four-color three-loop series through epsilon0 | [3loop.json](../../benchmarks/current/results/3loop.json) |

In the three-loop file, `channels.AA`, `channels.Af`, `channels.Ff` and `channels.ff` correspond to $C_FC_A^2$, $C_FC_AT$, $C_F^2T$ and $C_FT^2$, where $T=T_Fn_f$. Each channel stores `U`, `bare`, `counterterm`, `post_counterterm`, `logZ`, `F`, `f_at_L0`, `k3` and `h3`. The `native_common_z_HLOG_channels` field supplies the HLOG construction.

`F` is the finite coefficient of the renormalized connected logarithm, with its full $L$ dependence. `f_at_L0` gives its shape function at $L=0$. The article explains the relation of these functions to the bare Laurent coefficients.

Use the [current guide](../../docs/current/README.md) to inspect individual inputs and the [reproduction guide](../../docs/reproduction.md) to reconstruct the result. The [validation index](../../validation/README.md) collects exact assembly, function and RG checks. Source consumption is recorded in the [source/joint ledger](../2026-11-02-current/source-ledger.json) and [matter ledger](../2026-11-02-current/matter-ledger.json).

Earlier dated deliveries are indexed under [historical releases](../../docs/archive/README.md).
