"""Export checked compact matter coefficients without changing scientific inputs.

The Af compact bare series is a reconstruction from the selected compact finite
result, the selected poles and independently rebuilt coupling counterterms. Its
moment-to-Hlog back-substitution is checked against every selected raw/bare/post
coefficient. The other two raw sums are assembled from their current graph roots.
The existing RG replay supplies the separately checked 31-graph Af connection;
this script does not repeat its HyperInt change-of-argument calculation.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys
import time

sys.dont_write_bytecode = True
import sympy as s

HERE = Path(__file__).resolve().parent
POWERS = tuple(range(-6, 1))
CHECKER = "workflow/checks/matter-rg-20261017/check.py"
MOMENT_DIR = "workflow/checks/single-opposite/runs/20260916T090538777022Z"
COLORS = ("CF*CA*TF*nf", "CF^2*TF*nf", "CF*(TF*nf)^2")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("D:/wilson-agent"))
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument("--replay-dir", type=Path, default=HERE / "rg-replay")
    args = parser.parse_args()
    root, output, replay_dir = (x.resolve() for x in
                                (args.root, args.output_dir, args.replay_dir))
    started = time.time()
    hashes: dict[str, str] = {}
    replay_hashes: dict[str, str] = {}
    checks: list[dict] = []

    def digest(data: bytes) -> str:
        return hashlib.sha256(data).hexdigest()

    def text(rel: str) -> str:
        raw = (root / rel).read_bytes()
        h = digest(raw)
        if rel in hashes and hashes[rel] != h:
            raise AssertionError(("concurrent source modification", rel))
        hashes[rel] = h
        return raw.decode("utf-8-sig")

    def read(rel: str) -> dict:
        return json.loads(text(rel))

    def replay_read(name: str) -> dict:
        raw = (replay_dir / name).read_bytes()
        replay_hashes[name] = digest(raw)
        return json.loads(raw.decode("utf-8-sig"))

    text(CHECKER)
    spec = importlib.util.spec_from_file_location("matter_rg_helpers", root / CHECKER)
    assert spec is not None and spec.loader is not None
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    eps, L, rho, z, CA, CF, T = (getattr(helper, name) for name in
                               ("eps", "L", "rho", "z", "CA", "CF", "T"))
    parse, norm = helper.parse, helper.norm
    M01, M03 = s.symbols("M_0_1 M_0_3")

    def compact(expr):
        value = s.expand(s.sympify(expr).subs(M03, -M01**2))
        assert not value.atoms(s.Float), "inexact coefficient"
        return value

    def ck(name, residual, color=None, method="exact displayed-expression algebra"):
        value = norm(residual)
        row = dict(name=name, channel=color, method=method,
                   passed=value == 0, residual=str(value))
        checks.append(row)
        if value != 0:
            raise AssertionError((name, str(value)[:1800]))

    def nonzero(name, residual):
        value = norm(residual)
        checks.append(dict(name=name, negative_control=True,
                           passed=value != 0, residual=str(value)))
        if value == 0:
            raise AssertionError(name)

    def series(values, lo=-6, hi=0):
        actual = {int(q): parse(v) for q, v in values.items()}
        helper.require_budget(actual, lo, hi, "source series")
        return {q: compact(actual[q]) for q in range(lo, hi + 1)}

    def full_scale(values, n, lo=-6, hi=0):
        return {q: compact(sum(values[j] * helper.scale(n, q-j)
                               for j in range(lo, q+1)))
                for q in range(lo, hi+1)}

    def strip_scale(values, n, lo=-6, hi=0):
        unit = {}
        for q in range(lo, hi+1):
            value = compact((values[q] - sum(unit[j] * helper.scale(n, q-j)
                                             for j in range(lo, q))) / 2)
            # Preserve moment notation in the constant term. Terms depending
            # on the exterior scale must vanish using registered identities.
            for (power,), coefficient in s.Poly(value, L).terms():
                if power:
                    ck(f"strip S{n}: scale independence at eps{q}, L^{power}",
                       coefficient)
            unit[q] = compact(value.subs(L, 0))
        return unit

    def encode(values):
        return {str(q): str(compact(v)) for q, v in values.items()}

    replay_status = replay_read("checks.json")
    if replay_status.get("status") != "passed":
        raise AssertionError("current RG replay has not passed")
    replay_inputs = replay_read("inputs.json")
    replay_results = replay_read("current-coefficients.json")
    # Bind the reused current-graph bridge, including its dirty-workspace inputs.
    for rel, expected in replay_inputs["source_sha256"].items():
        text(rel)
        if hashes[rel] != expected:
            raise AssertionError(("replay source drift", rel))

    current = read("benchmarks/renormalization/three-loop/current.json")
    selected, selection = {}, {}
    for color in COLORS:
        pointer = current["channels"][color]
        rel = "benchmarks/renormalization/three-loop/" + pointer["result"]
        data = read(rel)
        if hashes[rel] != pointer["sha256"]:
            raise AssertionError(("selected result hash mismatch", color))
        selected[color] = data[pointer["field"]] if "field" in pointer else data
        selection[color] = dict(path=rel, sha256=hashes[rel],
                                field="/" + pointer.get("field", ""))

    moments, moment_records = {}, {}
    for name in ("M_0_1", "M_1_1", "M_2_1", "M_3_1", "M_1_2", "M_1_3"):
        rel = f"{MOMENT_DIR}/{name}.txt"
        expression = text(rel).strip().split("=", 1)[1].rstrip(";").strip()
        moments[s.Symbol(name)] = parse(expression)
        moment_records[name] = dict(source=rel, sha256=hashes[rel],
                                    maple_hlog=expression,
                                    sympy_hlog=str(moments[s.Symbol(name)]))
    moments[s.Symbol("N_1_2")] = (moments[s.Symbol("M_1_3")]
                                   - rho * moments[s.Symbol("M_1_2")])
    # Include the actual integral definition and generation convention in the
    # portable source closure, not just a textual pointer in the dictionary.
    text("benchmarks/2ne/families/single-opposite/three-loop.md")
    text(f"{MOMENT_DIR}/compute_moments.mpl")

    def hlog(expr):
        return norm(compact(expr).subs(moments, simultaneous=True))

    def bridge(name, left, right, color):
        ck(name, hlog(left) - hlog(right), color,
           "registered M-to-Hlog substitutions and exact shuffle normalization")

    one = read("benchmarks/renormalization/one-loop/result-r001.json")
    lower = read("benchmarks/renormalization/three-loop/lower-order-inputs-r001.json")
    two = read("benchmarks/renormalization/two-loop/result-r002.json")
    stored_ct = read("benchmarks/renormalization/three-loop/counterterms-r001.json")
    w1 = {q: compact(v.subs({s.Symbol("L_s"): L+rho/2,
                             s.Symbol("L_t"): L-rho/2}))
          for q, v in series(one["laurent_unit_coefficients"], -2, 2).items()}
    w2 = series(lower["hat_w2_over_CF"], -4, 1)
    w2_unit = series(lower["unscaled_group_over_CF"], -4, 1)
    w2_restored = full_scale(w2_unit, 2, -4, 1)
    for q in range(-4, 2):
        ck(f"lower-order stored unscaled sum times S2 eps{q}", w2_restored[q]-w2[q])
    b0 = (11*CA-4*T)/12
    b1 = (34*CA**2-20*CA*T-12*CF*T)/48
    ck("beta0 normalization", b0-parse(stored_ct["b0"]))
    ck("beta1 normalization", b1-parse(stored_ct["b1"]))
    ct_all = {q: compact(-2*b0*w2.get(q+1, 0)
                         + b0**2*w1.get(q+2, 0)-b1*w1.get(q+1, 0)/2)
              for q in POWERS}
    extract = {
        "CF*CA*TF*nf": lambda x: s.expand(x).coeff(CA).coeff(T),
        "CF^2*TF*nf": lambda x: s.expand(x).coeff(CF).coeff(T),
        "CF*(TF*nf)^2": lambda x: s.expand(x).coeff(T, 2),
    }
    ct = {color: {q: compact(extract[color](ct_all[q])) for q in POWERS}
          for color in COLORS}
    for color in COLORS:
        for q in POWERS:
            bridge(f"independently rebuilt CT vs stored eps{q}", ct[color][q],
                   parse(stored_ct["color_channels"][color][str(q)]), color)

    native = read("benchmarks/self-energy/current-r002.json")
    native_map = {entry["orbit_id"]: entry for entry in native["entries"]}
    channels = {}
    for color, key, short, ids, source_label in (
        ("CF^2*TF*nf", "cff", "Ff", ["2e-Pi2-01", "2e-Pi2-02"], "Ff"),
        ("CF*(TF*nf)^2", "nf2", "ff", ["2e-Pi1+Pi1-01", "2e-Pi1+Pi1-02"], "CF*(TF*nf)^2"),
    ):
        adopted = selected[color]
        unit = {q: s.S.Zero for q in POWERS}
        source_rows = []
        for oid in ids:
            pointer = native_map[oid]["result"]
            rel = pointer["path"]
            data = read(rel)
            if hashes[rel] != pointer["sha256"]:
                raise AssertionError(("native source hash mismatch", oid))
            if hashes[rel] != adopted["source_sha256"][rel]:
                raise AssertionError(("adopted source mismatch", oid))
            entry = next(x for x in data["color_channels"]
                         if x.get("monomial", x.get("channel")) == source_label)
            values = entry["coefficients"] if key == "nf2" else entry["laurent"]["coefficients"]
            helper.require_budget({int(q): v for q, v in values.items()}, -6, 0, oid)
            for q in POWERS:
                value = values[str(q)]
                moment = compact(parse(value["expression"]).subs(L, rho))
                if "hlog_terms" in value:
                    bridge(f"{oid} current Hlog/moment bridge eps{q}", moment,
                           parse(value), color)
                unit[q] += moment
            source_rows.append(dict(orbit_id=oid, path=rel, sha256=hashes[rel],
                                    field=source_label,
                                    field_sha256=helper.field_hash(values)))
        unit = {q: compact(v) for q, v in unit.items()}
        bare = full_scale(unit, 3)
        post = {q: compact(bare[q]+ct[color][q]) for q in POWERS}
        for q in POWERS:
            bridge(f"current raw-source sum vs selected bare eps{q}", bare[q],
                   parse(adopted[f"bare_w3_{key}"][str(q)]), color)
            bridge(f"current raw-source sum and CT vs selected post eps{q}", post[q],
                   parse(adopted[f"post_counterterm_{key}"][str(q)]), color)
        F, f = parse(adopted[f"F3_{key}"]), parse(adopted[f"f3_{key}"])
        bridge("finite coefficient vs selected F", post[0], F, color)
        bridge("finite shape vs selected f", F.subs(L, 0), f, color)
        channels[color] = dict(short_label=short, U_compact=unit,
                               bare_compact=bare, CT_compact=ct[color],
                               post_compact=post, F=compact(F), f=compact(f),
                               k=parse(adopted[f"k3_{key}"]),
                               h=parse(adopted[f"h3_{key}"]),
                               provenance=dict(method="direct current raw moment-series assembly",
                                               selection=selection[color], sources=source_rows))

    color = "CF*CA*TF*nf"
    adopted = selected[color]
    post = series(adopted["post_counterterm_coefficients"])
    post[0] = compact(parse(adopted["finite_MN"]))
    F, f = post[0], compact(parse(adopted["shape_MN"]))
    bare = {q: compact(post[q]-ct[color][q]) for q in POWERS}
    unit = strip_scale(bare, 3)
    for stage, values, field in (("raw", unit, "raw_unit_coefficients"),
                                  ("bare", bare, "bare_full_scale_coefficients"),
                                  ("post", post, "post_counterterm_coefficients")):
        for q in POWERS:
            bridge(f"Af reconstructed compact {stage} vs selected GPL eps{q}",
                   values[q], parse(adopted[field][str(q)]), color)
    bridge("Af finite MN vs finite GPL", F, parse(adopted["finite_GPL"]), color)
    bridge("Af shape MN vs shape GPL", f, parse(adopted["shape_GPL"]), color)
    bridge("Af F at L=0 vs selected shape MN", F.subs(L, 0), f, color)
    channels[color] = dict(short_label="Af", U_compact=unit, bare_compact=bare,
                           CT_compact=ct[color], post_compact=post,
                           F=F, f=f, k=parse(adopted["k3"]), h=parse(adopted["h3"]),
                           provenance=dict(
                               method="reconstruct compact post from selected poles and finite_MN; subtract independently rebuilt CT; divide S3 by exact coefficient recursion; verify every coefficient against original selected GPL",
                               reconstruction_is_not_a_new_independent_sum=True,
                               selection=selection[color],
                               sources=replay_results[color]["sources"],
                               current_graph_connection="reused passing rg-replay; original argument-change canonicalization in workflow/checks/w3-ca-nf/physical.mpl"))

    F1, F2 = parse(two["F1_over_CF"]), parse(two["F2_over_CF"])
    for color, record in channels.items():
        record["logZ"] = {q: record["post_compact"][q] if q < 0 else s.S.Zero
                          for q in POWERS}
        record["gamma4"] = 2*record["h"]
        for stage in ("U_compact", "bare_compact", "CT_compact", "post_compact", "logZ"):
            if set(record[stage]) != set(POWERS):
                raise AssertionError(("incomplete output series", color, stage))
            for q, value in record[stage].items():
                # From beta -> 1-beta: M_ij(1/z)=(-1)^(j+1) M_ij(z),
                # N_ij(1/z)=(-1)^j N_ij(z). L is fixed; rho changes sign.
                rules = {rho: -rho}
                for symbol in value.free_symbols:
                    match = re.fullmatch(r"([MN])_(\d+)_(\d+)", str(symbol))
                    if match:
                        parity = int(match[3]) + (1 if match[1] == "M" else 0)
                        rules[symbol] = (-1)**parity*symbol
                ck(f"compact crossing {stage} eps{q}", value-value.subs(rules, simultaneous=True),
                   color, "moment integral reflection beta -> 1-beta; rho -> -rho at fixed L")
        rescaled = full_scale(record["U_compact"], 3)
        for q in POWERS:
            ck(f"compact U times S3 equals bare eps{q}",
               rescaled[q]-record["bare_compact"][q], color)
            ck(f"compact bare plus CT equals post eps{q}",
               record["bare_compact"][q]+record["CT_compact"][q]-record["post_compact"][q], color)
            bridge(f"compact bare vs current replay eps{q}", record["bare_compact"][q],
                   parse(replay_results[color]["bare"][str(q)]), color)
            bridge(f"compact post vs current replay eps{q}", record["post_compact"][q],
                   parse(replay_results[color]["post"][str(q)]), color)
        ck("finite equals post epsilon zero", record["F"]-record["post_compact"][0], color)
        ck("shape equals finite at L=0", record["f"]-record["F"].subs(L, 0), color)
        ck("cusp extraction from single pole", -s.Rational(3, 2)*s.diff(record["logZ"][-1], L)-record["k"], color)
        ck("noncusp extraction from single pole", -3*record["logZ"][-1].subs(L, 0)-record["h"], color)
        lower_rg = compact(extract[color](2*b0*F2+b1*F1))
        ck("finite RG in compact representation", s.diff(record["F"], L)-lower_rg
           +2*record["k"]*L+record["h"], color)
        record["lower_RG"] = lower_rg
        for symbol in set().union(*(v.free_symbols for values in
                                    (record["U_compact"], record["CT_compact"], record["post_compact"])
                                    for v in values.values())):
            if symbol not in (L, rho) and not re.fullmatch(r"[MN]_\d+_\d+", str(symbol)):
                raise AssertionError(("unexpected symbol in color-stripped compact output", color, str(symbol)))
        record["moment_symbols"] = sorted({str(symbol) for stage in
            ("U_compact", "bare_compact", "CT_compact", "post_compact")
            for value in record[stage].values() for symbol in value.free_symbols
            if re.fullmatch(r"[MN]_\d+_\d+", str(symbol))})

    nonzero("negative: omit S3 constant factor two", channels[COLORS[0]]["U_compact"][-4])
    nonzero("negative: replace N12 by M13 without rho M12",
            moments[s.Symbol("N_1_2")]-moments[s.Symbol("M_1_3")])
    nonzero("negative: omit w2 epsilon-one from Af finite CT",
            extract[COLORS[0]](-2*b0*w2[1]))
    nonzero("negative: omit w1 epsilon-two from ff finite CT",
            extract[COLORS[2]](b0**2*w1[2]))

    dictionary = dict(
        expression_language="SymPy exact strings; ** powers; zeta(n); M_i_j and N_i_j are raw moments",
        domain="z>0; rho=log(z); L=(L_s+L_t)/2; L_s=L+rho/2 and L_t=L-rho/2",
        regulator="D=4-2*eps; coefficient [eps^q] carries no q!",
        moments=dict(
            r_z="z/(1+z)", W_beta="log(beta*(1-beta))",
            R_beta="log((1-beta)/beta)", L_beta="rho+R_beta",
            M_i_j="Integral(W_beta**i*L_beta**j/(beta-r_z), (beta,0,1))",
            N_i_j="Integral(W_beta**i*R_beta*L_beta**j/(beta-r_z), (beta,0,1))",
            index_domain="i>=0, j>=1; no factorial denominators; removable point beta=r_z",
            identities={"N_i_j": "M_i_(j+1)-rho*M_i_j", "M_0_1": "-(rho**2+pi**2)/2", "M_0_3": "-M_0_1**2"},
            reflection={"M_i_j(1/z)": "(-1)**(j+1)*M_i_j(z)", "N_i_j(1/z)": "(-1)**j*N_i_j(z)"},
            source="benchmarks/2ne/families/single-opposite/three-loop.md"),
        hlog=dict(definition="G(a1,...,an;z)=Integral(G(a2,...,an;t)/(t-a1),(t,0,z))",
                  word_order="first letter is outermost; empty word is 1",
                  zeros="G(0,...,0;z)=log(z)**n/n!", alphabet=[0, -1],
                  branch="real z>0; continuation outside this domain requires the source branch prescription",
                  symbol_encoding="H_z_m1_0_0 denotes Hlog(z,[-1,0,0]); all-zero words are replaced by rho**n/n!",
                  source="workflow/checks/matter-rg-20261017/check.py"),
        moment_to_hlog=moment_records,
        N_1_2_to_hlog=str(moments[s.Symbol("N_1_2")]),
        normalization=dict(T="TF*nf", B_f="-8*(1-eps)/(3-2*eps)",
                           S_n="exp(n*eps*(L+rho/2))+exp(n*eps*(L-rho/2))",
                           U_compact="scale-stripped bare exponent color coefficient",
                           bare_compact="S3*U_compact; lambda_b**3 and color stripped",
                           CT_compact="[-2*b0*w2/eps+(b0**2/eps**2-b1/(2*eps))*w1]_color",
                           post_compact="bare_compact+CT_compact",
                           logZ="negative-power coefficients of log(Z4), with epsilon-zero entry explicitly 0",
                           F="post_compact[0]", f="F.subs(L,0)",
                           k="-3/2*d(post_compact[-1])/dL", h="-3*post_compact[-1].subs(L,0)",
                           gamma4="2*h"),
        unused_in_these_totals=["tau_plus", "old weight-six E"],
        scope="Three matter color totals through eps^0 only. This is not a universal Hlog-to-M/N reduction of individual diagrams or an all-order closure claim.")

    lower_out = dict(
        normalization="All displayed lower-order quantities have the common outer CF removed; existing S1/S2 are retained in w1/w2.",
        b0=str(b0), b1=str(b1),
        w1_full_scale=encode(w1), w1_consumed_range=[-2, 2],
        w2_full_scale=encode(w2), w2_consumed_range=[-4, 1],
        w2_unscaled=encode(w2_unit),
        F1_over_CF=str(F1), F2_over_CF=str(F2),
        CT3_over_CF=encode(ct_all),
        sources={"w1": "benchmarks/renormalization/one-loop/result-r001.json#/laurent_unit_coefficients",
                 "w2": "benchmarks/renormalization/three-loop/lower-order-inputs-r001.json#/hat_w2_over_CF",
                 "F1_F2": "benchmarks/renormalization/two-loop/result-r002.json"})
    encoded_channels = {}
    for color in COLORS:
        encoded_channels[color] = {
            name: encode(value) if name in ("U_compact", "bare_compact", "CT_compact", "post_compact", "logZ")
            else str(compact(value)) if isinstance(value, s.Basic) else value
            for name, value in channels[color].items()}

    for rel, expected in hashes.items():
        if digest((root / rel).read_bytes()) != expected:
            raise AssertionError(("input changed during build", rel))
    for name, expected in replay_hashes.items():
        if digest((replay_dir / name).read_bytes()) != expected:
            raise AssertionError(("replay changed during build", name))
    summary = dict(status="passed", checks=len(checks),
                   negative_controls=sum(row.get("negative_control", False) for row in checks),
                   checks_detail=checks, elapsed_seconds=round(time.time()-started, 3),
                   sympy_version=s.__version__, integral_evaluation="not_run", formal="not_run",
                   scope="Exact algebra, selected-expression conversions and compact-stage reconstruction; not independent proofs of original integrals or external noncusp/finite-shape verification.")
    result = dict(
        format="wilson-three-matter-compact-coefficients/1", date="2026-10-23",
        status="passed", powers=list(POWERS), expression_language=dictionary["expression_language"],
        normalization="Coefficient of lambda**3 times the displayed color in log(W); no repeated exponentiation, C4 or Euler factor",
        channels=encoded_channels, lower_orders=lower_out, conversion_dictionary=dictionary,
        source_sha256=hashes,
        reused_rg_replay=dict(files_sha256=replay_hashes, checks=replay_status["checks"],
                              negative_controls=replay_status["negative_controls"],
                              role="current 31-graph Af source connection and existing RG checks, reused with all source hashes rechecked"),
        builder_sha256=digest(Path(__file__).read_bytes()),
        verification={k: v for k, v in summary.items() if k != "checks_detail"},
        limitations=["The pure CF*CA**2 channel is outside this export.",
                     "Af compact raw/bare series are reconstructed and checked, not a newly independent graph integration.",
                     "The external noncusp scale dictionary remains unresolved; internal gamma4=2*h is preserved.",
                     "A compact color total does not prove each individual diagram has an M/N representation."])
    output.mkdir(parents=True, exist_ok=True)
    for name, data in (("coefficients.json", result), ("coefficient-checks.json", summary),
                       ("conversion-dictionary.json", dictionary)):
        (output / name).write_text(json.dumps(data, ensure_ascii=True, indent=2)+"\n",
                                  encoding="utf-8", newline="\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "checks_detail"}), flush=True)


if __name__ == "__main__":
    main()
