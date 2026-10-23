"""Audit three selected matter channels against the October 19 UV bridge.

Read science only from --root; write new uv-*.json files only to --output-dir.
The external k/h dictionaries are read as restricted arithmetic ASTs from the
preserved bridge builder; that builder is never executed. Finite comparisons
reuse registered moment/GPL identities and are not new integral evaluations.
"""
from pathlib import Path
import argparse
import ast
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import sys
import time

import sympy as s

HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = Path("D:/wilson-agent")
BRIDGE_DIR = "workflow/checks/n4-w3-normalization-20261019/uv-bridge/"
HELPER = "workflow/checks/matter-rg-20261017/check.py"
MOMENT_DIR = "workflow/checks/single-opposite/runs/20260916T090538777022Z/"
CHANNELS = {
    "Af": ("CF*CA*TF*nf", None),
    "Ff": ("CF^2*TF*nf", "cff"),
    "ff": ("CF*(TF*nf)^2", "nf2"),
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def arithmetic(node, constants):
    """Read only exact constants, +,-,*,/,**, and s.Rational from a Python AST."""
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return s.Integer(node.value)
    if isinstance(node, ast.Name) and node.id in constants:
        return constants[node.id]
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
        v = arithmetic(node.operand, constants)
        return -v if isinstance(node.op, ast.USub) else v
    if isinstance(node, ast.BinOp):
        a, b = arithmetic(node.left, constants), arithmetic(node.right, constants)
        operations = {ast.Add: lambda: a+b, ast.Sub: lambda: a-b,
                      ast.Mult: lambda: a*b, ast.Div: lambda: a/b,
                      ast.Pow: lambda: a**b}
        if type(node.op) in operations:
            return operations[type(node.op)]()
    if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "s" and node.func.attr == "Rational"
            and len(node.args) == 2 and not node.keywords):
        return s.Rational(*(arithmetic(v, constants) for v in node.args))
    raise ValueError("Unsupported external reference arithmetic: " + ast.dump(node))


def reference_dictionaries(source):
    constants = {"z"+str(n): s.zeta(n) for n in (2, 3, 4, 5)}
    result = {}
    for node in ast.parse(source).body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name) or target.id not in ("kref", "href"):
            continue
        if target.id in result or not isinstance(node.value, ast.Dict):
            raise ValueError("Ambiguous external reference dictionary")
        result[target.id] = {ast.literal_eval(k): arithmetic(v, constants)
                             for k, v in zip(node.value.keys, node.value.values)}
    if set(result) != {"kref", "href"}:
        raise ValueError("Missing independently supplied kref/href dictionaries")
    for values in result.values():
        if not set(CHANNELS).issubset(values):
            raise ValueError("Missing matter channel reference")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                        help="Read-only scientific source root (or package/sources)")
    parser.add_argument("--output-dir", type=Path, default=HERE,
                        help="Directory without existing uv-checks.json or uv-source-hashes.json")
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output_dir.resolve()
    targets = [output / name for name in ("uv-checks.json", "uv-source-hashes.json")]
    if any(p.exists() for p in targets):
        raise FileExistsError("Refusing to overwrite UV evidence; use a fresh --output-dir")
    started = time.time()
    source_hashes, checks, field_bindings = {}, [], []

    def read_bytes(rel):
        path = (root / rel).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Source path escapes --root: " + rel)
        data = path.read_bytes()
        digest = sha(data)
        if rel in source_hashes and source_hashes[rel] != digest:
            raise RuntimeError("Source changed during audit: " + rel)
        source_hashes[rel] = digest
        return data

    def read_json(rel):
        return json.loads(read_bytes(rel).decode("utf-8-sig"))

    def fact(name, passed, category="source_binding", **details):
        row = dict(name=name, category=category, passed=bool(passed), **details)
        checks.append(row)
        if not passed:
            raise AssertionError(row)

    read_bytes(HELPER)
    spec = importlib.util.spec_from_file_location("matter_uv_helper", root / HELPER)
    helper = importlib.util.module_from_spec(spec)
    # Do not write __pycache__ into the read-only input tree.
    previous_bytecode = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec.loader.exec_module(helper)
    finally:
        sys.dont_write_bytecode = previous_bytecode
    helper.ROOT = root
    L, rho, CA, CF, T = helper.L, helper.rho, helper.CA, helper.CF, helper.T

    def eq(name, expression, category, channel=None):
        value = helper.norm(expression)
        fact(name, value == 0, category, channel=channel, residual=str(value))

    bridge = read_json(BRIDGE_DIR + "uv-coordinate-bridge-r002.json")
    qcd = read_json(BRIDGE_DIR + "qcd-total-up-to-C3-r001.json")
    builder = read_bytes(BRIDGE_DIR + "build_bridge.py").decode("utf-8-sig")
    references = reference_dictionaries(builder)
    fact("qcd-total binds the selected UV bridge bytes",
         qcd["bridge_sha256"] == source_hashes[BRIDGE_DIR + "uv-coordinate-bridge-r002.json"])
    fact("preserved UV bridge reports all checks passed", bridge["all_passed"] is True)
    for rel, expected in bridge["source_bindings"].items():
        fact("UV bridge input bytes retained: " + rel, sha(read_bytes(rel)) == expected)
    # This prose is provenance, not an executable mathematical premise.
    read_bytes("benchmarks/renormalization/theory/n4-w3-normalization-20261019/derivation-r002.md")

    # Repeat the exact scale algebra without fitting any anomalous dimension.
    gE = s.EulerGamma
    raw_logs = {"project": s.log(s.pi), "closed": s.log(2*s.pi),
                "FGM_local": s.log(4*s.pi)}
    scale_logs = {"project": 2*gE-s.log(4), "closed": 2*gE-s.log(2),
                  "FGM_local": 2*gE}
    for name in raw_logs:
        eq(name + " raw bare factor gives its native coordinate scale",
           s.expand_log(raw_logs[name]+gE-s.log(4*s.pi)+gE-scale_logs[name], force=True),
           "uv_scale_algebra")
    eq("mu_closed_squared equals twice mu_x_squared",
       scale_logs["closed"]-scale_logs["project"]-s.log(2), "uv_scale_algebra")
    eq("mu_FGM_local_squared equals twice mu_closed_squared",
       scale_logs["FGM_local"]-scale_logs["closed"]-s.log(2), "uv_scale_algebra")
    q, mu2 = s.symbols("q mu_MSbar_squared", positive=True)
    eq("s=2q makes the closed and project logarithm arguments equal",
       s.exp(2*gE)*mu2*q/2-s.exp(2*gE)*mu2*(2*q)/4, "uv_scale_algebra")

    moments = {}
    for name in ("M_0_1", "M_1_1", "M_2_1", "M_1_3", "M_1_2"):
        text = read_bytes(MOMENT_DIR + name + ".txt").decode("utf-8-sig").strip()
        moments[s.Symbol(name)] = helper.parse(text.split("=", 1)[1].rstrip(";"))
    moments[s.Symbol("N_1_2")] = moments[s.Symbol("M_1_3")]-rho*moments[s.Symbol("M_1_2")]

    def gpl(value):
        return helper.norm(value.subs(moments, simultaneous=True))

    current = read_json("benchmarks/renormalization/three-loop/current.json")
    two = read_json("benchmarks/renormalization/two-loop/result-r002.json")
    b0 = (11*CA-4*T)/12
    b1 = (34*CA**2-20*CA*T-12*CF*T)/48
    F1, F2 = helper.parse(two["F1_over_CF"]), helper.parse(two["F2_over_CF"])
    k2, h2 = helper.parse(two["k2"])/CF, helper.parse(two["h2"])/CF
    extract = {
        "Af": lambda x: s.expand(x).coeff(CA).coeff(T),
        "Ff": lambda x: s.expand(x).coeff(CF).coeff(T),
        "ff": lambda x: s.expand(x).coeff(T, 2),
    }
    encoded = {}
    for short, (color, key) in CHANNELS.items():
        pointer = current["channels"][color]
        rel = "benchmarks/renormalization/three-loop/" + pointer["result"]
        whole = read_json(rel)
        fact(short + " current selector source hash", source_hashes[rel] == pointer["sha256"])
        adopted = whole[pointer["field"]] if "field" in pointer else whole
        if key is None:
            post_field, finite_field, shape_field = "post_counterterm_coefficients", "finite_MN", "shape_MN"
            kfield, hfield = "k3", "h3"
        else:
            post_field, finite_field, shape_field = "post_counterterm_"+key, "F3_"+key, "f3_"+key
            kfield, hfield = "k3_"+key, "h3_"+key
        base_pointer = "/"+pointer["field"] if "field" in pointer else ""
        for field in (post_field, finite_field, shape_field, kfield, hfield):
            field_bindings.append(dict(channel=color, path=rel, json_pointer=base_pointer+"/"+field,
                                       field_sha256=helper.field_hash(adopted[field])))
        post = helper.cs(adopted[post_field])
        fact(short + " complete post-counterterm epsilon table",
             set(post) == set(range(-6, 1)), "input_budget")
        finite, shape = helper.parse(adopted[finite_field]), helper.parse(adopted[shape_field])
        k, h = helper.parse(adopted[kfield]), helper.parse(adopted[hfield])
        kref, href = references["kref"][short], references["href"][short]
        refrow = qcd["channels"][short]
        refpost = helper.cs(refrow["post_counterterm"])
        fact(short + " UV comparison epsilon table complete",
             set(refpost) == set(range(-6, 1)), "input_budget")

        eq(short + " selected k3 versus external cusp reference", k-kref, "external_uv_reference", color)
        eq(short + " selected h3 versus external Gamma_Box/2 reference", h-href, "external_uv_reference", color)
        eq(short + " October 19 k3 versus reference dictionary", helper.parse(refrow["k3"])-kref,
           "external_uv_reference", color)
        eq(short + " October 19 h3 versus reference dictionary", helper.parse(refrow["h3"])-href,
           "external_uv_reference", color)

        # The external reference constants, not the selected pole extraction,
        # enter this independently predicted log Z4 pole table.
        B, P, Q, D = (extract[short](v) for v in (b0**2, b0*k2, b0*h2, b1))
        expected = {
            -6: s.S.Zero, -5: s.S.Zero,
            -4: -s.Rational(11, 9)*B,
            -3: -s.Rational(2, 3)*B*L+s.Rational(5, 9)*P+s.Rational(8, 9)*D,
            -2: s.Rational(2, 3)*(P+D)*L+Q/3-2*kref/9,
            -1: -(2*kref*L+href)/3,
        }
        for power, predicted in expected.items():
            eq(f"{short} independently predicted log Z4 epsilon^{power}",
               post[power]-predicted, "independent_pole_prediction", color)
            eq(f"{short} selected pole versus October 19 epsilon^{power}",
               post[power]-refpost[power], "later_source_consistency", color)
        eq(short + " selected finite equals post-counterterm epsilon^0",
           gpl(finite)-gpl(post[0]), "finite_representation_identity", color)
        eq(short + " selected shape equals F3 at L=0", gpl(shape)-gpl(finite.subs(L, 0)),
           "finite_representation_identity", color)
        eq(short + " full finite equals October 19 QCD total", gpl(finite)-gpl(refpost[0]),
           "later_source_consistency", color)
        eq(short + " finite shape equals October 19 QCD total at L=0",
           gpl(shape)-gpl(refpost[0].subs(L, 0)), "later_source_consistency", color)
        lower_rg = extract[short](2*b0*F2+b1*F1)
        eq(short + " finite RG with external k3 and h3",
           gpl(s.diff(finite, L)-lower_rg+2*kref*L+href), "external_uv_finite_rg", color)
        if key is None:
            eq("Af compact finite_MN versus finite_GPL", gpl(finite)-helper.parse(adopted["finite_GPL"]),
               "finite_representation_identity", color)
            eq("Af compact shape_MN versus shape_GPL", gpl(shape)-helper.parse(adopted["shape_GPL"]),
               "finite_representation_identity", color)
            for power in range(-6, 0):
                eq(f"Af stored log_Z4 matches selected poles epsilon^{power}",
                   helper.parse(adopted["log_Z4_coefficients"][str(power)])-post[power],
                   "finite_representation_identity", color)
        allowed_moments = {s.Symbol(n) for n in ("M_0_1", "M_1_1", "M_2_1", "M_1_3", "N_1_2")}
        for label, value in (("finite", finite), ("shape", shape)):
            unexpected = value.free_symbols - allowed_moments - {L, rho}
            fact(short+" "+label+" requires only M/N moments and scale logs", not unexpected,
                 "function_budget", channel=color, unexpected_symbols=sorted(map(str, unexpected)))
        for label, value in (("selected finite", finite), ("UV total finite", refpost[0])):
            fact(short+" "+label+" has no C3_N4 dependence",
                 not any("C3" in str(symbol) for symbol in value.free_symbols),
                 "matter_completeness", channel=color)
        # A wrong factor of two and an unfixed finite shape must not be mistaken
        # for successful external/reference validation.
        fact(short+" negative control: h3=Gamma_Box is rejected", helper.norm(2*h-href) != 0,
             "negative_control", channel=color)
        fact(short+" negative control: finite constant shift is detected by comparison",
             gpl(finite+1-refpost[0]) != 0, "negative_control", channel=color)
        # This is intentionally a passing blindness demonstration, not a claim
        # that the finite shape has been independently established by RG.
        eq(short+" RG is blind to an added L-independent even shape",
           s.diff(finite+rho**2+1, L)-s.diff(finite, L), "validation_limit_demonstration", color)
        encoded[color] = dict(
            short_name=short, k3_external=str(helper.norm(kref)), h3_external=str(helper.norm(href)),
            independently_predicted_log_Z4={str(p): str(helper.norm(v)) for p, v in expected.items()},
            selected_finite_MN=str(finite), selected_shape_MN=str(shape),
            required_moments=sorted(str(v) for v in finite.free_symbols & allowed_moments),
            C3_N4_dependency=False, finite_comparison_independent_original_integral=False,
        )

    for rel, digest in source_hashes.items():
        if sha((root/rel).read_bytes()) != digest:
            raise RuntimeError("Source changed before publication: " + rel)
    timestamp = datetime.now(timezone.utc).isoformat()
    report = dict(
        format="wilson-three-matter-uv-delivery-audit/1", status="passed", run_at_utc=timestamp,
        source_root=str(root), check_count=len(checks),
        negative_controls=sum(v["category"] == "negative_control" for v in checks),
        checks=checks, channels=encoded, sympy_version=s.__version__,
        elapsed_seconds=round(time.time()-started, 3),
        scope=dict(fresh_exact_arithmetic=True, fresh_external_reference_pole_prediction=True,
                   external_reference_source="Preserved October 19 bridge kref/href and finite-line scale dictionary",
                   fresh_literature_retrieval=False, fresh_three_loop_integration=False,
                   fresh_moment_integration=False, finite_comparison_is_shared_input_consistency=True,
                   independent_Af_finite_shape=False, independent_primitive_three_loop_UV=False,
                   formal="not_run", full_four_color_result=False, science_selectors_changed=False),
    )
    manifest = dict(format="wilson-three-matter-uv-inputs/1", run_at_utc=timestamp,
                    source_root=str(root), source_sha256=dict(sorted(source_hashes.items())),
                    field_bindings=field_bindings, script_sha256=sha(Path(__file__).read_bytes()),
                    historical_downloaded_source_sha256=bridge["downloaded_source_sha256"],
                    historical_downloaded_sources_reread=False,
                    uv_checks_sha256=sha((json.dumps(report, ensure_ascii=False, indent=2)+"\n").encode()))
    output.mkdir(parents=True, exist_ok=True)
    for path, value in zip(targets, (report, manifest)):
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(value, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(dict(status="passed", checks=len(checks),
                          output_dir=str(output), source_files=len(source_hashes))))


if __name__ == "__main__":
    main()
