"""Check the new complex-6e prose convolutions against full generated HLOGs.

This is an independent implementation of the explicit finite convolution,
not independent integration: it shares the archived raw CAS primitives.
It uses a minimal scalar-AST decoder and a one-letter G(0) insertion shuffle,
and calls neither the runtime assembly nor a CAS integration process.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re

import sympy as sp

ROOT = Path(__file__).resolve().parents[4]
RUNTIME = ROOT / "integrals/families/three-loop-6e/runtime"
PRIMARY = ROOT / "workflow/checks/integral-documentation-6e-20261103/evidence/primary"
ORBITS = (10, 11, 12, 20, 21, 22, 23, 24, 25, 32, 33, 34, 37, 38)
Z = sp.Symbol("z")


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scalar_ast(node):
    if node["kind"] == "rational":
        return sp.Rational(node["numerator"], node["denominator"])
    if node["kind"] == "var":
        if node["name"] != "z":
            raise ValueError("Unexpected scalar variable")
        return Z
    args = [scalar_ast(arg) for arg in node.get("args", ())]
    operation = node["op"]
    if operation == "add":
        return sp.Add(*args)
    if operation == "mul":
        return sp.Mul(*args)
    if operation == "rpow":
        return args[0] ** args[1]
    if operation == "pi":
        return sp.pi
    if operation == "zeta":
        return sp.zeta(*args)
    raise ValueError("Unsupported scalar operation: " + operation)


def native_coordinates(path):
    """Parse the finite, linear HLOG output without the runtime parser."""
    words = []

    def replace_hlog(match):
        word = tuple(int(x.strip()) for x in match.group(1).split(",") if x.strip())
        if word not in words:
            words.append(word)
        return "g" + str(words.index(word))

    text = re.sub(r"Hlog\(z,\[([^]]*)\]\)", replace_hlog,
                  path.read_text(encoding="utf-8").strip())
    text = re.sub(r"zeta\[(\d+)\]", r"Z\1", text).replace("^", "**")
    local = {"z": Z, **{f"Z{k}": sp.zeta(k) for k in (2, 3, 4, 5, 6)},
             **{f"g{k}": sp.Symbol(f"g{k}") for k in range(len(words))}}
    expression = sp.sympify(text, locals=local).expand()
    coordinates = {word: expression.coeff(local[f"g{k}"])
                   for k, word in enumerate(words)}
    coordinates[()] = expression.subs({local[f"g{k}"]: 0 for k in range(len(words))})
    if sp.expand(expression - coordinates[()] - sum(
            coordinates[word] * local[f"g{k}"] for k, word in enumerate(words))) != 0:
        raise ValueError("Native output is not a finite linear HLOG expression")
    return coordinates


def add(*coordinates):
    result = defaultdict(lambda: sp.Integer(0))
    for item in coordinates:
        for word, coefficient in item.items():
            result[word] += coefficient
    return result


def scale(coordinates, coefficient):
    return {word: coefficient * value for word, value in coordinates.items()}


def multiply_by_log_z(coordinates):
    """G(0;z) G(w;z): insert 0 at every position, retaining multiplicity."""
    result = defaultdict(lambda: sp.Integer(0))
    for word, coefficient in coordinates.items():
        for position in range(len(word) + 1):
            result[word[:position] + (0,) + word[position:]] += coefficient
    return result


def generated_coordinates(generated, orbit, order):
    result = defaultdict(lambda: sp.Integer(0))
    payload = generated[orbit]["coefficients"].get(str(order), {"terms": []})
    for term in payload["terms"]:
        if term["periods"]:
            raise ValueError("Unexpected separate period object")
        word = tuple(int(scalar_ast(letter)) for letter in term["gpl"]["word"])
        result[word] += scalar_ast(term["coefficient"])
    return result


def check():
    generated_path = PRIMARY / "generation/generated.json"
    replay_path = PRIMARY / "generation/primitive-replay.json"
    generated, replay = load(generated_path), load(replay_path)
    rows, input_paths = [], {generated_path, replay_path}
    for number in ORBITS:
        orbit = f"6e-{number:02d}"
        sub_path = RUNTIME / f"integrals/three-loop/{orbit}/subintegrals.json"
        sub = load(sub_path)
        input_paths.add(sub_path)
        kappa = sp.Rational(sub["normalization"]["kappa"])
        kernel = sub["active_primitive_consumption"]["ordered_kernel"]
        event = replay["replay_events"][replay["events_by_orbit"][orbit][0]]
        raw = {}
        for item in event["raw_files"]:
            tail = item["path"].replace("\\", "/").split("/cas-runs/")[-1]
            path = RUNTIME / "primitives/objects/cas-runs" / tail
            input_paths.add(path)
            raw[item["order"], item["channel"]] = native_coordinates(path)

        if "gap" in kernel["operation"]:
            r = sp.Rational(kernel["prefactor_powers"]["z_epsilon_power"])
            a1 = scale(raw[1, "raw"], (1 + Z) ** 2)
            b1 = scale(raw[1, "reflected"], (1 + 1 / Z) ** 2)
            a2 = scale(raw[2, "raw"], (1 + Z) ** 2)
            b2 = scale(raw[2, "reflected"], (1 + 1 / Z) ** 2)
            predictions = {
                -1: scale(add(a1, b1), kappa / 2),
                0: scale(add(a2, b2, scale(multiply_by_log_z(add(a1, scale(b1, -1))),
                                           r - sp.Rational(3, 2))), kappa / 2),
            }
        else:
            predictions = {0: scale(add(raw[0, "raw"], raw[0, "reflected"]), kappa / 2)}

        for order, predicted in predictions.items():
            target = generated_coordinates(generated, orbit, order)
            difference = add(predicted, scale(target, -1))
            nonzero = {word: sp.cancel(sp.expand(value))
                       for word, value in difference.items()
                       if sp.cancel(sp.expand(value)) != 0}
            if nonzero:
                raise AssertionError((orbit, order, nonzero))
            rows.append({"orbit": orbit, "order": order, "passed": True,
                         "word_coordinates_compared": len(difference),
                         "exact_nonzero_residuals": 0})
        print(orbit, "explicit prose convolution equals complete generated HLOG coordinates")

    return {
        "status": "passed", "orbit_count": len(ORBITS), "comparisons": len(rows),
        "checks": rows,
        "scope": "Independent implementation of the newly written finite physical convolution; shared archived CAS primitives and generated comparison target. Not independent analytic integration or formal proof.",
        "operations": ["minimal scalar AST decoding", "native finite HLOG parsing",
                       "single-letter G(0) insertion shuffle", "exact rational coefficient subtraction"],
        "new_CAS_processes": 0, "runtime_assembly_called": False,
        "formal_state": "not_run",
        "inputs": [{"path": str(path.relative_to(ROOT)), "sha256": digest(path)}
                   for path in sorted(input_paths)],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = check()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("status", "orbit_count", "comparisons")}))


if __name__ == "__main__":
    main()
