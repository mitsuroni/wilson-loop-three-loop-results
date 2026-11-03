"""Portable P6 reassembly of 31 T1/T2/T3 orbits from evaluated primitives.

Only snapshot formula code and primitive HyperInt outputs are read.  In
particular, generate() never reads a selected physical-I coefficient record,
an old coefficients.json, a representation.json, or a basis containing I.
The formula modules are imported, never executed as scripts.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from fractions import Fraction as Q
from pathlib import Path

import sympy as sp

PACK = Path(__file__).resolve().parents[1]
SNAP = PACK / "snapshots"
RUNS = {
    "t1-empty": "workflow/checks/tgv-t1-empty/runs/20260917T115009990802Z",
    "t1-opposite": "workflow/checks/tgv-t1-opposite/runs/20260917T131211354092Z",
    "t3-empty": "workflow/checks/tgv-t3-empty/runs/20260917T142141597785Z",
    "t3-opposite": "workflow/checks/tgv-t3-opposite/runs/20260918T005641322714Z",
}
FAMILIES = {
    "t1-empty": ("01", "02", "03", "04", "17", "18", "19"),
    "t2-empty": ("05", "06", "07", "08", "20", "21", "31"),
    "t3-empty": ("09", "10", "26", "30", "38", "39"),
    "t1-opposite": ("11", "12", "32"),
    "t2-opposite": ("16", "36", "37"),
    "t3-opposite": ("22", "23", "24", "25", "44"),
}
T2_MAP = {"05": "04", "06": "03", "07": "02", "08": "01",
          "20": "19", "21": "18", "31": "17", "16": "32",
          "36": "12", "37": "11"}
MOTHER_RUN = "workflow/checks/single-opposite/runs/20260916T090538777022Z"


def _read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _formula(family):
    module = _module("p6_t123_" + family.replace("-", "_"),
                     SNAP / "workflow/checks" / ("tgv-" + family) / "formulas.py")
    # Bind the original mathematical producer to a fixed saved primitive run.
    # This avoids ambient current-run selectors while leaving its maths intact.
    module.ROOT = SNAP
    module.run = lambda: SNAP / RUNS[family]
    return module


def _source_contract(ids):
    """Check the root/normalization bridge before any coefficient generation."""
    native = {r["orbit_id"]: r for r in _read(
        SNAP / "benchmarks/current/prepared/3loop/5e-3g.json")["orbits"]}
    historical = _read(SNAP / "benchmarks/tgv/5e-3g/prepared/5e-3g-r001.json")
    old = {r["orbit_id"]: r for r in historical["orbits"]}
    expected_i = ("I=sym/32*Exp[-3*EulerGamma*eps]*Gamma[1-eps]*Gamma[2-2*eps]*"
                  "(J(z,eps)+z^(3*eps)*J(1/z,eps))/(1+z^(3*eps))")
    assert historical["I_definition"] == expected_i
    required = {
        "color": "I/4*CF", "coupling": "(4*Pi*alphaS)^3",
        "ft_prefactor": "Gamma[1-eps]^1*Gamma[1-2*eps]/((4*Pi^(2-eps))^1*64*Pi^(4-2*eps))",
        "derivative": "2*(1-2*eps)", "v14_scale": "(-2*v14+I0)^(3*eps)",
        "i_factor": "I", "twoscale": "(-2)^(-1)",
        "euler": "exp(-3*EulerGamma*eps)",
    }
    for oid in ids:
        row, prior = native[oid], old[oid]
        assert row["J_raw_s"] == prior["J_raw_s"], oid
        assert row["C"]["components"] == required, oid
        assert row["C"]["component_semantics"] == "multiplicative", oid
        assert row["normalization"]["assembly_sign"] == 1, oid
        assert row["normalization"]["euler_compensation_count"] == 1, oid
        assert row["c4_symmetry"]["orbit_multiplicity"] == 4, oid
    # Full-color conversion: CF = CA^2*CF - 2*CA*CF^2 under CA^2=1+2*CA*CF.
    ca, cf, eps = sp.symbols("CA CF eps")
    assert sp.expand((ca**2 * cf - 2 * ca * cf**2 - cf).subs(
        ca**2, 1 + 2 * ca * cf)) == 0
    # Including the native phase, derivative, Fourier constants and coupling
    # gives +CF*(1-2eps)/16, before the common Gamma/scale/pi factors.
    raw_prefactor = sp.I * 2*(1-2*eps) * (sp.I*cf/4) * 4**3/(4*64) * (-sp.Rational(1, 2))
    assert sp.expand(raw_prefactor - cf*(1-2*eps)/16) == 0
    return native


def _check_logs():
    expected = {"t1-empty": "COMPLETE_T1", "t1-opposite": "T1_OPPOSITE_COMPLETE",
                "t3-empty": "T3_EMPTY_COMPLETE", "t3-opposite": "T3_OPPOSITE_COMPLETE"}
    for family, token in expected.items():
        assert token in (SNAP / RUNS[family] / "maple.log").read_text(encoding="utf-8"), family
    assert "CONTROL_COMPLETE" in (SNAP / RUNS["t1-empty"] / "control.log").read_text(encoding="utf-8")


def _primitive_paths(family, oid):
    if family.startswith("t2-"):
        return _primitive_paths(family.replace("t2-", "t1-"), T2_MAP[oid])
    run = RUNS[family]
    if family == "t1-empty":
        # Gamma-only 01/17 do not consume F/Z; 18 inherits the 17/19 sum.
        if oid in ("01", "17"):
            return []
        if oid in ("02", "03", "04"):
            return [f"{run}/F{q}.txt" for q in range(6)]
        return [f"{run}/Z{q}.txt" for q in range(1, 5)]
    if family == "t1-opposite":
        prefixes = ("ET", "EK") if oid == "11" else (("ET",) if oid == "32" else ("EK",))
        return [f"{run}/{p}{q}.txt" for p in prefixes for q in range(1, 5)]
    if family == "t3-empty":
        base = {"09": 9, "39": 9, "10": 10, "38": 10, "26": 26, "30": 26}[oid]
        return [f"{run}/E{base}q{q}.txt" for q in (2, 3, 4)]
    if oid == "44":
        return [f"{MOTHER_RUN}/M_{q}_1.txt" for q in range(3)]
    base = 22 if oid in ("22", "25") else 23
    return [f"{run}/E{base}q{q}.txt" for q in (1, 2)]


def _budget(family, oid):
    if family.startswith("t2-"):
        return _budget(family.replace("t2-", "t1-"), T2_MAP[oid])
    if family == "t1-empty":
        return {"gamma_and_beta_through": 6, "F_through": 5 if oid in ("02", "03", "04") else None,
                "Z_through": 4 if oid in ("18", "19") else None,
                "explicit_denominator_powers": [6, 5, 4], "remainder": "O(eps) at fixed z>0"}
    if family == "t1-opposite":
        return {"ET_EK_through": 4, "ET_starts": 1, "EK_starts": 2,
                "analytic_multiplier_through": 3, "explicit_denominator_power": 4,
                "remainder": "O(eps) at fixed z>0"}
    if family == "t3-empty":
        return {"physical_E_orders": [2, 3, 4], "analytic_multiplier_through": 2,
                "explicit_denominator_power": 4, "remainder": "O(eps) at fixed z>0"}
    return {"physical_E_orders": [1, 2] if oid != "44" else None,
            "known_M_indices": [0, 1, 2] if oid == "44" else None,
            "analytic_multiplier_through": 0, "explicit_denominator_power": 2,
            "remainder": "O(eps) at fixed z>0"}


def generate():
    """Return every epsilon^-6,...,epsilon^0 coefficient for the 31 orbits."""
    ids = ["5e-3g-" + oid for members in FAMILIES.values() for oid in members]
    _source_contract(ids)
    _check_logs()
    # Explicitly replace any ambient hlog_algebra module with the frozen copy.
    ha = _module("hlog_algebra", SNAP / "workflow/checks/single-opposite/hlog_algebra.py")
    vectors = {}
    t1 = _formula("t1-empty")
    for oid, coefficients in t1.coefficients().items():
        vectors[oid] = {k: ha.from_sympy(v, {t1.ell: ha.LOG}) for k, v in coefficients.items()}
    for family in ("t1-opposite", "t3-empty"):
        vectors.update(_formula(family).coefficients())
    for oid, base, factor in (("22", 22, 32), ("25", 22, 32), ("23", 23, 64), ("24", 23, 64)):
        vectors[oid] = {k: {} for k in range(-6, 1)}
        vectors[oid][-1] = ha.scale(ha.parse(SNAP / RUNS["t3-opposite"] / f"E{base}q1.txt"), Q(1, factor))
        vectors[oid][0] = ha.scale(ha.parse(SNAP / RUNS["t3-opposite"] / f"E{base}q2.txt"), Q(1, factor))
    # 44 is a product of raw J factors.  Reconstruct its already-defined M
    # generators from their own primitive Hlog outputs, avoiding old I tables.
    m0, m1, m2 = [ha.parse(SNAP / MOTHER_RUN / f"M_{q}_1.txt") for q in range(3)]
    vectors["44"] = {k: {} for k in range(-6, 1)}
    vectors["44"][-1] = ha.scale(ha.mul(m0, ha.add(ha.scale(m0, 2), m1)), Q(-1, 32))
    vectors["44"][0] = ha.scale(ha.add(
        ha.mul(ha.add(ha.scale(ha.ONE, 4), ha.scale(ha.period(2), -1)), ha.mul(m0, m0)),
        ha.scale(ha.mul(m0, m1), 3),
        ha.scale(ha.add(ha.mul(m0, m2), ha.mul(m1, m1)), Q(1, 2))), Q(-1, 32))
    for target, source in T2_MAP.items():
        vectors[target] = vectors[source]
    assert set(vectors) == {oid[-2:] for oid in ids}
    result = {}
    for family, members in FAMILIES.items():
        for oid in members:
            values = vectors[oid]
            assert set(values) == set(range(-6, 1)), oid
            for v in values.values():
                assert all(set(word) <= {0, -1} and len(periods) == 3
                           for word, periods in v), oid
            entry = {
                "orbit_id": "5e-3g-" + oid,
                "family": family,
                "quantity": "I",
                "laurent": {"variable": "eps", "lowest_stored_power": -6,
                            "highest_stored_power": 0, "remainder": "O(eps)",
                            "coefficients": {str(k): {"hlog_terms": ha.serialize(values[k])}
                                             for k in range(-6, 1)}},
                "hlog_convention": {"endpoint": "z", "kernel": "dt/(t-a)",
                                    "alphabet": [0, -1], "period_generators": ["zeta(2)", "zeta(3)", "zeta(5)"],
                                    "domain": "z>0; real principal logarithms"},
                "primitive_inputs": ["snapshots/" + p for p in _primitive_paths(family, oid)],
                "epsilon_budget": _budget(family, oid),
                "normalization": {"I_bridge": "(1-2eps)*G/8*(J(z)+z^(3eps)*J(1/z))/(1+z^(3eps))",
                                  "G": "exp(-3*EulerGamma*eps)*Gamma(1-eps)*Gamma(1-2eps)",
                                  "c4_orbit_multiplicity": 4, "euler_compensation_count": 1,
                                  "native_assembly_sign": 1,
                                  "full_color_reduction": "CF = CA^2*CF - 2*CA*CF^2; select a1 after full reduction"},
                "reassembly": {"method": "exact Gamma/Beta jet and shuffle arithmetic from saved evaluated primitives",
                               "new_cold_integration": False, "formal": "not_run",
                               "final_I_reference_read": False},
            }
            if oid in T2_MAP:
                entry["full_epsilon_relation"] = {"equals": "5e-3g-" + T2_MAP[oid],
                    "scope": "full physical I from ordered-source reflection; z fixed",
                    "edge_map": {"1": 2, "2": 1, "3": 4, "4": 3},
                    "alpha_map": {"alpha1": "alpha3", "alpha2": "alpha2", "alpha3": "alpha1"}}
            if oid in ("30", "38", "39", "24", "25"):
                entry["full_epsilon_relation"] = {"equals": "5e-3g-" + {"30": "26", "38": "10", "39": "09", "24": "23", "25": "22"}[oid],
                    "scope": "full physical I; ordered J transforms as z^(3eps)*J(1/z)"}
            if oid == "44":
                entry["known_function_view"] = {
                    "-1": "-M_0_1*(2*M_0_1+M_1_1)/32",
                    "0": "-((4-zeta(2))*M_0_1^2+3*M_0_1*M_1_1+(M_0_1*M_2_1+M_1_1^2)/2)/32",
                    "status": "reassembled from separately evaluated existing M generators; no new family"}
            result[entry["orbit_id"]] = entry
    return result


if __name__ == "__main__":
    # Inspection only: the package integrator owns result files and freezing.
    data = generate()
    print(json.dumps({"orbits": sorted(data), "count": len(data),
                      "coefficient_count": sum(len(x["laurent"]["coefficients"]) for x in data.values()),
                      "formal": "not_run"}, indent=2))
