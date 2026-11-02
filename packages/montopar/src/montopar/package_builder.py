"""Build type-level initial-integral packages from montopar orbit outputs."""
from __future__ import annotations

import json
from collections import Counter
from copy import deepcopy
from pathlib import Path

from montopar.rotation import classify_structure, cyclic_dist, rotate_edge


_OPPOSITE = {(1, 3), (2, 4)}
_TGV_BASES = {
    "A": (3, 4, 4),
    "B": (3, 3, 4),
    "C": (1, 3, 3),
    "D": (1, 2, 3),
}


def _endpoint_data(raw: dict) -> tuple[dict[str, tuple[int, int]], dict[str, tuple[int, int]]]:
    """Return Lorentz-index -> (edge, slot) and loop -> eikonal endpoints."""
    slots: dict[int, int] = Counter()
    endpoints: dict[str, tuple[int, int]] = {}
    for vertex in raw["vertices"]:
        edge = vertex["line"]
        slots[edge] += 1
        endpoints[vertex["lorentz"]] = (edge, slots[edge])

    loop_ends: dict[str, tuple[int, int]] = {}
    for prop in raw["propagators"]:
        ends = tuple(endpoints[x] for x in prop["lorentz"] if x in endpoints)
        # A direct WL-WL propagator has two eikonal endpoints.  Internal legs
        # intentionally remain absent from this table.
        if len(ends) == 2:
            loop_ends[prop["loop"]] = ends  # type: ignore[assignment]
    return endpoints, loop_ends


def _rotated_fingerprint(raw: dict, rotation: int) -> tuple:
    """Fingerprint the complete attachment/connectivity structure under C4."""
    endpoints, _ = _endpoint_data(raw)

    def rotate_endpoint(point: tuple[int, int]) -> tuple[int, int]:
        return rotate_edge(point[0], rotation), point[1]

    direct = []
    for prop in raw["propagators"]:
        ends = [endpoints[x] for x in prop["lorentz"] if x in endpoints]
        if len(ends) == 2:
            direct.append(tuple(sorted(rotate_endpoint(x) for x in ends)))

    internal = []
    for vertex in raw.get("internal_vertices", []):
        legs = []
        vertex_loops = set(vertex.get("loops", []))
        for prop in raw["propagators"]:
            if prop["loop"] not in vertex_loops:
                continue
            ends = [endpoints[x] for x in prop["lorentz"] if x in endpoints]
            if len(ends) == 1:
                legs.append(rotate_endpoint(ends[0]))
        # Full four-gluon vertex is Bose symmetric: leg labels are dummy.
        if vertex.get("type") in ("4v", "fgv"):
            legs.sort()
        internal.append((vertex.get("type", "unknown"), tuple(legs)))
    return tuple(sorted(direct)), tuple(sorted(internal))


def c4_symmetry(raw: dict) -> dict:
    """Compute the C4 orbit size from a full topology fingerprint."""
    if raw.get("emission_type") == "3e-Gamma1":
        from montopar.gamma1_source import c4_symmetry as gamma1_symmetry
        return gamma1_symmetry(raw)
    if raw.get("emission_type") == "3e-3g+Pi1":
        from montopar.tgv_pi1 import c4_symmetry as tgv_pi1_symmetry
        return tgv_pi1_symmetry(raw)
    if raw.get("emission_type") == "2e-Pi2":
        from montopar.self_energy_pi2 import c4_symmetry as pi2_symmetry
        return pi2_symmetry(raw)
    if raw.get("emission_type") == "2e-Pi1+Pi1":
        from montopar.self_energy_chain import c4_symmetry as chain_symmetry
        return chain_symmetry(raw)
    if raw.get("emission_type") == "4e-3g+3g":
        from montopar.parser import parse_orbit
        from montopar.double_tgv_output import double_tgv_c4_symmetry
        return double_tgv_c4_symmetry(parse_orbit(raw))
    if any(v.get("type") == "se" for v in raw.get("internal_vertices", [])):
        from montopar.self_energy import c4_symmetry as se_symmetry
        return se_symmetry(raw)
    original = _rotated_fingerprint(raw, 0)
    stabilizer = [r for r in range(4) if _rotated_fingerprint(raw, r) == original]
    return {
        "orbit_multiplicity": 4 // len(stabilizer),
        "stabilizer_order": len(stabilizer),
        "stabilizer_rotations": stabilizer,
    }


def direct_opposite_edge_links(raw: dict) -> dict:
    """Project direct eikonal-to-eikonal propagators onto opposite channels."""
    endpoints, _ = _endpoint_data(raw)
    counts: Counter[tuple[int, int]] = Counter()
    for prop in raw["propagators"]:
        ends = [endpoints[x] for x in prop["lorentz"] if x in endpoints]
        if len(ends) == 2:
            pair = tuple(sorted(point[0] for point in ends))
            if pair in _OPPOSITE:
                counts[pair] += 1
    a, b = counts[(1, 3)], counts[(2, 4)]
    factors = []
    for name, power in (("L13", a), ("L24", b)):
        if power:
            factors.append(name if power == 1 else f"{name}^{power}")
    return {
        "definition": "direct_eikonal_to_eikonal_only",
        "channels": {"13": a, "24": b},
        "channel_partition": sorted(x for x in (a, b) if x),
        "structure": " × ".join(factors) if factors else "1",
        "total_order": a + b,
    }


def _tgv_core(raw: dict) -> dict | None:
    """Classify one-TGV cores using the same normalization as mixed montopar."""
    tgv = [v for v in raw.get("internal_vertices", []) if v.get("type") == "tgv"]
    if len(tgv) != 1 or len(tgv[0].get("loops", [])) != 3:
        return None
    endpoints, _ = _endpoint_data(raw)
    loops = set(tgv[0]["loops"])
    legs = []
    for prop in raw["propagators"]:
        if prop["loop"] in loops:
            ends = [endpoints[x] for x in prop["lorentz"] if x in endpoints]
            if len(ends) == 1:
                legs.append(ends[0][0])
    if len(legs) != 3:
        return None
    a, b, c = legs
    # This is the leg reordering used before montopar's canonical rotation.
    if a == b and a != c and cyclic_dist(a, c) == 2:
        legs = [c, a, b]
    elif a == c and a != b and cyclic_dist(a, b) == 2:
        legs = [b, a, c]
    elif len(set(legs)) == 3:
        for shift in (1, 2):
            candidate = legs[shift:] + legs[:shift]
            if cyclic_dist(candidate[0], candidate[1]) == 1 and cyclic_dist(candidate[1], candidate[2]) == 1:
                legs = candidate
                break
    structure = classify_structure(*legs)
    label = {
        "degenerate_bc_adjacent": "A",
        "degenerate_ab_adjacent": "B",
        "degenerate_bc_opposite": "C",
        "nondegenerate_consecutive": "D",
    }.get(structure)
    if label is None:
        raise ValueError(f"unsupported TGV core {tuple(legs)} ({structure})")
    return {
        "kind": "tgv",
        "topology_class": label,
        "canonical_attachment_edges": list(_TGV_BASES[label]),
    }


def integral_structure(raw: dict) -> dict:
    """Return the human- and machine-readable integral skeleton."""
    links = direct_opposite_edge_links(raw)
    if raw.get("emission_type") == "3e-Gamma1":
        from montopar.gamma1_source import validate_source
        data = validate_source(raw)
        return {"core": {"kind": "gamma1", "operator": "Gamma1", "subgraph_loop_order": 1,
                         "external_legs": 3, "independent_external_momenta": 2,
                         "endpoints": [[p["edge"], p["slot"]] for p in data["endpoints"]]},
                "direct_opposite_edge_links": links, "symbol": "central one-loop 1PI TGV"}
    if raw.get("emission_type") == "3e-3g+Pi1":
        from montopar.tgv_pi1 import validate_source
        data = validate_source(raw)
        return {"core": {"kind": "tgv_pi1", "operator": "Pi1", "subgraph_loop_order": 1,
                         "marked_leg": data["marked_leg"], "external_legs": 3,
                         "endpoints": [[p["edge"], p["slot"]] for p in data["endpoints"]]},
                "direct_opposite_edge_links": links, "symbol": "tree-TGV with Pi1 leg"}
    if raw.get("emission_type") == "2e-Pi2":
        from montopar.self_energy_pi2 import validate_source
        data = validate_source(raw)
        return {"core": {"kind": "self_energy_pi2", "operator": "Pi2", "insertion_count": 1,
                         "subgraph_loop_order": 2, "operator_kind": "1pi",
                         "endpoints": [[p["edge"], p["slot"]] for p in data["endpoints"]]},
                "direct_opposite_edge_links": links, "symbol": "Pi2"}
    if raw.get("emission_type") == "2e-Pi1+Pi1":
        from montopar.self_energy_chain import validate_source
        data = validate_source(raw)
        return {"core": {"kind": "self_energy_chain", "operator": "Pi1+Pi1", "insertion_count": 2,
                         "subgraph_loop_orders": [1, 1], "middle_propagators": 1,
                         "endpoints": [[p["edge"], p["slot"]] for p in data["endpoints"]]},
                "direct_opposite_edge_links": links, "symbol": "Pi1-D0-Pi1"}
    if any(v.get("type") == "se" for v in raw.get("internal_vertices", [])):
        from montopar.self_energy import validate_source
        data = validate_source(raw)
        se = data["lines"][0]
        return {"core": {"kind": "self_energy", "operator": "Pi1", "subgraph_loop_order": 1,
                         "endpoints": [[p["edge"], p["slot"]] for p in se["endpoints"]]},
                "direct_opposite_edge_links": links, "symbol": "Pi1 x " + links["structure"]}
    if raw.get("emission_type") == "4e-3g+3g":
        return {"core": {"kind": "double_tgv", "external_legs": 4,
                         "internal_vertices": 2, "internal_connectors": 1,
                         "oriented_pairing": "source vertex loop order"},
                "direct_opposite_edge_links": links, "symbol": "double-TGV"}
    tgv = _tgv_core(raw)
    if tgv:
        core = tgv
        prefix = tgv["topology_class"]
    elif not raw.get("internal_vertices"):
        core = {"kind": "identity"}
        prefix = ""
    else:
        kind = raw["internal_vertices"][0].get("type", "internal")
        core = {"kind": kind}
        prefix = kind
    link_symbol = links["structure"]
    symbol = " × ".join(x for x in (prefix, "" if link_symbol == "1" else link_symbol) if x)
    return {"core": core, "direct_opposite_edge_links": links, "symbol": symbol or "1"}


def build_initial_package(source: dict, workspace: Path, emission_type: str) -> tuple[dict, Path]:
    """Read existing orbit files and construct one type-level package in memory."""
    raws = [o for o in source["orbits"] if o["emission_type"] == emission_type]
    if not raws:
        raise ValueError(f"no orbits of emission type {emission_type}")
    nloops = raws[0]["nloops"]
    target_dir = workspace / f"{nloops}loop" / emission_type
    records = []
    geometry = None
    for raw in sorted(raws, key=lambda item: item["orbit_id"]):
        orbit_path = target_dir / raw["orbit_id"] / "montopar.json"
        if not orbit_path.is_file():
            raise FileNotFoundError(f"missing orbit result: {orbit_path}")
        result = json.loads(orbit_path.read_text(encoding="utf-8"))
        if "C" not in result or "J_raw" not in result:
            raise ValueError(f"invalid montopar result: {orbit_path}")
        if ((raw.get("emission_type") == "3e-Gamma1") != (result.get("dispatch") == "gamma1")):
            raise ValueError(f"gamma1_result_dispatch_mismatch: {orbit_path}")
        expected_se_dispatch = "tgv_pi1" if raw.get("emission_type") == "3e-3g+Pi1" else "self_energy_pi2" if raw.get("emission_type") == "2e-Pi2" else "self_energy_chain" if raw.get("emission_type") == "2e-Pi1+Pi1" else "self_energy"
        if (raw.get("emission_type") == "2e-Pi1+Pi1" or any(v.get("type") == "se" for v in raw.get("internal_vertices", []))) and result.get("dispatch") != expected_se_dispatch:
            raise ValueError(f"self_energy_result_dispatch_mismatch: {orbit_path}")
        geometry = geometry or deepcopy(result.get("geometry", {}))
        record = {
            "orbit_id": raw["orbit_id"],
            "dispatch": result.get("dispatch", "unknown"),
            "topology": {key: deepcopy(raw.get(key, [])) for key in (
                "vertices", "propagators", "internal_vertices", "phases", "bounds"
            )},
            "c4_symmetry": c4_symmetry(raw),
            "integral_structure": integral_structure(raw),
            "C": deepcopy(result["C"]),
            "J": deepcopy(result["J_raw"]),
        }
        if 'physical_convention' in raw:
            from montopar.sign_rules import validate_contract
            record['physical_convention']=validate_contract(raw,required=True)
            if result.get('physical_convention')!=record['physical_convention']:
                raise ValueError(f'regenerate_package_physical_convention_mismatch: {orbit_path}')
        if result.get("dispatch") == "gamma1":
            from montopar.gamma1 import validate_result
            validate_result(result)
            if result["representative"]["source_topology"] != raw:
                raise ValueError(f"gamma1_package_source_mismatch: {orbit_path}")
            record["source_topology"] = record.pop("topology")
            record["topology"] = deepcopy(raw)
            for key in ("schema_version", "representative", "source_sha256", "J_raw_s", "J_sector",
                        "normalization", "operator", "color_contraction", "branch_prescription",
                        "translation_status", "evaluation_status", "formal"):
                record[key] = deepcopy(result[key])
        if result.get("dispatch") in ("self_energy", "self_energy_chain", "self_energy_pi2", "tgv_pi1"):
            if result["dispatch"] == "tgv_pi1":
                from montopar.tgv_pi1 import validate_result
            elif result["dispatch"] == "self_energy_pi2":
                from montopar.self_energy_pi2 import validate_result
            elif result["dispatch"] == "self_energy_chain":
                from montopar.self_energy_chain import validate_result
            else:
                from montopar.self_energy import validate_result
            validate_result(result)
            if result["representative"]["source_topology"] != raw:
                raise ValueError(f"self_energy_package_source_mismatch: {orbit_path}")
            record["source_topology"] = record.pop("topology")
            record["topology"] = deepcopy(raw)
            for key in ("schema_version", "representative", "J_raw_s", "normalization", "operator",
                        "color_contraction", "branch_prescription", "translation_status", "evaluation_status", "formal"):
                record[key] = deepcopy(result[key])
            if result["dispatch"] == "tgv_pi1":
                record["source_sha256"] = result["source_sha256"]
                if "J_ibp" in result:
                    record["J_ibp"] = deepcopy(result["J_ibp"])
        if result.get("dispatch") in ("tgv", "mixed", "fourv", "double_tgv"):
            if result.get("dispatch") == "double_tgv":
                from montopar.double_tgv_output import validate_double_tgv_result
                validate_double_tgv_result(result)
            if result.get("dispatch") == "fourv":
                from montopar.fgv_output import validate_fgv_result
                validate_fgv_result(result)
            if result.get("schema_version") not in ("1.1", "1.2", "1.3") or "representative" not in result:
                raise ValueError(f"regenerate canonical self-contained TGV source: {orbit_path}")
            if "J_raw_s" not in result or "alpha_measure" not in result["J_raw_s"]:
                raise ValueError(f"missing ordered source semantics: {orbit_path}")
            record["source_topology"] = record.pop("topology")
            record["topology"] = deepcopy(result["representative"]["canonical_topology"])
            record["representative"] = deepcopy(result["representative"])
            record["J_raw_s"] = deepcopy(result["J_raw_s"])
            record["normalization"] = deepcopy(result["normalization"])
            if result.get("dispatch") in ("fourv", "double_tgv"):
                record["branch_prescription"] = deepcopy(result["branch_prescription"])
                record["evaluation_status"] = result["evaluation_status"]
                record["formal"] = result["formal"]
        if result.get("dispatch") == "double_tgv":
            if result["representative"]["source_topology"] != raw:
                raise ValueError(f"double_tgv_package_source_mismatch: {orbit_path}")
            record["color_contraction"] = deepcopy(result["color_contraction"])
            record["zero_reason"] = result["zero_reason"]
            record["zero_policy"] = result["zero_policy"]
            record["schema_version"] = result["schema_version"]
            for key in ("J_compact", "J_reduced", "J_operator"):
                if key in result:
                    record[key] = deepcopy(result[key])
        if "J_ibp" in result:
            record["J_ibp"] = deepcopy(result["J_ibp"])
        records.append(record)
    return {
        "schema_version": "1.8" if any(r["dispatch"] == "gamma1" for r in records) else "1.7" if any(r["dispatch"] == "tgv_pi1" for r in records) else "1.6" if any(r["dispatch"] == "self_energy_pi2" for r in records) else "1.5" if any(r["dispatch"] == "self_energy_chain" for r in records) else "1.4" if any(r["dispatch"] == "self_energy" for r in records) else "1.3" if any(r["dispatch"] == "double_tgv" for r in records) else "1.2" if any(r["dispatch"] == "fourv" for r in records) else ("1.1" if any("representative" in r for r in records) else "1.0"),
        "package_kind": "montopar-initial-integrals",
        "loop_order": nloops,
        "emission_type": emission_type,
        "geometry": geometry or {},
        "orbits": records,
    }, target_dir / "montopar-initial.json"


def write_initial_package(source: dict, workspace: Path, emission_type: str) -> Path:
    """Build and atomically write one package after all input validation succeeds."""
    package, path = build_initial_package(source, workspace, emission_type)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(package, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)
    return path
