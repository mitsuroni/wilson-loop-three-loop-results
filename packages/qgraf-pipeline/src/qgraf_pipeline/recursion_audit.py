"""Export pre-2PI Wilson diagrams for independent ordered-domain audits.

No colors or integral values are inferred here. Same-edge/scaleless rejection
is explicitly inventoried; an audit consumer must justify the corresponding
zero kernels. Default production classification remains unchanged.
"""
from collections import Counter
import hashlib
import json
from pathlib import Path

from .classifier import classify_orbits, _check_2pi
from .definition_exporter import _build_orbit_definition
from .generator import generate_diagrams
from .internal_types import enumerate_emission_types
from .json_emitter import _orbit_to_dict
from .topology import classify_topology

FAMILIES = {1: {'2e'}, 2: {'4e', '3e-3g'}, 3: {'6e', '5e-3g'}}


def build_audit(n):
    if n not in FAMILIES:
        raise ValueError('Recursion audit supports loop orders 1, 2, 3')
    types = enumerate_emission_types(n)
    diagrams = [d for d in generate_diagrams(n, types)
                if d.emission_type in FAMILIES[n]]
    orbits = classify_orbits(diagrams, types, include_reducible=True)
    membership = {idx: o.orbit_id for o in orbits for idx in o.diagram_indices}
    rows, counts = [], {}
    for d in diagrams:
        t = classify_topology(d)
        keep, reasons = _check_2pi(t)
        state = 'scaleless' if t.is_scaleless else ('production' if keep else 'auxiliary')
        counts.setdefault(d.emission_type, Counter())[state] += 1
        rows.append(dict(record_id=d.idx, emission_type=d.emission_type,
                         edge_counts=t.edge_counts, wl_wl_pairs=t.wl_wl_pairs,
                         internal_connections=t.internal_connections,
                         status=state, exclusion_reasons=([t.scaleless_reason] if
                         t.is_scaleless else reasons), orbit_id=membership.get(d.idx)))
    definitions = []
    for o in orbits:
        definition = _build_orbit_definition(o)
        if definition is None:
            raise ValueError(f'Missing definition: {o.orbit_id}')
        definitions.append(definition)
    sources = {}
    for name in ('generator.py', 'classifier.py', 'topology.py',
                 'definition_exporter.py', 'recursion_audit.py'):
        p = Path(__file__).with_name(name)
        sources[name] = hashlib.sha256(p.read_bytes()).hexdigest()
    return dict(schema_version='recursion-audit/1.0', n_loops=n,
                status='topology_export_only', color_status='not_computed',
                integral_status='not_computed', native_qgraf_binary_used=False,
                enumeration='existing generator; family-dependent edge-distribution preprocessing',
                source_hashes=sources, counts={k: dict(v) for k, v in counts.items()},
                records=rows, orbits=[_orbit_to_dict(o) for o in orbits]), dict(
                n_loops=n, n_orbits=len(definitions), orbits=definitions)


def write_audit(n, output_dir):
    audit, definitions = build_audit(n)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, data in [('recursion-audit.json', audit),
                       ('recursion-orbit-definitions.json', definitions)]:
        (output_dir / name).write_text(json.dumps(data, indent=2)+'\n', encoding='utf8')
    return audit
