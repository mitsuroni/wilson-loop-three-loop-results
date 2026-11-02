"""JSON emitter: serialize classified orbits to JSON."""
from __future__ import annotations

import json
from pathlib import Path
from typing import List

from .classifier import ClassifiedOrbit


def _orbit_to_dict(o: ClassifiedOrbit) -> dict:
    return {
        "orbit_id": o.orbit_id,
        "emission_type": o.emission_type,
        "category": o.category,
        "edge_counts": list(o.edge_counts),
        "eikonal_count": o.eikonal_count,
        "g_power": o.g_power,
        "n_diagrams": o.n_diagrams,
        "wl_wl_pairs": o.wl_wl_pairs,
        "mirror_relation": o.mirror_relation,
        "mirror_of": o.mirror_of,
        "is_2pi": o.is_2pi,
        "reducible_to": o.reducible_to,
        "diagram_indices": o.diagram_indices,
    }


def write_json(orbits: List[ClassifiedOrbit], n: int, path: Path) -> None:
    data = {
        "n_loops": n,
        "n_orbits": len(orbits),
        "orbits": [_orbit_to_dict(o) for o in orbits],
    }
    path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
