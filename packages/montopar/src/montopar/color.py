"""Color factor computation via Mathematica SUNSimplify.

Calls wolframscript to evaluate SUNSimplify[SUNT[...], SUNNToCACF -> True].
"""
import shutil
import subprocess
import tempfile
from pathlib import Path
from functools import lru_cache

from montopar.models import Vertex, Propagator
from montopar.sign_rules import transport_word
from montopar.sign_rules import transport_word

COLOR_TRACE_TEMPLATE = Path(__file__).parent / "templates" / "color_trace.wl"


def build_sunt_indices(
    vertices: list[Vertex], propagators: list[Propagator]
) -> list[str]:
    """Build the ordered list of color indices for SUNT[...].

    Indices are ordered by contour position: edge 1, then 2, then 3, then 4.
    Within each edge, vertices are ordered by their appearance in the input.
    Propagators handle the delta^{ab} contractions.
    """
    # Group vertices by edge, preserving order within each edge
    edge_order: dict[int, list[str]] = {}
    for v in vertices:
        edge_order.setdefault(v.line, []).append(v.color)

    # Flatten in contour order (edges 1,2,3,4)
    indices: list[str] = []
    for edge in [1, 2, 3, 4]:
        indices.extend(edge_order.get(edge, []))

    return transport_word(vertices)


def build_color_trace_script(
    indices: list[str],
    tgv_indices: list[str] | None = None,
    tgv_indices_2: list[str] | None = None,
) -> str:
    """Generate the Mathematica script for color trace computation.

    Args:
        indices: Ordered list of color index labels (e.g. ["a1","a2","a3","a4"]).
        tgv_indices: TGV color indices for first SUNF[a,b,c] factor.
        tgv_indices_2: TGV color indices for second SUNF[d,e,f] factor
            (used for double_tgv with two TGV vertices).

    Returns:
        Complete Mathematica script string.
    """
    template = COLOR_TRACE_TEMPLATE.read_text(encoding="utf-8")
    indices_str = ", ".join(indices)
    script = template.replace("{{INDICES}}", indices_str)

    tgv_parts = []
    if tgv_indices:
        tgv_str = ", ".join(tgv_indices)
        tgv_parts.append(f"SUNF[{tgv_str}]")
    if tgv_indices_2:
        tgv_str = ", ".join(tgv_indices_2)
        tgv_parts.append(f"SUNF[{tgv_str}]")

    if tgv_parts:
        script = script.replace("{{TGV_FACTOR}}", " * ".join(tgv_parts))
    else:
        script = script.replace("{{TGV_FACTOR}}", "1")
    return script


def parse_color_result(raw_output: str) -> str:
    """Parse the color factor from wolframscript output.

    Extracts the InputForm[...] expression line from the Mathematica output.
    Skips FeynCalc banner / comment / directive lines and returns the last
    meaningful line (the ``Print`` output appears last).

    Args:
        raw_output: stdout from wolframscript.

    Returns:
        Cleaned color factor expression string.
    """
    candidates: list[str] = []
    for line in raw_output.strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        if line.startswith("(*") or line.startswith("<<"):
            continue
        if line.startswith("Loading "):
            # FeynCalc banner line, e.g.
            # "Loading FeynCalc from C:\\...\\FeynCalc\\"
            continue
        candidates.append(line)
    if not candidates:
        raise ValueError(f"No color result found in output: {raw_output}")
    # The Print[InputForm[result]] line is emitted last.
    expr = candidates[-1]
    # Strip a literal InputForm[...] wrapper, which wolframscript prints
    # verbatim when there is no front end to format it.
    if expr.startswith("InputForm[") and expr.endswith("]"):
        expr = expr[len("InputForm["):-1].strip()
    return expr


def _contract_indices(
    indices: list[str], propagators: list[Propagator]
) -> list[str]:
    """Collapse color labels joined by propagator delta^{ab} contractions.

    Each ``Propagator.color`` pair ``(a, b)`` denotes a delta^{ab} that
    identifies the two labels. We union all such labels and, per group, use
    the member that appears first in ``indices`` as the display label so the
    resulting SUNT has repeated indices that SUNSimplify can reduce.

    Args:
        indices: Vertex color labels in contour order.
        propagators: Gluon propagators carrying color pairs.

    Returns:
        Contracted list of color labels (same length as ``indices``).
    """
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        root = x
        while parent[root] != root:
            root = parent[root]
        # path compression
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for p in propagators:
        if p.color and len(p.color) == 2:
            union(p.color[0], p.color[1])

    # Representative label = first-appearing member of each group.
    first_seen: dict[str, str] = {}
    for idx in indices:
        r = find(idx)
        first_seen.setdefault(r, idx)

    return [first_seen[find(idx)] for idx in indices]


def compute_color(
    vertices: list[Vertex],
    propagators: list[Propagator],
    tgv_color_indices: list[str] | None = None,
    tgv_color_indices_2: list[str] | None = None,
) -> str:
    """Compute the color factor via Mathematica SUNSimplify.

    Args:
        vertices: List of eikonal vertices.
        propagators: List of gluon propagators.
        tgv_color_indices: TGV color indices for SUNF[a,b,c].
            None for pure n-gluon (no TGV).
        tgv_color_indices_2: Second TGV color indices for double TGV.

    Returns:
        Color factor string in C_F, C_A basis (e.g. "-1/2*((CA - 2*CF)*CF)").
    """
    indices = build_sunt_indices(vertices, propagators)
    indices = _contract_indices(indices, propagators)
    factors=[tuple(v) for v in [tgv_color_indices,tgv_color_indices_2] if v]
    # Trace cyclicity and dummy-colour renaming preserve the same tensor.
    # The physical reversal has already happened in build_sunt_indices.
    candidates=[]
    for offset in range(len(indices)):
        word=indices[offset:]+indices[:offset]
        names={}
        def label(value):
            if value not in names:names[value]='c'+str(len(names)+1)
            return names[value]
        renamed=tuple(label(value) for value in word)
        ff=[]
        for factor in factors:
            mapped=tuple(label(value) for value in factor)
            ff.append(min(mapped[j:]+mapped[:j] for j in range(3)))
        candidates.append((renamed,tuple(sorted(ff))))
    if not candidates:raise ValueError('empty_Wilson_colour_word')
    return _evaluate_color_tensor(min(candidates))


@lru_cache(maxsize=512)
def _evaluate_color_tensor(key):
    indices,factors=key
    script=build_color_trace_script(list(indices),list(factors[0]) if factors else None,
        list(factors[1]) if len(factors)>1 else None)

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".wl", delete=False, encoding="utf-8"
    ) as f:
        f.write(script)
        tmp_path = f.name

    try:
        result = subprocess.run(
            ["wolframscript", "-file", tmp_path],
            capture_output=True,
            text=True,
            timeout=60,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"wolframscript failed (code {result.returncode}):\n"
                f"stderr: {result.stderr}\n"
                f"stdout: {result.stdout}"
            )
        return parse_color_result(result.stdout)
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def _wolframscript_available() -> bool:
    """Check whether wolframscript is on PATH."""
    return shutil.which("wolframscript") is not None
