"""Hypercube transformation: nested bounds -> [0,1]^N variables.

Algorithm:
  1. Build dependency graph from bounds (s_i depends on s_j if its
     upper bound is s_j; upper bound "0"/"1" means free).
  2. Topological sort (Kahn's algorithm, FIFO queue): free vars
     first, then dependent vars in breadth-first level order.
  3. Assign x variables in topological order:
       - free var   -> x_k            (just renamed)
       - dependent  -> parent_expr * x_k
     Because the parent is processed before its children, the parent's
     fully-expanded x-expression is already available, so no separate
     expansion pass is needed.
  4. Compute the Jacobian as the product of the parent-scale factors
     (one factor per dependent variable), then normalize it: collect
     like x-terms, sort by index, and write exponents as ``xN^c``.
"""
from collections import defaultdict

from montopar.models import Bound, HypercubeResult


def hypercube_transform(bounds: list[Bound]) -> HypercubeResult:
    """Transform nested bounds to [0,1] hypercube variables.

    Given bounds like ``s31 in [0, s32]``, ``s32 in [0, 1]`` this produces:
      - ``s32 -> x1``         (free, just renamed)
      - ``s31 -> x1*x3``       (dependent on s32, multiplied by new x3)
      - Jacobian: ``x1``       (the parent factor s32 = x1)

    Args:
        bounds: List of Bound objects with variable names and bounds.
            The lower bound is assumed to be ``0``; only the upper bound
            is used. An upper bound of ``"1"`` or ``"0"`` means the
            variable is free (already on [0, 1]); any other value is
            interpreted as the name of the parent variable.

    Returns:
        HypercubeResult with s->x mappings, Jacobian, and x variable list.

    Raises:
        ValueError: if a bound references an unknown parent variable or
            if the dependency graph contains a cycle.
    """
    # Build dependency graph: var -> parent_var (the var in its upper bound)
    deps: dict[str, str | None] = {}
    for b in bounds:
        hi = b.high
        if hi == "0" or hi == "1":
            deps[b.var] = None
        else:
            deps[b.var] = hi

    # Validate that every referenced parent is a known variable.
    var_set = set(deps)
    for v, parent in deps.items():
        if parent is not None and parent not in var_set:
            raise ValueError(
                f"Bound for {v} references unknown variable {parent!r}"
            )

    # Topological sort: free vars first, then dependent vars.
    sorted_vars = _topo_sort(list(deps.keys()), deps)

    # Assign x variables in topological order.  Because parents precede
    # their children, s_to_x[parent] is already the fully-expanded
    # x-expression when we build a child's expression.
    s_to_x: dict[str, str] = {}
    jacobian_parents: list[str] = []
    x_idx = 0

    for v in sorted_vars:
        parent = deps[v]
        x_idx += 1
        x_name = f"x{x_idx}"
        if parent is None:
            # Free variable: just rename.
            s_to_x[v] = x_name
        else:
            # Dependent: s_v = parent * x_k, parent already expanded.
            s_to_x[v] = f"{s_to_x[parent]}*{x_name}"
            jacobian_parents.append(parent)

    # Jacobian = product of the parent-scale factors (expanded), then
    # normalized (like terms collected, sorted by x-index).
    jacobian = _normalize_product([s_to_x[p] for p in jacobian_parents])

    x_vars = [f"x{i}" for i in range(1, x_idx + 1)]

    return HypercubeResult(
        s_to_x=s_to_x,
        jacobian=jacobian,
        x_vars=x_vars,
    )


def _topo_sort(variables: list[str], deps: dict[str, str | None]) -> list[str]:
    """Topological sort: depended-on variables before dependent ones.

    Uses Kahn's algorithm with a FIFO queue: the initial queue holds all
    in-degree-0 variables in input order, and children are appended to
    the back as they are released.  This yields a breadth-first level
    order (all free vars across edges, then all depth-1, then all
    depth-2, ...) when the bounds list is interleaved by level -- which
    is the format produced by the orbit parser.
    """
    # Adjacency: parent -> children, plus in-degree per variable.
    children: dict[str, list[str]] = defaultdict(list)
    in_degree: dict[str, int] = {v: 0 for v in variables}

    for v, parent in deps.items():
        if parent is not None:
            children[parent].append(v)
            in_degree[v] += 1

    # Kahn's algorithm.
    queue = [v for v in variables if in_degree[v] == 0]
    result: list[str] = []
    while queue:
        v = queue.pop(0)
        result.append(v)
        for child in children[v]:
            in_degree[child] -= 1
            if in_degree[child] == 0:
                queue.append(child)

    if len(result) != len(variables):
        cyclic = sorted(set(variables) - set(result))
        raise ValueError(
            f"Circular dependency detected in bounds (cycle involves: {cyclic})"
        )

    return result


def _normalize_product(factor_exprs: list[str]) -> str:
    """Normalize a list of x-factor expressions into a sorted product.

    Each input expression is a product of ``xN`` terms (e.g. ``"x1*x3"``).
    The result collects like terms, sorts by x-index, and writes an
    exponent as ``xN^c`` for ``c > 1`` (``xN`` for ``c == 1``).  An empty
    product (no dependent variables) returns ``"1"``.

    Example: ``["x1", "x2", "x1*x3", "x2*x4"]`` -> ``"x1^2*x2^2*x3*x4"``.
    """
    counts: dict[int, int] = {}
    for expr in factor_exprs:
        for term in expr.split("*"):
            term = term.strip()
            if not term:
                continue
            if not term.startswith("x") or not term[1:].isdigit():
                raise ValueError(
                    f"Unexpected Jacobian factor term: {term!r}"
                )
            idx = int(term[1:])
            counts[idx] = counts.get(idx, 0) + 1

    if not counts:
        return "1"

    parts = []
    for idx in sorted(counts):
        c = counts[idx]
        parts.append(f"x{idx}" if c == 1 else f"x{idx}^{c}")
    return "*".join(parts)
