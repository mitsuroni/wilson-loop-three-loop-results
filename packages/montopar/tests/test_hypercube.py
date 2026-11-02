"""Tests for montopar.hypercube."""
from montopar.models import Bound
from montopar.hypercube import hypercube_transform


def test_hypercube_4e01():
    """4e-01: edges 3 and 4 each have 2 attachment points."""
    bounds = [
        Bound(var="s31", low="0", high="s32"),
        Bound(var="s41", low="0", high="s42"),
        Bound(var="s32", low="0", high="1"),
        Bound(var="s42", low="0", high="1"),
    ]
    result = hypercube_transform(bounds)

    # s32 is free -> x1, s42 is free -> x2
    # s31 depends on s32 -> s32*x3 -> x1*x3
    # s41 depends on s42 -> s42*x4 -> x2*x4
    assert result.s_to_x["s32"] == "x1"
    assert result.s_to_x["s42"] == "x2"
    assert result.s_to_x["s31"] == "x1*x3"
    assert result.s_to_x["s41"] == "x2*x4"

    # Jacobian: s32 * s42 = x1 * x2
    assert result.jacobian == "x1*x2"

    assert result.x_vars == ["x1", "x2", "x3", "x4"]


def test_hypercube_6e01():
    """6e-01: edges 3 and 4 each have 3 attachment points."""
    bounds = [
        Bound(var="s31", low="0", high="s32"),
        Bound(var="s41", low="0", high="s42"),
        Bound(var="s32", low="0", high="s33"),
        Bound(var="s42", low="0", high="s43"),
        Bound(var="s33", low="0", high="1"),
        Bound(var="s43", low="0", high="1"),
    ]
    result = hypercube_transform(bounds)

    # Free: s33->x1, s43->x2
    assert result.s_to_x["s33"] == "x1"
    assert result.s_to_x["s43"] == "x2"
    # s32 = s33*x3 -> x1*x3
    assert result.s_to_x["s32"] == "x1*x3"
    # s42 = s43*x4 -> x2*x4
    assert result.s_to_x["s42"] == "x2*x4"
    # s31 = s32*x5 -> x1*x3*x5
    assert result.s_to_x["s31"] == "x1*x3*x5"
    # s41 = s42*x6 -> x2*x4*x6
    assert result.s_to_x["s41"] == "x2*x4*x6"

    # Jacobian: s33 * s43 * s32 * s42
    # = x1 * x2 * (x1*x3) * (x2*x4) = x1^2*x2^2*x3*x4
    assert result.jacobian == "x1^2*x2^2*x3*x4"

    assert len(result.x_vars) == 6


def test_hypercube_all_free():
    """All variables already [0,1] -- just rename."""
    bounds = [
        Bound(var="s31", low="0", high="1"),
        Bound(var="s41", low="0", high="1"),
    ]
    result = hypercube_transform(bounds)

    assert result.s_to_x["s31"] == "x1"
    assert result.s_to_x["s41"] == "x2"
    assert result.jacobian == "1"
    assert result.x_vars == ["x1", "x2"]


def test_hypercube_mixed_edges():
    """Edge 2 has 1 point, edge 3 has 2, edge 4 has 3.

    NOTE -- expected values corrected from the original brief.

    The original brief encoded a per-edge DFS labeling here
    (s31 == "x2*x3", s43 == "x4", s42 == "x4*x5", s41 == "x4*x5*x6",
    jacobian == "x2*x4^2*x5").  That labeling is inconsistent with
    test_hypercube_4e01 / test_hypercube_6e01 (and with the real
    4e-01 fixture data in tests/fixtures/), which use a BFS-by-level
    labeling produced by Kahn's algorithm with a FIFO queue.  The
    bounds list here is grouped per-edge rather than interleaved by
    level, but the *algorithm* must be the same across all orbits, so
    the expected values below are the BFS-consistent ones.
    """
    bounds = [
        Bound(var="s21", low="0", high="1"),
        Bound(var="s31", low="0", high="s32"),
        Bound(var="s32", low="0", high="1"),
        Bound(var="s41", low="0", high="s42"),
        Bound(var="s42", low="0", high="s43"),
        Bound(var="s43", low="0", high="1"),
    ]
    result = hypercube_transform(bounds)

    # Free (BFS, in input order among in-degree-0 vars):
    #   s21 -> x1, s32 -> x2, s43 -> x3
    assert result.s_to_x["s21"] == "x1"
    assert result.s_to_x["s32"] == "x2"
    assert result.s_to_x["s43"] == "x3"
    # Dependent: s31 = s32*x4 -> x2*x4
    assert result.s_to_x["s31"] == "x2*x4"
    # s42 = s43*x5 -> x3*x5
    assert result.s_to_x["s42"] == "x3*x5"
    # s41 = s42*x6 -> x3*x5*x6
    assert result.s_to_x["s41"] == "x3*x5*x6"

    # Jacobian: s32 * s43 * s42 = x2 * x3 * (x3*x5) = x2*x3^2*x5
    assert result.jacobian == "x2*x3^2*x5"

    assert len(result.x_vars) == 6


def test_hypercube_circular_dependency_raises():
    """Circular bounds should raise ValueError."""
    bounds = [
        Bound(var="s31", low="0", high="s32"),
        Bound(var="s32", low="0", high="s31"),
    ]
    try:
        hypercube_transform(bounds)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass
