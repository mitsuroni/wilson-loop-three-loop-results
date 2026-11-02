"""Tests for montopar.color."""
import pytest

from montopar.models import Vertex, Propagator
from montopar.color import (
    compute_color,
    build_sunt_indices,
    build_color_trace_script,
    parse_color_result,
    _wolframscript_available,
    COLOR_TRACE_TEMPLATE,
)

WOLFRAMSCRIPT_AVAILABLE = _wolframscript_available()


def test_build_sunt_indices_4e01():
    """4e-01: vertices on edges 3,3,4,4 -> SUNT[a1,a2,a3,a4]."""
    vertices = [
        Vertex(line=3, lorentz="mu1", color="a1"),
        Vertex(line=3, lorentz="mu2", color="a2"),
        Vertex(line=4, lorentz="mu3", color="a3"),
        Vertex(line=4, lorentz="mu4", color="a4"),
    ]
    propagators = [
        Propagator(loop="l1", lorentz=("mu1", "mu3"), color=("a1", "a3")),
        Propagator(loop="l2", lorentz=("mu2", "mu4"), color=("a2", "a4")),
    ]
    indices = build_sunt_indices(vertices, propagators)

    # Contour order: edge 1,2,3,4. Vertices on edge 3 first, then edge 4.
    # Indices are ordered by contour position.
    # All vertices on edge 3: a1, a2. Then edge 4: a3, a4.
    assert indices == ["a4", "a3", "a2", "a1"]


def test_build_sunt_indices_mixed_edges():
    """Vertices on edges 2,3,4 -> ordered by contour position."""
    vertices = [
        Vertex(line=4, lorentz="mu1", color="a1"),
        Vertex(line=2, lorentz="mu2", color="a2"),
        Vertex(line=3, lorentz="mu3", color="a3"),
    ]
    propagators = [
        Propagator(loop="l1", lorentz=("mu1", "mu2"), color=("a1", "a2")),
        Propagator(loop="l2", lorentz=("mu2", "mu3"), color=("a2", "a3")),
    ]
    indices = build_sunt_indices(vertices, propagators)

    # Contour order: edge 2 (a2), edge 3 (a3), edge 4 (a1)
    assert indices == ["a1", "a3", "a2"]


def test_build_color_trace_script():
    """Verify the generated Mathematica script structure."""
    script = build_color_trace_script(["a1", "a2", "a3", "a4"])
    assert "FeynCalc`" in script
    assert "SUNT[a1, a2, a3, a4]" in script
    assert "SUNNToCACF -> True" in script
    assert "InputForm" in script


def test_parse_color_result():
    """Parse typical SUNSimplify output."""
    # Example: SUNSimplify gives something like (-CA*CF)/2 + CF^2
    result = parse_color_result("(-CA*CF)/2 + CF^2")
    assert "CA" in result
    assert "CF" in result


def test_parse_color_result_with_spaces():
    result = parse_color_result("  -1/2*((CA - 2*CF)*CF)  ")
    assert result == "-1/2*((CA - 2*CF)*CF)"


def test_build_color_trace_script_with_tgv():
    """Template includes SUNF when tgv_indices provided."""
    script = build_color_trace_script(["a1", "a2", "a3", "a4", "a5"], ["a1", "a2", "a4"])
    assert "SUNF[a1, a2, a4]" in script
    assert "SUNT[a1, a2, a3, a4, a5]" in script


def test_build_color_trace_script_without_tgv():
    """Template uses 1 when no tgv_indices."""
    script = build_color_trace_script(["a1", "a2", "a3", "a4"])
    assert "SUNF" not in script
    assert "* 1" in script


def test_color_trace_template_exists():
    """Verify the template file exists and is non-empty."""
    import montopar.color
    assert COLOR_TRACE_TEMPLATE.exists()
    content = COLOR_TRACE_TEMPLATE.read_text()
    assert "SUNT" in content
    assert "SUNSimplify" in content


@pytest.mark.slow
@pytest.mark.skipif(not WOLFRAMSCRIPT_AVAILABLE, reason="wolframscript not available")
def test_compute_color_integration():
    """Integration test: actual Mathematica call for 4e-01."""
    vertices = [
        Vertex(line=3, lorentz="mu1", color="a1"),
        Vertex(line=3, lorentz="mu2", color="a2"),
        Vertex(line=4, lorentz="mu3", color="a3"),
        Vertex(line=4, lorentz="mu4", color="a4"),
    ]
    propagators = [
        Propagator(loop="l1", lorentz=("mu1", "mu3"), color=("a1", "a3")),
        Propagator(loop="l2", lorentz=("mu2", "mu4"), color=("a2", "a4")),
    ]
    result = compute_color(vertices, propagators)

    # Should contain CA and CF
    assert "CA" in result
    assert "CF" in result
    # Should be a valid expression (no error markers)
    assert "SUNSimplify" not in result
