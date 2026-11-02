"""Tests for montopar.c_factor."""
from montopar.c_factor import (
    compute_c, compute_c_tgv, compute_c_star, compute_c_fgv, compute_c_mixed,
)


def test_compute_c_n2():
    """C_n for n=2 (4e)."""
    c = compute_c(n=2, color="-1/2*((CA - 2*CF)*CF)")

    assert c.color == "-1/2*((CA - 2*CF)*CF)"
    assert c.coupling == "(4*Pi*alphaS)^2"
    assert c.ft_base == "Gamma[1-eps]^2/(4*Pi^(2-eps))^2"
    assert c.v14_scale == "(-2*v14+I0)^(2*eps)"
    assert c.twoscale == "(-2)^(-2)"


def test_compute_c_n3():
    """C_n for n=3 (6e)."""
    c = compute_c(n=3, color="CF*(CF - CA/2)^2")

    assert c.coupling == "(4*Pi*alphaS)^3"
    assert c.ft_base == "Gamma[1-eps]^3/(4*Pi^(2-eps))^3"
    assert c.v14_scale == "(-2*v14+I0)^(3*eps)"
    assert c.twoscale == "(-2)^(-3)"


def test_compute_c_n1():
    """C_n for n=1 (2e)."""
    c = compute_c(n=1, color="CF")

    assert c.coupling == "(4*Pi*alphaS)^1"
    assert c.v14_scale == "(-2*v14+I0)^(1*eps)"


def test_compute_c_custom_vref():
    """C_n with custom reference scale."""
    c = compute_c(n=2, color="CF", v_ref="vRef")

    assert "vRef" in c.v14_scale


def test_compute_c_tgv():
    """TGV C: 6 factors, all independent of edge distribution."""
    c = compute_c_tgv()
    assert c.color == "-I/2*CA*CF"
    assert c.coupling == "(4*Pi*alphaS)^2"
    assert "Gamma[1-2*eps]" in c.ft_base
    assert "64*Pi^(4-2*eps)" in c.ft_base
    assert c.v14_scale == "(-2*v14+I0)^(2*eps)"
    assert c.twoscale == "1"
    assert c.derivative == "2*(1-2*eps)"


def test_compute_c_star():
    c = compute_c_star()
    assert c.color == "CA*CF"
    assert c.coupling == "(4*Pi*alphaS)^2"
    assert c.ft_base == "Gamma[1-2*eps]/(256*Pi^(4-2*eps))"
    assert c.v14_scale == "(-2*v14+I0)^(2*eps)"
    assert c.derivative is None


def test_compute_c_fgv():
    """FGV C: 6 factors, lorentzZ is in J."""
    c = compute_c_fgv()
    assert c.color == "CA^2*CF"
    assert c.coupling == "(4*Pi*alphaS)^3"
    assert "Gamma[2-3*eps]" in c.ft_base
    assert "256*Pi^(6-3*eps)" in c.ft_base
    assert c.v14_scale == "(-2*v14+I0)^(3*eps)"
    assert c.twoscale == "1/4"
    assert c.derivative is None


def test_compute_c_mixed_n1():
    """Mixed C for n_pure=1 (5e-3g)."""
    c = compute_c_mixed(n_pure=1)
    assert c.color == "MIXED_COLOR_PLACEHOLDER"
    assert c.coupling == "(4*Pi*alphaS)^3"
    assert "Gamma[1-eps]^1*Gamma[1-2*eps]" in c.ft_base
    assert "64*Pi^(4-2*eps)" in c.ft_base
    assert c.v14_scale == "(-2*v14+I0)^(3*eps)"
    assert c.twoscale == "(-2)^(-1)"
    assert c.derivative == "2*(1-2*eps)"


def test_compute_c_mixed_n2():
    """Mixed C for n_pure=2 (7e-3g)."""
    c = compute_c_mixed(n_pure=2)
    assert c.coupling == "(4*Pi*alphaS)^4"
    assert "Gamma[1-eps]^2" in c.ft_base
    assert c.v14_scale == "(-2*v14+I0)^(4*eps)"
    assert c.twoscale == "(-2)^(-2)"


def test_compute_c_mixed_no_z():
    """C must not contain z1 (z belongs in J)."""
    c = compute_c_mixed(n_pure=1)
    full = f"{c.color} {c.coupling} {c.ft_base} {c.v14_scale} {c.twoscale} {c.derivative}"
    assert "z1" not in full
    assert "z" not in full.replace("alphaS", "").replace("eps", "")
