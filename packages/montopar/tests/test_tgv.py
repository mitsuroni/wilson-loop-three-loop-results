"""Tests for TGV Pnum and S computation."""
from montopar.models import TGVLeg
from montopar.tgv import compute_S, compute_Pnum


class TestComputeS:
    def test_3e3g_01_S(self):
        """3e-3g-01: legs on edges (3, 4, 4)."""
        legs = [
            TGVLeg(edge=3, s_var="s31", loop="l1"),
            TGVLeg(edge=4, s_var="s41", loop="l2"),
            TGVLeg(edge=4, s_var="s42", loop="l3"),
        ]
        S = compute_S(legs)
        # Δ'_31 = Δ'(4,3) = (1-s31)*s42*z1  (adjacent: 3→4, s31 is on earlier edge 3)
        assert "alpha1*alpha2*((1-s31)*s41*z1)" in S
        assert "alpha3*alpha1*((1-s31)*s42*z1)" in S
        # Δ'_23 should NOT appear (it's zero)
        assert "alpha2*alpha3" not in S

    def test_3e3g_04_S(self):
        """3e-3g-04: legs on edges (2, 3, 4) — non-degenerate."""
        legs = [
            TGVLeg(edge=2, s_var="s21", loop="l1"),
            TGVLeg(edge=3, s_var="s31", loop="l2"),
            TGVLeg(edge=4, s_var="s41", loop="l3"),
        ]
        S = compute_S(legs)
        # All three Δ'_ij are non-zero
        assert "alpha1*alpha2" in S
        assert "alpha2*alpha3" in S
        assert "alpha3*alpha1" in S


class TestComputePnum:
    def test_3e3g_01_Pnum_not_empty(self):
        """3e-3g-01 should produce non-zero Pnum."""
        legs = [
            TGVLeg(edge=3, s_var="s31", loop="l1"),
            TGVLeg(edge=4, s_var="s41", loop="l2"),
            TGVLeg(edge=4, s_var="s42", loop="l3"),
        ]
        Pnum = compute_Pnum(legs)
        assert Pnum != "0"
        # Should contain alpha factors
        assert "alpha1" in Pnum
        assert "alpha2" in Pnum or "alpha3" in Pnum
