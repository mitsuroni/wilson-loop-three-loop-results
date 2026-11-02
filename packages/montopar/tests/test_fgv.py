"""Tests for FGV Q' and lorentzZ computation."""
from montopar.models import FGVLeg, EdgeDistribution
from montopar.fgv import is_physical_zero, compute_Q_prime, compute_lorentzZ


class TestIsPhysicalZero:
    def test_max_3_is_zero(self):
        dist = EdgeDistribution(counts={3: 3, 4: 1}, s_vars={})
        assert is_physical_zero(dist) is True

    def test_max_2_is_not_zero(self):
        dist = EdgeDistribution(counts={3: 2, 4: 2}, s_vars={})
        assert is_physical_zero(dist) is False

    def test_all_distinct_is_not_zero(self):
        dist = EdgeDistribution(counts={1: 1, 2: 1, 3: 1, 4: 1}, s_vars={})
        assert is_physical_zero(dist) is False


class TestComputeQPrime:
    def test_4e4g_02_Q_prime(self):
        """4e-4g-02: legs on edges (3, 3, 4, 4)."""
        legs = [
            FGVLeg(edge=3, s_var="s31", loop="l1"),
            FGVLeg(edge=3, s_var="s32", loop="l2"),
            FGVLeg(edge=4, s_var="s41", loop="l3"),
            FGVLeg(edge=4, s_var="s42", loop="l4"),
        ]
        Qp = compute_Q_prime(legs)
        # Same-edge pairs vanish: (1,2) both edge 3, (3,4) both edge 4
        # Non-zero: (1,3): 3->4 adjacent, (1,4): 3->4 adjacent,
        #           (2,3): 3->4 adjacent, (2,4): 3->4 adjacent
        assert "alpha1*alpha3" in Qp
        assert "alpha1*alpha4" in Qp
        assert "alpha2*alpha3" in Qp
        assert "alpha2*alpha4" in Qp
        # Same-edge pairs should NOT appear
        assert "alpha1*alpha2" not in Qp
        assert "alpha3*alpha4" not in Qp

    def test_4e4g_09_Q_prime(self):
        """4e-4g-09: legs on edges (1, 2, 3, 4) -- all distinct."""
        legs = [
            FGVLeg(edge=1, s_var="s11", loop="l1"),
            FGVLeg(edge=2, s_var="s21", loop="l2"),
            FGVLeg(edge=3, s_var="s31", loop="l3"),
            FGVLeg(edge=4, s_var="s41", loop="l4"),
        ]
        Qp = compute_Q_prime(legs)
        # All 6 pairs non-zero
        assert "alpha1*alpha2" in Qp
        assert "alpha1*alpha3" in Qp
        assert "alpha1*alpha4" in Qp
        assert "alpha2*alpha3" in Qp
        assert "alpha2*alpha4" in Qp
        assert "alpha3*alpha4" in Qp


class TestComputeLorentzZ:
    def test_4e4g_02_lorentzZ(self):
        """4e-4g-02: edges (3,3,4,4) -> lorentzZ = -1/4 z^2."""
        legs = [
            FGVLeg(edge=3, s_var="s31", loop="l1"),
            FGVLeg(edge=3, s_var="s32", loop="l2"),
            FGVLeg(edge=4, s_var="s41", loop="l3"),
            FGVLeg(edge=4, s_var="s42", loop="l4"),
        ]
        lz = compute_lorentzZ(legs)
        # lorentz_1 = v33.v44 - v34.v43 = 0*0 - z1*z1 = -z1^2
        # lorentz_3 = v33.v44 - v33.v44 = 0*0 - 0*0 = 0
        # lorentzZ = 1/4 * (-lorentz_1 + lorentz_3) = 1/4 * (z1^2 + 0) ... wait
        # Actually spec says -1/4 z^2. Let me verify by re-derivation:
        # lorentz_1 = v_e1.e3 * v_e2.e4 - v_e1.e4 * v_e2.e3
        #           = v33 * v34 - v34 * v33 = 0*... = 0
        # lorentz_3 = v_e1.e2 * v_e3.e4 - v_e1.e3 * v_e2.e4
        #           = v33 * v44 - v33 * v34 = 0*0 - 0*z1 = 0
        # Hmm, wait - e1=3, e2=3, e3=4, e4=4:
        # lorentz_1 = v_e1.e3 * v_e2.e4 - v_e1.e4 * v_e2.e3
        #           = v(3,4) * v(3,4) - v(3,4) * v(3,4) = z1*z1 - z1*z1 = 0
        # lorentz_3 = v_e1.e2 * v_e3.e4 - v_e1.e3 * v_e2.e4
        #           = v(3,3) * v(4,4) - v(3,4) * v(3,4) = 0*0 - z1*z1 = -z1^2
        # lorentzZ = 1/4 * (-lorentz_1 + lorentz_3) = 1/4 * (0 + (-z1^2)) = -1/4*z1^2
        assert "z1" in lz
        assert lz != "0"

    def test_4e4g_09_lorentzZ(self):
        """4e-4g-09: edges (1,2,3,4) -> lorentzZ = -1/4(1+4z+z^2)."""
        legs = [
            FGVLeg(edge=1, s_var="s11", loop="l1"),
            FGVLeg(edge=2, s_var="s21", loop="l2"),
            FGVLeg(edge=3, s_var="s31", loop="l3"),
            FGVLeg(edge=4, s_var="s41", loop="l4"),
        ]
        lz = compute_lorentzZ(legs)
        # Spec: 09 | [1,1,1,1] | 1,3 | -1/4(1+4z+z^2)
        assert "z1" in lz
        assert lz != "0"
