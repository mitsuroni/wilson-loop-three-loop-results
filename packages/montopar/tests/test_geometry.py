"""Tests for N=4 geometry utilities."""
from montopar.geometry import (
    dot_product_factor,
    v_dot_V_diff,
    compute_Delta_prime,
)


class TestDotProductFactor:
    def test_adjacent_12(self):
        assert dot_product_factor(1, 2) == "z1"

    def test_adjacent_34(self):
        assert dot_product_factor(3, 4) == "z1"

    def test_adjacent_23(self):
        assert dot_product_factor(2, 3) == "1"

    def test_adjacent_41(self):
        assert dot_product_factor(4, 1) == "1"

    def test_opposite_13(self):
        assert dot_product_factor(1, 3) == "-(1+z1)"

    def test_opposite_24(self):
        assert dot_product_factor(2, 4) == "-(1+z1)"

    def test_symmetric(self):
        assert dot_product_factor(2, 1) == dot_product_factor(1, 2)
        assert dot_product_factor(4, 3) == dot_product_factor(3, 4)


class TestVDotVDiff:
    def test_same_point(self):
        assert v_dot_V_diff(1, 3, 3) == "0"

    def test_V3_minus_V1(self):
        # V_3 - V_1 = v_1 + v_2, v_1·(V_3-V_1) = v_1·v_1 + v_1·v_2 = 0 + z1 = z1
        assert v_dot_V_diff(1, 3, 1) == "z1"

    def test_V4_minus_V2(self):
        # V_4 - V_2 = v_2 + v_3, v_2·(V_4-V_2) = v_2·v_2 + v_2·v_3 = 0 + 1 = 1
        assert v_dot_V_diff(2, 4, 2) == "1"

    def test_V2_minus_V4(self):
        # V_2 - V_4 = -(v_2 + v_3), v_3·(V_2-V_4) = -v_3·v_2 - v_3·v_3 = -1 - 0 = -1
        assert v_dot_V_diff(3, 2, 4) == "-1"

    def test_V1_minus_V3(self):
        # V_1 - V_3 = -(v_1 + v_2), v_4·(V_1-V_3) = -v_4·v_1 - v_4·v_2 = -1 - (-(1+z1)) = -(1+z1) - 1... wait
        # v_4·v_1 = v_14 / v_ref = 1, v_4·v_2 = v_24 / v_ref = -(1+z1)
        # So v_4·(V_1-V_3) = -(1 + (-(1+z1))) = -(1 - 1 - z1) = z1
        assert v_dot_V_diff(4, 1, 3) == "z1"


class TestComputeDeltaPrime:
    def test_same_edge(self):
        assert compute_Delta_prime(3, 3, "s31", "s32") == "0"

    def test_adjacent_34(self):
        result = compute_Delta_prime(3, 4, "s31", "s41")
        assert result == "(1-s31)*s41*z1"

    def test_adjacent_43(self):
        result = compute_Delta_prime(4, 3, "s41", "s31")
        assert result == "(1-s31)*s41*z1"

    def test_opposite_13(self):
        result = compute_Delta_prime(1, 3, "s1", "s3")
        assert result == "z1*(1-s1)*(1-s3) + s1*s3"

    def test_opposite_24(self):
        result = compute_Delta_prime(2, 4, "s2", "s4")
        assert result == "(1-s2)*(1-s4) + z1*s2*s4"
