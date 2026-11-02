"""Tests for montopar CLI dispatch routing."""
import json
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

from montopar.cli import main
from montopar.models import (
    Orbit, EdgeDistribution, Bound, Vertex, Propagator,
    Pairing, CFactors, JFactors, HypercubeResult,
    TGVJFactors, FGVJFactors,
)
from montopar.parser import parse_all_orbits


FIXTURES = Path(__file__).parent / "fixtures"


def _make_tgv_orbit_json(orbit_id: str = "3e-3g-01") -> dict:
    """Minimal TGV orbit definition for testing."""
    return {
        "orbit_id": orbit_id,
        "emission_type": "3e-3g",
        "sector": "dg",
        "nloops": 2,
        "vertices": [
            {"line": 3, "lorentz": "mu1", "color": "a1"},
            {"line": 4, "lorentz": "mu2", "color": "a2"},
            {"line": 4, "lorentz": "mu3", "color": "a3"},
        ],
        "propagators": [
            {"loop": "l1", "lorentz": ["mu1", "nu1"], "color": ["a1", "c1"]},
            {"loop": "l2", "lorentz": ["mu2", "nu1"], "color": ["a2", "c1"]},
            {"loop": "l3", "lorentz": ["mu3", "nu1"], "color": ["a3", "c1"]},
        ],
        "internal_vertices": [
            {"type": "tgv", "loops": ["l1", "l2", "l3"], "lorentz": "nu1"},
        ],
        "phases": [
            {"loop": "l1", "expression": "-((1-s31)*v3)"},
            {"loop": "l2", "expression": "-((1-s41)*v4)"},
            {"loop": "l3", "expression": "-((1-s42)*v4)"},
        ],
        "bounds": [
            {"var": "s31", "low": "0", "high": "1"},
            {"var": "s41", "low": "0", "high": "s42"},
            {"var": "s42", "low": "0", "high": "1"},
        ],
    }


def _make_fgv_orbit_json(orbit_id: str = "4e-4g-02") -> dict:
    """Minimal FGV orbit definition for testing."""
    return {
        "orbit_id": orbit_id,
        "emission_type": "4e-4g",
        "sector": "dg",
        "nloops": 3,
        "vertices": [
            {"line": 3, "lorentz": "mu1", "color": "a1"},
            {"line": 3, "lorentz": "mu2", "color": "a2"},
            {"line": 4, "lorentz": "mu3", "color": "a3"},
            {"line": 4, "lorentz": "mu4", "color": "a4"},
        ],
        "propagators": [
            {"loop": "l1", "lorentz": ["mu1", "nu1"], "color": ["a1", "c1"]},
            {"loop": "l2", "lorentz": ["mu2", "nu2"], "color": ["a2", "c2"]},
            {"loop": "l3", "lorentz": ["mu3", "nu1"], "color": ["a3", "c1"]},
            {"loop": "l4", "lorentz": ["mu4", "nu2"], "color": ["a4", "c2"]},
        ],
        "internal_vertices": [
            {"type": "4v", "loops": ["l1", "l2", "l3", "l4"],
             "lorentz": ["nu1", "nu2"]},
        ],
        "phases": [
            {"loop": "l1", "expression": "-((1-s31)*v3)"},
            {"loop": "l2", "expression": "-((1-s32)*v3)"},
            {"loop": "l3", "expression": "-((1-s41)*v4)"},
            {"loop": "l4", "expression": "-((1-s42)*v4)"},
        ],
        "bounds": [
            {"var": "s31", "low": "0", "high": "s32"},
            {"var": "s32", "low": "0", "high": "1"},
            {"var": "s41", "low": "0", "high": "s42"},
            {"var": "s42", "low": "0", "high": "1"},
        ],
    }


def _make_source_json(orbits: list[dict]) -> dict:
    return {"orbits": orbits}


class TestCLIDispatchRouting:
    """Test that the CLI dispatches correctly to TGV/FGV/pure paths."""

    def test_parse_tgv_dispatch(self):
        """TGV orbit should have dispatch='tgv'."""
        orbit = parse_all_orbits(_make_source_json([_make_tgv_orbit_json()]))[0]
        assert orbit.dispatch == "tgv"
        assert len(orbit.tgv_legs) == 3

    def test_parse_fgv_dispatch(self):
        """FGV orbit should have dispatch='fourv'."""
        orbit = parse_all_orbits(_make_source_json([_make_fgv_orbit_json()]))[0]
        assert orbit.dispatch == "fourv"
        assert len(orbit.fgv_legs) == 4

    def test_type_filter_tgv(self, tmp_path):
        """--type 3e-3g should filter to TGV orbits."""
        source = _make_source_json([
            _make_tgv_orbit_json("3e-3g-01"),
            _make_tgv_orbit_json("3e-3g-02"),
        ])
        source_file = tmp_path / "source.json"
        source_file.write_text(json.dumps(source))

        out_dir = tmp_path / "out"

        with patch("sys.argv", [
            "montopar", "--source", str(source_file),
            "--type", "3e-3g", "--out-dir", str(out_dir),
        ]):
            with patch("montopar.cli.compute_color", return_value="CA*CF"):
                main()

        # Should have created output for both orbits
        assert (out_dir / "2loop" / "3e-3g" / "3e-3g-01" / "montopar.json").exists()
        assert (out_dir / "2loop" / "3e-3g" / "3e-3g-02" / "montopar.json").exists()

    def test_type_filter_fgv(self, tmp_path):
        """--type 4e-4g should filter to FGV orbits."""
        source = _make_source_json([
            _make_fgv_orbit_json("4e-4g-02"),
        ])
        source_file = tmp_path / "source.json"
        source_file.write_text(json.dumps(source))

        out_dir = tmp_path / "out"

        with patch("sys.argv", [
            "montopar", "--source", str(source_file),
            "--type", "4e-4g", "--out-dir", str(out_dir),
        ]):
            with patch("montopar.cli.compute_color", return_value="CA^2*CF"):
                main()

        # Should have created output
        assert (out_dir / "3loop" / "4e-4g" / "4e-4g-02" / "montopar.json").exists()

        # Verify dispatch type in output
        with open(out_dir / "3loop" / "4e-4g" / "4e-4g-02" / "montopar.json") as f:
            data = json.load(f)
        assert data["dispatch"] == "fourv"

    def test_type_filter_no_match(self, tmp_path, capsys):
        """--type with no matching orbits should print error."""
        source = _make_source_json([
            _make_tgv_orbit_json("3e-3g-01"),
        ])
        source_file = tmp_path / "source.json"
        source_file.write_text(json.dumps(source))

        with patch("sys.argv", [
            "montopar", "--source", str(source_file),
            "--type", "99e-99g",
        ]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 1

        captured = capsys.readouterr()
        assert "No orbits of type 99e-99g found" in captured.out

    def test_all_processes_all_types(self, tmp_path):
        """--all should process both TGV and FGV orbits."""
        source = _make_source_json([
            _make_tgv_orbit_json("3e-3g-01"),
            _make_fgv_orbit_json("4e-4g-02"),
        ])
        source_file = tmp_path / "source.json"
        source_file.write_text(json.dumps(source))

        out_dir = tmp_path / "out"

        with patch("sys.argv", [
            "montopar", "--source", str(source_file),
            "--all", "--out-dir", str(out_dir),
        ]):
            with patch("montopar.cli.compute_color", return_value="CA*CF"):
                main()

        assert (out_dir / "2loop" / "3e-3g" / "3e-3g-01" / "montopar.json").exists()
        assert (out_dir / "3loop" / "4e-4g" / "4e-4g-02" / "montopar.json").exists()

    def test_tgv_output_has_pnum(self, tmp_path):
        """TGV output should contain pnum and S_polynomial."""
        source = _make_source_json([_make_tgv_orbit_json()])
        source_file = tmp_path / "source.json"
        source_file.write_text(json.dumps(source))

        out_dir = tmp_path / "out"

        with patch("sys.argv", [
            "montopar", "--source", str(source_file),
            "--type", "3e-3g", "--out-dir", str(out_dir),
        ]):
            with patch("montopar.cli.compute_color", return_value="I/2*CA*CF"):
                main()

        with open(out_dir / "2loop" / "3e-3g" / "3e-3g-01" / "montopar.json") as f:
            data = json.load(f)

        assert data["dispatch"] == "tgv"
        assert data["J_raw"]["integrand"] is not None
        # TGV integrand contains alpha_monomial and S_polynomial
        assert "alpha" in data["J_raw"]["integrand"]
        assert data["C"]["components"]["derivative"] is not None

    def test_tgv_ibp_flag_adds_single_canonical_representative(self, tmp_path):
        source_file = tmp_path / "source.json"
        source_file.write_text(json.dumps(_make_source_json([_make_tgv_orbit_json()])))
        out_dir = tmp_path / "out"

        with patch("sys.argv", [
            "montopar", "--source", str(source_file), "--type", "3e-3g",
            "--out-dir", str(out_dir), "--tgv-ibp",
        ]):
            main()

        output = out_dir / "2loop" / "3e-3g" / "3e-3g-01" / "montopar.json"
        data = json.loads(output.read_text())
        assert data["J_ibp"]["status"] == "complete"
        assert not any(k.endswith("_rotated") for k in data)
        assert data["J_raw_s"]["alpha_measure"]["kind"] == "simplex"
        assert len(data["J_ibp"]["terms"]) == 3
        assert "J_raw" in data

    def test_fgv_output_has_lorentzZ(self, tmp_path):
        """FGV output should contain lorentzZ and Q_prime_polynomial."""
        source = _make_source_json([_make_fgv_orbit_json()])
        source_file = tmp_path / "source.json"
        source_file.write_text(json.dumps(source))

        out_dir = tmp_path / "out"

        with patch("sys.argv", [
            "montopar", "--source", str(source_file),
            "--type", "4e-4g", "--out-dir", str(out_dir),
        ]):
            with patch("montopar.cli.compute_color", return_value="CA^2*CF"):
                main()

        with open(out_dir / "3loop" / "4e-4g" / "4e-4g-02" / "montopar.json") as f:
            data = json.load(f)

        assert data["dispatch"] == "fourv"
        # FGV integrand is non-trivial and contains alpha variables (from
        # alpha_monomial / Q_prime_polynomial) and z1 (from lorentzZ).
        assert data["J_raw"]["integrand"] != "0"
        assert "alpha" in data["J_raw"]["integrand"]
        assert "z1" in data["J_raw"]["integrand"]

    def test_orbit_id_filter(self, tmp_path):
        """--orbit-id should filter to a specific orbit."""
        source = _make_source_json([
            _make_tgv_orbit_json("3e-3g-01"),
            _make_tgv_orbit_json("3e-3g-02"),
        ])
        source_file = tmp_path / "source.json"
        source_file.write_text(json.dumps(source))

        out_dir = tmp_path / "out"

        with patch("sys.argv", [
            "montopar", "--source", str(source_file),
            "--orbit-id", "3e-3g-01", "--out-dir", str(out_dir),
        ]):
            with patch("montopar.cli.compute_color", return_value="I/2*CA*CF"):
                main()

        assert (out_dir / "2loop" / "3e-3g" / "3e-3g-01" / "montopar.json").exists()
        assert not (out_dir / "2loop" / "3e-3g" / "3e-3g-02" / "montopar.json").exists()

    def test_pure_still_works(self, tmp_path):
        """--type 4e should still work for pure n-gluon orbits."""
        fixture = FIXTURES / "4e-01_orbit.json"
        with open(fixture) as f:
            orbit_data = json.load(f)

        source = _make_source_json([orbit_data])
        source_file = tmp_path / "source.json"
        source_file.write_text(json.dumps(source))

        out_dir = tmp_path / "out"

        with patch("sys.argv", [
            "montopar", "--source", str(source_file),
            "--type", "4e", "--out-dir", str(out_dir),
        ]):
            with patch("montopar.cli.compute_color", return_value="-1/2*((CA - 2*CF)*CF)"):
                main()

        assert (out_dir / "2loop" / "4e" / "4e-01" / "montopar.json").exists()

        with open(out_dir / "2loop" / "4e" / "4e-01" / "montopar.json") as f:
            data = json.load(f)
        assert data["dispatch"] == "pure"


def _declare_native_input(factory):
    from montopar.sign_rules import expected_contract
    def wrapped(*args, **kwargs):
        raw=factory(*args, **kwargs)
        raw['physical_convention']=expected_contract(raw['vertices'])
        return raw
    return wrapped

_make_tgv_orbit_json=_declare_native_input(_make_tgv_orbit_json)
_make_fgv_orbit_json=_declare_native_input(_make_fgv_orbit_json)
