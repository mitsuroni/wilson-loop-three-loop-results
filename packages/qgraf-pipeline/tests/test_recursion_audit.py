import json
from pathlib import Path
import pytest

from qgraf_pipeline.classifier import classify_orbits
from qgraf_pipeline.definition_exporter import _build_orbit_definition
from qgraf_pipeline.generator import generate_diagrams
from qgraf_pipeline.internal_types import enumerate_emission_types
from qgraf_pipeline.json_emitter import _orbit_to_dict
from qgraf_pipeline.recursion_audit import build_audit, write_audit


@pytest.mark.parametrize('n', [1, 2, 3])
def test_production_objects_unchanged(n):
    types = enumerate_emission_types(n)
    diagrams = generate_diagrams(n, types)
    standard = classify_orbits(diagrams, types)
    audit = classify_orbits(diagrams, types, include_reducible=True)
    retained = [o for o in audit if o.is_2pi]
    assert [_orbit_to_dict(o) for o in retained] == [_orbit_to_dict(o) for o in standard]
    assert [_build_orbit_definition(o) for o in retained] == [_build_orbit_definition(o) for o in standard]


def test_stable_auxiliary_ids_when_input_order_changes():
    types = enumerate_emission_types(3)
    diagrams = generate_diagrams(3, types)
    forward = classify_orbits(diagrams, types, include_reducible=True)
    backward = classify_orbits(list(reversed(diagrams)), types, include_reducible=True)
    assert {o.orbit_id for o in forward} == {o.orbit_id for o in backward}
    assert all('-aux-' in o.orbit_id for o in forward if not o.is_2pi)


def test_sextuple_and_scaleless_records():
    audit, definitions = build_audit(3)
    rows = [r for r in audit['records'] if r['emission_type']=='5e-3g'
            and tuple(r['edge_counts'])==(0,0,2,3)
            and sorted([r['wl_wl_pairs'][0][0],r['wl_wl_pairs'][0][2]])==[3,4]]
    assert len(rows)==6
    assert sum(r['status']=='production' for r in rows)==4
    assert sum(r['status']=='auxiliary' for r in rows)==2
    assert sum(r['status']=='scaleless' for r in audit['records'])==60
    assert len({o['orbit_id'] for o in definitions['orbits']})==148


def test_separate_outputs_and_invalid_loop_order(tmp_path):
    sentinel = tmp_path/'orbit_definitions.json'
    sentinel.write_bytes(b'production fixture')
    write_audit(2,tmp_path)
    assert sentinel.read_bytes()==b'production fixture'
    assert json.loads((tmp_path/'recursion-orbit-definitions.json').read_text())['n_orbits']==16
    with pytest.raises(ValueError):
        build_audit(4)


@pytest.mark.parametrize('tail', [['4','--audit-recursion'],
                                ['3','--audit-recursion','--out-file','existing.out']])
def test_cli_rejects_unsupported_audit_inputs(monkeypatch,tail):
    from qgraf_pipeline.cli import main
    monkeypatch.setattr('sys.argv',['qgraf-pipeline',*tail])
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code==2
