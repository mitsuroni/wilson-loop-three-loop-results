import pytest
from montopar.parser import parse_orbit

@pytest.mark.parametrize('raw',[
 {'emission_type':'3e-Gamma1','vertices':[],'propagators':[],'internal_vertices':[],'phases':[],'bounds':[]},
 {'emission_type':'3e','internal_vertices':[{'type':'gamma1'}]},
 {'emission_type':'3e','vertex_correction':{'operator':'Gamma1'}},
])
def test_central_vertex_cannot_fall_through_to_pure(raw):
    with pytest.raises(ValueError,match='gamma1_'):
        parse_orbit(raw)
