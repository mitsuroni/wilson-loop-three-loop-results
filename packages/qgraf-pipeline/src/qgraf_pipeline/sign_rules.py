"""Declare the action/path/Fourier meaning of source integrals without family signs."""
from collections import defaultdict
CONVENTION='wilson-qcd-action-path-fourier/1.0'

def attach_contract(raw):
    groups=defaultdict(list)
    for v in raw['vertices']:groups[v['line']].append(v['color'])
    chronological=[c for e in sorted(groups) for c in groups[e]]
    if len(set(chronological))!=len(chronological):
        raise ValueError('shared_signs_duplicate_Wilson_colour_slot')
    raw['physical_convention']={'id':CONVENTION,'wilson_transport':"U'=ig*A*U",
        'wilson_word':list(reversed(chronological)),'chronological_word':chronological,
        'vertex_labels_unchanged':True,'momentum_Fourier':'exp(-ipX)',
        'momentum_Wick':'I','coordinate_Wick':'-I','one_pi_includes_internal_loop_signs':True,
        'Dyson_external_coefficient':1,'raw_Euler_count':0}
    return raw
