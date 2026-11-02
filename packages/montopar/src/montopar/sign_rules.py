"""Action-derived Wilson ordering and exact i phases, independent of graph names."""
from collections import defaultdict
from copy import deepcopy
import re

CONVENTION = 'wilson-qcd-action-path-fourier/1.0'

def transport_word(vertices):
    groups=defaultdict(list)
    for vertex in vertices:
        edge=vertex['line'] if isinstance(vertex,dict) else vertex.line
        colour=vertex['color'] if isinstance(vertex,dict) else vertex.color
        groups[edge].append(colour)
    chronological=[colour for edge in sorted(groups) for colour in groups[edge]]
    if len(set(chronological))!=len(chronological):
        raise ValueError('shared_signs_duplicate_Wilson_colour_slot')
    return list(reversed(chronological))

def expected_contract(vertices):
    word=transport_word(vertices)
    return {'id':CONVENTION,'wilson_transport':"U'=ig*A*U",'wilson_word':word,
            'chronological_word':list(reversed(word)),'vertex_labels_unchanged':True,
            'momentum_Fourier':'exp(-ipX)','momentum_Wick':'I',
            'coordinate_Wick':'-I','one_pi_includes_internal_loop_signs':True,
            'Dyson_external_coefficient':1,'raw_Euler_count':0}

def expected_source_contract(raw):
    """Storage order is not contour order: recover slots from labelled endpoint phases."""
    if not raw.get('phases') or not raw.get('propagators'):
        return expected_contract(raw['vertices'])
    phases={p['loop']:p['expression'] for p in raw['phases']}
    positions=[]
    for v in raw['vertices']:
        matches=[]
        for prop in raw['propagators']:
            if v['lorentz'] in prop['lorentz']:
                for edge,slot in re.findall(r's([1-4])(\d+)',phases[prop['loop']]):
                    if int(edge)==v['line']:matches.append((int(edge),int(slot)))
        unique=set(matches)
        if len(unique)!=1:
            # Legacy topology fragments with ambiguous same-edge routing are not native CLI inputs.
            return expected_contract(raw['vertices'])
        positions.append((next(iter(unique)),v))
    chronological=[v for _,v in sorted(positions,key=lambda x:x[0])]
    return expected_contract(chronological)

def validate_contract(raw,required=False):
    actual=raw.get('physical_convention') if isinstance(raw,dict) else None
    if actual is None:
        if required:raise ValueError('shared_signs_source_requires_physical_convention')
        return expected_source_contract(raw) if isinstance(raw,dict) and 'vertices' in raw else None
    try:expected=expected_source_contract(raw)
    except (KeyError,TypeError,AttributeError):raise ValueError('shared_signs_invalid_source') from None
    if actual!=expected:
        raise ValueError('shared_signs_source_contract_mismatch')
    return expected

def phase(*exponents):
    """Exact power of i; -i has exponent 3 and -1 has exponent 2."""
    return ('1','I','-1','-I')[sum(exponents)%4]

def real_phase(*exponents):
    value=phase(*exponents)
    if value not in ('1','-1'):raise ValueError('shared_signs_expected_real_phase')
    return int(value)

def coordinate_phase(endpoints,propagators,centres,vertex4=0,derivatives=0,chain_negatives=0):
    return phase(endpoints,2*propagators,3*centres,3*vertex4,3*derivatives,2*chain_negatives)

def three_slot_trace_phase(raw):
    """Contract the actual transport word against a fixed ordered f tensor."""
    if len(raw['vertices'])!=3:raise ValueError('shared_signs_requires_three_slots')
    word=validate_contract(raw)['wilson_word']
    positions=[word.index(v['color']) for v in raw['vertices']]
    odd=sum(positions[i]>positions[j] for i in range(3) for j in range(i+1,3))%2
    return -1 if odd else 1

def annotate_result(result,raw):
    out=deepcopy(result)
    if 'vertices' not in raw:
        return out  # Low-level serialization unit fixtures are not native source integrals.
    out['physical_convention']=validate_contract(raw)
    out.setdefault('normalization',{})['physical_convention']=CONVENTION
    out['normalization']['assembly_sign']=1
    if 'physical_convention' not in raw:
        out['source_interpretation']='legacy topology only; physical integral rebuilt from the declared common action, not replay of a legacy amplitude'
    return out
