"""Prepare already signed native Montopar C/J without a legacy sign conversion.

This entrypoint is separate from historical self-energy-r002 replay and the
default 47 ordered roots. Internal flavour colours remain explicit root data.
"""
from copy import deepcopy
import re
import sympy as s

CONVENTION='wilson-qcd-action-path-fourier/1.0'

def _parse(text):
    text=re.sub(r'Gamma\[([^]]+)\]',r'gamma(\1)',text).replace('^','**')
    return s.sympify(text,locals={'I':s.I,'Pi':s.pi,'gamma':s.gamma})

def _gamma_euler_weight(expr):
    """Coefficient of EulerGamma*eps in integer-centred Gamma factors."""
    eps=s.Symbol('eps')
    if not expr.has(s.gamma):return s.S.Zero
    if expr.func==s.gamma:
        arg=expr.args[0];centre=arg.subs(eps,0);slope=s.diff(arg,eps)
        if centre.is_Integer is not True or centre<1 or slope.has(eps):
            raise ValueError('native_physical_unsupported_Gamma_argument')
        return -slope
    if expr.is_Mul:return sum(_gamma_euler_weight(a) for a in expr.args)
    if expr.is_Pow and expr.exp.is_Integer:
        return expr.exp*_gamma_euler_weight(expr.base)
    if expr.is_Add:
        weights={s.simplify(_gamma_euler_weight(a)) for a in expr.args}
        if len(weights)==1:return weights.pop()
    raise ValueError('native_physical_Gamma_weights_not_uniform')

def prepare_native_orbit(orbit,loop_order,*,connected_wilson_colour=None):
    if not isinstance(loop_order,int) or isinstance(loop_order,bool) or loop_order<1:
        raise ValueError('native_physical_invalid_loop_order')
    contract=orbit.get('physical_convention',{})
    if (contract.get('id')!=CONVENTION or contract.get('coordinate_Wick')!='-I'
        or contract.get('momentum_Wick')!='I' or contract.get('Dyson_external_coefficient')!=1
        or contract.get('raw_Euler_count')!=0):
        raise ValueError('native_physical_source_contract_required')
    normalization=orbit.get('normalization',{})
    if normalization.get('euler_compensation_count',0)!=0:
        raise ValueError('native_physical_already_prepared')
    if 'euler' in orbit['C']['components'] or 'euler_compensation' in orbit['C']['components']:
        raise ValueError('native_physical_duplicate_Euler')
    if any('EulerGamma' in str(v) for v in orbit['C']['components'].values()) or 'EulerGamma' in orbit['C']['combined']:
        raise ValueError('native_physical_duplicate_Euler')
    out=deepcopy(orbit)
    out['original_C']=deepcopy(out['C'])
    C=out['C']
    if C.get('component_semantics')!='multiplicative':
        raise ValueError('native_physical_requires_multiplicative_C')
    order=C['factor_order'];components=C['components']
    if len(order)!=len(set(order)) or set(order)!=set(components):
        raise ValueError('native_physical_C_factor_order_mismatch')
    parsed={k:_parse(v) for k,v in components.items()}
    if _parse(C['combined'])!=s.prod(parsed[k] for k in order):
        raise ValueError('native_physical_C_components_mismatch')
    weight=s.simplify(sum(_gamma_euler_weight(v) for v in parsed.values()))
    if weight!=loop_order:
        raise ValueError('native_physical_Gamma_weight_loop_mismatch')
    if connected_wilson_colour is not None:
        if 'wilson_color' not in C['components']:raise ValueError('native_physical_colour_slot_missing')
        C['components']['wilson_color']=connected_wilson_colour
    C['components']['euler']=f'exp(-{loop_order}*EulerGamma*eps)'
    C['factor_order'].append('euler')
    C['combined']='*'.join(f'({C["components"][k]})' for k in C['factor_order'])
    out.setdefault('normalization',{}).update(euler_compensation_count=1,assembly_sign=1,
        scheme='native-physical-CJ',physical_convention=CONVENTION)
    out['preparation']={'native_sign_conversion_count':0,'Euler_compensation_count':1,
        'J_changed':False,'colour_conversion':connected_wilson_colour is not None,
        'original_C_preserved':True,'Gamma_Euler_weight':str(weight),'formal':'not_run'}
    return out

def prepare_native_package(package,*,connected_wilson_colour=None):
    out=deepcopy(package)
    n=package['loop_order']
    out['orbits']=[prepare_native_orbit(o,n,connected_wilson_colour=connected_wilson_colour)
                   for o in package['orbits']]
    out['format']='wilson-native-physical-calculation-input/1'
    out['physical_convention']=CONVENTION
    out['default_47_inputs_modified']=False
    return out
