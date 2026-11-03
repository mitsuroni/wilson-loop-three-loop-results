"""Portable T4 reassembly from frozen original primitive/IBP outputs.

The eight empty orbits are assembled afresh from the raw Q jets, including an
exact rational inversion, absolute Laurent convolution and Gamma/cosh jet.
The opposite 27/28 hard integrations are explicitly reused from their original
CAS outputs; no selected result record is read. All individual wide coordinates
are kept. No Maple writer, CAS process, old selector, or legacy sign is invoked.
"""
from __future__ import annotations
import ast
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re
import sympy as sp
import t4_endpoints as ep

PACK = Path(__file__).resolve().parents[1]
SNAP = PACK / 'snapshots'
EMPTY_RUN = 'workflow/checks/tgv-t4-empty/runs/20260918T023134199721Z-assemble_13_14_15_40'
OPPOSITE_RUN = 'workflow/checks/tgv-t4-opposite/runs/20260918T034029006332Z-global_ibp'
OMEGA_RUN = 'workflow/checks/tgv-t4-opposite/runs/20260918T034133130688Z-omega_basis'
MOTHER_RUN = 'workflow/checks/single-opposite/runs/20260916T090538777022Z'
PERIODS = ('zeta_2', 'zeta_3', 'zeta_5', 'log_2', 'A4', 'A6')
WEIGHTS = (2, 3, 5, 1, 4, 6)
ZERO = (0,) * 6
ONE = {((), ZERO): Q(1)}
LOG = {((Q(0),), ZERO): Q(1)}
EMPTY_MAP = {'13': '13', '14': '14', '15': '15', '33': '15',
             '34': '14', '35': '13', '40': '40', '43': '40'}

def clean(v):
    return {key: Q(c) for key, c in v.items() if c}

def add(*vs):
    out = Counter()
    for v in vs:
        for key, c in v.items():
            out[key] += c
    return clean(out)

def scale(v, c):
    return clean({key: value * Q(c) for key, value in v.items()})

def mul(a, b):
    out = Counter()
    for (wa, pa), ca in a.items():
        for (wb, pb), cb in b.items():
            periods = tuple(x + y for x, y in zip(pa, pb))
            for w, n in ep.shuffles(wa, wb).items():
                out[w, periods] += ca * cb * n
    return clean(out)

def power(v, n):
    assert type(n) is int and n >= 0
    ans = ONE
    for _ in range(n):
        ans = mul(ans, v)
    return ans

def period(n):
    if n in (2, 3, 5):
        p = list(ZERO)
        p[(2, 3, 5).index(n)] = 1
        return {((), tuple(p)): Q(1)}
    if n == 4:
        return scale(power(period(2), 2), Q(2, 5))
    if n == 6:
        return scale(power(period(2), 3), Q(8, 35))
    raise ValueError(n)

def parse_text(text):
    """Parse the actual Maple Hlog output, retaining all six period powers."""
    text = text.strip().rstrip(':;')
    words = []
    def replace_hlog(match):
        w = tuple(Q(a.strip()) for a in match[1].split(',') if a.strip())
        assert set(w) <= {0, -1, Q(-1, 2), -2}, w
        words.append(w)
        return 'g' + str(len(words) - 1)
    text = re.sub(r'Hlog\(z,\[([^]]*)\]\)', replace_hlog, text)
    text = re.sub(r'zeta\[1,\s*-3\]', 'A4', text)
    text = re.sub(r'zeta\[1,\s*-5\]', 'A6', text).replace('ln(2)', 'L2')
    text = re.sub(r'zeta\[(\d+)\]', r'Z\1', text).replace('^', '**')
    names = {'g' + str(n): {(w, ZERO): Q(1)} for n, w in enumerate(words)}
    names.update({'Z' + str(n): period(n) for n in (2, 3, 4, 5, 6)})
    for name, pos in [('L2', 3), ('A4', 4), ('A6', 5)]:
        ps = list(ZERO)
        ps[pos] = 1
        names[name] = {((), tuple(ps)): Q(1)}
    def read(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return scale(ONE, node.value)
        if isinstance(node, ast.Name) and node.id in names:
            return names[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return scale(read(node.operand), -1)
        if isinstance(node, ast.BinOp):
            a, b = read(node.left), read(node.right)
            if isinstance(node.op, ast.Add):
                return add(a, b)
            if isinstance(node.op, ast.Sub):
                return add(a, scale(b, -1))
            if isinstance(node.op, ast.Mult):
                return mul(a, b)
            if isinstance(node.op, ast.Div):
                assert set(b) == set(ONE)
                return scale(a, 1 / b[(), ZERO])
            if isinstance(node.op, ast.Pow):
                assert set(b) <= set(ONE)
                exponent = b.get(((), ZERO), Q(0))
                assert exponent.denominator == 1
                return power(a, int(exponent))
        raise ValueError('Unsupported original Hlog expression: ' + ast.dump(node)[:160])
    return read(ast.parse(text, mode='eval').body)

def parse(path):
    return parse_text(Path(path).read_text(encoding='utf-8').split('=', 1)[1])

def vec(rows):
    return {(tuple(Q(a) for a in row['word']),
             tuple(row.get('period_powers', row.get('zeta_powers', []))) +
             ((0, 0, 0) if 'zeta_powers' in row else ())): Q(row['coefficient'])
            for row in rows}

def serialize(v):
    return [{'word': [int(a) if a.denominator == 1 else str(a) for a in w],
             'period_powers': list(p), 'coefficient': str(c)}
            for (w, p), c in sorted(v.items())]

def weight(key):
    return len(key[0]) + sum(n * p for n, p in zip(WEIGHTS, key[1]))

def _constant(expr):
    expr = sp.Poly(sp.expand(expr), ep.Z2, ep.Z3, ep.Z5)
    return {((), tuple(p) + (0, 0, 0)): Q(int(c.p), int(c.q))
            for p, c in expr.terms() if c}

@lru_cache(None)
def _inverse_word(w):
    if not w:
        return ONE
    assert set(w) <= {0, -1}
    tail = _inverse_word(w[1:])
    # dG(0,v;1/z)/dz=-G(v;1/z)/z;
    # dG(-1,v;1/z)/dz=(-1/z+1/(1+z))*G(v;1/z).
    primitive = {((Q(0),) + v, ps): -c for (v, ps), c in tail.items()}
    if w[0] == -1:
        primitive = add(primitive, {((Q(-1),) + v, ps): c for (v, ps), c in tail.items()})
    return add(primitive, _constant(ep.reginf_word(w)))

def invert_01(v):
    out = {}
    for (w, p), c in v.items():
        for (iw, ip), ic in _inverse_word(w).items():
            key = iw, tuple(x + y for x, y in zip(p, ip))
            out[key] = out.get(key, Q(0)) + c * ic
    return clean(out)

def _times_log(v, n, c=1):
    return scale(mul(power(LOG, n), v), c)

def empty_from_raw_q():
    """Return fresh physical I jets; absolute indices are never relative orders."""
    pref = {
        0: ONE, 1: {},
        2: add(scale(period(2), Q(5, 2)), _times_log(ONE, 2, Q(-9, 8))),
        3: scale(period(3), 3),
        4: add(scale(power(period(2), 2), Q(193, 40)),
               _times_log(period(2), 2, Q(-45, 16)),
               _times_log(ONE, 4, Q(135, 128))),
    }
    result, raw_q = {}, {}
    for oid in ('13', '14', '15', '40'):
        q = {k: parse(SNAP / EMPTY_RUN / f'Q{oid}q{k}.txt') for k in range(-2, 3)}
        # This involution is independent of stored physical I coefficients.
        assert all(invert_01(invert_01(v)) == v for v in q.values()), oid
        raw_q[oid] = q
        even = {}
        for k in range(-2, 3):
            even[k] = add(*[
                add(_times_log(q[k-j], j, Q((-1)**j, 2**(j+1) * sp.factorial(j))),
                    _times_log(invert_01(q[k-j]), j, Q(1, 2**(j+1) * sp.factorial(j))))
                for j in range(k + 3)])
        result[oid] = {k: ({} if k < -4 else scale(add(*[
            mul(pref[j], even[k+2-j]) for j in range(k + 5)]), Q(1, 32)))
            for k in range(-6, 1)}
        assert all(invert_01(v) == v for v in result[oid].values()), oid
    return result, raw_q

def factorized_41():
    """Source-factorized T4 x H0 with its own reflection and orbit weight."""
    m0, m1, m2 = [parse(SNAP / MOTHER_RUN / f'M_{q}_1.txt') for q in range(3)]
    bare = {-2: scale(period(2), Q(1, 8)),
            -1: scale(add(period(3), scale(m1, -1)), Q(1, 8)),
             0: add(scale(m2, Q(-1, 8)), scale(mul(m0, m0), Q(1, 48)),
                    scale(power(period(2), 2), Q(9, 10)))}
    h0 = {0: m0, 1: scale(m1, Q(1, 2)),
          2: add(scale(m2, Q(1, 8)), scale(mul(m0, m0), Q(-1, 24)))}
    pref = {0: ONE, 2: add(scale(m0, Q(5, 4)), scale(period(2), Q(17, 4)))}
    return {k: scale(add(*[mul(mul(a, b), c)
                    for ia, a in bare.items() for ib, b in h0.items()
                    for ic, c in pref.items() if ia + ib + ic == k]), Q(1, 8))
            for k in range(-6, 1)}

def omega_data():
    omega = parse(SNAP / OMEGA_RUN / 'Omega6.txt')
    phi = parse(SNAP / OMEGA_RUN / 'Phi5.txt')
    primitive = parse(SNAP / OMEGA_RUN / 'Tprimitive.txt')
    assert all(w and w[0] == Q(-1, 2) for w, _ in primitive)
    assert {(w[1:], p): c for (w, p), c in primitive.items()} == phi
    assert all(weight(key) == 5 for key in phi)
    assert all(weight(key) == 6 for key in omega)
    assert omega.get(((Q(-1, 2), 0, 0, 0, 0, 0), ZERO)) == 8
    return omega, phi, primitive

def _source_contract(ids):
    native = json.loads((SNAP / 'benchmarks/current/prepared/3loop/5e-3g.json').read_text(encoding='utf-8'))
    historical = json.loads((SNAP / 'benchmarks/tgv/5e-3g/prepared/5e-3g-r001.json').read_text(encoding='utf-8'))
    old = {r['orbit_id']: r for r in historical['orbits']}
    rows = {r['orbit_id']: r for r in native['orbits']}
    for oid in ids:
        r = rows[oid]
        assert r['J_raw_s'] == old[oid]['J_raw_s'], oid
        assert r['C']['component_semantics'] == 'multiplicative', oid
        assert r['C']['components']['euler'] == 'exp(-3*EulerGamma*eps)', oid
        assert r['normalization']['assembly_sign'] == 1, oid
        assert r['normalization']['euler_compensation_count'] == 1, oid
        assert r['c4_symmetry']['orbit_multiplicity'] == 4, oid
    return rows

def generate():
    ids = ['5e-3g-' + oid for oid in (*EMPTY_MAP, '27', '28', '29', '41', '42')]
    _source_contract(ids)
    assert 'ASSEMBLY_COMPLETE' in (SNAP / EMPTY_RUN / 'maple.log').read_text(encoding='utf-8')
    assert 'GLOBAL_COMPLETE' in (SNAP / OPPOSITE_RUN / 'maple.log').read_text(encoding='utf-8')
    assert 'OMEGA_COMPLETE' in (SNAP / OMEGA_RUN / 'maple.log').read_text(encoding='utf-8')
    empty, raw_q = empty_from_raw_q()
    omega, phi, primitive = omega_data()
    vectors = {oid: empty[representative] for oid, representative in EMPTY_MAP.items()}
    for oid in ('27', '28'):
        vectors[oid] = {k: {} for k in range(-6, 1)}
        vectors[oid][-1] = parse(SNAP / OPPOSITE_RUN / f'I{oid}m1.txt')
        vectors[oid][0] = parse(SNAP / OPPOSITE_RUN / f'I{oid}zero.txt')
    vectors['29'] = vectors['27']
    vectors['41'] = factorized_41()
    vectors['42'] = vectors['41']
    omega_coefficient = {'27': Q(-1, 64), '28': Q(1, 32), '29': Q(-1, 64), '41': Q(0), '42': Q(0)}
    cancellation = add(scale(vectors['27'][0], 2), vectors['28'][0])
    assert all(set(w) <= {0, -1} and not any(p[3:]) for w, p in cancellation)
    result = {}
    for oid, values in vectors.items():
        is_empty = oid in EMPTY_MAP
        assert all(weight(key) == 6+k for k, v in values.items() for key in v), oid
        primitive_inputs = ([f'snapshots/{EMPTY_RUN}/Q{EMPTY_MAP[oid]}q{k}.txt' for k in range(-2, 3)]
                            if is_empty else
                            [f'snapshots/{OPPOSITE_RUN}/I{("27" if oid == "29" else oid)}{q}.txt' for q in ('m1', 'zero')]
                            if oid in ('27', '28', '29') else
                            [f'snapshots/{MOTHER_RUN}/M_{q}_1.txt' for q in range(3)])
        entry = {
            'orbit_id': '5e-3g-' + oid, 'family': 't4-empty' if is_empty else 't4-opposite',
            'quantity': 'I',
            'laurent': {'variable': 'eps', 'lowest_stored_power': -6, 'highest_stored_power': 0,
                        'remainder': 'O(eps)', 'coefficients': {str(k): {'hlog_terms': serialize(values[k])} for k in range(-6, 1)}},
            'hlog_convention': {'endpoint': 'z', 'kernel': 'dt/(t-a)',
                'alphabet': [0, -1] if is_empty or oid in ('41', '42') else [0, -1, '-1/2', -2],
                'period_generators': list(PERIODS), 'period_weights': list(WEIGHTS),
                'domain': 'z>0; real principal logarithms; no new complex continuation claimed'},
            'primitive_inputs': primitive_inputs,
            'epsilon_budget': ({'Q_absolute_min': -2, 'Q_absolute_max': 2,
                'prefactor_absolute_max': 4, 'even_Q_absolute_max': 2,
                'relative_order_used': False, 'explicit_denominator_power': 2,
                'F_G_through': 4, 'Y_faces_through': 3, 'hard_recoil_TW_through': 2,
                'remainder': 'O(eps) at fixed z>0'} if is_empty else
                {'stored_hard_jet': [-1, 0], 'new_cold_integration': False,
                 'remainder': 'O(eps) at fixed z>0'} if oid in ('27', '28', '29') else
                {'bare_T4_orders': [-2, -1, 0], 'H0_through': 2, 'analytic_multiplier_through': 2,
                 'remainder': 'O(eps) at fixed z>0'}),
            'normalization': {'I_bridge': '(1-2eps)*G/8*(J(z)+z^(3eps)*J(1/z))/(1+z^(3eps))',
                'G': 'exp(-3*EulerGamma*eps)*Gamma(1-eps)*Gamma(1-2eps)',
                'c4_orbit_multiplicity': 4, 'euler_compensation_count': 1, 'native_assembly_sign': 1,
                'full_color_reduction': 'CF = CA^2*CF - 2*CA*CF^2; select a1 after full reduction'},
            'reassembly': {'method': 'fresh raw-Q inversion and absolute Gamma/cosh convolution' if is_empty else
                'reparse original completed global-IBP HyperInt coefficients, preserving wide coordinates' if oid in ('27', '28', '29') else
                'fresh product of admitted bare-T4 and H0 jets with its own physical prefactor',
                'new_cold_integration': False, 'final_I_reference_read': False, 'formal': 'not_run'},
        }
        dependency_paths = primitive_inputs + [
            'snapshots/benchmarks/current/prepared/3loop/5e-3g.json',
            'snapshots/benchmarks/tgv/5e-3g/prepared/5e-3g-r001.json',
            f'snapshots/{EMPTY_RUN if is_empty else OPPOSITE_RUN if oid in ("27", "28", "29") else MOTHER_RUN}/'
            + ('maple.log' if is_empty or oid in ('27', '28', '29') else 'compute_moments.mpl'),
            'scripts/t4.py', 'scripts/t4_endpoints.py']
        if not is_empty and oid in ('27', '28', '29'):
            dependency_paths += [f'snapshots/{OMEGA_RUN}/{name}' for name in
                                 ('Omega6.txt', 'Phi5.txt', 'Tprimitive.txt', 'maple.log')]
        entry['dependency_sha256'] = {p: hashlib.sha256((PACK / p).read_bytes()).hexdigest()
                                      for p in sorted(set(dependency_paths))}
        if is_empty:
            entry['budget_correction'] = {'affected_by_historical_relative_series_bug': oid in ('13', '14', '15', '33', '34', '35'),
                'correction': 'retain absolute even-Q powers through eps^2 before dividing by eps^2',
                'historical_implementation': 'Maple series(evenQ,eps=0,3) followed a relative Laurent order',
                'current_evaluation': 'all absolute-power convolutions computed from raw Q, not corrected I outputs'}
            entry['primitive_expansion'] = {'Q': {str(k): serialize(v) for k, v in raw_q[EMPTY_MAP[oid]].items()},
                'Q_definition': 'z*cuspX-z*hardY+(1+z)*(Y1-Y0)/eps-hardW+cuspW+derivative+4*(1-2eps)*recoil'}
            if oid != EMPTY_MAP[oid]:
                entry['full_epsilon_relation'] = {'equals': '5e-3g-' + EMPTY_MAP[oid],
                    'scope': 'full physical I at fixed z; ordered J has its own reflected domain and z^(3eps) transport'}
        else:
            oc = omega_coefficient[oid]
            if oc:
                residual = add(values[0], scale(omega, -oc))
                assert all(set(w) <= {0, -1} and not any(p[3:]) for w, p in residual), oid
                entry['known_function_view'] = {'finite': 'known old-alphabet part + Omega6_coefficient*Omega6',
                    'Omega6_coefficient': str(oc), 'known_part_hlog_terms': serialize(residual),
                    'scope': 'individual finite coefficient; cancellation in the sum is a separate identity'}
                entry['function_definitions'] = {'Omega6': {
                    'definition': 'T6(z)+T6(1/z)',
                    'T6_definition': 'integral_0^z Phi5(t)/(t+1/2) dt; T6(0)=0',
                    'hlog_terms': serialize(omega), 'Phi5_hlog_terms': serialize(phi),
                    'Tprimitive_hlog_terms': serialize(primitive),
                    'A4_definition': 'sum(n>=1,(-1)^n*H_(n-1)/n^3)',
                    'A6_definition': 'sum(n>=1,(-1)^n*H_(n-1)/n^5)',
                    'branch': 'z>0, ordinary convergent T6 at zero and real periods'}}
            if oid == '29':
                entry['full_epsilon_relation'] = {'equals': '5e-3g-27',
                    'scope': 'full physical I; ordered J transforms with reflection and z^(3eps)'}
            if oid in ('41', '42'):
                entry['function_definitions'] = {'M_' + str(q) + '_1': {
                    'definition': 'existing evaluated M-family generator in the fixed z-Hlog convention',
                    'hlog_terms': serialize(parse(SNAP / MOTHER_RUN / f'M_{q}_1.txt')),
                    'source': f'snapshots/{MOTHER_RUN}/M_{q}_1.txt'} for q in range(3)}
                entry['known_function_view'] = {
                    '-2': 'zeta(2)*M_0_1/64',
                    '-1': '(M_0_1*(zeta(3)-M_1_1)+zeta(2)*M_1_1/2)/64',
                    '0': '(M_0_1*(-M_2_1/8+M_0_1^2/48+9*zeta(2)^2/10)+M_1_1*(zeta(3)-M_1_1)/16+zeta(2)*(M_2_1/8-M_0_1^2/24)/8+M_0_1*zeta(2)*(5*M_0_1/4+17*zeta(2)/4)/8)/8',
                    'factorized_I': 'exp(-EulerGamma*eps)*Gamma(1-eps)*cosh(eps*log(z))/(8*cosh(3*eps*log(z)/2))*I_T4(z,eps)*H0(z,eps)',
                    'bare_T4': 'zeta(2)/(8eps^2)+(zeta(3)-M_1_1)/(8eps)-M_2_1/8+M_0_1^2/48+9*zeta(2)^2/10+O(eps)',
                    'H0': 'M_0_1+eps*M_1_1/2+eps^2*(M_2_1/8-M_0_1^2/24)+O(eps^3)',
                    'relative_prefactor': '1+eps^2*(5*M_0_1/4+17*zeta(2)/4)+O(eps^3)',
                    'source_weight': '1/8 includes each orbit own reflection and symmetry weight'}
                if oid == '42':
                    entry['full_epsilon_relation'] = {'equals': '5e-3g-41',
                        'scope': 'full physical I at fixed z; retain complementary ordered exchange domain'}
        result[entry['orbit_id']] = entry
    return result

if __name__ == '__main__':
    data = generate()
    print(json.dumps({'count': len(data), 'coefficient_count': sum(len(r['laurent']['coefficients']) for r in data.values()),
                      'orbits': sorted(data), 'formal': 'not_run'}, indent=2))
