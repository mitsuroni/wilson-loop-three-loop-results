"""Independent exact replay of the tau/E dictionary and indexed-source conventions.

No fitting, integration, or writes to old evidence. The scalar-period shuffle
algebra below is independent of the archived hlog7/hlog8/exact_span adapters.
"""
from __future__ import annotations

import ast
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import re

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OLD = ROOT / 'workflow/checks/gluonic-beta-sources'
SOURCES: dict[str, str] = {}
CHECKS: list[dict] = []
EMPTY = (0, 0, 0, 0, 0)
ONE = {((), EMPTY): Q(1)}
LOG = {((0,), EMPTY): Q(1)}
PERIOD_WEIGHTS = (2, 3, 5, 7, 8)


def read(path: Path) -> str:
    SOURCES[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path.read_text(encoding='utf-8')


def data(path: Path):
    return json.loads(read(path))


def check(name: str, condition: bool, **extra):
    CHECKS.append(dict(name=name, passed=bool(condition), **extra))
    assert condition, name


def scale(v, q):
    return {k: c * q for k, c in v.items() if c * q}


def add(a, b):
    v = dict(a)
    for k, c in b.items():
        v[k] = v.get(k, Q(0)) + c
        if not v[k]:
            del v[k]
    return v


@lru_cache(None)
def shuffle(a, b):
    if not a:
        return ((b, 1),)
    if not b:
        return ((a, 1),)
    out = defaultdict(int)
    for w, n in shuffle(a[1:], b):
        out[(a[0],) + w] += n
    for w, n in shuffle(a, b[1:]):
        out[(b[0],) + w] += n
    return tuple(sorted(out.items()))


def mul(a, b):
    out = defaultdict(Q)
    for (wa, pa), ca in a.items():
        for (wb, pb), cb in b.items():
            p = tuple(x + y for x, y in zip(pa, pb))
            for w, n in shuffle(wa, wb):
                out[w, p] += ca * cb * n
    return {k: c for k, c in out.items() if c}


def power(v, n):
    assert isinstance(n, int) and n >= 0
    out = ONE
    for _ in range(n):
        out = mul(out, v)
    return out


def scalar(v):
    assert not (set(v) - {((), EMPTY)})
    return v.get(((), EMPTY), Q(0))


def expression(text, names):
    def visit(n):
        if isinstance(n, ast.Constant) and type(n.value) is int:
            return scale(ONE, Q(n.value))
        if isinstance(n, ast.Name):
            return names[n.id]
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            return scale(visit(n.operand), -1)
        if isinstance(n, ast.BinOp):
            a, b = visit(n.left), visit(n.right)
            if isinstance(n.op, ast.Add):
                return add(a, b)
            if isinstance(n.op, ast.Sub):
                return add(a, scale(b, -1))
            if isinstance(n.op, ast.Mult):
                return mul(a, b)
            if isinstance(n.op, ast.Div):
                return scale(a, 1 / scalar(b))
            if isinstance(n.op, ast.Pow):
                q = scalar(b)
                assert q.denominator == 1
                return power(a, int(q))
        raise ValueError(ast.dump(n))
    return visit(ast.parse(text.replace('^', '**'), mode='eval').body)


def parse_raw(path):
    text = read(path).split('=', 1)[1].strip()
    names = {}
    for i, w in enumerate(PERIOD_WEIGHTS):
        p = tuple(int(k == i) for k in range(5))
        names['Z' + str(w)] = {((), p): Q(1)}

    def word(m):
        w = tuple(int(x.strip()) for x in m[1].split(',') if x.strip())
        assert set(w) <= {0, -1}
        name = 'H' + str(len(names))
        names[name] = {(w, EMPTY): Q(1)}
        return name

    text = re.sub(r'Hlog\(z,\[([^]]*)\]\)', word, text)
    text = re.sub(r'zeta\[(2|3|5|7)\]', lambda m: 'Z' + m[1], text)
    assert 'zeta[' not in text
    return expression(text, names)


def lift(rows):
    return {(tuple(r['word']), tuple(r['zeta_powers']) + (0,) * (5 - len(r['zeta_powers']))):
            Q(r['coefficient']) for r in rows}


def serialize(v):
    return [dict(word=w, zeta_powers=p, coefficient=str(c)) for (w, p), c in sorted(v.items())]


# Read the raw independent integration outputs, not the saved fitted G vectors.
raw = {}
for s in (1, 2):
    for prefix in ('J', 'Jb'):
        path = OLD / f'runs/free-beta-w6/{prefix}_0_1_1_{s}.txt'
        raw[prefix, s] = parse_raw(path)
        repeated = parse_raw(OLD / f'runs/free-beta-w8/{path.name}')
        check('raw low/high job agreement ' + path.name, raw[prefix, s] == repeated)

T12 = add(raw['J', 2], raw['Jb', 2])
S11 = add(raw['J', 1], scale(raw['Jb', 1], -1))
LS11 = mul(LOG, S11)
report6 = data(OLD / 'comparison.json')
gmap = {g['name']: lift(g['hlog_terms']) for g in report6['generators']}
check('T12 is G_0_1_1_2_0', T12 == gmap['G_0_1_1_2_0'])
check('ell S11 is G_0_1_1_1_1', LS11 == gmap['G_0_1_1_1_1'])
names = {n[0] + '_' + n.split('_', 1)[1] if n.startswith(('M0_', 'N0_')) else n: v
         for n, v in gmap.items()}
names.update(T12=T12, S11=S11, ell=LOG)

# These readable expressions are written independently of the saved fit list.
RTAU = ('-11*N_2_2/8 + M_1_1**2/4 + 5*M_2_3/64 + 29*M_0_1**3/288 '
        '+31*M_4_1/128 -12*M_1_1*Z3 -459*M_0_1*Z2**2/40 '
        '-23*M_0_1**2*Z2/48 +3*M_0_1*M_2_1/8 +39*M_2_1*Z2/8')
RE = ('-M_2_3/4 +N_2_2/4 +3*M_0_1**3/20 -79*M_0_1**2*Z2/60 '
      '-53*M_0_1*Z2**2/10 -3*M_2_1*Z2/4 -M_0_1*M_2_1/4')
tau = lift(data(ROOT / 'benchmarks/2ne/families/single-opposite/basis-functions.json')['tauDq4']['hlog_terms'])
E = lift(data(ROOT / 'benchmarks/2ne/families/same-channel-two/structural-basis-functions.json')['functions']['chi_id_0_1_1']['hlog_terms'])
tau_back = expression('-3*T12-2*ell*S11+(' + RTAU + ')', names)
E_back = expression('-ell*S11+(' + RE + ')', names)
check('full tau dictionary exact', tau_back == tau, target_coordinates=len(tau))
check('full E dictionary exact', E_back == E, target_coordinates=len(E))
check('isolated tau minus twice E exact',
      add(tau, scale(E, -2)) == expression('-3*T12+(' + RTAU + ')-2*(' + RE + ')', names))
check('wrong T coefficient rejected', add(tau_back, T12) != tau)
check('wrong S sign rejected', add(E_back, scale(LS11, 2)) != E)

def weight(key):
    return len(key[0]) + sum(a * b for a, b in zip(key[1], PERIOD_WEIGHTS))

for name, v, w in [('T12', T12, 6), ('S11', S11, 5), ('ell S11', LS11, 6),
                   ('tau', tau, 6), ('E', E, 6)]:
    check('homogeneous weight ' + name, all(weight(k) == w for k in v))

# Diagnose the smaller odd seed after explicitly allowing ell and all products.
# This is a declared constant-Q span; it does not allow arbitrary rational z.
base_generators = [(g['name'], g['weight'], g['parity'], lift(g['hlog_terms']))
                   for g in report6['generators'] if not g['name'].startswith('G_')]
base_generators.append(('ell', 1, -1, LOG))


def candidates(total_weight, wanted_parity):
    cols = []
    def rec(start, left, parity, value, factors):
        if left == 0:
            if parity == wanted_parity and value:
                cols.append((factors, value))
            return
        for i in range(start, len(base_generators)):
            name, w, pi, v = base_generators[i]
            if w <= left:
                rec(i, left-w, parity*pi, mul(value, v), factors+[name])
    rec(0, total_weight, 1, ONE, [])
    return cols


def echelon(columns):
    pivots = {}
    for col in columns:
        value = dict(col)
        for pivot, basis in sorted(pivots.items()):
            if pivot in value:
                value = add(value, scale(basis, -value[pivot]))
        if value:
            pivot = min(value)
            pivots[pivot] = scale(value, 1/value[pivot])
    return pivots


def residual(target, basis):
    value = dict(target)
    for pivot, row in sorted(basis.items()):
        if pivot in value:
            value = add(value, scale(row, -value[pivot]))
    return value


def pairing(functional, vector):
    return sum((c*vector.get(k, Q(0)) for k, c in functional.items()), Q(0))


def nonmember_witness(target, basis):
    rem = residual(target, basis)
    if not rem:
        return {}
    f = {min(rem): Q(1)}
    for pivot, row in reversed(sorted(basis.items())):
        factor = pairing(f, row)
        f[pivot] = f.get(pivot, Q(0)) - factor
        if not f[pivot]:
            del f[pivot]
    return f


span_diagnostics = []
for w, parity, targets in [(5, -1, [('S11', S11)]), (6, 1, [('T12', T12), ('ell S11', LS11)])]:
    options = candidates(w, parity)
    basis = echelon([v for _, v in options])
    target_rows = []
    for name, v in targets:
        rem = residual(v, basis)
        witness = nonmember_witness(v, basis)
        if rem:
            check('nonmember witness annihilates all candidates ' + name,
                  all(pairing(witness, col) == 0 for _, col in options), candidate_count=len(options))
            check('nonmember witness separates target ' + name, pairing(witness, v) != 0)
        target_rows.append(dict(name=name, member=not rem, residual=serialize(rem),
                                witness=serialize(witness), obstruction=str(pairing(witness, v))))
    augmented = echelon([v for _, v in options] + [v for _, v in targets])
    span_diagnostics.append(dict(weight=w, parity=parity, candidate_count=len(options),
                                 rank=len(basis), augmented_rank=len(augmented), targets=target_rows))
    check(f'constant-Q span rank consistency weight {w}',
          len(basis) <= len(augmented) <= len(basis)+len(targets))
    if w == 6:
        with_T = echelon([v for _, v in options] + [T12])
        witness = nonmember_witness(LS11, with_T)
        check('two weight-six directions independent modulo MN+ell',
              all(pairing(witness, col) == 0 for _, col in options)
              and pairing(witness, T12) == 0 and pairing(witness, LS11) != 0)
        span_diagnostics[-1]['joint_witness'] = serialize(witness)
        span_diagnostics[-1]['joint_obstruction'] = str(pairing(witness, LS11))

# Factorial, reflection and removable-pole conventions are independent symbolic checks.
p, c, S, D = sp.symbols('p c S D')
for m in range(5):
    for s in range(5 - m):
        value = sp.diff(sp.exp(p * S + c * D), p, m, c, s).subs({p: 0, c: 0})
        check(f'derivative-index convention {m},{s}', sp.expand(value - S**m * D**s) == 0,
              Taylor_denominator=math.factorial(m) * math.factorial(s))

z, beta, u = sp.symbols('z beta u', positive=True)
ell, R = sp.symbols('ell R', real=True)
r = z / (1 + z)
K = (ell + R) / (beta - r)
Kref = (-ell - R) / ((1 - beta) - (1 / z) / (1 + 1 / z))
check('reflected outer kernel invariant', sp.cancel(Kref - K) == 0)
check('Dplus maps to Dminus', sp.expand((1-u*beta).subs(beta, 1-beta) - (1-u+u*beta)) == 0)
check('pole location maps to complement', sp.cancel((1/z)/(1+1/z) - (1-r)) == 0)
check('removable pole value', sp.simplify((-1/(beta*(1-beta))).subs(beta, r) + (1+z)**2/z) == 0)

# The m=0 slice integrates in u by elementary FTC for every integer s>=0.
for s in range(7):
    primitive = -sp.log(1-u*beta)**(s+1) / (s+1)
    check(f'zero Beta index inner FTC s={s}',
          sp.simplify(sp.diff(primitive, u) - beta*sp.log(1-u*beta)**s/(1-u*beta)) == 0)

for s in range(4):
    P = parse_raw(OLD / f'runs/free-beta-w6/J_0_1_0_{s}.txt')
    Pmirror = parse_raw(OLD / f'runs/free-beta-w6/Jb_0_1_0_{s}.txt')
    for parity, label, target in [(0, 'T', add(P, Pmirror)), (1, 'S', add(P, scale(Pmirror, -1)))]:
        built = {}
        for v in range(parity, s+2, 2):
            for h in range(v+1):
                factor = -Q(math.comb(s+1, v)*math.comb(v, h)*(-1)**(v-h), (s+1)*2**s)
                term = mul(power(LOG, v-h), gmap[f'M0_{s+1-v}_{h+1}'])
                built = add(built, scale(term, factor))
        check(f'zero Beta index explicit M reduction {label}0,{s}', built == target)

# Match every used archived byte with the previous validation source ledger.
validation = data(OLD / 'validation.json')
old_bindings = validation['source_sha256']
bound = 0
for rel, digest in list(SOURCES.items()):
    if rel in old_bindings:
        check('archived byte binding ' + rel, digest == old_bindings[rel])
        bound += 1
check('all eight reused raw outputs were historically bound',
      all((OLD / f'runs/free-beta-w{w}/{prefix}_0_1_1_{s}.txt').relative_to(ROOT).as_posix()
          in old_bindings for w in (6, 8) for prefix in ('J', 'Jb') for s in (1, 2)))

report = dict(
    passed=True, checks=CHECKS, check_count=len(CHECKS), formal='not_run',
    new_integration_runs=0,
    reused_raw_files=sum('/runs/free-beta-' in rel and rel.endswith('.txt') for rel in SOURCES),
    historical_bound_inputs=bound,
    convention='T and S are raw derivatives, not factorial-normalized Taylor coefficients',
    formulas={'tau': '-3*T12 -2*ell*S11 + (' + RTAU + ')',
              'E': '-ell*S11 + (' + RE + ')'},
    source_vectors={'T12': serialize(T12), 'S11': serialize(S11)},
    MN_log_span=span_diagnostics,
    source_sha256=SOURCES,
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    scope='Exact archived-vector replay and analytic source conventions; no new all-epsilon C/D reduction, completeness or minimality claim.'
)
(HERE / 'verification.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps(dict(passed=True, checks=len(CHECKS), T12_coordinates=len(T12),
                     S11_coordinates=len(S11), historical_bound_inputs=bound,
                     MN_log_span=[dict(weight=r['weight'], rank=r['rank'], augmented=r['augmented_rank'],
                                       targets=[dict(name=t['name'], member=t['member']) for t in r['targets']])
                                  for r in span_diagnostics], formal='not_run')))
