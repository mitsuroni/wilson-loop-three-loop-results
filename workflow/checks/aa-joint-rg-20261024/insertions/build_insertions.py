"""Fresh native-input AA insertion sum through eps^-2, with exact lineage.

Only this directory is written. No RG target is read or used. Existing integral
values and moment identities are reused, never reclassified as fresh integration.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import sympy as s

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
POWERS = tuple(range(-6, -1))
SOURCES = {}
CHECKS = []


def sha(data):
    return hashlib.sha256(data).hexdigest()


def raw(path, expected=None):
    path = str(path).replace('\\', '/')
    data = (ROOT / path).read_bytes()
    digest = sha(data)
    assert path not in SOURCES or SOURCES[path]['sha256'] == digest, ('changed input', path)
    SOURCES[path] = {'sha256': digest, 'bytes': len(data)}
    if expected is not None:
        assert digest == expected, ('hash mismatch', path, expected, digest)
    return data


def read(path, expected=None):
    return json.loads(raw(path, expected).decode('utf-8-sig'))


def ref(path, pointer, value):
    raw(path)
    return {'path': path, 'sha256': SOURCES[path]['sha256'], 'json_pointer': pointer,
            'field_sha256': sha(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                           separators=(',', ':')).encode('utf-8'))}


PARSER = 'workflow/checks/matter-rg-20261017/check.py'
raw(PARSER)
spec = importlib.util.spec_from_file_location('existing_exact_parser', ROOT / PARSER)
math = importlib.util.module_from_spec(spec)
spec.loader.exec_module(math)
eps, rho = math.eps, math.rho
a = (10 - 6*eps)/(3 - 2*eps)


def parse(value):
    if isinstance(value, str):
        value = re.sub(r'M\((\d+),(\d+),z\)', r'M_\1_\2', value)
    return math.parse(value).subs(math.L, rho)


def simple(value):
    # Keep M symbols in the compact ledger; shuffle only z Hlog products.
    return s.cancel(math.shuffle_normalize(s.expand(value)))


def check(name, residual, negative=False):
    value = simple(residual)
    ok = value != 0 if negative else value == 0
    CHECKS.append({'name': name, 'passed': ok, 'negative_control': negative,
                   'residual': str(value)})
    assert ok, (name, str(value)[:1000])


def complete(data, values, pointer, path, channel=None):
    output, origins = {}, {}
    meta = data if channel is None else channel
    for p in POWERS:
        k = str(p)
        if k in values:
            output[p] = parse(values[k])
            origins[k] = ref(path, pointer+'/'+k, values[k])
        elif k in data.get('vanishing_higher_poles', {}):
            output[p] = parse(data['vanishing_higher_poles'][k])
            assert output[p] == 0
            origins[k] = ref(path, '/vanishing_higher_poles/'+k,
                             data['vanishing_higher_poles'][k])
        elif p in data.get('zero_powers', []):
            output[p] = s.S.Zero
            origins[k] = ref(path, '/zero_powers', data['zero_powers'])
        elif 'highest_pole_order' in meta and p < -meta['highest_pole_order']:
            output[p] = s.S.Zero
            origins[k] = ref(path, ('/channels/g' if channel else '')+'/highest_pole_order',
                             meta['highest_pole_order'])
        else:
            raise ValueError(('unknown, not zero', path, pointer, k))
    return output, origins


def convolve(values, multiplier):
    series = s.series(multiplier, eps, 0, 5).removeO().expand()
    result = {p: simple(sum(v*series.coeff(eps, p-q)
                           for q, v in values.items() if q <= p)) for p in POWERS}
    polynomial = s.expand(sum(v*eps**q for q, v in values.items()) * series)
    for p in POWERS:
        check('convolution polynomial agrees at '+str(p), result[p]-polynomial.coeff(eps,p))
    return result


def serial(values):
    return {str(p): str(simple(values[p])) for p in POWERS}


def source_contracts(data, selector, family):
    out = []
    if family in ('Pi1', 'chain', 'Pi2'):
        for field in ('raw', 'prepared'):
            item = selector[field]
            read(item['path'], item['sha256'])
            out.append(dict(item, kind=field))
    else:
        item = data['input']
        read(item['native_package'], item['native_package_sha256'])
        assert item['assembly_sign'] == 1
        out.append({'path': item['native_package'], 'sha256': item['native_package_sha256'],
                    'kind': 'native C/J', 'assembly_sign': 1})
    return out


def main():
    outer = read('benchmarks/self-energy/current.json')
    assert outer['manifest'] == 'current-r002.json'
    read('benchmarks/self-energy/'+outer['manifest'],outer['sha256'])
    selectors = {
        'se': read('benchmarks/self-energy/current-r002.json'),
        'legs': read('benchmarks/self-energy/tgv-pi1/current.json'),
        'Gamma': read('benchmarks/self-energy/gamma1/current.json'),
    }
    assert selectors['se']['external_source_sign'] == 1
    for path in ('benchmarks/self-energy/conventions.md',
                 'benchmarks/theory/montopar/shared-sign-conventions.md',
                 'workflow/checks/shared-sign-production/native-assembly-check.json',
                 'workflow/checks/shared-sign-production/publication-check.json',
                 'workflow/checks/tgv-pi1-representation-r003/verification.json'):
        raw(path)
    dictionaries = {'z': {}, 'r': {}}
    for arg, base in [('z', 'workflow/checks/single-opposite/runs/20260916T090538777022Z'),
                      ('r', 'workflow/checks/tgv-pi1-07-mn')]:
        for name in ('M_0_1', 'M_1_1', 'M_2_1'):
            text = raw(f'{base}/{name}.txt').decode('utf-8-sig').strip()
            dictionaries[arg][s.Symbol(name)] = parse(text.split('=', 1)[-1])

    rows = []
    for sel in selectors['se']['entries']:
        family = {'4e-Pi1':'Pi1','2e-Pi1+Pi1':'chain','2e-Pi2':'Pi2'}.get(sel['family'])
        if family is None:
            continue
        path = sel['result']['path']
        data = read(path, sel['result']['sha256'])
        assert data['normalization']['assembly_sign'] == sel['assembly_sign'] == 1
        assert data['normalization']['Euler_count'] == 1
        if family == 'Pi1':
            pointer = '/laurent/coefficients'
            fields = data['laurent']['coefficients']
            multiplier = a
        else:
            index = next(i for i,c in enumerate(data['color_channels'])
                         if c.get('monomial',c.get('color_monomial')) == 'CA^2*CF')
            channel = data['color_channels'][index]
            pointer = f'/color_channels/{index}/'+('coefficients' if family=='chain' else 'laurent/coefficients')
            fields = channel['coefficients'] if family == 'chain' else channel['laurent']['coefficients']
            multiplier = s.S.One
        values, origins = complete(data, fields, pointer, path)
        weighted = convolve(values,multiplier)
        if family == 'chain':
            scalar, _ = complete(data,data['laurent']['coefficients'],'/laurent/coefficients',path)
            independent = convolve(scalar,a*a)
            for p in POWERS:
                check(sel['orbit_id']+' stored AA equals a^2 scalar '+str(p), weighted[p]-independent[p])
        rows.append({'orbit_id':sel['orbit_id'],'family':family,
                     'source':ref(path,pointer,fields),'coefficient_sources':origins,
                     'native_coefficients':serial(values),'weight':str(multiplier),
                     'weighted_coefficients':serial(weighted),
                     'normalization':data['normalization'],
                     'contracts':source_contracts(data,sel,family)})

    for family,base in [('legs','benchmarks/self-energy/tgv-pi1'),
                         ('Gamma','benchmarks/self-energy/gamma1')]:
        for oid,sel in selectors[family]['orbits'].items():
            chosen = sel.get('preferred_representation',sel)
            path = base+'/'+chosen['result']
            data = read(path,chosen['sha256'])
            assert data['normalization']['assembly_sign'] == 1
            assert data['normalization']['euler_compensation_count'] == 1
            assert data['normalization']['c4_multiplicity'] == 4
            if family == 'Gamma':
                assert data['channels']['g']['color_factor'] == 'CF*CA^2'
                assert not data['normalization']['external_leg_self_energy_included']
                pointer = '/channels/g/I_coefficients'
                fields = data['channels']['g']['I_coefficients']
                values,origins = complete(data,fields,pointer,path,data['channels']['g'])
                if chosen is not sel:
                    read(base+'/'+sel['result'],sel['sha256'])
                    assert data['mathematical_result_sha256'] == sel['sha256']
                    for leaf, digest in data['legacy_evidence']['evidence_sha256'].items():
                        raw('workflow/checks/gamma1-03-z/'+leaf,digest)
                multiplier = s.S.One
            else:
                pointer = '/laurent_coefficients'
                fields = data['laurent_coefficients']
                values,origins = complete(data,fields,pointer,path)
                if 'second_pole_in_existing_M_basis' in data:
                    compact = parse(data['second_pole_in_existing_M_basis'])
                    check(oid+' double pole: r-moment dictionary to primary Hlog',
                          compact.subs(dictionaries['r'],simultaneous=True)-values[-2])
                    values[-2] = compact
                    origins['-2'] = ref(path,'/second_pole_in_existing_M_basis',
                                        data['second_pole_in_existing_M_basis'])
                if 'laurent_Hlog_coefficients' in data:
                    for p in POWERS:
                        if str(p) not in data['laurent_Hlog_coefficients']:
                            continue
                        alternate = parse(data['laurent_Hlog_coefficients'][str(p)])
                        arg = 'r' if any(str(x).startswith('H_r_') for x in alternate.free_symbols) else 'z'
                        check(oid+' compact to '+arg+' primary Hlog '+str(p),
                              values[p].subs(dictionaries[arg],simultaneous=True)-alternate)
                multiplier = a
            weighted = convolve(values,multiplier)
            if family == 'legs' and 'color_channel_coefficients' in data:
                aa = data['color_channel_coefficients'].get('CF*CA^2')
                if aa:
                    for p in POWERS:
                        if str(p) in aa:
                            recorded = parse(aa[str(p)])
                            candidate = weighted[p]
                            if any(str(x).startswith('H_r_') for x in recorded.free_symbols):
                                candidate = candidate.subs(dictionaries['r'],simultaneous=True)
                            elif any(str(x).startswith('H_z_') for x in recorded.free_symbols):
                                candidate = candidate.subs(dictionaries['z'],simultaneous=True)
                            check(oid+' native stored AA convolution '+str(p),candidate-recorded)
            # Test the published native sign against its bound legacy result.
            migration = data.get('migration',{})
            if family == 'Gamma' and migration.get('factor_from_native_CJ') == -1:
                oldref = migration['legacy_result']
                old = read(oldref['path'],oldref['sha256'])
                for p in POWERS:
                    oldv = old['channels']['g']['I_coefficients'].get(str(p),'0')
                    check(oid+' native g sign bridge '+str(p),values[p]+parse(oldv))
            rows.append({'orbit_id':oid,'family':family,'source':ref(path,pointer,fields),
                         'coefficient_sources':origins,'native_coefficients':serial(values),
                         'weight':str(multiplier),'weighted_coefficients':serial(weighted),
                         'normalization':data['normalization'],
                         'contracts':source_contracts(data,sel,family)})

    def compactify(value, label):
        """Eliminate existing moment pivots, then use rho^2=-2 M01-pi^2.

        All coefficients come from the computed expression, not an RG target.
        The complete candidate is expanded back to the original z expression.
        """
        value = simple(value)
        residual, moments = s.expand(value), s.S.Zero
        for name,pivot in [('M_2_1','H_z_0_m1_0_0'),('M_1_1','H_z_m1_0_0')]:
            symbol,pivot = s.Symbol(name),s.Symbol(pivot)
            coefficient = s.expand(residual).coeff(pivot)/dictionaries['z'][symbol].expand().coeff(pivot)
            residual = s.expand(residual-coefficient*dictionaries['z'][symbol])
            moments += coefficient*symbol
        assert not any(str(x).startswith('H_') for x in residual.free_symbols), (label,residual)
        polynomial = s.Poly(residual,rho)
        assert all(monomial[0]%2==0 for monomial in polynomial.monoms()), (label,residual)
        compact = moments + sum(coefficient*(-2*s.Symbol('M_0_1')-s.pi**2)**(monomial[0]//2)
                                for monomial,coefficient in polynomial.terms())
        compact = s.expand(compact)
        check(label+' compact M representation back-substitution',
              compact.subs(dictionaries['z'],simultaneous=True)-value)
        return compact

    expected = {'Pi1':11,'chain':2,'Pi2':2,'legs':12,'Gamma':4}
    families = {}
    for family,count in expected.items():
        group = [row for row in rows if row['family']==family]
        assert len(group)==count
        compact = {p:simple(sum(parse(row['weighted_coefficients'][str(p)]) for row in group)) for p in POWERS}
        zvalues = {p:simple(v.subs(dictionaries['z'],simultaneous=True)) for p,v in compact.items()}
        assert not any(str(x).startswith('H_r_') for v in zvalues.values() for x in v.free_symbols)
        unified_compact = {p:compactify(v,family+' '+str(p)) for p,v in zvalues.items()}
        families[family] = {'diagram_count':count,'diagram_ids':[row['orbit_id'] for row in group],
                            'coefficients':serial(zvalues),'compact_coefficients':serial(unified_compact)}
    total = {p:simple(sum(parse(row['coefficients'][str(p)]) for row in families.values())) for p in POWERS}
    compact_total = {p:compactify(v,'total '+str(p)) for p,v in total.items()}
    # Controls target actual hazards, without using expected RG coefficients.
    p1 = next(row for row in rows if row['family']=='Pi1')
    v1 = {p:parse(p1['native_coefficients'][str(p)]) for p in POWERS}
    full = convolve(v1,a)
    truncated = convolve(v1,a.subs(eps,0))
    check('reject dropping epsilon dependence of a at eps^-4',full[-4]-truncated[-4],True)
    g1 = next(row for row in rows if row['orbit_id']=='3e-Gamma1-01')
    check('reject a second legacy sign on native Gamma g',2*parse(g1['weighted_coefficients']['-6']),True)
    chain = next(row for row in rows if row['family']=='chain')
    cv = {p:parse(chain['weighted_coefficients'][str(p)]) for p in POWERS}
    check('reject double color weight on already-colored chain',convolve(cv,a*a)[-4]-cv[-4],True)
    # End-of-run source stability, including selectors.
    for path, info in list(SOURCES.items()):
        raw(path,info['sha256'])
    output = {
        'format':'wilson-aa-native-insertion-sum/1','date':'2026-10-24',
        'producer':{'path':Path(__file__).relative_to(ROOT).as_posix(),
                    'sha256':sha(Path(__file__).read_bytes())},
        'status':'exact_existing-input_assembly_through_eps_minus2',
        'color':'CF*CA^2','quantity':'U_3,AA: coupling and S3 stripped connected exponent',
        'epsilon_range':[-6,-2], 'diagram_count':len(rows),
        'normalization':'native physical signs; Euler and C4 already inside I; a only for Pi1/legs; stored colored chain/Pi2/Gamma fields once',
        'families':families,'total_coefficients':serial(total),
        'total_compact_coefficients':serial(compact_total),'diagrams':rows,
        'function_conventions':{'rho':'log(z)','M_k_j':'original raw moment, no factorial division',
                                'H_z_word':'G(word;z), kernels dt/(t-a); H_z_m1=log(1+z)',
                                'r':'1/(1+z); primary r fields compared in r algebra before replacing with the existing M dictionary'},
        'rg_targets_used':False,'formal':'not_run','new_integration_performed':False,
        'limitations':['Reuses selected integral results and pre-existing exact function dictionaries.',
                       'Gamma03 preferred z conversion is reused with original inverse/crossing evidence, not freshly rerun in Maple.',
                       'Only coefficients eps^-6 through eps^-2 are assembled.'],
        'checks':CHECKS,'source_sha256':SOURCES,
    }
    (HERE/'insertions.json').write_text(json.dumps(output,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'status':'passed','diagrams':len(rows),'checks':len(CHECKS),
                      'negative_controls':sum(c['negative_control'] for c in CHECKS),
                      'sources':len(SOURCES),'families':families,'total_coefficients':serial(total)},indent=2))


if __name__=='__main__':
    main()
