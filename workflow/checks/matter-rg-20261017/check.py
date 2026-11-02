"""Read current matter inputs; write only this new evidence directory.

No single-graph integration, selector update, or historical evidence rewrite.
SymPy treats mixed Hlog words as independent symbols: equality is exact in
the displayed expression algebra, without assuming unregistered identities.
"""
from pathlib import Path
import argparse
import collections
import hashlib
import json
import re
import time
from functools import lru_cache
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TASK_INPUTS = ROOT / '.wilson/task-inputs.json'
if TASK_INPUTS.is_file():
    ROOT = Path(json.loads(TASK_INPUTS.read_text(encoding='utf-8'))['read_only_input_root'])
eps, L, rho, z, r, CA, CF, T = s.symbols('eps L rho z r CA CF T')
CHECKS, SOURCES, LINEAGE, CHANGES, GAPS = [], {}, [], [], []
START = time.time()

def digest(data):
    return hashlib.sha256(data).hexdigest()

def read(rel):
    data = (ROOT / rel).read_bytes()
    h = digest(data)
    if rel in SOURCES:
        assert SOURCES[rel] == h, ('concurrent input change', rel)
    SOURCES[rel] = h
    return json.loads(data.decode('utf-8-sig'))

def read_text(rel):
    data = (ROOT / rel).read_bytes()
    h = digest(data)
    if rel in SOURCES:
        assert SOURCES[rel] == h, ('concurrent input change', rel)
    SOURCES[rel] = h
    return data.decode('utf-8-sig')

def field_hash(value):
    return digest(json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(',', ':')).encode())

def hl(arg, word):
    if not word:
        return s.S.One
    if arg == 'z' and all(x == 0 for x in word):
        return rho ** len(word) / s.factorial(len(word))
    return s.Symbol('H_' + arg + '_' + '_'.join('m1' if x == -1 else str(x) for x in word))

def parse(expr):
    if isinstance(expr, dict):
        if 'hlog_terms' in expr:
            return s.expand(sum(parse(t['coefficient']) * hl('z', t['word']) *
                                s.prod(s.zeta(n)**p for n,p in zip([2,3,5],t['zeta_powers']))
                                for t in expr['hlog_terms']))
        expr = expr.get('sympy', expr.get('expression'))
    if not isinstance(expr, str):
        return s.sympify(expr)
    expr = re.sub(r'Hlog\(([zr]),\s*\[([^]]*)\]\)',
                  lambda m: '(' + str(hl(m[1], [int(x) for x in m[2].split(',') if x.strip()])) + ')', expr)
    expr = re.sub(r'\b(?:Z|zeta_)([235])\b', r'zeta(\1)', expr)
    expr = re.sub(r'[Zz]eta\[(\d+)\]', r'zeta(\1)', expr)
    expr = re.sub(r'Mom\[(\d+),(\d+)\]', r'M_\1_\2', expr)
    expr = re.sub(r'Nom\[(\d+),(\d+)\]', r'N_\1_\2', expr)
    expr = re.sub(r'\bM([0-9])([0-9])\b', r'M_\1_\2', expr)
    expr = re.sub(r'\bN([0-9])([0-9])\b', r'N_\1_\2', expr)
    expr = expr.replace('Pi', 'pi').replace('Zeta(', 'zeta(').replace('^','**').replace('ln(', 'log(')
    value = s.sympify(expr)
    assert not value.atoms(s.Float), ('inexact parsed input', expr)
    return value.subs(s.log(z), rho)

ELEMENTARY = {s.Symbol('M_0_1'): -(rho**2+s.pi**2)/2,
              s.Symbol('M_0_3'): -(rho**2+s.pi**2)**2/4}

@lru_cache(None)
def shuffle(a,b):
    if not a: return {b:1}
    if not b: return {a:1}
    out=collections.Counter()
    for word,n in shuffle(a[1:],b).items(): out[(a[0],)+word]+=n
    for word,n in shuffle(a,b[1:]).items(): out[(b[0],)+word]+=n
    return dict(out)

def shuffle_normalize(expr):
    """Standard iterated-integral shuffle; rho=Hlog(z,[0])."""
    output=0
    for term in s.Add.make_args(s.expand(expr)):
        powers=term.as_powers_dict(); words=[]; coefficient=term
        has_h=any(str(f).startswith('H_z_') for f in powers)
        if not has_h:
            output+=term;continue
        for factor,power in powers.items():
            if str(factor).startswith('H_z_') and power.is_Integer and power>0:
                word=tuple(-1 if t=='m1' else int(t) for t in str(factor)[4:].split('_'))
                words.extend([word]*int(power));coefficient/=factor**power
        rp=powers.get(rho,s.S.Zero)
        if rp.is_Integer and rp>0:
            words.extend([(0,)]*int(rp));coefficient/=rho**rp
        if len(words)<2:
            output+=term;continue
        merged={():1}
        for word in words:
            nxt=collections.Counter()
            for a,m in merged.items():
                for b,n in shuffle(a,word).items(): nxt[b]+=m*n
            merged=nxt
        output+=coefficient*sum(n*hl('z',word) for word,n in merged.items())
    return s.expand(output)

def norm(x):
    x = shuffle_normalize(s.expand(x.subs(ELEMENTARY)))
    if x != 0:
        x = s.cancel(x)
    return x

def ck(name, value, channel=None):
    v = norm(value)
    CHECKS.append(dict(name=name, channel=channel, passed=v == 0, residual=str(v)))
    if v != 0:
        raise AssertionError((name, str(v)[:1500]))

def nonzero(name, value):
    v = norm(value)
    CHECKS.append(dict(name=name, passed=v != 0, residual=str(v), negative_control=True))
    assert v != 0, name

def cs(values):
    return {int(q): parse(v) for q,v in values.items()}

def scale(n, j):
    return s.expand(n**j*((L+rho/2)**j+(L-rho/2)**j)/s.factorial(j))

def fullscale(values, n, lo=-6, hi=0):
    return {q: norm(sum(values.get(j,0)*scale(n,q-j) for j in range(lo,q+1))) for q in range(lo,hi+1)}

def strip_scale(values,n,lo,hi):
    unit={}
    for q in range(lo,hi+1):
        unit[q]=norm((values[q]-sum(unit[j]*scale(n,q-j) for j in range(lo,q)))/2)
    return unit

def require_budget(values,lo,hi,label):
    missing=sorted(set(range(lo,hi+1))-set(values))
    if missing:
        raise ValueError(f'{label}: missing epsilon powers {missing}')

BF = -8*(1-eps)/(3-2*eps)
BFC = {j: s.expand(s.series(BF,eps,0,7).removeO()).coeff(eps,j) for j in range(7)}

def multiply_bf(values):
    return {q: norm(sum(v*BFC.get(q-j,0) for j,v in values.items())) for q in range(-6,1)}

def bound_read(rel, sha):
    data = read(rel)
    assert SOURCES[rel] == sha, ('selector hash mismatch',rel,sha,SOURCES[rel])
    return data

def changed_paths(a,b,prefix=''):
    if isinstance(a,dict) and isinstance(b,dict):
        out=[]
        for key in sorted(set(a)|set(b)):
            path=prefix+'/'+key
            if key not in a or key not in b: out.append(path)
            else: out.extend(changed_paths(a[key],b[key],path))
        return out
    if a != b:
        return [prefix]
    return []

def audit_original_bindings(label, data):
    for rel, expected in data.get('source_sha256',{}).items():
        actual=digest((ROOT/rel).read_bytes()) if (ROOT/rel).is_file() else None
        LINEAGE.append(dict(owner=label,path=rel,original_sha256=expected,current_sha256=actual,
                            original_bytes_present=actual==expected))

def check_current_31():
    old=read('workflow/checks/w3-ca-nf/inputs.json')
    se=read('benchmarks/self-energy/current-r002.json')
    native={x['orbit_id']: x for x in se['entries']}
    legs=read('benchmarks/self-energy/tgv-pi1/current.json')
    gamma=read('benchmarks/self-energy/gamma1/current.json')
    assert legs['physical_convention']==gamma['physical_convention']=='wilson-qcd-action-path-fourier/1.0'
    bridge=read('workflow/checks/w3-ca-nf/physical-bridge.json')
    assert bridge['passed'] and not bridge['RG_used_to_derive_sign'] and len(bridge['checks'])==31
    assert bridge['weights']==dict(Pi1=1,chain=1,Pi2=1,legs=-1,Gamma=-1)
    frozen='workflow/checks/shared-sign-production/publication/'
    frozen_gamma=read(frozen+'benchmarks/self-energy/gamma1/current.json')
    current_sum={q:s.S.Zero for q in range(-6,1)}
    out=[]
    assert collections.Counter(x['part'] for x in old['parts']) == dict(Pi1=11,chain=2,Pi2=2,legs=12,Gamma=4)
    for row in old['parts']:
        oid,part=row['orbit_id'],row['part']
        if part=='legs':
            pointer=legs['orbits'][oid];rel='benchmarks/self-energy/tgv-pi1/'+pointer['result']
            data=bound_read(rel,pointer['sha256'])
            field='laurent_Hlog_coefficients' if 'laurent_Hlog_coefficients' in data else 'laurent_coefficients'
            raw=data[field]; actual=multiply_bf(cs(raw))
            old_native=read(frozen+f'benchmarks/self-energy/tgv-pi1/results/{oid}/r002.json')
            assert data['normalization']['euler_compensation_count']==1
            assert data['normalization']['c4_multiplicity']==4
            frozen_values=old_native.get('laurent_Hlog_coefficients',old_native['laurent_coefficients'])
            assert data['normalization']==old_native['normalization']
            for q in set(cs(raw))|set(cs(frozen_values)):
                ck(f'{oid} current consumed raw coefficient {q}',cs(raw).get(q,0)-cs(frozen_values).get(q,0),'CF*CA*TF*nf')
            CHANGES.append(dict(orbit_id=oid,current=rel,consumed_field='/'+field,
                                consumed_field_hash=field_hash(raw),frozen_field_hash=field_hash(frozen_values),
                                all_changed_paths=changed_paths(old_native,data)))
        elif part=='Gamma':
            pointer=gamma['orbits'][oid];selected=pointer.get('preferred_representation',pointer)
            rel='benchmarks/self-energy/gamma1/'+selected['result'];data=bound_read(rel,selected['sha256'])
            field='channels/q/I_coefficients';raw=data['channels']['q']['I_coefficients'];actual=cs(raw)
            prior_pointer=frozen_gamma['orbits'][oid];prior=prior_pointer.get('preferred_representation',prior_pointer)
            prior_rel=frozen+'benchmarks/self-energy/gamma1/'+prior['result'];old_native=read(prior_rel)
            frozen_values=old_native['channels']['q']['I_coefficients']
            assert data['normalization']==old_native['normalization']
            assert data['normalization']['euler_compensation_count']==1
            assert not data['normalization']['external_leg_self_energy_included']
            for q in set(cs(raw))|set(cs(frozen_values)):
                ck(f'{oid} current consumed q coefficient {q}',cs(raw).get(q,0)-cs(frozen_values).get(q,0),'CF*CA*TF*nf')
            CHANGES.append(dict(orbit_id=oid,current=rel,consumed_field='/'+field,
                                consumed_field_hash=field_hash(raw),frozen_field_hash=field_hash(frozen_values),
                                all_changed_paths=changed_paths(old_native,data)))
            # Also inspect default r003: preferred z coefficients are a separate identity bridge.
            default_rel='benchmarks/self-energy/gamma1/'+pointer['result'];default=bound_read(default_rel,pointer['sha256'])
            if pointer['revision']=='r003':
                older=read(default_rel.replace('/r003.json','/r002.json'))
                CHANGES.append(dict(orbit_id=oid,representation='default',current=default_rel,
                                    q_I_same=default['channels']['q']['I_coefficients']==older['channels']['q']['I_coefficients'],
                                    all_changed_paths=changed_paths(older,default)))
        else:
            pointer=native[oid]['result'];rel=pointer['path'];data=bound_read(rel,pointer['sha256'])
            if part=='Pi1':
                field='laurent/coefficients';raw=data['laurent']['coefficients'];actual=multiply_bf(cs(raw))
            elif part=='chain':
                channel=next(x for x in data['color_channels'] if x['monomial']=='CA*CF*TF*nf')
                field='color_channels[CA*CF*TF*nf]/coefficients';raw=channel['coefficients'];actual=cs(raw)
            else:
                channel=next(x for x in data['color_channels'] if x['channel']=='Af')
                field='color_channels[Af]/laurent/coefficients';raw=channel['laurent']['coefficients'];actual=cs(raw)
            assert rel==row['result'] and SOURCES[rel]==row['sha256']
        historical_sign=-1 if part in ['legs','Gamma'] else 1
        for q in range(-6,1):
            ck(f'{oid} current native versus physical legacy {q}',actual.get(q,0)-historical_sign*parse(row['coefficients'][str(q)]),'CF*CA*TF*nf')
            current_sum[q]+=actual.get(q,0)
        out.append(dict(orbit_id=oid,part=part,path=rel,sha256=SOURCES[rel],field='/'+field,
                        field_sha256=field_hash(raw),scalar_multiplier=str(BF) if part in ['Pi1','legs'] else '1',
                        native_physical_sign=1,legacy_physical_sign=historical_sign,
                        consumed_powers=list(range(-6,1))))
    # No r->z Hlog bridge is invented here: historical exact HyperInt normalization is explicit.
    report=read('workflow/checks/shared-sign-production/native-assembly-check.json')
    assert report['status']=='passed' and report['exact_graph_coefficient_checks']==217
    return out, current_sum

def main():
    global ROOT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT,help='Read-only Wilson repository root')
    parser.add_argument('--output-dir',type=Path,default=HERE,help='New output directory; use a fresh directory to preserve published evidence')
    args=parser.parse_args()
    ROOT=args.root.resolve();output=args.output_dir.resolve()
    output.mkdir(parents=True,exist_ok=True)
    current=read('benchmarks/renormalization/three-loop/current.json')
    selected={}
    for color in ['CF*CA*TF*nf','CF^2*TF*nf','CF*(TF*nf)^2']:
        p=current['channels'][color];data=bound_read('benchmarks/renormalization/three-loop/'+p['result'],p['sha256'])
        selected[color]=data[p['field']] if 'field' in p else data
        audit_original_bindings(color,selected[color])
    published=selected['CF*CA*TF*nf']
    for name,expected in published['evidence_sha256'].items():
        rel='workflow/checks/w3-ca-nf/'+name
        assert digest(read_text(rel).encode('utf-8'))==expected or digest((ROOT/rel).read_bytes())==expected, ('historical evidence changed',rel)
    rows,current_sum=check_current_31()
    one=read('benchmarks/renormalization/one-loop/result-r001.json')
    lower=read('benchmarks/renormalization/three-loop/lower-order-inputs-r001.json')
    two=read('benchmarks/renormalization/two-loop/result-r002.json')
    stored_ct=read('benchmarks/renormalization/three-loop/counterterms-r001.json')
    audit_original_bindings('lower w2 through eps1',lower)
    audit_original_bindings('coupling counterterms',stored_ct)
    moments={}
    moment_dir='workflow/checks/single-opposite/runs/20260916T090538777022Z/'
    for name in ['M_0_1','M_1_1','M_2_1','M_3_1','M_1_3','M_1_2']:
        text=read_text(moment_dir+name+'.txt').strip().split('=',1)[1].rstrip(';')
        moments[s.Symbol(name)]=parse(text)
    moments[s.Symbol('N_1_2')]=moments[s.Symbol('M_1_3')]-rho*moments[s.Symbol('M_1_2')]
    w1={q:norm(v.subs({s.Symbol('L_s'):L+rho/2,s.Symbol('L_t'):L-rho/2})) for q,v in cs(one['laurent_unit_coefficients']).items()}
    w2=cs(lower['hat_w2_over_CF'])
    require_budget(w1,-2,2,'w1')
    require_budget(w2,-4,1,'w2')
    b0=(11*CA-4*T)/12;b1=(34*CA**2-20*CA*T-12*CF*T)/48
    ck('beta0 normalization',b0-parse(stored_ct['b0']))
    ck('beta1 normalization',b1-parse(stored_ct['b1']))
    ct={q:norm(-2*b0*w2.get(q+1,0)+b0**2*w1.get(q+2,0)-b1*w1.get(q+1,0)/2) for q in range(-6,1)}
    color_coeff={'CF*CA*TF*nf':lambda x:s.expand(x).coeff(CA).coeff(T),
                 'CF^2*TF*nf':lambda x:s.expand(x).coeff(CF).coeff(T),
                 'CF*(TF*nf)^2':lambda x:s.expand(x).coeff(T,2)}
    for color,extract in color_coeff.items():
        for q in range(-6,1): ck(f'{color} independent coupling CT {q}',extract(ct[q])-parse(stored_ct['color_channels'][color][str(q)]),color)
    native=read('benchmarks/self-energy/current-r002.json');native_map={x['orbit_id']:x for x in native['entries']}
    results={}
    for key,color,family,ids,label in [
        ('nf2','CF*(TF*nf)^2','pi1-chain',['2e-Pi1+Pi1-01','2e-Pi1+Pi1-02'],'CF*(TF*nf)^2'),
        ('cff','CF^2*TF*nf','pi2',['2e-Pi2-01','2e-Pi2-02'],'Ff')]:
        adopted=selected[color]; geom={q:s.S.Zero for q in range(-6,1)};deps=[]
        for oid in ids:
            pointer=native_map[oid]['result'];data=bound_read(pointer['path'],pointer['sha256'])
            if key=='nf2':
                ch=next(x for x in data['color_channels'] if x['monomial']==label);values=ch['coefficients']
            else:
                ch=next(x for x in data['color_channels'] if x['channel']==label);values=ch['laurent']['coefficients']
            assert adopted['source_sha256'][pointer['path']]==pointer['sha256']
            for q,v in values.items():
                # These channels publish both M and Hlog views; consume the
                # same M expression as their original RG and check its bridge.
                expression=parse(v['expression']).subs(L,rho)
                if 'hlog_terms' in v:
                    ck(f'{oid} {label} Hlog-to-M bridge {q}',(parse(v)-expression).subs(moments),color)
                geom[int(q)]+=expression.subs(L,rho)
            deps.append(dict(orbit_id=oid,path=pointer['path'],sha256=pointer['sha256'],field=label,field_sha256=field_hash(values)))
        bare=fullscale(geom,3)
        post={q:norm(bare[q]+color_coeff[color](ct[q])) for q in range(-6,1)}
        for q in range(-6,1):
            ck(f'{key} current-source S3 bare coefficient {q}',bare[q]-parse(adopted[f'bare_w3_{key}'][str(q)]),color)
            ck(f'{key} independently reconstructed post CT {q}',post[q]-parse(adopted[f'post_counterterm_{key}'][str(q)]),color)
        k3=parse(adopted[f'k3_{key}']);h3=parse(adopted[f'h3_{key}'])
        results[color]=dict(bare=bare,post=post,k3=k3,h3=h3,finite=post[0],sources=deps)
    color='CF*CA*TF*nf';adopted=selected[color]
    raw=cs(adopted['raw_unit_coefficients']);bare=fullscale(raw,3)
    # Canonical moment bridge used by the original exact GPL export, re-read here.
    post=cs(adopted['post_counterterm_coefficients'])
    for q in range(-6,1):
        ck(f'Af exact S3 convolution {q}',bare[q]-parse(adopted['bare_full_scale_coefficients'][str(q)]),color)
        ck(f'Af independent coupling CT to published post {q}',bare[q]+color_coeff[color](ct[q]).subs(moments)-post[q],color)
    # Current per-graph coefficients equal the physical historical input above;
    # the r->z sum uses the preserved HyperInt witness, not a new integration.
    results[color]=dict(bare=bare,post=post,k3=parse(adopted['k3']),h3=parse(adopted['h3']),finite=parse(adopted['finite_MN']),sources=rows,
                        raw_sum_bridge='workflow/checks/w3-ca-nf/physical.mpl + physical-coefficients.txt; historical exact HyperInt canonicalization')
    F1=parse(two['F1_over_CF']);F2=parse(two['F2_over_CF']);k2=parse(two['k2'])/CF;h2=parse(two['h2'])/CF
    for color,rec in results.items():
        extract=color_coeff[color];post=rec['post'];k3,h3=rec['k3'],rec['h3']
        for q in [-6,-5]: ck(f'{color} superleading pole {q}',post[q],color)
        b02=extract(b0**2);b0k2=extract(s.expand(b0*k2));b0h2=extract(s.expand(b0*h2));bc1=extract(b1)
        expected={-4:-s.Rational(11,9)*b02,
                  -3:-s.Rational(2,3)*b02*L+s.Rational(5,9)*b0k2+s.Rational(8,9)*bc1,
                  -2:s.Rational(2,3)*(b0k2+bc1)*L+b0h2/3-2*k3/9,
                  -1:-(2*k3*L+h3)/3}
        for q,v in expected.items():ck(f'{color} RG pole template {q}',post[q]-v,color)
        ck(f'{color} cusp pole extraction',-s.Rational(3,2)*s.diff(post[-1],L)-k3,color)
        ck(f'{color} noncusp pole extraction',-3*post[-1].subs(L,0)-h3,color)
        lower_rg=extract(s.expand(2*b0*F2+b1*F1))
        ck(f'{color} finite RG',s.diff(rec['finite'],L)-lower_rg+2*k3*L+h3,color)
        ck(f'{color} cusp independent geometry',s.diff(k3,rho),color)
        ck(f'{color} noncusp independent geometry',s.diff(h3,rho),color)
    ck('nf2 registered cusp reference',results['CF*(TF*nf)^2']['k3']+s.Rational(1,27))
    ck('cff registered cusp reference',results['CF^2*TF*nf']['k3']-s.zeta(3)+s.Rational(55,48))
    ck('Af registered cusp reference',results['CF*CA*TF*nf']['k3']+s.Rational(209,216)-5*s.zeta(2)/9+7*s.zeta(3)/6)
    # Concrete controls: wrong source sign, wrong common scale, positive-order omission.
    leg=sum(parse(x['coefficients']['-5']) for x in read('workflow/checks/w3-ca-nf/inputs.json')['parts'] if x['part'] in ['legs','Gamma'])
    nonzero('negative: apply historical minus sign again to native I',2*leg)
    unit1=strip_scale(w1,1,-2,2);unit2=strip_scale(w2,2,-4,1)
    wrong_w1=fullscale(unit1,3,-2,2);wrong_w2=fullscale(unit2,3,-4,1)
    nonzero('negative: replace S1 by S3 in one-loop finite counterterm',color_coeff['CF*(TF*nf)^2'](b0**2*(wrong_w1[2]-w1[2])))
    nonzero('negative: replace S2 by S3 in two-loop finite counterterm',color_coeff['CF*CA*TF*nf'](-2*b0*(wrong_w2[1]-w2[1])))
    nonzero('negative: omit w1 eps2 finite counterterm',color_coeff['CF*(TF*nf)^2'](b0**2*w1[2]))
    nonzero('negative: omit w2 eps1 finite counterterm',color_coeff['CF*CA*TF*nf'](-2*b0*w2[1]))
    for source,lo,hi,label in [(w1,-2,2,'w1'),(w2,-4,1,'w2')]:
        missing={q:v for q,v in source.items() if q!=hi}
        try:
            require_budget(missing,lo,hi,label)
        except ValueError as exc:
            CHECKS.append(dict(name=f'negative: {label} missing positive order rejected',passed=True,negative_control=True,message=str(exc)))
        else:
            raise AssertionError('missing positive-order budget was accepted')
    extra1={**w1,3:s.Symbol('unused_high_w1')};extra2={**w2,2:s.Symbol('unused_high_w2')}
    for q in range(-6,1):
        unchanged=-2*b0*extra2.get(q+1,0)+b0**2*extra1.get(q+2,0)-b1*extra1.get(q+1,0)/2
        ck(f'unconsumed higher orders do not change CT {q}',unchanged-ct[q])
    # A display/numeric-only change must leave consumed coefficient identity intact.
    for item in CHANGES:
        if item.get('representation')=='default': assert item['q_I_same']
    # Freeze every actually consumed file and reject concurrent modifications.
    for rel,h in SOURCES.items(): assert digest((ROOT/rel).read_bytes())==h, ('input changed before publication',rel)
    save=lambda name,obj:(output/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
    save('inputs.json',dict(format='wilson-matter-rg-current-inputs/1',root=str(ROOT),source_sha256=SOURCES,graph_consumption=rows,historical_bindings=LINEAGE))
    save('revision-fields.json',dict(revisions=CHANGES))
    encoded={color:{k:({str(q):str(v) for q,v in value.items()} if k in ['bare','post'] else str(value) if isinstance(value,s.Basic) else value) for k,value in rec.items()} for color,rec in results.items()}
    save('current-coefficients.json',encoded)
    summary=dict(status='passed',date='2026-10-17',graphs=31,checks=len(CHECKS),negative_controls=sum(bool(x.get('negative_control')) for x in CHECKS),
                 colors=list(results),w1_maximum_consumed_epsilon_power=2,w2_maximum_consumed_epsilon_power=1,w3_maximum_consumed_epsilon_power=0,
                 current_analytic_reuse=True,historical_sum_canonicalization_reused=True,formal='not_run',
                 complete_w3=False,independent_h3_verification=False,independent_finite_shape_verification=False,
                 elapsed_seconds=round(time.time()-START,3),python_sympy_version=s.__version__,checks_detail=CHECKS)
    save('checks.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k!='checks_detail'},ensure_ascii=True),flush=True)

if __name__=='__main__':
    main()
