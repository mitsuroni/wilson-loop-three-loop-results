"""WILSON-4.46.6: fresh accepted AA known subtotal, no A/B finite import.

All output is restricted to this directory. Historical builders are imported
only for pure parser/arithmetic functions, never executed via main().
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys
import time
from fractions import Fraction as Q
sys.dont_write_bytecode=True
import sympy as s

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]
POWERS=tuple(range(-6,1));NEGATIVE=tuple(range(-6,0))
HASHES={};FIELDS=[];CHECKS=[]
BLOCKED=('4e-4g-05/','4e-4g-08/','4e-4g-09/','4e-3g+3g-13/',
         '4e-3g+3g-15/','4e-3g+3g-22/','4e-3g+3g-24/','4e-3g+3g-25/')


def raw(path,expected=None):
    path=str(path).replace('\\','/')
    data=(ROOT/path).read_bytes();h=hashlib.sha256(data).hexdigest()
    assert path not in HASHES or HASHES[path]==h,('input changed',path)
    assert expected is None or h==expected,('hash mismatch',path)
    HASHES[path]=h;return data


def read(path,expected=None):
    assert not any(part in str(path) for part in BLOCKED),('whole unknown source read forbidden',path)
    return json.loads(raw(path,expected).decode('utf-8-sig'))


def bind_field(path,pointer,value):
    row=dict(path=path,sha256=HASHES[path],json_pointer=pointer,
             field_sha256=hashlib.sha256(json.dumps(value,sort_keys=True,ensure_ascii=False,separators=(',',':')).encode()).hexdigest())
    FIELDS.append(row);return row


def selective(path,pointer,expected=None):
    """Navigate JSON without deserializing unselected values, notably eps^0.

    Full bytes are hash-bound. Only key strings and the requested negative
    coefficient value are decoded; other values are skipped lexically.
    """
    if any(part in str(path) for part in BLOCKED):
        assert re.search(r'/(?:eps\^)?-[1-6](?:/|$)',pointer),('nonnegative unknown-field read forbidden',path,pointer)
    text=raw(path,expected).decode('utf-8-sig');decoder=json.JSONDecoder()
    def ws(i):
        while i<len(text) and text[i].isspace():i+=1
        return i
    def skip(i):
        i=ws(i)
        if text[i]=='"':
            i+=1
            while True:
                if text[i]=='\\':i+=2
                elif text[i]=='"':return i+1
                else:i+=1
        if text[i] in '[{':
            end=']' if text[i]=='[' else '}';i+=1
            while True:
                i=ws(i)
                if text[i]==end:return i+1
                if text[i] in ',:':i+=1
                else:i=skip(i)
        while i<len(text) and text[i] not in ',]}\r\n\t ':i+=1
        return i
    i=0
    for key in pointer.strip('/').split('/'):
        i=ws(i)
        if text[i]!='{':raise KeyError(pointer)
        i+=1;found=False
        while True:
            i=ws(i)
            if text[i]=='}':break
            name,j=decoder.raw_decode(text,i);j=ws(j)
            assert text[j]==':';j=ws(j+1)
            if name==key:i=j;found=True;break
            i=ws(skip(j))
            if text[i]==',':i+=1
        if not found:raise KeyError(pointer)
    value,end=decoder.raw_decode(text,i)
    bind_field(path,pointer,value);return value


def load_module(name,path):
    raw(path);spec=importlib.util.spec_from_file_location(name,ROOT/path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return module


base=load_module('accepted_pole_insertion_helpers','workflow/checks/aa-joint-rg-20261024/insertions/build_insertions.py')
math=base.math
raw('workflow/checks/matter-rg-20261017/check.py')
assert math.ROOT.resolve()==ROOT.resolve(),('unexpected parser input root',math.ROOT)
eps,rho,z=math.eps,math.rho,math.z
a=(10-6*eps)/(3-2*eps)
ha=load_module('known_hlog_algebra','workflow/checks/single-opposite/hlog_algebra.py')


def dump(name,value):
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/name).write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')


def check(name,left,right=0,negative=False):
    value=s.cancel(math.shuffle_normalize(s.expand(left-right)))
    ok=value!=0 if negative else value==0
    assert ok,(name,str(value)[:1200]);CHECKS.append(dict(name=name,negative=negative,residual=str(value)))


def parse(value):
    if isinstance(value,list):
        return s.expand(sum(parse(t['coefficient'])*math.hl('z',t['word'])*
                            s.prod(s.zeta(n)**p for n,p in zip((2,3,5),t.get('zeta_powers',(0,0,0)))) for t in value))
    if isinstance(value,dict) and 'hlog_terms' in value:return parse(value['hlog_terms'])
    return base.parse(value)


def vec_expr(rows):return parse(rows)


def normalize(value):
    value=parse(value).subs({s.Symbol('ell'):rho,s.Symbol('chi'):(1-z)/(1+z)})
    value=value.subs(MOMENTS,simultaneous=True)
    value=s.cancel(math.shuffle_normalize(s.expand(value)))
    assert not any(str(x).startswith('M_') or str(x).startswith('N_') or str(x).startswith('H_r_') for x in value.free_symbols),str(value)[:1500]
    assert value.free_symbols <= {z,rho,*[v for v in value.free_symbols if str(v).startswith('H_z_')]},str(value.free_symbols)
    return value


def series_values(data,fields,pointer,path,highest=None):
    out={};origins={}
    for p in POWERS:
        if str(p) in fields:
            val=fields[str(p)];origins[str(p)]=bind_field(path,pointer+'/'+str(p),val);out[p]=parse(val)
        elif highest is not None and p < -highest:
            out[p]=s.S.Zero;origins[str(p)]=dict(zero_reason='below explicit highest_pole_order',highest_pole_order=highest)
        elif str(p) in data.get('vanishing_higher_poles',{}):
            val=data['vanishing_higher_poles'][str(p)];out[p]=parse(val);assert out[p]==0
            origins[str(p)]=bind_field(path,'/vanishing_higher_poles/'+str(p),val)
        elif p in data.get('zero_powers',[]):out[p]=s.S.Zero;origins[str(p)]=dict(zero_reason='explicit zero_powers')
        else:raise ValueError(('missing coefficient, not zero',path,p))
    return out,origins


def convolve(values,multiplier):
    # Absolute epsilon budget: native -6..0 with analytic multiplier to +6.
    series=s.series(multiplier,eps,0,7).removeO().expand()
    result={p:s.expand(sum(v*series.coeff(eps,p-q) for q,v in values.items() if q<=p)) for p in POWERS}
    # Shift the Laurent polynomial before coefficient extraction; otherwise
    # SymPy may combine eps*(1+z) into a denominator invisible to .coeff.
    poly=s.Poly(s.expand(sum(v*eps**(q+6) for q,v in values.items())*series),eps)
    for p in POWERS:check('absolute-power convolution '+str(p),result[p],poly.coeff_monomial(eps**(p+6)))
    return result


def extract_coefficient(path,p,expected=None,unknown=False):
    pointers=[f'/whole_I_laurent/coefficients/{p}',f'/whole_I_coefficients/eps^{p}',f'/whole_I_coefficients/{p}',
              f'/laurent/coefficients/{p}',f'/laurent_coefficients/{p}']
    for pointer in pointers:
        try:return parse(selective(path,pointer,expected)),pointer
        except KeyError:pass
    raise KeyError(('no accepted coefficient field',path,p))


def expression_vector(value):
    """Common z HLOG rows with exact rational-z coefficients and periods."""
    value=normalize(value);symbols=sorted((v for v in value.free_symbols if str(v).startswith('H_z_')),key=str)
    rows=[]
    for atom in symbols:
        coefficient=s.expand(value).coeff(atom)
        assert rho not in coefficient.free_symbols
        word=[-1 if t=='m1' else int(t) for t in str(atom)[4:].split('_')]
        rows.append(dict(word=word,coefficient=str(coefficient)))
    rest=s.expand(value.subs(dict.fromkeys(symbols,0)))
    for (n,),coefficient in s.Poly(rest,rho).terms():
        if coefficient:rows.append(dict(word=[0]*n,coefficient=str(coefficient*s.factorial(n))))
    out={}
    Z2,Z3,Z5=s.symbols('Z2 Z3 Z5')
    for row in rows:
        coefficient=parse(row['coefficient']).xreplace({s.pi**6:216*Z2**3,s.pi**4:36*Z2**2,s.pi**2:6*Z2,
                         s.zeta(3):Z3,s.zeta(5):Z5})
        coefficient=s.expand(coefficient)
        for powers,c in s.Poly(coefficient,Z2,Z3,Z5).terms():
            c=s.cancel(c)
            if c:out[(tuple(row['word']),powers)]=out.get((tuple(row['word']),powers),0)+c
    return {k:s.cancel(c) for k,c in out.items() if c}


def serialize_vec(vec):
    return [dict(word=list(w),zeta_powers=list(ps),coefficient=str(c)) for (w,ps),c in sorted(vec.items()) if c]


def main():
    global MOMENTS
    started=time.monotonic();raw(str(Path(__file__).relative_to(ROOT)).replace('\\','/'))
    dictionary=read('workflow/checks/four-leg-bootstrap-20261101/basis/dictionary.json')
    MOMENTS={s.Symbol(row['name']):vec_expr(row['hlog_terms']) for row in dictionary['atoms']}
    MOMENTS[s.Symbol('ell')]=rho
    K4=read('workflow/checks/four-leg-bootstrap-20261101/basis/accepted-k4.json')
    effective={name:read(f'benchmarks/recursion-audit/results/{name}-effective-sum/r001.json') for name in ('6e','5e')}
    native=[];families={}
    seouter=read('benchmarks/self-energy/current.json')
    se=read('benchmarks/self-energy/'+seouter['manifest'],seouter['sha256'])
    legs=read('benchmarks/self-energy/tgv-pi1/current.json')
    gamma=read('benchmarks/self-energy/gamma1/current.json')
    for sel in se['entries']:
        family={'4e-Pi1':'Pi1','2e-Pi1+Pi1':'chain','2e-Pi2':'Pi2'}.get(sel['family'])
        if family is None:continue
        path=sel['result']['path'];data=read(path,sel['result']['sha256'])
        assert data['normalization']['assembly_sign']==sel['assembly_sign']==1
        assert data['normalization']['Euler_count']==1
        if family=='Pi1':pointer='/laurent/coefficients';fields=data['laurent']['coefficients'];multiplier=a
        else:
            idx=next(i for i,c in enumerate(data['color_channels']) if c.get('monomial',c.get('color_monomial'))=='CA^2*CF')
            ch=data['color_channels'][idx];pointer=f'/color_channels/{idx}/'+('coefficients' if family=='chain' else 'laurent/coefficients')
            fields=ch['coefficients'] if family=='chain' else ch['laurent']['coefficients'];multiplier=s.S.One
        vals,origins=series_values(data,fields,pointer,path)
        weighted=convolve(vals,multiplier)
        if family=='chain':
            scalars,_=series_values(data,data['laurent']['coefficients'],'/laurent/coefficients',path)
            expected=convolve(scalars,a*a)
            for p in POWERS:check(sel['orbit_id']+' stored AA a^2 '+str(p),normalize(weighted[p]),normalize(expected[p]))
        for contract in ('raw','prepared'):raw(sel[contract]['path'],sel[contract]['sha256'])
        native.append(dict(orbit_id=sel['orbit_id'],family=family,path=path,revision=data.get('revision'),
              sha256=HASHES[path],coefficient_sources=origins,weight=str(multiplier),
              coefficients={str(p):str(normalize(v)) for p,v in weighted.items()}))
    for family,selector,folder in [('legs',legs,'tgv-pi1'),('Gamma',gamma,'gamma1')]:
        for oid,sel in selector['orbits'].items():
            chosen=sel.get('preferred_representation',sel);path=f'benchmarks/self-energy/{folder}/'+chosen['result']
            data=read(path,chosen['sha256']);norm=data['normalization'];assert norm['assembly_sign']==1 and norm['euler_compensation_count']==1
            raw(data.get('input',{}).get('native_package',f'benchmarks/self-energy/{folder}/current.json'),data.get('input',{}).get('native_package_sha256'))
            if family=='Gamma':
                fields=data['channels']['g']['I_coefficients'];pointer='/channels/g/I_coefficients'
                vals,origins=series_values(data,fields,pointer,path,data['channels']['g']['highest_pole_order']);multiplier=s.S.One
                if chosen is not sel:
                    assert data['mathematical_result_sha256']==sel['sha256']
                    raw(f'benchmarks/self-energy/{folder}/'+sel['result'],sel['sha256'])
                    for leaf,h in data['legacy_evidence']['evidence_sha256'].items():raw('workflow/checks/gamma1-03-z/'+leaf,h)
            else:
                vals,origins=series_values(data,data['laurent_coefficients'],'/laurent_coefficients',path,data['highest_pole_order']);multiplier=a
                if oid.endswith(('07','08')):
                    vals[-2]=parse(data['second_pole_in_existing_M_basis']);origins['-2']=bind_field(path,'/second_pole_in_existing_M_basis',data['second_pole_in_existing_M_basis'])
                    if oid.endswith('07'):
                        alternate='benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-07/mn-r002.json';mn=read(alternate)
                        assert mn['target_result_sha256']==HASHES[path]
                        assert mn['target_source_sha256']==data['input']['native_package_sha256']
                        for p in (-1,0):vals[p]=parse(mn['coefficients'][str(p)]);origins[str(p)]=bind_field(alternate,'/coefficients/'+str(p),mn['coefficients'][str(p)])
                        raw('workflow/checks/tgv-pi1-representation-r003/verification.json')
                    else:
                        for p in (-1,0):
                            val=data['M_representation']['coefficients'][str(p)];vals[p]=parse(val);origins[str(p)]=bind_field(path,'/M_representation/coefficients/'+str(p),val)
            weighted=convolve(vals,multiplier)
            if family=='legs' and 'color_channel_coefficients' in data:
                aa=data['color_channel_coefficients'].get('CF*CA^2',{})
                for p in POWERS:
                    if str(p) in aa and 'Hlog(r' not in str(aa[str(p)]):
                        check(oid+' accepted AA convolution '+str(p),normalize(weighted[p]),normalize(aa[str(p)]))
            native.append(dict(orbit_id=oid,family=family,path=path,revision=chosen['revision'],sha256=HASHES[path],
                  coefficient_sources=origins,weight=str(multiplier),coefficients={str(p):str(normalize(v)) for p,v in weighted.items()}))
    expected={'Pi1':11,'chain':2,'Pi2':2,'legs':12,'Gamma':4}
    for family,count in expected.items():
        group=[r for r in native if r['family']==family];assert len(group)==count
        families[family]={p:normalize(sum(parse(row['coefficients'][str(p)]) for row in group)) for p in POWERS}
        print('assembled',family,count,'finite HLOG coordinates',len(expression_vector(families[family][0])),flush=True)
    assert len(native)==31
    insertions={p:normalize(sum(v[p] for v in families.values())) for p in POWERS}
    old=read('workflow/checks/aa-joint-rg-20261024/insertions/insertions.json')
    for p in range(-6,-1):
        check('accepted historical insertion pole '+str(p),insertions[p],normalize(old['total_coefficients'][str(p)]))
        for family in families:
            check('accepted historical '+family+' pole '+str(p),families[family][p],normalize(old['families'][family]['coefficients'][str(p)]))
    dump('insertions.json',dict(status='exact31_source_assembly_through_finite',formal='not_run',
         diagrams=native,families={f:{str(p):str(v) for p,v in vals.items()} for f,vals in families.items()},
         total={str(p):str(v) for p,v in insertions.items()},input_sha256=HASHES,coefficient_sources=FIELDS,
         normalization='a(eps) only on scalar Pi1/legs; stored chain/Pi2/Gamma AA fields once; absolute analytic budget through eps6'))
    k4values={0:normalize(K4['hlog_coefficients'])}
    for p in NEGATIVE:
        terms=[]
        for item in K4['sources']:
            val,pointer=extract_coefficient(item['path'],p,item['sha256'])
            terms.append((2 if item['group']=='I3444' else 1)*val)
        k4values[p]=normalize(sum(terms))
    # Consume the explicit HLOG fields, avoiding any accidental interpretation
    # of the project function E as SymPy's Euler-number constant.
    trees={name:{p:normalize(data['hlog_coefficients'][str(p)]) for p in POWERS} for name,data in effective.items()}
    known={p:normalize(k4values[p]+trees['6e'][p]+trees['5e'][p]+insertions[p]) for p in POWERS}
    dump('known-finite.json',dict(status='exact_accepted_finite_known_subtotal',formal='not_run',
         definition='Kknown[0]=K4[0]+6e[0]+5e[0]+31acceptedInsertions[0]; excludes whole2A+B finite objects',
         Kknown_finite=str(known[0]),Kknown_finite_hlog=serialize_vec(expression_vector(known[0])),
         unknown_finite_fields_read=False,input_sha256=HASHES))
    # The five uncertain-source files are traversed only by selective negative
    # JSON pointers. Their finite fields are never deserialized or evaluated.
    selectors={'FGV':read('benchmarks/fgv/4e-4g/current.json'),'DT':read('benchmarks/double-tgv/4e-3g+3g/current.json')}
    unknown_groups={'A':[('FGV','05'),('DT','13'),('DT','15')],'B':[('FGV','09'),('DT','25')]}
    unknown_poles={}
    for group,ids in unknown_groups.items():
        vv={p:[] for p in NEGATIVE}
        for family,idx in ids:
            oid=('4e-4g-' if family=='FGV' else '4e-3g+3g-')+idx;sel=selectors[family]['orbits'][oid]
            folder='benchmarks/fgv/4e-4g/' if family=='FGV' else 'benchmarks/double-tgv/4e-3g+3g/'
            path=folder+sel['result']
            for p in NEGATIVE:
                val,pointer=extract_coefficient(path,p,sel['sha256'],unknown=True);vv[p].append(val)
        unknown_poles[group]={p:normalize(sum(v)) for p,v in vv.items()}
    full_poles={p:normalize(known[p]+2*unknown_poles['A'][p]+unknown_poles['B'][p]) for p in NEGATIVE}
    control=read('workflow/checks/aa-joint-rg-20261024/coefficients.json')
    # Only the historical independent known-pole fields are used as controls.
    for p in range(-6,-1):check('full accepted historical U3 pole '+str(p),full_poles[p],normalize(control['scale_stripped_bare'][str(p)]))
    pirow=next(row for row in native if row['family']=='Pi1')
    check('negative: duplicate Pi1 a factor changes finite',normalize(a.subs(eps,0)*parse(pirow['coefficients']['0'])),normalize(pirow['coefficients']['0']),negative=True)
    grow=next(row for row in native if row['family']=='Gamma')
    check('negative: extra legacy Gamma sign changes finite',normalize(grow['coefficients']['0']),-normalize(grow['coefficients']['0']),negative=True)
    blocked_path='benchmarks/double-tgv/4e-3g+3g/results/4e-3g+3g-25/ibp-progress-r005.json'
    try:selective(blocked_path,'/whole_I_coefficients/0')
    except AssertionError as exc:
        assert exc.args[0][0]=='nonnegative unknown-field read forbidden'
        CHECKS.append(dict(name='negative: unknown finite pointer rejected before IO',negative=True,residual='rejected'))
    else:raise AssertionError('finite firewall failed')
    for path,h in HASHES.items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h,('end input change',path)
    dump('assembly.json',dict(status='assembled_exact_accepted_known_subtotal',formal='not_run',color='CF*CA^2',
        quantity='scale-stripped bare connected U3',definition='U3[0]=Kknown[0]+H; H=2 A+B; A=I2344[0], B=I1234[0]',
        normalization='All I include native sign, C4 and one Euler factor. a(eps) applied only to scalar Pi1/legs; chain/Pi2/Gamma colored fields once.',
        insertion_diagrams=native,insertion_family_coefficients={f:{str(p):str(v) for p,v in vals.items()} for f,vals in families.items()},
        insertion_total={str(p):str(v) for p,v in insertions.items()},
        Kknown_coefficients={str(p):str(v) for p,v in known.items()},Kknown_finite_hlog=serialize_vec(expression_vector(known[0])),
        full_negative_U3={str(p):str(v) for p,v in full_poles.items()},
        full_negative_U3_hlog={str(p):serialize_vec(expression_vector(v)) for p,v in full_poles.items()},
        accepted_unknown_object_poles={g:{str(p):str(v) for p,v in vv.items()} for g,vv in unknown_poles.items()},
        scale_transport='[eps0] S3 U3 = 2 H+2 Kknown[0]+sum_{p<0} S3[-p] U3[p]. Counterterm CT0 is separate.',
        unknown_finite_fields_read=False,new_original_integral_evaluation=False,N4_predictions_used=False,
        coefficient_sources=FIELDS,input_sha256=HASHES,checks=CHECKS,end_input_stability=True,elapsed_seconds=time.monotonic()-started))
    print('passed known assembly; checks',len(CHECKS),'inputs',len(HASHES),'elapsed',round(time.monotonic()-started,2),flush=True)


if __name__=='__main__':main()
