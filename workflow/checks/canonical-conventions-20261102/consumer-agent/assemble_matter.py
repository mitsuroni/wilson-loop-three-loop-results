"""Direct native-current matter assembly, with no legacy sign or total inputs.

Existing exact scalar-view and moment dictionaries are reused. All old totals
are read only after direct native U, S3, CT and post series have been formed.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,sys,time
sys.dont_write_bytecode=True
import sympy as s
H=Path(__file__).resolve().parent;R=H.parents[3]
INPUTS={};CHECKS=[];ROWS=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(rel):
    raw=(R/rel).read_bytes();digest=hashlib.sha256(raw).hexdigest()
    assert rel not in INPUTS or INPUTS[rel]==digest,('input changed',rel)
    INPUTS[rel]=digest;return json.loads(raw.decode('utf-8-sig'))
def module(name,rel):
    path=R/rel;INPUTS[rel]=sha(path)
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
k=module('native_current_matter_helpers','workflow/checks/four-leg-bootstrap-20261101/known-total/assemble_known.py')
assert k.ROOT.resolve()==R.resolve() and k.math.ROOT.resolve()==R.resolve()
eps,rho,L,CA,CF,T=k.eps,k.rho,k.math.L,k.math.CA,k.math.CF,k.math.T
BF=-8*(1-eps)/(3-2*eps);POWERS=range(-6,1)
COLORS={'Af':CA*T,'Ff':CF*T,'ff':T*T}

def physical(text):return k.math.parse(text)
def normL(value):
    value=s.expand(value);poly=s.Poly(value,L)
    return s.expand(sum(k.normalize(c)*L**monomial[0] for monomial,c in poly.terms()))
def check(name,value,negative=False):
    value=normL(value);assert value!=0 if negative else value==0,(name,str(value)[:1500])
    CHECKS.append({'name':name,'negative_control':negative,'residual':str(value)})
def series(fields,highest,oid):
    out={}
    for p in POWERS:
        if str(p) in fields:out[p]=k.parse(fields[str(p)])
        elif highest is not None and p < -highest:out[p]=s.S.Zero
        else:raise ValueError(('Unknown coefficient, not zero',oid,p))
    return out
def consume(label,oid,path,fields,pointer,weight,highest,view_bridges=None):
    values=series(fields,highest,oid);ws=s.series(weight,eps,0,7).removeO().expand()
    weighted={p:k.normalize(sum(v*ws.coeff(eps,p-q) for q,v in values.items() if q<=p)) for p in POWERS}
    origins={}
    for p in POWERS:
        bridge=next((x for x in (view_bridges or []) if x['power']==p),None)
        target=bridge['path'] if bridge else path
        field=bridge['pointer'] if bridge else pointer+'/'+str(p)
        if str(p) in fields:
            actual=read(target)
            for item in field.strip('/').split('/'):
                actual=actual[int(item)] if isinstance(actual,list) else actual[item]
            assert actual==fields[str(p)],(oid,p,field)
            origins[str(p)]={'path':target,'sha256':INPUTS[target],'json_pointer':field,
              'field_sha256':hashlib.sha256(json.dumps(actual,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()).hexdigest()}
        else:
            origins[str(p)]={'path':path,'sha256':INPUTS[path],'zero_reason':'below explicit highest_pole_order','highest_pole_order':highest}
    ROWS.append({'color':label,'orbit_id':oid,'path':path,'sha256':INPUTS[path],
      'source_field':pointer,'selected_field_sha256':hashlib.sha256(json.dumps(fields,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
      'physical_convention':'wilson-qcd-action-path-fourier/1.0','assembly_sign':1,
      'weight':str(weight),'view_bridges':view_bridges or [],'source_fields_by_power':origins,
      'coefficients':{str(p):str(v) for p,v in weighted.items()}})
    print('read',label,oid,flush=True)
    return weighted
def encode(values):return {str(p):str(v) for p,v in values.items()}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--lower-input',type=Path,default=Path('benchmarks/current/low-loops/coefficients.json'))
    a=ap.parse_args();out=a.output_dir.resolve()
    allowed=(not out.is_relative_to(R.resolve()) or out.is_relative_to(H.resolve()) or
      any(out.is_relative_to((R/name).resolve()) for name in ['.wilson','replay-output']))
    assert allowed and not out.exists(),'Use a new output directory in replay-output, .wilson, this evidence directory, or outside the checkout.'
    out.mkdir(parents=True);tick=time.monotonic();INPUTS[Path(__file__).relative_to(R).as_posix()]=sha(Path(__file__))
    dictionary=read('workflow/checks/four-leg-bootstrap-20261101/basis/dictionary.json')
    k.MOMENTS={s.Symbol(x['name']):k.vec_expr(x['hlog_terms']) for x in dictionary['atoms']};k.MOMENTS[s.Symbol('ell')]=rho
    native=read('benchmarks/self-energy/current-r002.json')
    for row in native['entries']:
        fam=row['family']
        if fam not in ['4e-Pi1','2e-Pi1+Pi1','2e-Pi2']:continue
        path=row['result']['path'];data=read(path);assert INPUTS[path]==row['result']['sha256']
        assert data['normalization']['assembly_sign']==row['assembly_sign']==1
        if fam=='4e-Pi1':
            consume('Af',row['orbit_id'],path,data['laurent']['coefficients'],'/laurent/coefficients',BF,data.get('highest_pole_order'))
        elif fam=='2e-Pi1+Pi1':
            for label,monomial in [('Af','CA*CF*TF*nf'),('ff','CF*(TF*nf)^2')]:
                i=next(i for i,x in enumerate(data['color_channels']) if x['monomial']==monomial)
                consume(label,row['orbit_id'],path,data['color_channels'][i]['coefficients'],f'/color_channels/{i}/coefficients',s.S.One,data.get('highest_pole_order'))
        else:
            for label in ['Af','Ff']:
                i=next(i for i,x in enumerate(data['color_channels']) if x['channel']==label)
                consume(label,row['orbit_id'],path,data['color_channels'][i]['laurent']['coefficients'],f'/color_channels/{i}/laurent/coefficients',s.S.One,data.get('highest_pole_order'))
    for fam,folder in [('legs','tgv-pi1'),('Gamma','gamma1')]:
        selectors=read(f'benchmarks/self-energy/{folder}/current.json')
        assert selectors['physical_convention']=='wilson-qcd-action-path-fourier/1.0'
        for oid,row in selectors['orbits'].items():
            chosen=row.get('preferred_representation',row);path=f'benchmarks/self-energy/{folder}/'+chosen['result'];data=read(path)
            assert INPUTS[path]==chosen['sha256'] and data['normalization']['assembly_sign']==1
            bridges=[]
            if fam=='Gamma':
                fields=data['channels']['q']['I_coefficients'];ptr='/channels/q/I_coefficients';weight=s.S.One;highest=data['channels']['q']['highest_pole_order']
            else:
                fields=dict(data['laurent_coefficients']);ptr='/laurent_coefficients';weight=BF;highest=data['highest_pole_order']
                if oid.endswith(('07','08')):
                    fields['-2']=data['second_pole_in_existing_M_basis']
                    bridges.append({'power':-2,'path':path,'pointer':'/second_pole_in_existing_M_basis'})
                    if oid.endswith('07'):
                        alt='benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-07/mn-r002.json';mn=read(alt)
                        assert mn['target_result_sha256']==INPUTS[path] and mn['target_source_sha256']==data['input']['native_package_sha256']
                        for p in [-1,0]:fields[str(p)]=mn['coefficients'][str(p)];bridges.append({'power':p,'path':alt,'pointer':f'/coefficients/{p}'})
                    else:
                        for p in [-1,0]:fields[str(p)]=data['M_representation']['coefficients'][str(p)];bridges.append({'power':p,'path':path,'pointer':f'/M_representation/coefficients/{p}'})
            consume('Af',oid,path,fields,ptr,weight,highest,bridges)
    assert len({x['orbit_id'] for x in ROWS})==31 and len(ROWS)==35
    assert {c:sum(x['color']==c for x in ROWS) for c in COLORS}=={'Af':31,'Ff':2,'ff':2}
    U={c:{p:k.normalize(sum(physical(x['coefficients'][str(p)]) for x in ROWS if x['color']==c)) for p in POWERS} for c in COLORS}
    lowerpath=(R/a.lower_input).resolve()
    lower=read(lowerpath.relative_to(R.resolve()).as_posix() if lowerpath.is_relative_to(R.resolve()) else str(lowerpath))
    w1={int(p):physical(v) for p,v in lower['one_loop']['bare'].items()}
    w2={int(p):physical(v) for p,v in lower['two_loop']['bare'].items()}
    assert set(range(-2,3))<=set(w1) and set(range(-4,2))<=set(w2)
    b0=(11*CA-4*T)/12;b1=(34*CA**2-20*CA*T-12*CF*T)/48
    ct={p:s.expand(-2*b0*w2.get(p+1,0)+b0**2*w1.get(p+2,0)-b1*w1.get(p+1,0)/2) for p in POWERS}
    coeff=lambda v,c:s.Poly(s.expand(v),CA,CF,T).coeff_monomial(COLORS[c])
    scale={p:s.expand(3**p*((L+rho/2)**p+(L-rho/2)**p)/s.factorial(p)) for p in range(7)}
    results={}
    for c in COLORS:
        bare={p:normL(sum(scale[p-q]*v for q,v in U[c].items() if q<=p)) for p in POWERS}
        CT={p:normL(coeff(v,c)) for p,v in ct.items()};post={p:normL(bare[p]+CT[p]) for p in POWERS}
        cusp=normL(-s.Rational(3,2)*s.diff(post[-1],L));noncusp=normL(-3*post[-1].subs(L,0))
        assert not cusp.free_symbols and not noncusp.free_symbols,(c,cusp,noncusp)
        rg=normL(s.diff(post[0],L)-coeff(2*b0*physical(lower['two_loop']['F'])+b1*physical(lower['one_loop']['F']),c)+2*cusp*L+noncusp)
        check(c+' direct finite RG',rg)
        for p in [-6,-5]:check(c+' superleading cancellation '+str(p),post[p])
        results[c]={'U':encode(U[c]),'bare':encode(bare),'counterterm':encode(CT),'post_counterterm':encode(post),'F':str(post[0]),'k3':str(cusp),'h3':str(noncusp),'finite_RG_residual':str(rg)}
    # Freeze the direct construction before opening any old three-loop target.
    candidate={'channels':results,'source_rows':ROWS,'inputs_at_candidate_freeze':dict(INPUTS),'old_three_loop_targets_read':False}
    (out/'constructed-before-holdout.json').write_text(json.dumps(candidate,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    target=read('benchmarks/renormalization/three-loop/first-complete-20261102/coefficients.json')
    for c,rec in results.items():
        old=target['channels'][c]
        for stage in ['U','bare','counterterm','post_counterterm']:
            for p,value in rec[stage].items():check(c+' final-stage '+stage+'['+p+']',physical(value)-physical(old[stage][p]))
        for field in ['k3','h3']:check(c+' extracted '+field+' versus accepted coefficient',physical(rec[field])-physical(old[field]))
    check('negative doubled native leg sign',2*physical(next(x for x in ROWS if x['orbit_id']=='3e-3g+Pi1-01')['coefficients']['-4']),True)
    pi=next(x for x in ROWS if x['orbit_id']=='4e-Pi1-01')
    bseries=s.series(BF,eps,0,7).removeO().expand()
    duplicated=sum(physical(v)*bseries.coeff(eps,-int(p)) for p,v in pi['coefficients'].items())
    check('negative duplicate scalar Bf on Pi1',duplicated-physical(pi['coefficients']['0']),True)
    for rel,digest in INPUTS.items():assert sha(R/rel)==digest,('input changed during run',rel)
    record={'format':'wilson-current-native-matter-assembly/1','status':'passed','task':'WILSON-6.3','scope':'31 original native orbit IDs,35 color contributions; direct current scalar fields into all seven U/bare/CT/post coefficients. Exact accepted representations and moment identities reused; no fresh integral/CAS claim.',
        'convention':'wilson-qcd-action-path-fourier/1.0','native_assembly_sign':1,'legacy_family_sign_applied':False,
        'color_monomials_over_CF':{c:str(v) for c,v in COLORS.items()},'channels':results,'source_rows':ROWS,'checks':CHECKS,
        'inputs':INPUTS,'helper_inputs':k.HASHES,'candidate_sha256':sha(out/'constructed-before-holdout.json'),
        'seconds':time.monotonic()-tick,'formal':'not_run','selectors_changed':False,'historical_outputs_changed':False}
    (out/'matter.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'passed','orbits':31,'color_contributions':35,'checks':len(CHECKS),'seconds':record['seconds']}))

if __name__=='__main__':main()
