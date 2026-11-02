"""Fresh AA sum from current tree/insertion sources and complete joint objects.

Historical pure parser/HLOG functions are reused; no historical writer main is
called. No sign conversion is applied. Accepted family/AA totals are holdouts.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
import argparse,hashlib,importlib.util,json,re,sys,time
sys.dont_write_bytecode=True
import sympy as s

CODE=Path(__file__).resolve().parent;ROOT=CODE.parents[3];H=CODE
P=tuple(range(-6,1));INPUTS={};CHECKS=[];FIELDS=[]
def raw(rel,expected=None):
    p=ROOT/rel;b=p.read_bytes();v=hashlib.sha256(b).hexdigest()
    assert rel not in INPUTS or INPUTS[rel]==v,rel
    assert expected is None or v==expected,(rel,v,expected)
    INPUTS[rel]=v;return b
def read(rel,expected=None):return json.loads(raw(rel,expected).decode('utf-8-sig'))
def load(name,rel):
    raw(rel);spec=importlib.util.spec_from_file_location(name,ROOT/rel);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
def dump(name,v):
    p=H/name;assert not p.exists(),p;p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def field(path,pointer,value):
    FIELDS.append({'path':path,'sha256':INPUTS[path],'json_pointer':pointer,'field_sha256':hashlib.sha256(json.dumps(value,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()).hexdigest()})
    return value

h=load('canonical_aa_helpers','workflow/checks/four-leg-bootstrap-20261101/known-total/assemble_known.py')
INPUTS.update(h.HASHES)
h.HASHES=INPUTS
t=load('canonical_tree_helpers','workflow/checks/aa-joint-rg-20261024/tree-sums/build.py')
math=h.math;eps,rho,z=h.eps,h.rho,h.z
dictionary=read('workflow/checks/four-leg-bootstrap-20261101/basis/dictionary.json')
h.MOMENTS={s.Symbol(x['name']):h.vec_expr(x['hlog_terms']) for x in dictionary['atoms']};h.MOMENTS[s.Symbol('ell')]=rho
def normal(x):return h.normalize(x)
def same(name,a,b=0,negative=False):
    v=s.cancel(normal(a)-normal(b));ok=v!=0 if negative else v==0
    assert ok,(name,str(v)[:1500]);CHECKS.append({'name':name,'passed':True,'negative_control':negative,'residual':str(v)})
def vec(rows,labels=None):
    out=defaultdict(Q)
    for x in rows:
        ps=x.get('period_powers',x.get('zeta_powers'))
        factor=Q(1)
        if ps is None and 'pi_power' in x:
            assert x['pi_power']%2==0
            ps=(x['pi_power']//2,x['zeta3_power'],x['zeta5_power']);factor=Q(6)**(x['pi_power']//2)
        assert ps is not None,x
        ps=tuple(ps)+(0,)*(6-len(ps));assert len(ps)==6
        w=tuple(Q(a) for a in x.get('word',x.get('Hlog_word')));out[(w,ps)]+=factor*Q(x['coefficient'])
    return {k:v for k,v in out.items() if v}
def standard(v):
    assert all(all(a in (0,-1) for a in w) and not any(ps[3:]) for w,ps in v), 'Extra letters/periods survive the complete tree sum'
    return {(tuple(map(int,w)),ps[:3]):c for (w,ps),c in v.items()}
def expr(v):return sum(s.Rational(c.numerator,c.denominator)*math.hl('z',w)*s.prod(s.zeta(k)**n for k,n in zip((2,3,5),ps)) for (w,ps),c in v.items())
def coeff(data,q,path):
    a=data['laurent']['coefficients'][str(q)];ptr=f'/laurent/coefficients/{q}'
    if isinstance(a,dict) and 'hlog_terms' in a:rows=a['hlog_terms'];field(path,ptr+'/hlog_terms',rows);return vec(rows)
    if str(q) in data.get('hlog_coefficients',{}):rows=data['hlog_coefficients'][str(q)];field(path,f'/hlog_coefficients/{q}',rows);return vec(rows)
    field(path,ptr,a);vv,pointer=t.parse_coeff(data,q)
    return {(w,ps+(0,0,0)):c for (w,ps),c in vv.items()}
def trees():
    rows=[]
    for p in sorted((ROOT/'benchmarks/2ne/families').glob('*/catalog.json')):
        rel=p.relative_to(ROOT).as_posix();d=read(rel)
        for x in d['items']:
            if x['orbit_id'].startswith('6e-'):rows.append({'id':x['orbit_id'],'path':x.get('latest_result',x.get('result')),'expected':x.get('result_sha256'),'family':'6e','selector':rel})
    rel='benchmarks/tgv/5e-3g/current.json';c=read(rel)
    rows.extend({'id':k,'path':'benchmarks/tgv/5e-3g/'+v['result'],'expected':v['sha256'],'family':'5e','selector':rel} for k,v in sorted(c['orbits'].items()))
    assert len({x['id'] for x in rows})==len(rows)==83
    groups={f:{p:{} for p in P} for f in ['6e','5e']};ledger=[];special={}
    for x in rows:
        d=read(x['path'],x['expected']);start=len(FIELDS)
        for p in P:
            v=coeff(d,p,x['path']);groups[x['family']][p]=h.ha.add(groups[x['family']][p],v)
            if x['id'] in ['5e-3g-27','5e-3g-28','5e-3g-29']:special.setdefault(p,[]).append(v)
        ledger.append(x|{'sha256':INPUTS[x['path']],'assembly_weight':'1','fields':FIELDS[start:]})
    for p,vv in special.items():standard(h.ha.add(*vv))
    for family in groups:
        for p in P:groups[family][p]=expr(standard(groups[family][p]))
    return groups,ledger

def insertions():
    a=(10-6*eps)/(3-2*eps);out=[]
    p='benchmarks/self-energy/current.json';outer=read(p);manifest=read('benchmarks/self-energy/'+outer['manifest'],outer['sha256'])
    for sel in manifest['entries']:
        family={'4e-Pi1':'Pi1','2e-Pi1+Pi1':'chain','2e-Pi2':'Pi2'}.get(sel['family'])
        if family is None:continue
        path=sel['result']['path'];d=read(path,sel['result']['sha256']);assert sel['assembly_sign']==d['normalization']['assembly_sign']==1
        assert d['normalization']['Euler_count']==1
        if family=='Pi1':fields=d['laurent']['coefficients'];pointer='/laurent/coefficients';weight=a
        else:
            i=next(i for i,v in enumerate(d['color_channels']) if v.get('monomial',v.get('color_monomial'))=='CA^2*CF');c=d['color_channels'][i]
            fields=c['coefficients'] if family=='chain' else c['laurent']['coefficients'];pointer=f'/color_channels/{i}/'+('coefficients' if family=='chain' else 'laurent/coefficients');weight=s.S.One
        vals,orig=h.series_values(d,fields,pointer,path)
        vals=h.convolve(vals,weight);out.append({'id':sel['orbit_id'],'family':family,'path':path,'weight':str(weight),'coefficients':{p:normal(v) for p,v in vals.items()},'fields':orig})
    for family,folder in [('legs','tgv-pi1'),('Gamma','gamma1')]:
        c=read(f'benchmarks/self-energy/{folder}/current.json')
        for oid,sel in c['orbits'].items():
            chosen=sel.get('preferred_representation',sel);path=f'benchmarks/self-energy/{folder}/'+chosen['result'];d=read(path,chosen['sha256']);assert d['normalization']['assembly_sign']==1 and d['normalization']['euler_compensation_count']==1
            if family=='Gamma':vals,orig=h.series_values(d,d['channels']['g']['I_coefficients'],'/channels/g/I_coefficients',path,d['channels']['g']['highest_pole_order']);weight=s.S.One
            else:
                vals,orig=h.series_values(d,d['laurent_coefficients'],'/laurent_coefficients',path,d['highest_pole_order']);weight=a
                if oid.endswith(('07','08')):
                    vals[-2]=h.parse(field(path,'/second_pole_in_existing_M_basis',d['second_pole_in_existing_M_basis']))
                    if oid.endswith('07'):
                        pp='benchmarks/self-energy/tgv-pi1/results/3e-3g+Pi1-07/mn-r002.json';mn=read(pp);assert mn['target_result_sha256']==INPUTS[path]
                        for p in [-1,0]:vals[p]=h.parse(field(pp,f'/coefficients/{p}',mn['coefficients'][str(p)]))
                    else:
                        for p in [-1,0]:vals[p]=h.parse(field(path,f'/M_representation/coefficients/{p}',d['M_representation']['coefficients'][str(p)]))
            vals=h.convolve(vals,weight);out.append({'id':oid,'family':family,'path':path,'weight':str(weight),'coefficients':{p:normal(v) for p,v in vals.items()},'fields':orig})
    assert len(out)==31
    return {p:normal(sum(x['coefficients'][p] for x in out)) for p in P},out

def selected(oid,fg,dt):
    base='benchmarks/fgv/4e-4g/' if oid.startswith('4e-4g-') else 'benchmarks/double-tgv/4e-3g+3g/'
    x=(fg if oid.startswith('4e-4g-') else dt)['orbits'][oid];return base+x['result'],x['sha256']
def fourleg():
    k4=read('workflow/checks/four-leg-bootstrap-20261101/basis/accepted-k4.json')
    fg=read('benchmarks/fgv/4e-4g/current.json');dt=read('benchmarks/double-tgv/4e-3g+3g/current.json')
    sums={p:0 for p in P};ledger=[]
    for x in k4['sources']:
        path,digest=selected(x['orbit_id'],fg,dt);read(path,digest);factor=2 if x['group']=='I3444' else 1
        for p in P:
            value,ptr=h.extract_coefficient(path,p,digest);sums[p]+=factor*normal(value)
        ledger.append({'orbit_id':x['orbit_id'],'joint_group':x['group'],'path':path,'sha256':digest,'reflection_weight':factor})
    A={p:0 for p in P}
    for oid in ['4e-4g-05','4e-3g+3g-13','4e-3g+3g-15']:
        path,digest=selected(oid,fg,dt);raw(path,digest)
        for p in range(-6,0):v,ptr=h.extract_coefficient(path,p,digest,unknown=True);A[p]+=normal(v)
        ledger.append({'orbit_id':oid,'joint_group':'I2344','path':path,'sha256':digest,'reflection_weight':2,'finite_consumed_via':'complete I2344 HLOG below'})
    pp='workflow/checks/four-leg-2344-multiroute-20261024/weight-physical-hlog.json';d=read(pp);field(pp,'/coordinates',d['coordinates'])
    for x in d['coordinates']:
        v=x['coefficient'].replace('^','**').replace('z1','z');v=re.sub(r'zeta\[(\d+)\]',r'zeta(\1)',v)
        A[0]+=s.sympify(v,locals={'z':z,'Pi':s.pi,'zeta':s.zeta})*math.hl('z',x['word'])
    pp='workflow/checks/four-leg-1234-b-h1-20261101/physical-c0-candidate/I1234-Laurent-candidate.json';d=read(pp);field(pp,'/coefficients',d['coefficients']);B={p:s.S.Zero for p in P}
    for p,v in d['coefficients'].items():
        for x in v['terms']:
            B[int(p)]+=s.sympify(x['coefficient'],locals={'Z2':s.pi**2/6,'Z3':s.zeta(3),'Z5':s.zeta(5)})*math.hl('z',x['word'])
    ledger.append({'joint_group':'I1234','source_orbits':['4e-4g-09','4e-3g+3g-25','4e-3g+3g-26'],'path':pp,'sha256':INPUTS[pp],'weight':1,'DT25_rotation_counted_again':False})
    # This old file supplies only the attachment/source-ID grouping. No stored
    # coefficient is used. Current selectors bind all35 source records below.
    groups=read('workflow/checks/aa-joint-rg-20261024/four-leg/four-leg.json')['groups']
    source_ids=[];joint_map=[]
    for group,row in groups.items():
        owner={'3334':'3444','2234':'2344'}.get(group,group);members=[]
        for oid in row['source_orbits']:
            source_ids.append(oid);sel=(fg if oid.startswith('4e-4g-') else dt)['orbits'][oid]
            member={'orbit_id':oid,'selector_row':sel,'added_separately_after_joint':False}
            if sel.get('result'):
                path,digest=selected(oid,fg,dt);raw(path,digest);member.update(result=path,sha256=digest)
            members.append(member)
        joint_map.append({'joint_object':group,'represented_by':owner,'joint_consumption_count':1,'members':members,'reflection_input':'Established common-normalization reflection identity reused' if owner!=group else None})
    assert len(source_ids)==len(set(source_ids))==35 and len(joint_map)==9
    raw('benchmarks/renormalization/three-loop/four-leg-grouping.md')
    return {p:normal(sums[p]+2*A[p]+B[p]) for p in P},ledger,joint_map

def constant_mn(values):
    """Exact HLOG-coordinate reduction in the pre-existing raw moment family.

    Enumerate all MN/zeta monomials of weight<=6 and even inversion parity.
    No accepted AA coefficient or its support is read to select this space.
    """
    atoms=[x for x in dictionary['atoms'] if x['name'] not in ['ell','tauDq4','E']]
    # Prefer smaller indices/raw low-weight moments in underdetermined spaces.
    atoms.sort(key=lambda x:(x['weight'],x['name']))
    monomials=[]
    def rec(index,remaining,parity,value):
        if index==len(atoms):
            if parity==1:monomials.append(s.expand(value))
            return
        a=atoms[index]
        for n in range(remaining//a['weight']+1):rec(index+1,remaining-n*a['weight'],parity*a['inversion_parity']**n,value*s.Symbol(a['name'])**n)
    rec(0,6,1,s.S.One)
    monomials.sort(key=lambda x:(sum(int(v)*next(a['weight'] for a in atoms if str(k)==a['name']) for k,v in x.as_powers_dict().items() if k!=1),str(x)))
    vectors=[h.expression_vector(normal(v)) for v in monomials]
    targets=[h.expression_vector(values[p]) for p in P];keys=sorted(set().union(*(set(v) for v in vectors+targets)))
    matrix=s.Matrix([[v.get(k,0) for v in vectors] for k in keys]);rhs=s.Matrix([[v.get(k,0) for v in targets] for k in keys])
    solutions,free=matrix.gauss_jordan_solve(rhs)
    solutions=solutions.subs(dict.fromkeys(free.free_symbols,0))
    result={p:s.expand(sum(v*c for v,c in zip(monomials,solutions[:,j]))) for j,p in enumerate(P)}
    for p in P:same(f'fresh raw-MN reduction exact back substitution[{p}]',values[p],result[p])
    return result,{'candidate_monomials':len(monomials),'HLOG_coordinates':len(keys),'free_coefficients_fixed_to_zero':len(free.free_symbols),'selection':'All even MN/zeta monomials of weight<=6 before reading accepted AA target; exact rational coordinate solve, not numerical fitting'}

def main():
    global H
    ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args();H=args.output_dir.resolve()
    assert not H.exists(),'Use a new output directory'
    assert not H.is_relative_to(ROOT) or H.is_relative_to(CODE) or any(H.is_relative_to(ROOT/x) for x in ['.wilson','replay-output']),'Output would enter a source tree'
    H.mkdir(parents=True)
    started=time.monotonic();raw(Path(__file__).relative_to(ROOT).as_posix());treesums,treeledger=trees();print('trees built',flush=True)
    ins,insledger=insertions();print('insertions built',flush=True)
    four,fourledger,jointmap=fourleg();print('fourlegs built',flush=True)
    U={p:normal(treesums['6e'][p]+treesums['5e'][p]+ins[p]+four[p]) for p in P}
    mn,reduction=constant_mn(U);print('MN reduction built',flush=True)
    dump('predicted-U.json',{'quantity':'U3 AA; CF CA^2 removed','coefficients':{str(p):h.serialize_vec(h.expression_vector(v)) for p,v in U.items()},'U_constant_MN':{str(p):str(v) for p,v in mn.items()},'MN_reduction':reduction,'construction':'83 individually selected tree coefficients +31 native insertion coefficients + complete fourleg objects; no total/aggregate used as construction input','formal':'not_run'})
    reference=read('benchmarks/renormalization/three-loop/first-complete-20261102/coefficients.json')
    compact=reference['channels']['AA']['U']
    for p in P:same(f'fresh AA epsilon{p} versus accepted full AA',U[p],h.parse(compact[str(p)]))
    for family in treesums:
        control=read(f'benchmarks/recursion-audit/results/{family}-effective-sum/r001.json')
        for p in P:same(f'fresh {family}[{p}] versus accepted aggregate HOLDOUT',treesums[family][p],h.parse(control['hlog_coefficients'][str(p)]))
    control=next(x for x in treeledger if x['id']=='6e-01');d=read(control['path'],control['sha256'])
    omitted=expr(standard(coeff(d,0,control['path'])))
    same('negative: omitting actual6e01 finite contribution is detected',U[0],U[0]-omitted,True)
    control=next(x for x in insledger if x['id']=='3e-Gamma1-01')
    same('negative: reapplying old Gamma family minus sign is detected',U[0],U[0]-2*control['coefficients'][0],True)
    for path,digest in t.SOURCES.items():raw(path,digest)
    for path,digest in h.HASHES.items():raw(path,digest)
    for path,digest in INPUTS.items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    dump('assembly.json',{'status':'passed_fresh_current_AA_reassembly','quantity':'scale-stripped connected U3 AA; common CF CA^2 removed','normalization':'Current physical I, assembly_sign+1; genuine a(eps) on scalar Pi1/legs, stored colored channels once; fourleg joint/reflection weights only','powers':list(P),'U':{str(p):str(v) for p,v in U.items()},'U_HLOG':{str(p):h.serialize_vec(h.expression_vector(v)) for p,v in U.items()},'U_constant_MN':{str(p):str(v) for p,v in mn.items()},'MN_reduction':reduction,'family_HLOG':{f:{str(p):h.serialize_vec(h.expression_vector(v)) for p,v in d.items()} for f,d in (treesums|{'insertions':ins,'four-leg':four}).items()},'tree_sources':treeledger,'insertion_sources':[{k:({str(p):str(v) for p,v in v.items()} if k=='coefficients' else v) for k,v in x.items()} for x in insledger],'fourleg_sources':fourledger,'fourleg_nine_joint_ownership':jointmap,'source_orbit_coverage':{'tree':83,'insertions':31,'fourleg':35,'total':149},'checks':CHECKS,'absolute_epsilon_convolution_checks':h.CHECKS,'tree_nonstandard_letter_period_cancellation':'5e27+28+29 have zero nonstandard letter/period coordinates at every power, retaining all rational letters and six period powers before summation; ordinary0/-1 coordinates may remain','construction_did_not_consume_old_aggregate_or_AA_total':True,'original_integrations_reused':True,'I2344_source':'accepted existing full research solution; no new independent general-z integration','I1234_source':'new admitted original-kernel solution reused','input_sha256':INPUTS,'field_bindings':FIELDS+h.FIELDS,'seconds':time.monotonic()-started,'formal':'not_run'})
    print(json.dumps({'status':'passed','checks':len(CHECKS),'source_inputs':len(INPUTS),'seconds':time.monotonic()-started}))
if __name__=='__main__':main()
