"""Rebuild low-loop color sums, scales, counterterms, RG and finite checks.

Input coefficients come from the curated individual integral packets. Frozen
historical totals and the independently translated author finite result are
read only after assembly, as comparisons. No CAS/network fallback is used.
"""
from pathlib import Path
import argparse, collections, hashlib, json, re, subprocess, sys, time
from itertools import permutations

sys.dont_write_bytecode=True
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=Path(__file__).resolve().parents[2]
eps,L,rho,CA,CF,T=s.symbols('eps L rho CA CF T')
M=lambda k,j:s.Symbol(f'M_{k}_{j}')
N=lambda k,j:s.Symbol(f'N_{k}_{j}')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):
    if isinstance(v,dict):v=v.get('sympy',v.get('expression'))
    text=re.sub(r'\b([MN])([0-9])([0-9])\b',r'\1_\2_\3',str(v))
    text=re.sub(r'\b([MN])([0-9])_([0-9])\b',r'\1_\2_\3',text)
    loc={f'zeta_{n}':s.zeta(n) for n in range(2,8)}|{'ell':rho,'tau_D4':s.Symbol('tau_D4')}
    value=s.sympify(text,locals=loc)
    assert not value.atoms(s.Float),text
    return value

def clean(v):
    v=s.sympify(v)
    v=v.subs({M(0,3):-M(0,1)**2,N(0,2):M(0,1)**2/3+4*s.zeta(2)*M(0,1)},simultaneous=True)
    return s.expand(v.subs(M(0,1),-(rho**2+s.pi**2)/2))

def canon(v):
    v=clean(v)
    for a in list(v.free_symbols):
        m=re.fullmatch(r'N_(\d+)_(\d+)',str(a))
        if m:
            k,j=map(int,m.groups());v=v.subs(a,M(k,j+1)-rho*M(k,j))
    return s.expand(v)

def scale(n,hi):
    return {q:s.expand(n**q*((L+rho/2)**q+(L-rho/2)**q)/s.factorial(q)) for q in range(hi+1)}

def scaled(values,n,lo,hi):
    cs=scale(n,hi-min(values))
    return {q:clean(sum(cs[q-p]*v for p,v in values.items() if p<=q)) for q in range(lo,hi+1)}

def series_product(coefficients,weight,lo,hi):
    ws=s.series(weight,eps,0,hi-min(coefficients)+1).removeO().expand()
    return {q:clean(sum(v*ws.coeff(eps,q-p) for p,v in coefficients.items() if p<=q)) for q in range(lo,hi+1)}

def encode(d):return {str(q):str(clean(v)) for q,v in sorted(d.items())}
def dump(path,value):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--regenerate-integrals',action='store_true',help='Run the actual curated mother/series/HLOG engines first, then consume their fresh coefficients')
    ap.add_argument('--reference-dir',type=Path,help='Optional original author ancillary directory; exact hashes checked; no source downloaded or executed')
    args=ap.parse_args();out=args.output_dir.resolve()
    if any(out.is_relative_to((ROOT/name).resolve()) for name in ['data','packages','results','integrals']):
        raise SystemExit('Output must be outside the frozen source/data/result/article trees.')
    if out.exists():raise SystemExit('Use a new output directory; frozen input and prior output are never overwritten.')
    out.mkdir(parents=True);tick=time.monotonic();inputs={};checks=[]
    inputs[Path(__file__).resolve().relative_to(ROOT).as_posix()]=sha(Path(__file__).resolve())
    def read(p):
        p=Path(p);inputs[p.relative_to(ROOT).as_posix() if p.is_relative_to(ROOT) else str(p)]=sha(p)
        return json.loads(p.read_text(encoding='utf-8-sig'))
    def check(name,value,negative=False):
        residual=canon(value)
        assert residual!=0 if negative else residual==0,(name,str(residual)[:2500])
        checks.append({'name':name,'passed':True,'negative_control':negative,'residual':str(residual)})
    onebase=ROOT/'integrals/one-loop';twobase=ROOT/'integrals/two-loop'
    child_runs=[]
    if args.regenerate_integrals:
        for name,base in [('one-loop',onebase),('two-loop',twobase)]:
            dest=out/(name+'-recomputed')
            cmd=[sys.executable,'-X','utf8','-B',str(base/'reproduce.py'),'--output-dir',str(dest)]
            result=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,encoding='utf-8')
            child_runs.append({'family':name,'command':cmd,'exit_code':result.returncode,'stdout':result.stdout,'stderr':result.stderr})
            assert result.returncode==0,child_runs[-1]
        onebase=out/'one-loop-recomputed';twobase=out/'two-loop-recomputed'
    one={}
    for oid in ['2e-01','2e-02']:
        d=read(onebase/oid/'result.json')
        cs={int(q):parse(v) for q,v in d['laurent']['coefficients'].items()}
        assert set(range(-2,4))<=set(cs),('one-loop budget',oid)
        one[oid]=cs
    U1={q:clean(sum(v[q] for v in one.values())) for q in range(-2,4)}
    bare1=scaled(U1,1,-2,3)
    catalog=read(twobase/'catalog.json');two={};ledger=[]
    expected=collections.Counter({'pure':6,'TGV':4,'Pi1':2})
    assert collections.Counter(r['family'] for r in catalog['integrals'])==expected
    for row in catalog['integrals']:
        oid=row['orbit_id'];assert oid not in two and row['assembly_sign']==1
        d=read(twobase/row['result'])
        cs={int(q):parse(v) for q,v in d['coefficients'].items()}
        assert set(range(-4,3))<=set(cs),(oid,'two-loop budget')
        weight=CA*parse(row['color_weights_over_CF']['CA'])+T*parse(row['color_weights_over_CF']['T'])
        two[oid]={'cs':cs,'weight':weight,'family':row['family']}
        ledger.append({'orbit_id':oid,'family':row['family'],'result_path':str((twobase/row['result']).relative_to(ROOT)) if (twobase/row['result']).is_relative_to(ROOT) else str(twobase/row['result']),'canonical_convention':row['canonical_convention'],'assembly_weight_over_CF':str(weight),'assembly_sign_applied':1,'legacy_to_canonical_multiplier_already_absorbed':row['legacy_to_canonical_multiplier'],'orbit_factor_reapplied':False})
    weighted={oid:series_product(d['cs'],d['weight'],-4,2) for oid,d in two.items()}
    U2={q:clean(sum(v[q] for v in weighted.values())) for q in range(-4,3)}
    family_scalar={fam:{q:clean(sum(d['cs'][q] for d in two.values() if d['family']==fam)) for q in range(-4,3)} for fam in expected}
    B1=((10-6*eps)*CA-8*(1-eps)*T)/(3-2*eps)
    grouped={q:clean(CA*(family_scalar['pure'][q]+family_scalar['TGV'][q])+series_product(family_scalar['Pi1'],B1,-4,2)[q]) for q in range(-4,3)}
    for q in U2:check(f'per-orbit versus family color convolution[{q}]',U2[q]-grouped[q])
    bare2=scaled(U2,2,-4,2);b0=(11*CA-4*T)/12
    ct2={q:clean(-b0*bare1.get(q+1,0)) for q in range(-4,3)}
    assert 3 in bare1,'two-loop post epsilon2 requires one-loop epsilon3'
    post2={q:clean(bare2[q]+ct2[q]) for q in range(-4,3)}
    F1,F2=bare1[0],post2[0]
    # A separate ordered-domain/exponent bridge, specific to two loops.
    pair_counts=collections.Counter()
    for p in permutations(range(1,5)):
        pair=tuple(sorted((tuple(sorted(p[:2])),tuple(sorted(p[2:])))))
        pair_counts[pair]+=1
    assert sorted(pair_counts.values())==[8,8,8] and len(pair_counts)==3
    check('Abelian ordered-domain shuffle coefficient for each pairing',s.Rational(8,8)-1)
    check('crossed raw trace minus Abelian trace gives exponent color',CF*(CF-CA/2)-CF**2+CF*CA/2)
    # The author result is an independent post-assembly comparison.
    reference=read(HERE/'reference/fixture.json')
    if args.reference_dir:
        from reference_extract import extract
        fresh=extract(args.reference_dir)
        for key in ['F1_over_CF','F2_over_CF','WR2_over_CF','k2_over_CF','h2_over_CF']:
            check('fresh author-file translation equals derived fixture: '+key,parse(fresh[key])-parse(reference[key]))
    kref,href=parse(reference['k2_over_CF']),parse(reference['h2_over_CF'])
    B,Li2,Li3=s.symbols('B Li2 Li3')
    M11=2*s.zeta(3)+(rho**2+s.pi**2)*B-s.pi**2*rho/2+2*rho*Li2-2*Li3-rho**3/6
    check('full one-loop finite versus author W1',F1-parse(reference['F1_over_CF']))
    check('full two-loop log-W finite versus author W2-W1^2/2',canon(F2).subs(M(1,1),M11)-parse(reference['F2_over_CF']))
    WR2=clean(F2+CF*F1**2/2)
    check('ordinary renormalized W2 versus author ordinary W2',canon(WR2).subs(M(1,1),M11)-parse(reference['WR2_over_CF']))
    expected_poles={-4:s.S.Zero,-3:3*b0/2,-2:b0*L-kref/2,-1:-kref*L-href/2}
    for q,v in expected_poles.items():
        check(f'full post-counterterm epsilon^{q} versus independent RG prediction',post2[q]-v)
        for color in [CA,T]:check(f'{color} pole epsilon^{q}',s.expand(post2[q]-v).coeff(color))
    k2=-s.diff(post2[-1],L);h2=-2*post2[-1].subs(L,0)
    check('graph-extracted cusp versus external reference',k2-kref)
    check('graph-extracted noncusp versus external threshold reference',h2-href)
    check('one-loop cusp pole',bare1[-2]+2)
    check('one-loop simple pole',bare1[-1]+2*L)
    check('one-loop finite RG',s.diff(F1,L)+2*L)
    check('one-loop finite crossing',F1-F1.subs(rho,-rho))
    check('two-loop finite crossing with even M11',F2-F2.subs(rho,-rho))
    check('full two-loop finite RG',s.diff(F2,L)-b0*F1+2*kref*L+href)
    for color in [CA,T]:check(f'{color} finite RG',s.expand(s.diff(F2,L)-b0*F1+2*kref*L+href).coeff(color))
    Z1=sum(bare1[q]*eps**q for q in [-2,-1]);Z2=sum(post2[q]*eps**q for q in range(-4,0))
    dr=s.expand(s.diff(Z2,L)-2*eps*Z2-b0*Z1)
    for q in range(-4,0):check(f'D-dimensional RG negative power {q}',dr.coeff(eps,q))
    check('pole RG finite versus finite-function RG',dr-2*kref*L-href)
    # Existing totals are holdouts, not coefficient construction inputs.
    old1=read(HERE/'sources/benchmarks/renormalization/one-loop/result-r001.json')
    for q in range(-2,3):
        target=parse(old1['laurent_unit_coefficients'][str(q)]).subs({s.Symbol('L_s'):L+rho/2,s.Symbol('L_t'):L-rho/2})
        check(f'fresh one-loop assembly versus accepted total[{q}]',bare1[q]-target)
    old2=read(HERE/'sources/benchmarks/renormalization/two-loop/result-r002.json')
    for q in range(-4,1):
        check(f'fresh two-loop bare versus accepted native total[{q}]',bare2[q]-parse(old2['bare_w2_coefficients_over_CF'][str(q)]))
        check(f'fresh two-loop CT versus accepted total[{q}]',ct2[q]-parse(old2['ordinary_counterterm_over_CF'][str(q)]))
    for key,v in [('F1_over_CF',F1),('F2_over_CF',F2),('f2_over_CF',F2.subs(L,0))]:check('current accepted '+key,v-parse(old2[key]))
    lower=read(HERE/'sources/benchmarks/renormalization/three-loop/lower-order-inputs-r001.json')
    for q in range(-3,2):check(f'computed two-loop bare budget versus three-loop consumed data[{q}]',bare2[q]-parse(lower['hat_w2_over_CF'][str(q)]))
    # Meaningful errors and the finite-shape blind direction are explicit.
    check('negative control duplicated 4e-01 leading pole',two['4e-01']['cs'][-4]*2,True)
    check('negative control missing finite CT',-s.diff(ct2[0],L),True)
    wrong_se=clean(-2*weighted['2e-Pi1-01'][-3])
    check('negative control reapplying legacy Pi1 sign',wrong_se,True)
    check('negative control treating S2[0] as one',-U2[-3],True)
    check('negative control forgetting positive one-loop term in CT',b0*bare1[1],True)
    check('RG is blind to an arbitrary L-independent shape',s.diff(F2+s.Symbol('Phi_of_z'),L)-s.diff(F2,L))
    for rel,digest in inputs.items():
        path=Path(rel) if Path(rel).is_absolute() else ROOT/rel
        assert sha(path)==digest,('Input changed during calculation',rel)
    result={'format':'wilson-low-loop-full-flow/1','quantity':'Connected exponent in lambda=alpha_s/pi, with external CF divided out. T=TF*nf.','domain':'fixed real z>0, rho=log(z), inherited source boundary values','normalization':'All individual I already contain their recorded orbit/crossing and one Euler compensation; restore S1 or S2 once. Native Pi1 already contains the physical sign.','one_loop':{'input_ids':list(one),'U':encode(U1),'bare':encode(bare1),'post_counterterm':encode(bare1),'logZ':encode({q:bare1[q] for q in [-2,-1]}),'F':str(clean(F1)),'k_over_CF':'1','h_over_CF':'0','ordinary_WR1_over_CF':str(clean(F1))},'two_loop':{'input_ids':list(two),'family_scalar_I':{f:encode(v) for f,v in family_scalar.items()},'U':encode(U2),'bare':encode(bare2),'counterterm':encode(ct2),'post_counterterm':encode(post2),'logZ':encode({q:post2[q] for q in range(-4,0)}),'F':str(clean(F2)),'f_at_L0':str(clean(F2.subs(L,0))),'k_over_CF':str(clean(k2)),'h_over_CF':str(clean(h2)),'ordinary_WR2_over_CF':str(clean(WR2))},'three_loop_input_budget':{'one_loop_through':2,'two_loop_through':1,'independently_reassembled_from_curated_integrals':True,'one_loop_extra_order3_for_two_loop_post_order2':True,'post_counterterm_positive_epsilon_terms_are_not_the_physical_epsilon0_finite_term':True},'ordered_shuffle_pairing_counts':{str(k):v for k,v in pair_counts.items()},'all_checks_passed':True,'formal':'not_run','original_integral_evaluation_scope':'See individual packets. This driver assembles and verifies coefficients; optional regenerate-integrals executes their declared mother/HLOG routes, not a new14-integral CAS derivation.'}
    dump(out/'coefficients.json',result)
    dump(out/'assembly-ledger.json',{'one_loop':list(one),'two_loop':ledger,'coefficient_range':{'one':[-2,3],'two':[-4,2]},'source_input_hashes':inputs})
    dump(out/'checks.json',{'status':'passed','checks':len(checks),'negative_controls':sum(r['negative_control'] for r in checks),'checks_detail':checks,'child_runs':child_runs,'seconds':time.monotonic()-tick,'formal':'not_run'})
    dump(out/'three-loop-budget.json',{'status':'reassembled_and_matched','one_loop_bare':encode({q:bare1[q] for q in range(-2,3)}),'two_loop_bare':encode({q:bare2[q] for q in range(-3,2)}),'note':'Actual source-based sums match the accepted three-loop input series. Higher one-loop order3 is additionally supplied for two-loop post epsilon2.'})
    print(json.dumps({'status':'passed','checks':len(checks),'negative_controls':sum(r['negative_control'] for r in checks),'integrals':14,'seconds':time.monotonic()-tick,'three_loop_budget_matched':True,'formal':'not_run'}))

if __name__=='__main__':main()
