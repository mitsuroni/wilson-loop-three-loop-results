"""Declared rational weight-nine spaces, exact fits and dual obstructions."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
import argparse
import hashlib
import importlib.util
import json
import math
import sys
import hlog9 as h

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
spec=importlib.util.spec_from_file_location('_pure3_exact_span',ROOT/'workflow/checks/three-kernel-eps2/exact_span.py')
span=importlib.util.module_from_spec(spec);spec.loader.exec_module(span)
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def signature(v):return tuple(sorted(v.items()))

def unique(options):
    seen=set();out=[]
    for option in options:
        sig=signature(option['v'])
        if sig and sig not in seen:seen.add(sig);out.append(option)
    return out

def options(gens,total=9):
    result=[]
    def visit(start,left,parity,chosen):
        if left==0:
            if parity==1:
                value=h.ONE
                for i in chosen:value=h.mul(value,gens[i]['v'])
                result.append({'factors':[gens[i]['name'] for i in chosen],'v':value})
            return
        for i in range(start,len(gens)):
            g=gens[i]
            if g['weight']>left or (g['name']=='ell' and i in chosen):continue
            visit(i,left-g['weight'],parity*g['parity'],chosen+[i])
    visit(0,total,1,[])
    return unique(sorted(result,key=lambda r:(len(r['factors']),r['factors'])))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('physical_run');ap.add_argument('m9_run')
    ap.add_argument('--revision',default='r002');args=ap.parse_args()
    assert args.revision.startswith('r') and args.revision[1:].isdigit()
    prefix='' if args.revision=='r001' else args.revision+'-'
    run=ROOT/args.physical_run;m9=ROOT/args.m9_run
    assert read(run/'assembly-verification.json')['passed'] and read(m9/'status.json')['passed']
    assert span.check()
    dependencies={Path(__file__).resolve(),HERE/'hlog9.py',run/'coefficients.json',run/'assembly-verification.json',m9/'status.json'}
    def source(path):dependencies.add(path);return read(path)
    def parse(path):dependencies.add(path);return h.parse(path)
    pointers={
        6:ROOT/'workflow/checks/single-opposite/current-moments.json',
        7:ROOT/'workflow/checks/weight-seven/active-compute_single.json',
        8:ROOT/'workflow/checks/4e04-high-order/current-run.json'}
    locations={w:ROOT/source(p)['run'] for w,p in pointers.items()}
    gens=[{'name':'zeta_'+str(n),'weight':n,'parity':1,'v':h.period(n)} for n in (2,3,5,7,9)]
    gens += [{'name':'Z53','weight':8,'parity':1,'v':h.MZV53}]
    for weight in range(2,10):
        folder=locations[6] if weight<=6 else locations[weight] if weight<=8 else m9
        for j in range(1,weight):
            k=weight-j-1;vec=parse(folder/f'M_{k}_{j}.txt')
            assert all(h.weight(key)==weight for key in vec)
            gens.append({'name':f'M_{k}_{j}','weight':weight,'parity':(-1)**(j+1),'v':vec})
    gens.append({'name':'ell','weight':1,'parity':-1,'v':h.LOG})
    lowptr=ROOT/'workflow/checks/unified-weight-six/current-run.json'
    low=ROOT/source(lowptr)['run']
    oldgens=source(low/'basis-and-conversions.json')['generators']
    lowextra=[{'name':g['name'],'weight':g['weight'],'parity':1,'v':h.lift(g['hlog_terms'])}
              for g in oldgens if g['name']=='E' or 'tau' in g['name'].lower()]
    assert len(lowextra)==2,[g['name'] for g in lowextra]
    rho_log=h.add(h.LOG,h.scale({((-1,),h.ZERO):Q(1)},-2))
    normal={'name':'log_rho','weight':1,'parity':1,'v':rho_log}
    data=source(run/'coefficients.json')
    previous=[]
    for name in ('W','Lsum'):
        for q in range(3):previous.append({'name':name+str(q),'weight':6+q,'parity':1,'v':h.lift(data['bare'][name][q])})
    tiers={
        'MN_periods':gens,
        'MN_old_weight6_products':gens+lowextra,
        'known_normalization_closure':gens+lowextra+[normal],
        'previous_three_kernel_coefficients':gens+lowextra+[normal]+previous}
    # Save the finite candidate rules and complete columns before any fitting.
    declaration={'coefficient_field':'Q','weight':9,'crossing_parity':'even',
                 'period_weights':list(h.ZW),'M_range':'all k>=0,j>=1 with k+j+1<=9',
                 'N_coverage':'N_k,j=M_k,j+1-ell*M_k,j; ell exponent at most one via ell^2=-2M01-6*zeta2.',
                 'tiers':{name:[{k:v for k,v in g.items() if k!='v'} for g in gs] for name,gs in tiers.items()},
                 'limits':'Finite registered source space. No arbitrary rational-z coefficients, no complete enumeration of all decorated lower-kernel sources.'}
    dest=run/(prefix+'space-declaration.json');assert not dest.exists()
    dest.write_text(json.dumps(declaration,indent=2)+'\n',encoding='utf-8')
    option_sets={name:options(gs) for name,gs in tiers.items()}
    cols={name:[{'factors':o['factors'],'terms':h.serialize(o['v'])} for o in opts] for name,opts in option_sets.items()}
    (run/(prefix+'candidate-columns.json')).write_text(json.dumps(cols,indent=2)+'\n',encoding='utf-8')
    targets={'I20_3':h.lift(data['physical']['6e-20'][3]),'I21_3':h.lift(data['physical']['6e-21'][3]),
             'W3':h.lift(data['bare']['W'][3]),'Lambda3':h.lift(data['bare']['Lsum'][3])}
    for name in ('A','B'):
        seq=[h.lift(v) for v in data['bare'][name]]
        stripped=h.add(*(h.scale(h.mul(h.power(rho_log,r),seq[3-r]),Q((-3)**r,2**r*math.factorial(r))) for r in range(4)))
        targets[name+'3_without_rho']=stripped
    reports={}
    for tier,opts in option_sets.items():
        columns=[o['v'] for o in opts]
        rank=len(span.build(columns))
        pairs={label:len(span.build(columns+[targets[a],targets[b]]))-rank
               for label,a,b in (('physical','I20_3','I21_3'),('without_rho','A3_without_rho','B3_without_rho'),
                                 ('natural_sources','W3','Lambda3'))}
        fitted={name:span.fit(v,opts,h.serialize) for name,v in targets.items()}
        reports[tier]={'candidate_count':len(opts),'rank':rank,'pair_quotient_dimensions':pairs,'targets':fitted}
        print(tier,len(opts),rank,pairs,flush=True)
    last=reports['previous_three_kernel_coefficients']['pair_quotient_dimensions']
    assert last['physical']==last['without_rho']==last['natural_sources']
    # This equivalence is required only in the final tier, which explicitly
    # includes all lower source coefficients and known normalization factors.
    snapshot=run/('compare-weight9-'+args.revision+'-source.py')
    assert not snapshot.exists();snapshot.write_bytes(Path(__file__).read_bytes())
    report={'passed':True,'formal':'not_run','weight':9,'declaration_sha256':sha(dest),
            'candidate_columns_sha256':sha(run/(prefix+'candidate-columns.json')),'spaces':reports,
            'normalization_equivalence_checked':True,'executed_source_snapshot':snapshot.name,
            'inputs_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in sorted(dependencies)},
            'limits':declaration['limits']}
    out=run/(prefix+'comparison-weight9.json');assert not out.exists()
    out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Exact declared-space comparisons complete.',flush=True)

if __name__=='__main__':main()
