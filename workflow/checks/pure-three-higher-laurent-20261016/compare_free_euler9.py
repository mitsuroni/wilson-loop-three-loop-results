"""Full graded comparison with three independent Euler indices and Mellin."""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import json
import math
import compare_weight9 as c
from compare_gauss9 import independent
import hlog9 as h

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('physical');ap.add_argument('m9');ap.add_argument('euler');args=ap.parse_args()
    run,m9,euler=[ROOT/x for x in (args.physical,args.m9,args.euler)]
    assert read(run/'assembly-verification.json')['passed'] and read(m9/'status.json')['passed'] and read(euler/'status.json')['passed']
    deps={Path(__file__).resolve(),HERE/'hlog9.py',HERE/'compare_weight9.py',HERE/'compare_gauss9.py'}
    cache={}
    def source(p):deps.add(p);return read(p)
    def parse(p):
        deps.add(p)
        if p not in cache:cache[p]=h.parse(p)
        return cache[p]
    gens=[{'name':'zeta_'+str(n),'weight':n,'parity':1,'v':h.period(n)} for n in (2,3,5,7,9)]
    gens += [{'name':'Z53','weight':8,'parity':1,'v':h.MZV53},
             {'name':'ell','weight':1,'parity':-1,'v':h.LOG},
             {'name':'log_rho','weight':1,'parity':1,'v':h.add(h.LOG,h.scale({((-1,),h.ZERO):Q(1)},-2))}]
    paths={6:ROOT/source(ROOT/'workflow/checks/single-opposite/current-moments.json')['run'],
           7:ROOT/source(ROOT/'workflow/checks/weight-seven/active-compute_single.json')['run'],
           8:ROOT/source(ROOT/'workflow/checks/4e04-high-order/current-run.json')['run'],9:m9}
    for weight in range(2,10):
        for j in range(1,weight):
            k=weight-j-1
            gens.append({'name':f'M_{k}_{j}','weight':weight,'parity':(-1)**(j+1),'v':parse(paths[max(6,weight)]/f'M_{k}_{j}.txt')})
    low=ROOT/source(ROOT/'workflow/checks/unified-weight-six/current-run.json')['run']
    for g in source(low/'basis-and-conversions.json')['generators']:
        if g['name']=='E' or 'tau' in g['name'].lower():gens.append({'name':g['name'],'weight':g['weight'],'parity':1,'v':h.lift(g['hlog_terms'])})
    old=source(ROOT/'benchmarks/2ne/theory/generation/three-kernel-eps2.json')
    for name in ('W','Lsum'):
        for q in range(3):gens.append({'name':name+str(q),'weight':6+q,'parity':1,'v':h.lift(old['bare'][name][q])})
    family=source(euler/'free-euler-families.json')
    for name,row in family['records'].items():
        plus=parse(euler/(name+'_plus.txt'));minus=parse(euler/(name+'_minus.txt'))
        assert h.add(plus,minus)==h.lift(row['terms']),name
        for parity,value,channel in ((1,plus,'plus'),(-1,minus,'minus')):
            if value:
                assert all(h.weight(key)==row['weight'] for key in value)
                gens.append({'name':name+'_'+channel,'weight':row['weight'],'parity':parity,'v':value})
    # Check the old diagonal source embedding in both crossing parts, not just raw values.
    diagonal=HERE/'gauss-runs/20261016T023349031430Z'
    diagonal_data=source(diagonal/'gauss-families.json');diagonal_checks=0
    for name,row in diagonal_data['records'].items():
        numbers=list(map(int,name.split('_')[1:]));k,j=numbers[:2]
        for channel in ('plus','minus'):
            if name.startswith('Phat'):
                m=numbers[2];names=[(math.comb(k,r),f'Ehat_{r}_{k-r}_{j}_{m}_{channel}.txt') for r in range(1,k+1)]
            else:names=[(math.comb(k,r),f'E_{r}_{k-r}_{j}_{channel}.txt') for r in range(1,k+1)]
            rebuilt=h.add(*(h.scale(parse(euler/file),co) for co,file in names))
            assert rebuilt==parse(diagonal/(name+'_'+channel+'.txt')),name
            diagonal_checks+=1
    declaration={'weight':9,'coefficient_field':'Q','closure':'All homogeneous products, both crossing parities retained at intermediate weights.',
                 'source':'E(a,b,c;z) with independent a,b,c, plus one independent Mellin exponent; includes prior P family by 240 exact crossing regressions.',
                 'generators':[{k:v for k,v in g.items() if k!='v'} for g in gens],
                 'limits':'Does not enumerate arbitrary decorated two-kernel sources or rational-z coefficient modules.'}
    dest=run/'free-euler-space-declaration.json';assert not dest.exists()
    dest.write_text(json.dumps(declaration,indent=2)+'\n',encoding='utf-8')
    bases={};ranks=[]
    for weight in range(1,10):
        for parity in (1,-1):
            opts=[{'factors':[g['name']],'v':g['v']} for g in gens if g['weight']==weight and g['parity']==parity]
            for i in range(1,weight//2+1):
                for p in (1,-1):
                    for a in bases[(i,p)]:
                        for b in bases[(weight-i,parity*p)]:
                            opts.append({'factors':sorted(a['factors']+b['factors']),'v':h.mul(a['v'],b['v'])})
            bases[(weight,parity)]=independent(opts)
        row={'weight':weight,'even':len(bases[(weight,1)]),'odd':len(bases[(weight,-1)])}
        ranks.append(row);print('free Euler closure',row,flush=True)
    opts=bases[(9,1)]
    columns=[{'factors':o['factors'],'terms':h.serialize(o['v'])} for o in opts]
    (run/'free-euler-weight9-columns.json').write_text(json.dumps(columns,indent=2)+'\n',encoding='utf-8')
    data=source(run/'coefficients.json')
    targets={name:h.lift(data['physical'][oid][3]) for name,oid in (('I20_3','6e-20'),('I21_3','6e-21'))}
    targets.update({name:h.lift(data['bare'][key][3]) for name,key in (('W3','W'),('Lambda3','Lsum'))})
    reports={name:c.span.fit(value,opts,h.serialize) for name,value in targets.items()}
    cols=[o['v'] for o in opts];rank=len(cols)
    augmented=len(c.span.build(cols+[targets['I20_3'],targets['I21_3']]))
    source_augmented=len(c.span.build(cols+[targets['W3'],targets['Lambda3']]))
    assert augmented==source_augmented and rank>=373
    out={'passed':True,'formal':'not_run','rank':rank,'augmented_rank':augmented,'pair_quotient_dimension':augmented-rank,
         'graded_ranks':ranks,'targets':reports,'diagonal_crossing_regressions':diagonal_checks,
         'declaration_sha256':sha(dest),'columns_sha256':sha(run/'free-euler-weight9-columns.json'),
         'inputs_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in sorted(deps)},'limits':declaration['limits']}
    (run/'comparison-free-euler-weight9.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'passed':True,'rank':rank,'augmented_rank':augmented,'pair_quotient_dimension':augmented-rank}))

if __name__=='__main__':main()
