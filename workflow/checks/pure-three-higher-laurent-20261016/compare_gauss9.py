"""Exact graded algebra closure after adding the existing Gauss/Mellin source."""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import json
import compare_weight9 as c
import hlog9 as h

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def independent(opts):
    opts=c.unique(opts)
    basis=c.span.build([o['v'] for o in opts])
    chosen=sorted({max(rep) for vec,rep in basis.values()})
    assert len(chosen)==len(basis)
    return [opts[i] for i in chosen]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('physical');ap.add_argument('m9');ap.add_argument('gauss');args=ap.parse_args()
    run,m9,gauss=[ROOT/x for x in (args.physical,args.m9,args.gauss)]
    assert read(run/'assembly-verification.json')['passed'] and read(m9/'status.json')['passed'] and read(gauss/'status.json')['passed']
    deps={Path(__file__).resolve(),HERE/'hlog9.py',HERE/'compare_weight9.py',gauss/'gauss-families.json'}
    def source(p):deps.add(p);return read(p)
    def parse(p):deps.add(p);return h.parse(p)
    gens=[{'name':'zeta_'+str(n),'weight':n,'parity':1,'v':h.period(n)} for n in (2,3,5,7,9)]
    gens += [{'name':'Z53','weight':8,'parity':1,'v':h.MZV53},
             {'name':'ell','weight':1,'parity':-1,'v':h.LOG},
             {'name':'log_rho','weight':1,'parity':1,'v':h.add(h.LOG,h.scale({((-1,),h.ZERO):Q(1)},-2))}]
    paths={6:ROOT/source(ROOT/'workflow/checks/single-opposite/current-moments.json')['run'],
           7:ROOT/source(ROOT/'workflow/checks/weight-seven/active-compute_single.json')['run'],
           8:ROOT/source(ROOT/'workflow/checks/4e04-high-order/current-run.json')['run'],9:m9}
    for weight in range(2,10):
        for j in range(1,weight):
            k=weight-j-1;v=parse(paths[max(6,weight)]/f'M_{k}_{j}.txt')
            gens.append({'name':f'M_{k}_{j}','weight':weight,'parity':(-1)**(j+1),'v':v})
    low=ROOT/source(ROOT/'workflow/checks/unified-weight-six/current-run.json')['run']
    for g in source(low/'basis-and-conversions.json')['generators']:
        if g['name']=='E' or 'tau' in g['name'].lower():gens.append({'name':g['name'],'weight':g['weight'],'parity':1,'v':h.lift(g['hlog_terms'])})
    baseline=source(ROOT/'benchmarks/2ne/theory/generation/three-kernel-eps2.json')
    for name in ('W','Lsum'):
        for q in range(3):gens.append({'name':name+str(q),'weight':6+q,'parity':1,'v':h.lift(baseline['bare'][name][q])})
    family=source(gauss/'gauss-families.json')
    for name,row in family['records'].items():
        plus=parse(gauss/(name+'_plus.txt'));minus=parse(gauss/(name+'_minus.txt'))
        assert h.add(plus,minus)==h.lift(row['terms']),name
        for parity,value,channel in ((1,plus,'plus'),(-1,minus,'minus')):
            if value:
                assert all(h.weight(key)==row['weight'] for key in value)
                gens.append({'name':name+'_'+channel,'weight':row['weight'],'parity':parity,'v':value})
    declaration={'weight':9,'coefficient_field':'Q','closure':'All homogeneous products through weight nine, both crossing parities retained at intermediate weights.',
                 'sources':'M/N plus periods, known normalization logs, old weight-six extras, old W/Lambda through order two, existing Gauss P and Mellin-P.',
                 'generators':[{k:v for k,v in g.items() if k!='v'} for g in gens],
                 'limits':'Does not include every independently weighted two-kernel source or arbitrary rational-z prefactors.'}
    dest=run/'gauss-space-declaration.json';assert not dest.exists()
    dest.write_text(json.dumps(declaration,indent=2)+'\n',encoding='utf-8')
    bases={};ranks=[]
    for weight in range(1,10):
        for parity in (1,-1):
            opts=[{'factors':[g['name']],'v':g['v']} for g in gens if g['weight']==weight and g['parity']==parity]
            for i in range(1,weight//2+1):
                j=weight-i
                for p in (1,-1):
                    for a in bases[(i,p)]:
                        for b in bases[(j,parity*p)]:
                            opts.append({'factors':sorted(a['factors']+b['factors']),'v':h.mul(a['v'],b['v'])})
            bases[(weight,parity)]=independent(opts)
        ranks.append({'weight':weight,'even':len(bases[(weight,1)]),'odd':len(bases[(weight,-1)])})
        print('graded closure',ranks[-1],flush=True)
    opts=bases[(9,1)]
    columns=[{'factors':o['factors'],'terms':h.serialize(o['v'])} for o in opts]
    (run/'gauss-weight9-columns.json').write_text(json.dumps(columns,indent=2)+'\n',encoding='utf-8')
    # New targets are first used after the declared graded closure is complete.
    target=source(run/'coefficients.json')
    targets={name:h.lift(target['physical'][oid][3]) for name,oid in (('I20_3','6e-20'),('I21_3','6e-21'))}
    targets.update({name:h.lift(target['bare'][key][3]) for name,key in (('W3','W'),('Lambda3','Lsum'))})
    reports={name:c.span.fit(value,opts,h.serialize) for name,value in targets.items()}
    cols=[o['v'] for o in opts];rank=len(opts)
    augmented=len(c.span.build(cols+[targets['I20_3'],targets['I21_3']]))
    srcaug=len(c.span.build(cols+[targets['W3'],targets['Lambda3']]))
    assert augmented==srcaug
    out={'passed':True,'formal':'not_run','rank':rank,'augmented_rank':augmented,
         'pair_quotient_dimension':augmented-rank,'graded_ranks':ranks,'targets':reports,
         'declaration_sha256':sha(dest),'columns_sha256':sha(run/'gauss-weight9-columns.json'),
         'inputs_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in sorted(deps)},'limits':declaration['limits']}
    (run/'comparison-gauss-weight9.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'passed':True,'rank':rank,'augmented_rank':augmented,'pair_quotient_dimension':augmented-rank}))

if __name__=='__main__':main()
