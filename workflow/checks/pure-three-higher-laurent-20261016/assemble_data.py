"""Parse complete fresh physical outputs, regress old orders, export new data."""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import json
import sys
import hlog9 as h

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('run');args=ap.parse_args()
    run=ROOT/args.run
    status=read(run/'status.json');assert status['passed'] and status['order']==3
    for name,expected in status['outputs_sha256'].items():assert sha(run/name)==expected
    deps={run/'status.json',run/'provenance.json',Path(__file__).resolve(),HERE/'hlog9.py'}
    vectors={}
    for name in ('W','Lsum','A','B','H'):
        vectors[name]=[]
        for q in range(4):
            p=run/f'{name}{q}.txt';deps.add(p);v=h.parse(p)
            assert all(h.weight(k)==(2 if name=='H' else 6)+q for k in v)
            vectors[name].append(v)
    for name in ('I20','I21'):
        vectors[name]=[]
        for q in range(4):
            p=run/f'{name}_{q}.txt';deps.add(p);v=h.parse(p)
            assert all(h.weight(k)==6+q for k in v)
            vectors[name].append(v)
    # Baseline is only read after the new source expressions and physical
    # normalization have been computed and saved in an independent CAS run.
    base=ROOT/'benchmarks/2ne/theory/generation/three-kernel-eps2.json'
    deps.add(base);old=read(base)
    expected=read(HERE/'runs/r001/baseline.json')['sha256']
    assert sha(base)==expected
    exact=[]
    for target,oldname in (('I20','6e-20'),('I21','6e-21')):
        for q in range(3):
            assert vectors[target][q]==h.lift(old['physical'][oldname][q]),(target,q)
            exact.append({'target':target,'epsilon_power':q,'difference':'0','old_coordinates':len(old['physical'][oldname][q])})
    for name in ('W','Lsum','A','B'):
        for q in range(3):
            assert vectors[name][q]==h.lift(old['bare'][name][q]),(name,q)
            exact.append({'target':name,'epsilon_power':q,'difference':'0'})
    oldrun=ROOT/'workflow/checks/three-kernel-eps2/runs/20260928T031243018668Z'
    for n in (3,4,5):
        for kind in ('W','L'):
            p=oldrun/f'{kind}_num{n}.txt';deps.add(p)
            assert h.parse(run/p.name)==h.parse(p)
            exact.append({'target':f'{kind}_num{n}','source_difference':'0'})
    p2=h.add(h.scale(h.period(2),Q(3,2)),h.scale(h.power(h.LOG,2),Q(-9,8)))
    p3=h.period(3)
    negative=[]
    for name,bare in (('I20','A'),('I21','B')):
        factor=Q(-1,16 if name=='I20' else 32)
        for q in (2,3):
            wrong=h.scale(vectors[bare][q],factor)
            residual=h.add(vectors[name][q],h.scale(wrong,-1))
            assert residual
            predicted=h.scale(h.mul(p2,vectors[bare][q-2]),factor)
            if q==3:predicted=h.add(predicted,h.scale(h.mul(p3,vectors[bare][0]),factor))
            assert residual==predicted
            negative.append({'target':name,'epsilon_power':q,'mutation':'omit Gamma/cosh normalization',
                             'rejected':True,'nonzero_coordinates':len(residual)})
    coefficients={'domain':'z>0; epsilon analytic near zero','orders':[0,1,2,3],
                  'remainder':'O(epsilon^4)','negative_powers_zero':'reused analytic source bound',
                  'period_convention':{'weights':list(h.ZW),'names':['zeta2','zeta3','zeta5','zeta7','zeta(5,3)','zeta9']},
                  'physical':{'6e-20':[h.serialize(v) for v in vectors['I20']],
                              '6e-21':[h.serialize(v) for v in vectors['I21']],
                              '6e-22':[h.serialize(v) for v in vectors['I21']]},
                  'bare':{name:[h.serialize(v) for v in vectors[name]] for name in ('A','B','W','Lsum','H')},
                  'formal':'not_run','source_run':run.relative_to(ROOT).as_posix()}
    path=run/'coefficients.json';assert not path.exists()
    path.write_text(json.dumps(coefficients,indent=2)+'\n',encoding='utf-8')
    report={'passed':True,'formal':'not_run','run':run.relative_to(ROOT).as_posix(),
            'exact_regressions':exact,'negative_controls':negative,
            'new_weight9_coordinates':{name:len(vectors[name][3]) for name in vectors},
            'new_full_physical_epsilon3':True,'independent_new_numeric':'pending',
            'coefficients_sha256':sha(path),
            'inputs_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in sorted(deps)},
            'scope':'Complete new physical coefficient calculation and exact old-order regression; weight-nine numerical and function-space checks separate.'}
    (run/'assembly-verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'passed':True,'exact_regressions':len(exact),'negative_controls':len(negative),
                      'new_weight9_coordinates':report['new_weight9_coordinates']}))

if __name__=='__main__':main()
