"""Freeze the existing physical baseline and export exact new source words.

No old report is rerun or overwritten. The mathematical source constructors are
imported read-only; this run writes only a fresh numbered directory below HERE.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import importlib
import json
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OLD=ROOT/'workflow/checks/pure-three-generators-20261015'
sys.dont_write_bytecode=True
sys.path.insert(0,str(OLD))
w=importlib.import_module('w_source_check')
l=importlib.import_module('w_lambda_gpl_check')

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text(encoding='utf-8'))
def write_new(path,data):
    with path.open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(data,stream,ensure_ascii=False,indent=2)
        stream.write('\n')

def baseline():
    p=ROOT/'benchmarks/2ne/theory/generation/three-kernel-eps2.json'
    data=read(p)
    assert data['orders']==[0,1,2] and data['negative_powers_zero'] is True
    assert data['physical']['6e-21']==data['physical']['6e-22']
    counts={};coordinates=0
    for group in ('physical','bare'):
        counts[group]={}
        for name,series in data[group].items():
            sizes=[]
            for q,terms in enumerate(series):
                assert q<=2
                seen=set()
                for term in terms:
                    key=(tuple(term['word']),tuple(term['zeta_powers']))
                    assert key not in seen and Fraction(term['coefficient'])!=0
                    seen.add(key)
                    assert set(key[0])<={0,-1} and len(key[1])==5
                    # Original hlog8 convention: zeta2,zeta3,zeta5,zeta7,Z53.
                    assert len(key[0])+sum(a*b for a,b in zip(key[1],(2,3,5,7,8)))==6+q
                    coordinates+=1
                sizes.append(len(terms))
            counts[group][name]=sizes
    return {'passed':True,'scope':'Read-only schema, exact rational/weight and I21=I22 audit of existing results; not a new source-to-physical regression.',
            'source':p.relative_to(ROOT).as_posix(),'sha256':sha(p),
            'known_orders':[0,1,2],'known_physical_epsilon3':False,
            'negative_powers_zero_reused':True,'period_weights':[2,3,5,7,8],
            'counts':counts,'coordinates_checked':coordinates,'formal':'not_run'}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',default='r001')
    parser.add_argument('--max-source-order',type=int,default=6)
    parser.add_argument('--skip-new-numeric',action='store_true')
    args=parser.parse_args()
    assert args.run.startswith('r') and args.run[1:].isdigit()
    assert 3<=args.max_source_order<=7
    out=HERE/'runs'/args.run
    out.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter()
    paths=[Path(__file__).resolve(),OLD/'w_source_check.py',OLD/'w_lambda_gpl_check.py',
           OLD/'w_source_report.json',OLD/'w_lambda_gpl_report.json',
           ROOT/'benchmarks/2ne/theory/generation/three-kernel-eps2.json']
    inputs={p.relative_to(ROOT).as_posix():sha(p) for p in paths}
    write_new(out/'baseline.json',baseline())
    old_w={tuple(row['powers']):row['terms'] for row in read(OLD/'w_source_report.json')['finite_rule_checks']}
    old_l={row['source_degree']:row for row in read(OLD/'w_lambda_gpl_report.json')['checks']}
    f_records=[]
    for total in range(1,args.max_source_order+1):
        for i in range(total):
            for j in range(total-i):
                k=total-i-j
                terms=w.coefficient_words(i,j,k)
                assert all(len(cw)+len(vw)==total+3 for cw,vw in terms)
                if (i,j,k) in old_w:assert len(terms)==old_w[(i,j,k)]
                rows=[{'coefficient':str(c),'constant_word':list(cw),'outer_word':list(vw)}
                      for (cw,vw),c in sorted(terms.items())]
                # Serialization is part of the deliverable, so check exact roundtrip.
                restored={(tuple(r['constant_word']),tuple(r['outer_word'])):Fraction(r['coefficient']) for r in rows}
                assert restored=={key:Fraction(str(v)) for key,v in terms.items()}
                f_records.append({'powers':[i,j,k],'source_degree':total,'weight':total+3,'terms':rows})
    write_new(out/'W-parent-words.json',{
        'definition':'[a^i b^j c^k] F(a,b,c;X,Y,Z), k>=1; all G arguments are 1',
        'letter_map':{'0':'0','A':'1/Y','PX':'1/X','PXY':'Y/X','PZ':'1/Z'},
        'scope':'Generic parent source; physical endpoint projection and boundary assembly not performed.',
        'coefficients':f_records})
    print('W parent exact exports:',len(f_records),flush=True)
    l_records=[];l_counts=[];new_words=None
    for total in range(3,args.max_source_order+1):
        terms=l.lambda_words(total)
        endpoint_counts=l.endpoint_log_audit(terms)
        assert all(len(cw)+len(vw)+(cd+1 if cd else 0)==total+3 for cw,cd,vw in terms)
        if total in old_l:assert len(terms)==old_l[total]['gpl_terms']
        rows=[{'coefficient':str(c),'constant_word':list(cw),'boundary_J_degree':cd,'outer_word':list(vw)}
              for (cw,cd,vw),c in sorted(terms.items())]
        restored={(tuple(r['constant_word']),r['boundary_J_degree'],tuple(r['outer_word'])):Fraction(r['coefficient']) for r in rows}
        assert restored=={key:Fraction(str(v)) for key,v in terms.items()}
        l_records.append({'source_degree':total,'weight':total+3,
                          'endpoint_log_cancellation':endpoint_counts,'terms':rows})
        l_counts.append({'source_degree':total,'weight':total+3,'terms':len(rows)})
        new_words=terms
        print('Lambda exact export:',total,len(rows),flush=True)
    write_new(out/'Lambda-sector-words.json',{
        'definition':'[epsilon^N] complete S_epsilon(A,B,C;d); all G arguments are 1',
        'letter_map':{'0':'0','A':'1/A_phys','PX':'1/C_phys','PXY':'A_phys/C_phys','PZ':'1/B_phys'},
        'boundary_factor':'If boundary_J_degree=k>0, multiply by G(0,(1/d)^k;1); otherwise multiply by 1.',
        'scope':'Generic complete sector. Three physical sectors, mirror and rho prefactor not assembled.',
        'coefficients':l_records})
    numeric=None
    if args.max_source_order>=6 and not args.skip_new_numeric:
        mp=l.mp;mp.mp.dps=60
        A,B,C,d=[-mp.mpf(1)/n for n in (3,4,5,6)]
        ev=w.gpl_numeric_factory(C,A,B,cutoff=140)
        ed=w.gpl_numeric_factory(d,d,d,cutoff=140)
        value=mp.fsum(mp.mpf(str(co))*ev(cw)*ev(vw)*(ed(('0',)+('A',)*cd) if cd else 1)
                      for (cw,cd,vw),co in new_words.items())
        direct=l.coefficient_quad(args.max_source_order,A,B,C,d,nodes=30,terms=100)
        delta=abs(value-direct)
        assert delta<mp.mpf('1e-40'),str(delta)
        numeric={'passed':True,'source_degree':args.max_source_order,'weight':args.max_source_order+3,
                 'arguments':{'A':str(A),'B':str(B),'C':str(C),'d':str(d)},
                 'gpl_value':str(value),'independent_expanded_kernel_quadrature':str(direct),
                 'absolute_difference':str(delta),'digits':60,'gpl_cutoff':140,
                 'quadrature_nodes':30,'J_coefficient_cutoff':100,'rigorous_error_bound':False}
        write_new(out/'new-source-numeric.json',numeric)
    assert all(sha(ROOT/rel)==expected for rel,expected in inputs.items())
    outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()}
    report={'passed':True,'formal':'not_run','run':out.relative_to(ROOT).as_posix(),
            'inputs_sha256':inputs,'outputs_sha256':outputs,'W_parent_coefficients':len(f_records),
            'Lambda_coefficients':l_counts,'new_generic_source_check':numeric,
            'new_physical_epsilon3':False,'new_physical_baseline_regression':False,
            'elapsed_seconds':time.perf_counter()-started,
            'scope':'Explicit exact source exports and a new generic source coefficient; no new physical I coefficient or function-space conclusion.'}
    write_new(out/'verification.json',report)
    print(json.dumps({k:report[k] for k in ('passed','run','W_parent_coefficients','Lambda_coefficients','elapsed_seconds')}))

if __name__=='__main__':main()
