"""Independent finite-epsilon and coefficient checks of the intrinsic source."""
from datetime import datetime, timezone
from pathlib import Path
import json
import time
import numpy as np
import intrinsic
from original_adapter import HERE, ROOT, full_epsilon, sha, provenance


def main():
    start=time.perf_counter()
    reference_path=HERE/'original-reference.json'
    reference=json.loads(reference_path.read_text(encoding='utf-8'))
    assert reference['passed'] and not reference['candidate_expression_read']
    for path,expected in reference['inputs_sha256'].items():
        assert sha(ROOT/path)==expected,(path,'frozen dependency changed')
    cases=[];negative=[];coefficient_cases=[];rectangle=[];integer_controls=[]
    full_targets=list(reference['full_epsilon_cases'])
    for e,z in ((.13,.7),(-.05,.7),(0.,.7)):
        vals={str(n):full_epsilon(e,z,n) for n in (28,48)}
        full_targets.append({'epsilon':e,'z':z,'original_values_by_nodes':vals,
                             'observed_refinement':abs(vals['48']-vals['28'])})
    for row in full_targets:
        e,z=row['epsilon'],row['z']
        vals={str(n):intrinsic.full(e,z,n) for n in (96,176)}
        delta=abs(vals['176']-row['original_values_by_nodes']['48'])
        ref=abs(vals['176']-vals['96'])
        tol=max(2e-8,10*ref,10*row['observed_refinement'])
        assert delta<tol,(e,z,delta,tol)
        cases.append({**row,'intrinsic_values_by_nodes':vals,
                      'absolute_difference':delta,'intrinsic_refinement':ref,
                      'tolerance':tol,'passed':True})
        if e:
            wrong=intrinsic.full(e,z,96,wrong_scale=True)
            diff=abs(wrong-row['original_values_by_nodes']['48'])
            assert diff>max(100*tol,1e-5),(e,z,'undetected wrong scale',diff)
            negative.append({'epsilon':e,'z':z,'mutation':'omit half log((1-t)*s) from U',
                             'wrong_value':wrong,'absolute_difference':diff,'rejected':True})
        print('full',e,z,delta,flush=True)
    for z,old in reference['reused_frozen_coefficients'].items():
        levels={}
        for n in (96,160,256):
            levels[str(n)]=intrinsic.coefficients(z,n)
            print('coefficients',z,n,levels[str(n)],flush=True)
        for q in range(5):
            delta=abs(levels['256'][q]-old['D0_to_D4'][q])
            ref=abs(levels['256'][q]-levels['160'][q])
            tol=max(2e-8,10*ref,10*old['observed_refinement'][q])
            assert delta<tol,(z,q,delta,tol)
            coefficient_cases.append({'z':z,'epsilon_power':q,
                    'intrinsic_by_nodes':{n:v[q] for n,v in levels.items()},
                    'original_512':old['D0_to_D4'][q],
                    'original_refinement':old['observed_refinement'][q],
                    'absolute_difference':delta,'intrinsic_refinement':ref,
                    'tolerance':tol,'passed':True})
        wrong=intrinsic.coefficients(z,64,wrong_scale=True)
        for q in (1,2,3,4):
            delta=abs(wrong[q]-old['D0_to_D4'][q])
            assert delta>1e-3
            negative.append({'z':z,'epsilon_power':q,
                'mutation':'omit half log((1-t)*s) from U',
                'wrong_value':wrong[q],'absolute_difference':delta,'rejected':True})
    for e in (-.05,0.,.13,.31):
        for eta,kappa in ((.3,1.2),(-.6,1.5),(.1,1.8)):
            low=intrinsic.rectangle_controls(e,eta,kappa,96)
            high=intrinsic.rectangle_controls(e,eta,kappa,176)
            mirror=float(intrinsic.rectangle_j(e,np.asarray(-eta),np.asarray(kappa)))
            tol=max(2e-10,10*abs(high['quadrature']-low['quadrature']))
            assert high['absolute_difference']<tol,(e,eta,kappa,high,tol)
            assert abs(mirror-high['J_closed'])<1e-12
            rectangle.append({'epsilon':e,'eta':eta,'kappa':kappa,**high,
                  'quadrature_refinement':abs(high['quadrature']-low['quadrature']),
                  'mirror_difference':abs(mirror-high['J_closed']),'tolerance':tol,'passed':True})
    for eta,kappa in ((.3,1.2),(-.6,1.5),(.1,1.8)):
        y,w=intrinsic.grid(96)
        v,c,d=intrinsic.rectangle_interval(np.asarray(eta),np.asarray(kappa),y)
        expected={1:-kappa,2:-kappa*((np.exp(eta)+np.exp(-eta))/2-kappa/4)}
        for e in (1,2):
            value=float(-kappa*np.sum(w*v*np.exp(e*c)*intrinsic.sinhc(e*d/2)))
            diff=abs(value-expected[e]);assert diff<2e-12
            integer_controls.append({'epsilon':e,'eta':eta,'kappa':kappa,
                 'quadrature':value,'exact_polynomial_value':expected[e],
                 'absolute_difference':diff,'passed':True})
    bindings=provenance()
    for path in (Path(__file__),HERE/'intrinsic.py',reference_path):
        bindings[str(path.relative_to(ROOT))]=sha(path)
    report={'passed':True,'formal':'not_run','recorded_at_utc':datetime.now(timezone.utc).isoformat(),
            'full_epsilon_cases':cases,'coefficient_cases':coefficient_cases,
            'rectangle_controls':rectangle,'rectangle_integer_controls':integer_controls,
            'negative_controls':negative,'inputs_sha256':bindings,
            'wall_seconds':time.perf_counter()-start,
            'normalization':'Bare geometric D=R231. No Gamma, color, cosh, or external rho factor.',
            'independence':'Full source is assembled from root intrinsic cosh difference and universal rectangle J formula. Coefficients are extracted before quadrature from an independently integrated rectangle coordinate; the full source and coefficient route do not use original-D code. Comparison is to old original ordered-coordinate kernels and provenance-bound old D0..D4 values.',
            'limits':'Finite real-parameter float tests and observed refinement; no interval certification or all-parameter/formal proof. Coefficients q0..q4 only. Laurent comparison targets are explicitly reused, not freshly reintegrated.'}
    output=HERE/'verification.json';assert not output.exists(),output
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'passed':True,'full_cases':len(cases),'coefficients':len(coefficient_cases),
                     'rectangle_controls':len(rectangle),'integer_controls':len(integer_controls),
                     'negatives':len(negative),'max_full_abs':max(x['absolute_difference'] for x in cases),
                     'max_coefficient_abs':max(x['absolute_difference'] for x in coefficient_cases),
                     'seconds':report['wall_seconds']}),flush=True)


if __name__=='__main__':
    main()
