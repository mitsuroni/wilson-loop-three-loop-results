"""Individual intrinsic X moments, crossing, and lower-loop degenerations."""
from datetime import datetime, timezone
from pathlib import Path
import json
import math
import time
import numpy as np
import intrinsic
from original_adapter import HERE, ROOT, original, sha


INDICES = [(m,n,j) for q in range(5) for j in range(1,(q+2)//2+1)
           for m in range(q+3-2*j) for n in [q+2-2*j-m]]


def key(index):
    return ','.join(str(v) for v in index)


def individual_moments(z,nodes,block=8):
    z=np.longdouble(str(z));x,w=intrinsic.grid(nodes)
    s=x[None,:,None];y=x[None,None,:]
    answer={index:np.longdouble(0) for index in INDICES}
    for lo in range(0,nodes,block):
        t=x[lo:lo+block,None,None]
        U,ws,wt,factor,eta,kappa=intrinsic.geometry(t,s,z)
        value,center,width=intrinsic.rectangle_interval(eta,kappa,y)
        raw=(np.ones_like(center),center,center**2+width**2/12,
             center**3+center*width**2/4,
             center**4+center**2*width**2/2+width**4/80)
        # B_n = partial_e^n B_e at zero, without a 1/n! convention.
        Bn=[np.sum(-kappa*value*p*w[None,None,:],axis=2) for p in raw]
        u=U[:,:,0];f=factor[:,:,0]
        L=((ws+wt)/2)[:,:,0];V=((ws-wt)/2)[:,:,0]
        diff_polys={j:sum(L**(2*(j-1-k))*V**(2*k) for k in range(j)) for j in (1,2,3)}
        measure=w[lo:lo+block,None]*w[None,:]
        for m,n,j in INDICES:
            answer[(m,n,j)]+=np.sum(measure*f*u**m*diff_polys[j]*Bn[n])
    return {key(index):float(value) for index,value in answer.items()}


def assemble(moments):
    return [sum(2*moments[key((m,n,j))]/(math.factorial(m)*math.factorial(n)*math.factorial(2*j))
                for m,n,j in INDICES if m+n+2*j==q+2) for q in range(5)]


def vm_integral(e,z,nodes):
    z=np.longdouble(str(z));beta,w=intrinsic.grid(nodes)
    r=z/(1+z);delta=beta-r
    W=np.log(beta*(1-beta))
    ratio=-(1+z)*delta/beta
    L=np.where(abs(ratio)<.3,np.log1p(np.maximum(ratio,-.5)),np.log(z*(1-beta)/beta))
    quotient=np.divide(L,delta,out=np.asarray(-(1+z)/beta).copy(),where=abs(delta)>1e-12)
    return float(np.sum(w*np.exp(e*W/2)*quotient*intrinsic.sinhc(e*L/2)))


def intrinsic_C(e,z,nodes,omit_density=False):
    z=np.longdouble(str(z));x,w=intrinsic.grid(nodes)
    t=x[:,None];s=x[None,:]
    U,ws,wt,factor,eta,kappa=intrinsic.geometry(t,s,z)
    bs=intrinsic.rectangle_j(e,eta,kappa)
    _,_,_,_,other_eta,other_kappa=intrinsic.geometry(s,t,z)
    bt=intrinsic.rectangle_j(e,other_eta,other_kappa)
    h=(1-t)*(1-s)+z*t*s
    density=1 if omit_density else (1+z)/h
    return float(-np.sum(w[:,None]*w[None,:]*density*np.exp(e*U)*bs*bt))


def main():
    start=time.perf_counter()
    reference=HERE/'original-reference.json'
    old=json.loads(reference.read_text(encoding='utf-8'))
    indices=[{'m':m,'n':n,'j':j,'epsilon_power':m+n+2*j-2,
              'reconstruction_factor':f'2/({m}! {n}! {2*j}!)'} for m,n,j in INDICES]
    values={};reconstructions={};controls=[];crossing=[];boundaries=[];c_cases=[];negatives=[]
    z0=np.longdouble('0.7');points=(z0,1/z0)
    for z in points:
        label=str(z);levels={}
        for nodes in (96,160,256):
            levels[str(nodes)]=individual_moments(z,nodes)
            print('moments',label,nodes,assemble(levels[str(nodes)]),flush=True)
        values[label]=levels
        reconstructions[label]={n:assemble(table) for n,table in levels.items()}
    label0,label1=map(str,points)
    for index in INDICES:
        k=key(index)
        a=values[label0]['256'][k];b=values[label1]['256'][k]
        delta=abs(a-b)
        refinement=max(abs(a-values[label0]['160'][k]),abs(b-values[label1]['160'][k]))
        tolerance=max(2e-10*max(1,abs(a),abs(b)),10*refinement)
        assert delta<tolerance,(index,a,b,delta,tolerance)
        crossing.append({'indices':list(index),'z':label0,'reciprocal_z':label1,
                         'X_at_z':a,'X_at_reciprocal':b,'absolute_difference':delta,
                         'observed_refinement':refinement,'tolerance':tolerance,'passed':True})
    target=old['reused_frozen_coefficients']['0.7']
    for label in (label0,label1):
        for q in range(5):
            high=reconstructions[label]['256'][q]
            low=reconstructions[label]['160'][q]
            delta=abs(high-target['D0_to_D4'][q])
            tolerance=max(2e-8,10*abs(high-low),10*target['observed_refinement'][q])
            assert delta<tolerance,(label,q,delta,tolerance)
            controls.append({'z':label,'epsilon_power':q,'X_reconstruction':high,
                             'original_at_z_0_7':target['D0_to_D4'][q],
                             'original_refinement':target['observed_refinement'][q],
                             'X_refinement':abs(high-low),'absolute_difference':delta,
                             'tolerance':tolerance,'passed':True})
    for e in (-.05,.07,.13):
        for z in (.7,1.4):
            eta=-np.log(z)/2;kappa=(1+z)/np.sqrt(z)
            # At this exact boundary the product of the two J arguments is
            # exactly one. Do not reconstruct that endpoint by rounded eta
            # and kappa: a value just above one leaves the real J API domain.
            B=float(-kappa**(-e)/e*(intrinsic.j(e,np.asarray(-z))
                   +intrinsic.j(e,np.asarray(-1/z))
                   -2*intrinsic.j(e,np.asarray(1.))))
            low=vm_integral(e,z,96);high=vm_integral(e,z,176)
            delta=abs(B-high);tolerance=max(2e-10,10*abs(high-low))
            assert delta<tolerance,(e,z,B,high,delta,tolerance)
            boundaries.append({'epsilon':e,'z':z,'rectangle_eta':eta,'rectangle_kappa':kappa,
                      'B_boundary':B,'VM_96':low,'VM_176':high,
                      'absolute_difference':delta,'refinement':abs(high-low),
                      'tolerance':tolerance,'passed':True})
    for e in (0.,.13):
        z=.7;new={str(n):intrinsic_C(e,z,n) for n in (96,176)}
        # The unchanged original C implementation rounds its nested upper
        # endpoint to exactly one at larger grids when epsilon=0. Smaller
        # nodes avoid that representational log(0); the known C0 is a second
        # independent control of the candidate at its full 176-node grid.
        old_nodes=(20,24) if e==0 else (28,48)
        oldC={str(n):original.original_C(e,z,n) for n in old_nodes}
        oldlow,oldhigh=(oldC[str(n)] for n in old_nodes)
        delta=abs(new['176']-oldhigh)
        tolerance=max(2e-8,10*abs(new['176']-new['96']),10*abs(oldhigh-oldlow))
        assert delta<tolerance,(e,delta,tolerance)
        row={'epsilon':e,'z':z,'intrinsic_by_nodes':new,'original_by_nodes':oldC,
             'absolute_difference':delta,'tolerance':tolerance,'passed':True}
        if e==0:
            known=(-((np.log(z)**2+np.pi**2)/2))**3/18
            assert abs(new['176']-known)<2e-9
            row['known_H0_cubed_over_18']=known
            row['known_value_difference']=abs(new['176']-known)
        c_cases.append(row)
        wrong=intrinsic_C(e,z,96,omit_density=True)
        assert abs(wrong-oldhigh)>1e-3
        negatives.append({'target':'C','epsilon':e,'z':z,
                          'mutation':'omit (1+z)/h density','wrong_value':wrong,
                          'absolute_difference':abs(wrong-oldhigh),'rejected':True})
    paths=[Path(__file__),HERE/'intrinsic.py',HERE/'original_adapter.py',reference,
           ROOT/'workflow/checks/pure-three-cd-20261016/numeric_direct.py',
           ROOT/'workflow/checks/pure-three-higher-laurent-20261016/numeric_physical.py',
           ROOT/'workflow/checks/pure-three-higher-laurent-20261016/numeric_shared_physical_generator.py']
    report={'passed':True,'formal':'not_run','recorded_at_utc':datetime.now(timezone.utc).isoformat(),
            'index_convention':'X_mnj = integral U^m (L^(2j)-V^(2j)) partial_e^n B_e|0 / ((t-alpha)(s-alpha)); raw n-th derivative, no factorial inside X.',
            'indices':indices,'individual_X_by_z_and_nodes':values,
            'D0_to_D4_reconstructed_by_z_and_nodes':reconstructions,
            'reconstruction_checks':controls,'individual_crossing_checks':crossing,
            'M_boundary_checks':boundaries,'C_shared_rectangle_checks':c_cases,
            'negative_controls':negatives,'wall_seconds':time.perf_counter()-start,
            'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
            'limits':'Float numerical moments through q4, with empirical refinement and no interval certificate. Independent VM is the original one-dimensional M kernel. C comparison uses unchanged original ordered-coordinate code. Index transport and all-weight identities need analytic reasoning; these finite controls do not establish basis minimality.'}
    out=HERE/'moments-verification.json';assert not out.exists(),out
    out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'passed':True,'moments':len(INDICES),'reconstructions':len(controls),
                     'crossing':len(crossing),'M_boundary':len(boundaries),'C':len(c_cases),
                     'max_coefficient_difference':max(r['absolute_difference'] for r in controls),
                     'max_crossing_difference':max(r['absolute_difference'] for r in crossing),
                     'max_M_boundary_difference':max(r['absolute_difference'] for r in boundaries),
                     'seconds':report['wall_seconds']}),flush=True)


if __name__=='__main__':
    main()
