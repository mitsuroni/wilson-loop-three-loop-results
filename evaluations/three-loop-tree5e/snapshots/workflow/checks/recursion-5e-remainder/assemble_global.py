"""Assemble the independently computed P_low - A25 before reading 5e targets.

All products use absolute coefficient convolution. Old 5e results enter
only the final comparison or the explicitly separated revision publication.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
import importlib.util,sys,json,re
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE));from compare_t3 import read,sha,vec
sys.path.insert(0,str(ROOT/'workflow/checks/single-opposite'));import hlog_algebra as ha
POWERS=range(-6,1)
BINDINGS={}
def bind(p):
    BINDINGS[p.relative_to(ROOT).as_posix()]=sha(p);return p
def run(case):
    ptr=read(bind(HERE/f'current-{case}.json'));p=ROOT/ptr['run']
    assert read(bind(p/'process.json'))['success'];return p
def table(path,prefix,lo,hi):return {k:ha.parse(bind(path/f'{prefix}{k}.txt')) for k in range(lo,hi+1)}
def add(*tables):return {k:ha.add(*(t.get(k,{}) for t in tables)) for k in set().union(*tables)}
def scale(t,c):return {k:ha.scale(v,c) for k,v in t.items()}
def shift(t,n):return {k+n:v for k,v in t.items()}
def mul(a,b,hi):
    out=defaultdict(dict)
    for i,x in a.items():
        for j,y in b.items():
            if i+j<=hi:out[i+j]=ha.add(out[i+j],ha.mul(x,y))
    return dict(out)
def coefficients(f,e,ell,lo,hi):
    expanded=s.expand(f)
    return {k:ha.from_sympy(expanded.coeff(e,k),{ell:ha.LOG}) for k in range(lo,hi+1)}
def serialize(t):return {str(k):ha.serialize(t.get(k,{})) for k in sorted(t)}
def projected(t,inv,n,e,ell,hi,left=0,right=None):
    # Default P_n; optional z^left*z^eps and z^right*eps weighted projection.
    if right is None:right=n
    wl=s.series(s.exp(left*e*ell)/(1+s.exp(n*e*ell)),e,0,hi-min(t)+1).removeO()
    wr=s.series(s.exp(right*e*ell)/(1+s.exp(n*e*ell)),e,0,hi-min(t)+1).removeO()
    return add(mul(coefficients(wl,e,ell,0,hi-min(t)),t,hi),
               mul(coefficients(wr,e,ell,0,hi-min(t)),inv,hi))

def construct():
    # Reuse only the six already checked primitive F moments and Gamma algebra.
    libpath=bind(ROOT/'workflow/checks/tgv-t1-empty/formulas.py')
    spec=importlib.util.spec_from_file_location('primitive_F_and_Gamma',libpath)
    lib=importlib.util.module_from_spec(spec);spec.loader.exec_module(lib)
    e,ell=lib.e,lib.ell;G3,A,B,V,_,Rt=lib.gamma_series();U=2*A-B
    F=sum(lib.moment('F',q)*e**q for q in range(6))
    for q in range(6):bind(lib.run()/f'F{q}.txt')
    G1=lib.exp_series({k:s.zeta(k)/k for k in range(2,7)})
    G2=lib.exp_series({k:2**k*s.zeta(k)/k for k in range(2,7)})
    R=s.series(s.cosh(e*ell/2)*s.cosh(e*ell)/s.cosh(3*e*ell/2),e,0,7).removeO()
    r=coefficients(R,e,ell,0,6)
    parent1=coefficients(G2*U/(32*e**4),e,ell,-4,2)
    bareadj=coefficients(-G1/e**2,e,ell,-2,4)
    ingredient=run('global-ingredients')
    bareJ=table(ingredient,'Jopp',0,4);bareJinv=table(ingredient,'JoppInv',0,4)
    bareopp=scale(mul(coefficients(G1,e,ell,0,4),projected(bareJ,bareJinv,1,e,ell,4),4),Q(-1,2))
    c3=table(ingredient,'T3C',1,4);c3inv=table(ingredient,'T3CInv',1,4)
    parent3=shift(scale(mul(coefficients(G2,e,ell,0,3),projected(c3,c3inv,2,e,ell,4),4),Q(-1,8)),-2)
    parent4=table(run('t4-parent'),'ParentI',-2,2)
    lower=scale(mul(r,mul(add(scale(parent1,2),parent3,parent4),add(bareadj,bareopp),0),0),-1)
    # Direct auxiliary values from gap integrals; none are defined from old I.
    a1adj=coefficients(G3*(5*U+9*U*V*Rt+(1-2*e)*B-4*(1-2*e)*e*F)/(288*e**6),e,ell,-6,0)
    L=table(ingredient,'L',0,4);Linv=table(ingredient,'LInv',0,4)
    a1opp=shift(scale(mul(coefficients(G3*U,e,ell,0,4),projected(L,Linv,3,e,ell,4,left=2,right=1),4),Q(1,32)),-4)
    r3=run('t3-assemble');r4a=run('t4-assemble');r4o=run('t4-opp')
    auxparts={'T1_T2_adjacent':a1adj,'T1_T2_opposite':a1opp,
        'T3_adjacent':table(r3,'Aadj',-6,0),'T3_opposite':table(r3,'Aopp',-6,0),
        'T4_adjacent':table(r4a,'AdjacentAuxiliary',-6,0),'T4_opposite':table(r4o,'OppAuxiliary',-6,0)}
    aux=add(*auxparts.values());rhs=add(lower,scale(aux,-1))
    s1adj=coefficients(G3*(4*U+9*Rt*U*(1-V)-(1-2*e)*B+4*(1-2*e)*e*F)/(288*e**6),e,ell,-6,0)
    groups={'T1_T2_adjacent':s1adj,'T1_T2_opposite':table(run('t1-opp'),'Sopp',-6,0),
        'T3_adjacent':table(r3,'Sadj',-6,0),'T3_opposite':table(r3,'Sopp',-6,0),
        'T4_adjacent':table(run('t4-final'),'AdjacentEffective',-6,0),'T4_opposite':table(r4o,'OppEffective',-6,0)}
    direct=add(*groups.values())
    checks=[]
    def check(name,k,a,b):
        delta=ha.add(a,ha.scale(b,-1));checks.append(dict(check=name,power=k,residual=ha.serialize(delta),passed=not delta))
    for k in POWERS:check('global_P_low_minus_A25_equals_six_group_sum',k,rhs.get(k,{}),direct.get(k,{}))
    for prefix,new,parentprefix,newparent in [('T3',r3,'IT',parent3)]:
        for k,name in [(-1,'ITm1'),(0,'ITzero')]:check('fresh_parent3_old_budget_overlap',k,newparent[k],ha.parse(bind(new/(name+'.txt'))))
    for k,name in [(0,'Igzero'),(1,'Igone')]:check('bare_overlap',k,bareopp[k],ha.parse(bind(r3/(name+'.txt'))))
    # The controls ensure a truncated lower parent and the historical bug fail.
    truncated={k:v for k,v in parent4.items() if k<=0}
    badlower=scale(mul(r,mul(add(scale(parent1,2),parent3,truncated),add(bareadj,bareopp),0),0),-1)
    neg={str(k):ha.serialize(ha.add(badlower.get(k,{}),ha.scale(lower.get(k,{}),-1))) for k in POWERS}
    assert neg['-1'] and neg['0']
    out=dict(format='wilson-recursion-5e-construction/1',formal='not_run',old_5e_targets_read=False,
        passed=all(x['passed'] for x in checks),powers=list(POWERS),checks=checks,
        lower_product=serialize({k:lower.get(k,{}) for k in POWERS}),
        auxiliary_25=serialize({k:aux.get(k,{}) for k in POWERS}),
        auxiliary_groups={n:serialize(t) for n,t in auxparts.items()},
        effective_groups={n:serialize(t) for n,t in groups.items()},
        S5e=serialize({k:rhs.get(k,{}) for k in POWERS}),
        independent_lower_T3=serialize(parent3),independent_lower_T4=serialize(parent4),
        parent_truncation_negative=neg,source_sha256=dict(BINDINGS),
        dependence='Original C/J, Wilson conventions, HyperInt backend and six checked F primitives shared; old two-loop T4 -2..0 used only in opposite subgroup and verified against fresh parent. No old 5e target defines any auxiliary value.')
    (HERE/'construction.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    assert out['passed'],[(r['check'],r['power'],len(r['residual'])) for r in checks if not r['passed']]
    return rhs

def old_target_sum(replacements=None):
    replacements=replacements or {};wide={k:defaultdict(Q) for k in POWERS}
    for i in range(1,45):
        d=read(bind(ROOT/f'benchmarks/tgv/5e-3g/results/5e-3g-{i:02}/r001.json'))
        for k in POWERS:
            if i in replacements:terms=ha.serialize(replacements[i][k])
            elif all('word' in t for t in d['hlog_coefficients'][str(k)]):terms=d['hlog_coefficients'][str(k)]
            else:
                st=d['laurent']['coefficients'][str(k)]
                for n in (2,3,5):st=re.sub(rf'\bzeta_{n}\b',f'zeta({n})',st)
                terms=ha.serialize(ha.from_sympy(s.sympify(st,locals={'E':s.Symbol('E')}),{s.Symbol('ell'):ha.LOG}))
            for t in terms:
                p=t.get('period_powers',t.get('zeta_powers'));p=tuple(p)+(0,)*(6-len(p))
                w=tuple(Q(x) for x in t['word']);wide[k][(w,p)]+=Q(t['coefficient'])
    out={}
    for k,t in wide.items():
        t={key:c for key,c in t.items() if c}
        assert all(not any(p[3:]) and all(a in (0,-1) for a in w) for w,p in t),'Uncancelled additional function'
        out[k]={(tuple(int(a) for a in w),p[:3]):c for (w,p),c in t.items()}
    return out

def compare(rhs):
    runb=run('t4-budget');correct={}
    for representative,pair in [(13,35),(14,34),(15,33)]:
        correct[representative]=table(runb,f'Correct{representative}q',-6,0)
        correct[pair]=correct[representative]
    old=old_target_sum();revised=old_target_sum(correct);rows=[]
    for k in POWERS:
        before=ha.add(rhs[k],ha.scale(old[k],-1));after=ha.add(rhs[k],ha.scale(revised[k],-1))
        rows.append(dict(power=k,old_r001_residual=ha.serialize(before),corrected_residual=ha.serialize(after),passed=not after))
    assert all(r['passed'] for r in rows)
    assert [r['power'] for r in rows if r['old_r001_residual']]==[-1,0]
    out=dict(passed=True,formal='not_run',target_integrals=44,target_positions=7,
        construction_sha256=sha(HERE/'construction.json'),rows=rows,source_sha256=dict(BINDINGS),
        old_results_passed=False,corrected_results_passed=True,
        corrected_individual_ids=[13,14,15,33,34,35],
        scope='Physical effective sum only; cancellation directions and Omega6 are not separately validated.')
    (HERE/'global-comparison.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print('Global S5e: seven corrected Laurent positions match P_low - A25 exactly.')
    print('Archived r001 sum fails only at -1 and 0; six individual revisions required.')

if __name__=='__main__':compare(construct())
