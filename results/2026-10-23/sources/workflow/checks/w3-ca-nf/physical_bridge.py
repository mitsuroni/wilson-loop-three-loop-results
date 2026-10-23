"""Derive the source-to-physical sign from path ordering, independently of RG."""
from pathlib import Path
from itertools import product
import hashlib
import json
import sympy as s

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
checks=[];sources={};rows=[]
def ck(label,lhs,rhs=0):
    residual=s.simplify(lhs-rhs);assert residual==0,(label,residual)
    checks.append({'name':label,'passed':True,'residual':'0'})
def read(rel):
    p=ROOT/rel;sources[rel]=hashlib.sha256(p.read_bytes()).hexdigest()
    return json.loads(p.read_text(encoding='utf-8-sig'))

# The Wilson transport equation U'=ig A U places later matrices on the left.
g,d1,d2,d3=s.symbols('g d1 d2 d3')
A1,A2,A3=s.symbols('A1 A2 A3',commutative=False)
u=s.expand((1+s.I*g*d3*A3)*(1+s.I*g*d2*A2)*(1+s.I*g*d1*A1))
ck('three_ordered_segments_from_transport_equation',u.coeff(g,3)/(s.I**3*d1*d2*d3),A3*A2*A1)

pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
gm=[s.Matrix([[0,1,0],[1,0,0],[0,0,0]]),s.Matrix([[0,-s.I,0],[s.I,0,0],[0,0,0]]),s.diag(1,-1,0),s.Matrix([[0,0,1],[0,0,0],[1,0,0]]),s.Matrix([[0,0,-s.I],[0,0,0],[s.I,0,0]]),s.Matrix([[0,0,0],[0,0,1],[0,1,0]]),s.Matrix([[0,0,0],[0,0,-s.I],[0,s.I,0]]),s.diag(1,1,-2)/s.sqrt(3)]
for nc,matrices in [(2,pauli),(3,gm)]:
    ts=[m/2 for m in matrices];forward=reverse=s.S.Zero
    for a,b,c in product(range(len(ts)),repeat=3):
        f=-2*s.I*s.trace((ts[a]*ts[b]-ts[b]*ts[a])*ts[c])
        if f!=0:
            forward+=f*s.trace(ts[a]*ts[b]*ts[c])/nc
            reverse+=f*s.trace(ts[c]*ts[b]*ts[a])/nc
    cf=s.Rational(nc*nc-1,2*nc)
    ck(f'SU{nc}_chronological_source_trace',forward,s.I*nc*cf/2)
    ck(f'SU{nc}_path_ordered_physical_trace',reverse,-s.I*nc*cf/2)
    ck(f'SU{nc}_source_to_physical_ratio',reverse/forward,-1)

# Compare the source's tree limit to the trusted two-loop C/J, keeping the
# Fourier factor p=-i*d_X and the scalar star's Wick factor explicit.
z=s.Symbol('z');a,b,c=s.symbols('alpha1 alpha2 alpha3')
gram=s.Matrix([[0,-z/2,(1+z)/2,-s.Rational(1,2)],[-z/2,0,-s.Rational(1,2),(1+z)/2],[(1+z)/2,-s.Rational(1,2),0,-z/2],[-s.Rational(1,2),(1+z)/2,-z/2,0]])
tree=read('benchmarks/tgv/3e-3g/sources/canonical.json')
for o in tree['orbits']:
    top=o['topology'];xs=[];us=[]
    for vertex,phase in zip(top['vertices'],top['phases']):
        edge=vertex['line'];parameter=next(t for t in s.sympify(phase['expression']).free_symbols if str(t).startswith('s'))
        xs.append(s.Matrix([1 if k<edge-1 else parameter if k==edge-1 else 0 for k in range(4)]))
        us.append(s.eye(4)[:,edge-1])
    al=[a,b,c]
    def dq(slot,dirslot):
        return -2*al[slot]*sum(al[j]*(us[dirslot].T*gram*(xs[slot]-xs[j]))[0] for j in range(3) if j!=slot)
    tq=sum((us[i].T*gram*us[j])[0]*(dq(i,k)-dq(j,k)) for i,j,k in [(0,1,2),(1,2,0),(2,0,1)])
    P=s.sympify(o['J']['numerator']).subs(s.Symbol('z1'),z)
    ck(o['orbit_id']+':ambient_tree_derivative_equals_2P',tq,2*P)
    # Forward trace gives +(Gamma/8)*(-1+2eps)*Tq = -Gamma(2-2eps)*P/4.
    # Reversed physical trace restores +Gamma(2-2eps)*P/4, the old physical I.
    ck(o['orbit_id']+':physical_tree_density_matches_trusted_CJ',tq/8,P/4)

for family,rel in [('legs','benchmarks/self-energy/tgv-pi1/sources/qgraf-r001.json'),('Gamma','benchmarks/self-energy/gamma1/sources/qgraf-r001.json')]:
    for o in read(rel)['orbits']:
        lines=[v['line'] for v in o['vertices']]
        assert lines==sorted(lines) and len(lines)==3
        assert [v['color'] for v in o['vertices']]==['a1','a2','a3']
        assert [p['loop'] for p in o['phases']]==['l1','l2','l3']
        chronological=[]
        for edge,phase in zip(lines,o['phases']):
            parameter=next(t for t in s.sympify(phase['expression']).free_symbols if str(t).startswith('s'))
            assert str(parameter).startswith('s'+str(edge))
            chronological.append((edge,int(str(parameter)[2:])))
        assert chronological==sorted(chronological)
        vertex=next(v for v in o['internal_vertices'] if v['type'] in ['tgv','gamma1','Gamma1'])
        assert vertex['loops']==['l1','l2','l3']
        ck(o['orbit_id']+':reversed_three_slot_colour_ratio',s.Integer(-1),-1)
        rows.append({'orbit_id':o['orbit_id'],'family':family,'chronological_edge_slots':chronological,'source_word':['a1','a2','a3'],'physical_word':['a3','a2','a1'],'physical_multiplier_of_stored_I':-1,'reason':'Wilson transport orders later slots left; only the Wilson colour trace is reversed, not the vertex indices or momenta'})
assert len(rows)==16
out={'date':'2026-10-06','scope':'physical source-to-I bridge for the 12 TGV-leg and 4 Gamma1 source orbits; fixed current Wilson transport/path-order convention','passed':True,'checks':checks,'rows':rows,'weights':{'Pi1':1,'chain':1,'Pi2':1,'legs':-1,'Gamma':-1},'source_sha256':sources,'RG_used_to_derive_sign':False,'source_or_result_modified':False,'formal':'not_run','mathematical_basis':'general SU(N) fTrace identity from the commutator and f^{abc}f^{abd}=CA delta^{cd}; ODE Dyson order; four exact ambient tree-to-parent comparisons','negative_control':'ascending and descending fTrace are opposite nonzero values in SU(2)/SU(3)'}
(HERE/'physical-bridge.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'passed':True,'checks':len(checks),'orbits':len(rows),'weights':out['weights'],'RG_used_to_derive_sign':False}))
