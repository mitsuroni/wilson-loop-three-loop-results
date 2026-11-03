"""Physical assembly of directly integrated ordered gap coefficients."""
import time
import sympy as sp
from .algebra import gamma_series
from .assembly import scalar_hlog
from .hlog import HLogSum
from .exchange import EPS,Z1
from .backends import find_maple_config
from .ordered_backend import compute_ordered_moments
from .cas_budget import CASBudget


def evaluate_ordered(analysis,data,request,normalization,*,run_dir,cas_budget=None):
    budget=request['through']+2
    if not 0<=budget<=2:
        return {'status':'unsupported','reason':{'code':'ordered_jet_scope','message':'This ordered gap adapter currently covers I/J through epsilon zero.'}}
    config=find_maple_config(max_words=request['limits']['max_terms'],memory_mb=request['limits']['memory_mb'],
                             timeout_seconds=min(120,request['limits']['seconds']))
    if config is None:return {'status':'failed','reason':{'code':'backend_unavailable','message':'Maple/HyperInt not configured.'}}
    cas_budget=cas_budget or CASBudget(request['limits']['seconds'],request['limits'].get('cas_seconds',request['limits']['seconds']))
    orbit=data['orbit']
    answer=compute_ordered_moments(orbit['J'],max_order=budget,run_dir=run_dir,loop_order=data['loop_order'],
                endpoint_edges=orbit.get('endpoint_edges'),topology=orbit.get('topology'),config=config,cas_budget=cas_budget)
    complete=-1
    for q in range(budget+1):
        if str(q) in answer.get('coefficients',{}):complete=q
        else:break
    values={p:HLogSum([]) for p in range(-2,complete-1)}
    n=data['loop_order'];r=analysis['descriptor']['opposite_count'];adj=analysis['descriptor']['adjacent_z_count']
    if request['quantity']=='I':
        g=sp.exp(-sp.EulerGamma*EPS)*sp.gamma(1-EPS)
        common=sp.Rational(normalization['kappa'])*g**n/(1+Z1**(n*EPS))
        channels=(('raw',common*(-1-Z1)**r*Z1**(adj*EPS)),
                  ('reflected',common*(-1-1/Z1)**r*Z1**((n-adj)*EPS)))
    else:channels=(('raw',(-1-Z1)**r*Z1**(adj*EPS)),)
    for channel,prefactor in channels:
        jet=gamma_series(prefactor,EPS,max(0,complete))
        for q in range(complete+1):
            h=HLogSum.from_json(answer['coefficients'][str(q)][channel])
            for k in range(complete-q+1):
                scalar=HLogSum.from_json(scalar_hlog(jet.coeff(EPS,k),Z1))
                power=q+k-2
                values[power]=values[power].add(h.multiply(scalar,max_terms=request['limits']['max_terms']))
    return {'status':answer['status'],'reason':{'code':'complete' if answer['status']=='computed' else answer.get('reason_code','ordered_backend'),
                'message':answer.get('diagnostic','Original ordered gap kernel and physical convolution evaluated.')},
            'coefficients':{str(p):v.to_json() for p,v in values.items()},'complete_powers':list(values),
            'lower_bound':-2,'backend':answer.get('provenance',{}),'ordered_run':answer.get('run_directory'),
            'cas_elapsed_seconds':cas_budget.elapsed}
