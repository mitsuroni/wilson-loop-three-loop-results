"""Finite coefficient generation from actual one-opposite endpoint weights."""
from pathlib import Path
import time
import sympy as sp
from .algebra import gamma_series
from .assembly import scalar_hlog
from .exchange import single_opposite_weights,reflect_single_weight,single_weight_terms,EPS,Z1
from .hlog import HLogSum,HLogError
from .backends import compute_single_opposite_moments,find_maple_config
from interfaces.validate import digest
from .cas_budget import CASBudget


def evaluate_single(analysis,request,normalization,*,run_dir,previous=None,kernel_cache=None,cas_budget=None):
    weight=single_opposite_weights(analysis)
    n=analysis['descriptor']['loop_order']
    phase=sp.Rational(analysis['descriptor']['adjacent_z_count'])+sp.Rational(1,2)
    i_symmetrized=False
    if (request['quantity']=='I' and weight['status']=='partial' and 2*phase==n
            and weight.get('residuals') and all(r.get('kind')=='ordered_BA_insertion' for r in weight['residuals'])):
        # t->1-t,s->1-s sends H_W(1/z) to H_R(W)(z). Equal
        # physical phase weights permit W -> (W+R(W))/2 only at I level.
        symmetric=sp.expand((weight['normalized_weight']+reflect_single_weight(weight['normalized_weight']))/2)
        try:
            terms=single_weight_terms(symmetric)
        except ValueError:
            terms=None
        if terms is not None:
            weight['terms']=terms
            weight['status']='reduced';i_symmetrized=True
    if weight['status']!='reduced':
        return {'status':'unsupported','reason':{'code':'ordered_insertion_remaining','message':'An ordered insertion or kernel transport remains.'},
                'residuals':weight.get('residuals',[]),'coefficients':{},'complete_powers':[],'moments':{}}
    transport=weight.get('transport','identity')
    if transport not in ('identity','z-inversion'):raise ValueError('unknown single-kernel transport')
    raw,reflected=('reflected','raw') if transport=='z-inversion' else ('raw','reflected')
    config=find_maple_config(max_words=int(request['limits']['max_terms']),memory_mb=request['limits']['memory_mb'])
    if config is None:
        return {'status':'failed','reason':{'code':'backend_unavailable','message':'Configure WILSON_MAPLE and WILSON_HYPERINT.'},
                'coefficients':{},'complete_powers':[],'moments':{}}
    pole=2*(n-1);budget=request['through']+pole
    if budget<0:
        return {'status':'computed','reason':{'code':'complete','message':'Target is below an established pole bound.'},
                'coefficients':{},'complete_powers':[],'lower_bound':-pole,'moments':{}}
    cas_budget=cas_budget or CASBudget(request['limits']['seconds'],request['limits'].get('cas_seconds',request['limits']['seconds']))
    stop=None;runs=[];deadline=cas_budget.deadline
    fingerprint=config.fingerprint()
    previous=dict(previous or {});moments={};groups=[]
    allowed={','.join(map(str,t['powers'])) for t in weight['terms']}
    if set(previous)-allowed:raise ValueError('checkpoint has a different endpoint-weight family')

    def validate_moments(known,maximum=None):
        if not isinstance(known,dict):raise ValueError('invalid checkpoint moment mapping')
        for power,pair in known.items():
            if (not isinstance(power,str) or not power.isdecimal() or str(int(power))!=power
                    or (maximum is not None and int(power)>maximum)
                    or not isinstance(pair,dict) or set(pair)!={'raw','reflected'}):
                raise ValueError('invalid checkpoint moment precision')
            for channel in pair.values():
                parsed=HLogSum.from_json(channel)
                if set(parsed.parameters)-{'z'}:raise ValueError('checkpoint moment parameter mismatch')
                if any(g.to_json()['branch']!='positive-real' for t in parsed.terms for g in (t.gpl,*t.periods)):
                    raise ValueError('checkpoint moment branch mismatch')
        return dict(known)

    # Validate *all* resume/session entries before starting any new CAS work.
    for term in weight['terms']:
        parameters=term['powers'];key=','.join(map(str,parameters))
        known=validate_moments(previous.get(key,{}),budget)
        session_key=digest({'definition':'centered-single-opposite/0.1','parameters':parameters,'backend':fingerprint})
        if kernel_cache is not None:
            cached=validate_moments(kernel_cache.get(session_key,{}))
            for q,pair in cached.items():
                if int(q)<=budget:
                    if q in known and known[q]!=pair:raise ValueError('resume and session moments disagree')
                    known[q]=pair
        moments[key]=known
        groups.append((term,parameters,key,session_key))
    values={};complete=-1;scalar_cache={}
    if request['quantity']=='I':
        g=sp.exp(-sp.EulerGamma*EPS)*sp.gamma(1-EPS)
        common=sp.Rational(normalization['kappa'])*g**n/(1+Z1**(n*EPS))
        channels=((raw,phase),(reflected,n-phase))
    else:common=sp.Integer(1);channels=((raw,phase),)
    # For positive z, D_tilde(z)=z*D(1/z).  The centered normalization
    # H_W(1/z)=-(1+z)*z^(eps/2-1)*integral D(1/z)^(eps-1)*W
    # absorbs that complete z^(eps-1) factor.  Only the H channels exchange;
    # the phase adj+1/2 and the physical reflection weight stay unchanged.

    def assemble(order):
        value=HLogSum([])
        for term,_,key,_ in groups:
            for channel,exponent in channels:
                if time.monotonic()>deadline:raise TimeoutError('finite convolution time budget exhausted')
                scalar_key=(term['coefficient'],exponent)
                scalar=scalar_cache.get(scalar_key,[])
                if len(scalar)<=order:
                    prefactor=sp.expand(gamma_series(common*term['coefficient']*Z1**(exponent*EPS),EPS,order))
                    scalar=[HLogSum.from_json(scalar_hlog(prefactor.coeff(EPS,k),Z1)) for k in range(order+1)]
                    scalar_cache[scalar_key]=scalar
                for q in range(order+1):
                    if time.monotonic()>deadline:raise TimeoutError('finite convolution time budget exhausted')
                    hlog=HLogSum.from_json(moments[key][str(q)][channel])
                    value=value.add(hlog.multiply(scalar[order-q],max_terms=request['limits']['max_terms']))
                    if len(value.terms)>request['limits']['max_terms']:raise HLogError('assembled sum exceeds term budget')
        return value

    # Complete one order across every endpoint weight before asking for the
    # next.  Commit its physical convolution immediately, so a later CAS
    # timeout cannot erase coefficients already assembled from a common prefix.
    for q in range(budget+1):
        for _,parameters,key,session_key in groups:
            if str(q) in moments[key]:continue
            remaining=deadline-time.monotonic()
            if remaining<=0:
                stop={'status':'resource_exhausted','diagnostic':'total coefficient budget exhausted'};break
            cas_remaining=cas_budget.cas_remaining
            if cas_remaining<=0:
                stop={'status':'resource_exhausted','reason_code':'cas_budget',
                      'diagnostic':'CAS budget exhausted; cached moments and assembled coefficients retained'};break
            started_cas=cas_budget.elapsed
            generated=compute_single_opposite_moments(parameters,budget,run_dir=Path(run_dir),config=config,
                       orders=[q],total_seconds=min(remaining,cas_remaining),cas_budget=cas_budget)
            elapsed=cas_budget.elapsed-started_cas
            runs.append({'parameters':list(parameters),'directory':generated.get('run_directory'),
                         'orders':[q],'elapsed_seconds':elapsed,
                         'complete_orders':generated.get('complete_orders',[]),'status':generated['status']})
            produced=validate_moments(generated['coefficients'],budget)
            if set(produced)-{str(q)}:raise ValueError('backend returned an unrequested moment order')
            moments[key].update(produced)
            if kernel_cache is not None:kernel_cache.setdefault(session_key,{}).update(moments[key])
            if generated['status']!='computed':stop=generated;break
            if str(q) not in produced:
                stop={'status':'failed','diagnostic':'backend marked an absent moment computed'};break
        if all(str(q) in moments[key] for _,_,key,_ in groups):
            try:
                coefficient=assemble(q)
            except TimeoutError as exc:
                stop={'status':'resource_exhausted','reason_code':'assembly_time_budget','diagnostic':str(exc)}
            except HLogError as exc:
                if 'budget' not in str(exc).lower():raise
                stop={'status':'resource_exhausted','reason_code':'word_budget','diagnostic':str(exc)}
            else:
                values[q-pole]=coefficient;complete=q
        if stop:break
    status='computed' if complete==budget else (stop or {}).get('status','partial')
    if status=='unsupported':reason='unsupported_backend_output'
    elif status=='failed':reason='backend_failed'
    elif status=='resource_exhausted':reason='timeout_or_word_budget'
    else:reason='complete' if status=='computed' else 'missing_moments'
    return {'status':status,'reason':{'code':(stop or {}).get('reason_code',reason),'message':(stop or {}).get('diagnostic','Centered moments and physical finite convolution evaluated.')},
            'coefficients':{str(p):v.to_json() for p,v in values.items()},'complete_powers':list(values),
            'lower_bound':-pole,'moments':moments,'backend':fingerprint,
            'i_only_symmetrization':i_symmetrized,'kernel_transport':transport,
            'cas_runs':runs,'cas_elapsed_seconds':cas_budget.elapsed}
