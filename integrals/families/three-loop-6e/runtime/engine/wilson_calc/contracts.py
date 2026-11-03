"""Versioned calculation contracts. These do not admit mathematical evidence."""
import copy
import re

REQUEST_FORMAT = 'wilson-calculation-request/0.1'
RESULT_FORMAT = 'wilson-calculation-result/0.1'
CONTEXT = {'domain':'z > 0','branch':'positive-real','ordered_measure':'iterated_real_interval_integral'}


def validate_request(value):
    if not isinstance(value,dict): raise ValueError('request object required')
    allowed = {'format','source','quantity','epsilon','through','context','representation','limits','cache_policy','resume_from','assembly'}
    if set(value)-allowed: raise ValueError('unknown request fields: '+','.join(sorted(set(value)-allowed)))
    if value.get('format') != REQUEST_FORMAT: raise ValueError('unsupported request version')
    if type(value.get('through')) is not int: raise ValueError('through must be an integer')
    if value.get('quantity','I') not in ('I','J'): raise ValueError('quantity must be I or J')
    if value.get('epsilon','eps') != 'eps': raise ValueError('only the declared eps regulator is supported')
    if value.get('representation','hlog') != 'hlog': raise ValueError('unsupported output representation')
    context = value.get('context',CONTEXT)
    if context != CONTEXT: raise ValueError('unsupported domain or branch; positive-real context required')
    source = value.get('source')
    if not isinstance(source,dict) or source.get('kind') not in ('prepared','explicit'): raise ValueError('source kind required')
    limit = value.get('limits',{})
    if not isinstance(limit,dict) or set(limit)-{'seconds','cas_seconds','max_order','max_terms','memory_mb'}: raise ValueError('invalid limits')
    for key,x in limit.items():
        if type(x) not in (int,float) or not (0 <= x < float('inf') if key=='cas_seconds' else 0 < x < float('inf')): raise ValueError('positive finite limits required')
        if key not in ('seconds','cas_seconds') and type(x) is not int:raise ValueError('integer resource limit required')
    if value.get('cache_policy','reuse') not in ('reuse','recompute'): raise ValueError('invalid cache policy')
    if 'resume_from' in value and (not isinstance(value['resume_from'],str) or not re.fullmatch('[0-9a-f]{64}',value['resume_from'])): raise ValueError('invalid checkpoint digest')
    if 'resume_from' in value and value.get('cache_policy','reuse')=='recompute':raise ValueError('recompute cannot consume a historical checkpoint')
    result = copy.deepcopy(value)
    result.update(quantity=value.get('quantity','I'),epsilon='eps',context=copy.deepcopy(CONTEXT),representation='hlog',
                  limits={'seconds':120,'cas_seconds':limit.get('seconds',120),'max_order':12,'max_terms':100000,'memory_mb':2048,**limit},cache_policy=value.get('cache_policy','reuse'))
    return result


def result_shell(request):
    return {'format':RESULT_FORMAT,'status':'unsupported','reason':{'code':'unsupported_structure','message':''},
            'quantity':request.get('quantity','I'),'through':request.get('through'), 'coefficients':{},
            'complete_powers':[],'lower_bound':None,'residuals':[], 'checks':[],
            'formal_state':'not_run','provenance':{}}


def validate_result(value,*,semantic_key,through,quantity):
    from .hlog import HLogSum
    if not isinstance(value,dict) or value.get('format')!=RESULT_FORMAT:raise ValueError('unknown cached result version')
    if value.get('status') not in ('computed','partial','unsupported','resource_exhausted','failed'):raise ValueError('invalid cached status')
    if value.get('quantity')!=quantity or type(value.get('through')) is not int or value['through']!=through:raise ValueError('cached target mismatch')
    if value.get('provenance',{}).get('semantic_key')!=semantic_key:raise ValueError('cached context mismatch')
    if value.get('formal_state')!='not_run':raise ValueError('calculation cache cannot grant formal evidence')
    powers=value.get('complete_powers');coefficients=value.get('coefficients');lower=value.get('lower_bound')
    if not isinstance(powers,list) or any(type(p) is not int for p in powers) or len(set(powers))!=len(powers):raise ValueError('invalid completed powers')
    if not isinstance(coefficients,dict) or set(coefficients)!={str(p) for p in powers}:raise ValueError('coefficient inventory mismatch')
    if lower is not None and type(lower) is not int:raise ValueError('invalid lower bound')
    if any(p>through or lower is None or p<lower for p in powers):raise ValueError('cached coefficient outside declared precision')
    if value['status']=='computed' and (lower is None or sorted(powers)!=list(range(lower,through+1)) or value.get('residuals')):
        raise ValueError('incomplete result labeled computed')
    for coefficient in coefficients.values():
        expression=HLogSum.from_json(coefficient)
        if set(expression.parameters)-{'z'}:raise ValueError('cached HLOG parameter differs from the target context')
        for term in expression.terms:
            for gpl in (term.gpl,*term.periods):
                if gpl.to_json()['branch']!='positive-real':raise ValueError('cached HLOG branch differs from the target context')
    return value
