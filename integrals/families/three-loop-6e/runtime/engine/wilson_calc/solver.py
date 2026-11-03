"""Deterministic dispatch from original ordered roots to supported calculations."""
import copy
import time
from dataclasses import replace
from fractions import Fraction
from interfaces.validate import digest
from .contracts import validate_request,result_shell,validate_result
from .inputs import load_input,InputBudgetError
from .assembly import prepared_normalization,expand_scalar
from .runtime import Runtime,engine_fingerprint,scalar_worker
from .algebra import from_sympy,freeze_math
from .cas_budget import CASBudget,CASBudgetExhausted


def solve(request,*,store=None):
    result=result_shell(request if isinstance(request,dict) else {})
    runtime=None
    try:
        request=validate_request(request);result=result_shell(request)
        data=load_input(request['source'])
        normalization=prepared_normalization(data) if request['quantity']=='I' and data['kind']=='prepared' else None
        if 'assembly' in request and request['assembly']!=normalization:raise ValueError('caller assembly differs from the prepared C recipe')
        from .backends import find_maple_config
        configuration=find_maple_config()
        backend=configuration.fingerprint() if configuration else None
        engine=engine_fingerprint()
        semantic=digest({'J':data['J_ast'],'C':data['orbit'].get('C') if request['quantity']=='I' else None,
                         'loop_order':data['loop_order'],'assembly':normalization,'quantity':request['quantity'],
                         'context':request['context'],'through':request['through'],'engine':engine,'backend':backend,'root_provenance':data['provenance']})
        result['provenance'].update(data['provenance'],semantic_key=semantic,engine=engine,
                                    context=request['context'],normalization=normalization)
        runtime=Runtime(store)
        previous=None
        if request.get('resume_from'):
            previous=runtime.read(request['resume_from'])
            validate_result(previous,semantic_key=semantic,through=request['through'],quantity=request['quantity'])
            if previous['status']=='computed':
                previous=copy.deepcopy(previous);previous['provenance']['resumed_from']=request['resume_from']
                previous['provenance']['source_cas_elapsed_seconds']=previous['provenance'].get('source_cas_elapsed_seconds',previous['provenance'].get('cas_elapsed_seconds',0.0))
                previous['provenance'].update(cas_elapsed_seconds=0.0,cas_process_runs=[])
                return runtime.save(previous)
        if request['cache_policy']=='reuse':
            cached=runtime.lookup(semantic,request['through'],request['quantity'])
            if cached is not None:
                cached=copy.deepcopy(cached);cached['provenance']['cache_hit']=True
                cached['provenance']['source_cas_elapsed_seconds']=cached['provenance'].get('source_cas_elapsed_seconds',cached['provenance'].get('cas_elapsed_seconds',0.0))
                cached['provenance'].update(cas_elapsed_seconds=0.0,cas_process_runs=[])
                return runtime.save(cached)
        result['residuals']=[{'J_ast':data['J_ast'],'source':data['provenance'],'through':request['through'],
                              'context':request['context'],'combination':'original task, no terms dropped'}]
        if request['quantity']=='I' and normalization is None:
            result['reason']={'code':'unsupported_assembly','message':'Explicit research I assembly requires a supported versioned bridge.'}
            return runtime.save(result)
        if request['through']+2*data['loop_order']>request['limits']['max_order']:
            result.update(status='resource_exhausted',reason={'code':'order_budget','message':'Required sufficient jet exceeds max_order.'})
            return runtime.save(result)
        start=time.monotonic()
        cas_budget=CASBudget(request['limits']['seconds'],request['limits']['cas_seconds'],provenance=result['provenance'])
        from .exchange import analyze_exchange,evaluate_empty
        orbit=data['orbit']
        analysis=analyze_exchange(orbit['J'],loop_order=data['loop_order'],
                                  topology=orbit.get('topology'),endpoint_edges=orbit.get('endpoint_edges'),
                                  max_word_terms=request['limits']['max_terms'])
        result['provenance']['structure']=freeze_math(analysis.get('descriptor',{}))
        if analysis['status']!='reduced':
            result['reason']={'code':'unsupported_structure','message':str(analysis.get('reason','Unrecognized ordered structure.'))}
            return runtime.save(result)
        evaluated=evaluate_empty(analysis['descriptor'])
        if evaluated['status']=='computed':
            from .exchange import EPS,Z1
            coefficients,lower=scalar_worker({'J':from_sympy(evaluated['J']),'normalization':normalization,
                'quantity':request['quantity'],'through':request['through'],'loop_order':data['loop_order']},
                max(0.01,request['limits']['seconds']-(time.monotonic()-start)),request['limits']['memory_mb'])
            result.update(status='computed',reason={'code':'complete','message':'Original ordered elementary chains evaluated.'},
                          coefficients=coefficients,complete_powers=list(range(lower,request['through']+1)),lower_bound=lower,residuals=[])
            result['checks'].append({'kind':'source_connection','scope':'ordered adjacent-chain primitive, Gamma recurrence and physical assembly',
                                     'formal':'not_run'})
        elif analysis['descriptor']['opposite_count']==1:
            from .single import evaluate_single
            limited=copy.deepcopy(request)
            limited['limits']['seconds']=max(0.01,request['limits']['seconds']-(time.monotonic()-start))
            if not hasattr(runtime.store,'_session_moments'):runtime.store._session_moments={}
            computed=evaluate_single(analysis,limited,normalization,run_dir=runtime.store.path/'cas-runs',
                                    previous=(previous or {}).get('provenance',{}).get('moments'),kernel_cache=runtime.store._session_moments,cas_budget=cas_budget)
            for field in ('status','reason','coefficients','complete_powers','lower_bound'):
                if field in computed:result[field]=computed[field]
            result['provenance']['moments']=computed.get('moments',{})
            result['provenance']['cas_runs']=computed.get('cas_runs',[])
            result['provenance']['cas_elapsed_seconds']=computed.get('cas_elapsed_seconds',0.0)
            if 'kernel_transport' in computed:result['provenance']['kernel_transport']=computed['kernel_transport']
            if computed.get('i_only_symmetrization'):result['provenance']['i_only_symmetrization']=True
            if 'backend' in computed:result['provenance']['backend']=computed['backend']
            if computed.get('residuals'):result['residuals'][0]['remaining_insertions']=freeze_math(computed['residuals'])
            if result['status']=='computed':result['residuals']=[]
        elif analysis['descriptor']['opposite_count']==data['loop_order'] and request['through']<=0:
            from .ordered_backend import compute_ordered_constant
            answer=compute_ordered_constant(orbit['J'],run_dir=runtime.store.path/'cas-runs',loop_order=data['loop_order'],
                       endpoint_edges=orbit.get('endpoint_edges'),topology=orbit.get('topology'),
                       cas_budget=cas_budget,
                       config=replace(configuration,timeout_seconds=min(120,request['limits']['seconds']),
                           memory_mb=request['limits']['memory_mb'],max_words=request['limits']['max_terms']) if configuration else None)
            if answer['status']=='computed':
                from .hlog import HLogSum
                from .algebra import rational
                pair=answer['coefficients']['0'];raw=HLogSum.from_json(pair['raw'])
                value=raw if request['quantity']=='J' else raw.add(HLogSum.from_json(pair['reflected'])).scale(rational(Fraction(normalization['kappa'])/2))
                result.update(status='computed',reason={'code':'complete','message':'Original ordered opposite kernels integrated at epsilon zero.'},
                              coefficients={'0':value.to_json()} if request['through']==0 else {},
                              complete_powers=[0] if request['through']==0 else [],lower_bound=0,residuals=[])
                result['provenance']['ordered_kernel']=answer.get('provenance',{})
            else:
                result.update(status=answer['status'],reason={'code':answer.get('reason_code','ordered_backend'),'message':answer.get('diagnostic','Ordered integration incomplete.')})
        elif analysis['descriptor']['opposite_count']==2 and analysis['descriptor']['adjacent_count']==1:
            from .ordered import evaluate_ordered
            limited=copy.deepcopy(request)
            limited['limits']['seconds']=max(0.01,cas_budget.seconds_remaining)
            computed=evaluate_ordered(analysis,data,limited,normalization,run_dir=runtime.store.path/'cas-runs',cas_budget=cas_budget)
            for field in ('status','reason','coefficients','complete_powers','lower_bound'):
                if field in computed:result[field]=computed[field]
            result['provenance']['ordered_kernel']=computed.get('backend',{})
            if computed.get('ordered_run'):result['provenance']['ordered_run']=computed['ordered_run']
            if result['status']=='computed':result['residuals']=[]
        else:
            result['reason']={'code':'unsupported_structure','message':'Structure recognized; its coefficient integration route is not yet implemented.'}
        result['provenance']['elapsed_seconds']=time.monotonic()-start
        validate_result(result,semantic_key=semantic,through=request['through'],quantity=request['quantity'])
        return runtime.save(result)
    except InputBudgetError as exc:
        result.update(status='resource_exhausted',reason={'code':'input_budget','message':str(exc)})
        return runtime.save(result) if runtime is not None else result
    except CASBudgetExhausted as exc:
        result.update(status='resource_exhausted',reason={'code':exc.code,'message':str(exc)})
        return runtime.save(result) if runtime is not None else result
    except TimeoutError as exc:
        result.update(status='resource_exhausted',reason={'code':'timeout','message':str(exc)})
        return runtime.save(result) if runtime is not None else result
    except RuntimeError as exc:
        result.update(status='resource_exhausted',reason={'code':'recognizer_or_process_limit','message':str(exc)})
        return runtime.save(result) if runtime is not None else result
    except (ValueError,TypeError,KeyError,OSError,ArithmeticError,SyntaxError,RecursionError) as exc:
        budget='budget' in str(exc).lower()
        result.update(status='resource_exhausted' if budget else 'failed',reason={'code':'algebra_budget' if budget else 'invalid_input_or_backend','message':str(exc)})
        # A failed parse or invalid resume must not replace a valid cache entry.
        return result
