"""Describe actual per-entry primitive consumption, with portable exact refs."""
from pathlib import Path
from prepare import HERE,ROOT,RUNTIME,read,put,sha

def bound(rel,object,role,**other):
    p=RUNTIME/rel;assert p.is_file(),p
    return {'path':rel,'sha256':sha(p),'object':object,'role':role,**other}

def main():
    replay=read(HERE/'runs/primary-r001/generation/primitive-replay.json');generated=read(HERE/'runs/primary-r001/generation/generated.json');count=0;unique_raw=set();owners={event_id:oid for oid,ids in replay['events_by_orbit'].items() for event_id in ids};single={}
    for event_id,event in enumerate(replay['replay_events']):
        if event['operation']=='centered-single-opposite-moments':
            source=Path(event['source_directory']);request=read(source.parent/'request.json')
            for q in event['orders']:
                key=(','.join(map(str,request['parameters'])),str(q));assert key not in single,key;single[key]=event_id
    for entry in read(HERE/'entries.json'):
        oid=entry['orbit_id'];own=ROOT/entry['path'];old=read(RUNTIME/f'integrals/three-loop/{oid}/subintegrals.json')
        units=[];primitives=[];codes=[]
        codes.append(bound('scripts/evaluate_batch.py','batch entry','Calls the original-root recognizers and actual assembly, with only archived CAS and output adapters.'))
        codes.append(bound('engine/wilson_calc/assembly.py','physical_expression / expand_scalar','Actual kappa, Gamma, crossing and finite-order assembly.'))
        if entry['family']=='empty':codes.append(bound('engine/wilson_calc/exchange.py','whole_word_value / evaluate_empty','Recognize this original source and evaluate its ordered endpoint chains using Beta/Gamma identities.'))
        elif entry['family']=='single-opposite':
            codes.append(bound('engine/wilson_calc/exchange.py','single_opposite_weights','Derive this source-specific endpoint weight.'))
            codes.append(bound('engine/wilson_calc/single.py','evaluate_single','Use required weighted-moment jets and assemble this original physical I.'))
        else:codes.append(bound('engine/wilson_calc/ordered.py','evaluate_ordered','Derive the ordered gap kernel or keep the ordinary all-opposite J0, then assemble this I.'))
        new_event_ids=replay['events_by_orbit'][oid]
        if entry['family']=='single-opposite':
            event_ids=sorted({single[(parameters,q)] for parameters,orders in generated[oid]['provenance']['moments'].items() for q in orders})
        else:event_ids=new_event_ids
        for event_id in event_ids:
            event=replay['replay_events'][event_id];source=Path(event['source_directory']);rel=source.relative_to(RUNTIME).as_posix();operation=event['operation']
            if operation=='centered-single-opposite-moments':kind='centered weighted opposite moment'
            elif operation=='ordered-opposite-gap-moments':kind='regular ordered gap jet; pole and external powers excluded'
            else:kind='original ordinary ordered J at epsilon0; not normalized I'
            request=source/'request.json' if (source/'request.json').is_file() else source.parent/'request.json'
            result=source/'result.json' if (source/'result.json').is_file() else source.parent/'result.json'
            units.append({'operation':operation,'orders':event['orders'],'first_replayed_by_orbit':owners[event_id],'from_shared_batch_cache':event_id not in new_event_ids,'request':bound(request.relative_to(RUNTIME).as_posix(),kind+' mathematical request','Exact definition, domain, variable order and backend binding used by the original integration.'),'execution_record':bound(result.relative_to(RUNTIME).as_posix(),kind+' historical execution record','Records original CAS execution; success is not a mathematical proof. No new CAS is run now.')})
            codes.append(bound(rel+'/compute.mpl',kind+' original HyperInt integration','Historical actual integrator. Its archived absolute CAS/library paths are provenance and are not executed by current Python replay.'))
            for raw in event['raw_files']:
                p=Path(raw['path']);assert sha(p)==raw['sha256']
                primitives.append(bound(p.relative_to(RUNTIME).as_posix(),kind,'Adopted original integrated raw output, parsed and connected to this current request before assembly.',order=raw['order'],channel=raw['channel'],operation=operation))
                unique_raw.add(p.relative_to(RUNTIME).as_posix())
        if units:codes.append(bound('scripts/primitive_replay.py','PrimitiveReplay','Checks mathematical request, regenerated script, original execution and raw parser, then returns adopted primitive jets.'))
        index={'orbit_id':oid,'runtime_base':'../../../families/three-loop-6e/runtime','paths_relative_to':'runtime_base; the original_ordered_root field is explicitly relative to this data/index.json file','original_ordered_root':'../input.json','archived_subintegral_budget':bound(f'integrals/three-loop/{oid}/subintegrals.json','original root, reductions and epsilon budget','Historical exact root/consumption record. Absolute historical paths inside it are not current output paths.'),'required_budget':{'target_powers':old['all_requested_powers'],'rule':old['budget_rule'],'normalization':old['normalization']},'new_primitive_events_in_this_orbit':len(new_event_ids),'primitive_units':units,'adopted_integrated_primitives':primitives,'actual_code':codes,'comparison_targets':[bound(f'references/results/{oid}.json','historical complete I','Only opened after generation for comparison; not a primitive input.'),bound('references/basis-and-conversions.json','moment dictionary plus target conversion table','Used by reading/verification; complete target rows are not opened by generation.')],'scope':'All required primitive dependencies are listed, including shared batch-cache inputs first replayed by another orbit. A lack of new events does not mean no mathematical dependency. Historical root, subintegrals, whole I, representation and checks remain distinct; this index is not a proof receipt.'}
        put(own/'data/index.json',index);count+=len(primitives)
    put(HERE/'data-index-check.json',{'status':'bound','entries':39,'per_entry_raw_primitive_references':count,'distinct_raw_primitive_files':len(unique_raw),'scope':'Exact object-level dependencies from actual primary replay, including transitive batch-cache inputs. Repeated shared raw objects are not independent integrals.'})
    put(HERE/'task-binding.json',{'task':'WILSON-5.25','path':'workflow/backlog/tasks/wilson-5.25 - Complete-continuous-derivations-and-audit-all39-three-loop6e-integrals.md','note':'Created with pinned Backlog after initial sandbox lock-directory failure; no commit/staging/remote operation.'})
    print({'entries':39,'raw_primitive_references':count})
if __name__=='__main__':main()
