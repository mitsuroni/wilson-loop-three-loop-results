"""Index precise scalar, operator and dictionary roles for each of15 sources."""
import re
from prepare import HERE,ROOT,RUNTIME,read,put,sha
def ref(rel,obj,role,**more):
    p=RUNTIME/rel;assert p.is_file(),p
    return {'path':rel,'sha256':sha(p),'object':obj,'role':role,**more}
def main():
    pointer=read(RUNTIME/'snapshots/workflow/checks/single-opposite/current-moments.json')['run']
    run='snapshots/'+pointer;global_M=[(k,j) for k in range(4) for j in range(1,5-k)]
    count=0
    for e in read(HERE/'entries.json'):
        oid=e['orbit_id'];own=ROOT/e['path'];old=read(RUNTIME/f'integrals/three-loop/{oid}/subintegrals.json');result=read(own/'result.json');pi1=oid.startswith('4e-')
        mathematical=[];operators=[];code=[ref('scripts/evaluate_batch.py','15-root generator','Reads declared source and admitted function data, then invokes the actual scalar/colour assembly. Does not read whole-I references.')]
        if pi1:
            expression=old['source_paired_J_jet_expression'];needed=sorted({tuple(map(int,t)) for t in re.findall(r'M_(\d+)_(\d+)',expression)})
            for k,j in needed:mathematical.append(ref(run+f'/M_{k}_{j}.txt',f'M_{{{k},{j}}}','Evaluated original no-factorial moment HLOG coordinates required by this route, not a complete-I table.'))
            for atom in ('rho_id_0_0_0','rho_id_0_0_1'):
                if atom in expression:mathematical.append(ref('snapshots/benchmarks/2ne/families/same-channel-two/structural-basis-functions.json',atom,'Evaluated ordered T coefficient, not log(z) and not a complete I.',json_pointer='/functions/'+atom))
            global_reads=[ref(run+f'/M_{k}_{j}.txt',f'M_{{{k},{j}}}','Globally loaded for the eleven-entry batch; does not mean this graph mathematically requires it.') for k,j in global_M]
            global_reads.append(ref('snapshots/benchmarks/2ne/families/same-channel-two/structural-basis-functions.json','rho0/rho1 dictionary','Globally loaded by pi1.primitives; route-specific need is listed separately.'))
            code.append(ref('scripts/pi1.py','source_check, paired_kernels, generate','Checks original tensor/C/J, uses six analytically established route reductions, expands Gamma/rational jets and external B1. The Python program does not derive those six reductions from scratch.'))
            operators=[{'object':'B1(eps)','definition':'[(10-6eps)CA-8(1-eps)TF nf]/(3-2eps)','role':'Regular external connected colour operator; its Laurent/Taylor convolution through the actual scalar pole is required. Original complete Wilson trace is preserved separately.'}]
            budget={k:old[k] for k in ('required_jet_orders','actual_conservative_Gamma_jet','external_operator_Taylor_through')}
        else:
            needed=[] if oid.endswith('-01') else [(0,1),(1,1),(2,1),(0,3)]
            for k,j in needed:mathematical.append(ref(run+f'/M_{k}_{j}.txt',f'M_{{{k},{j}}}','Evaluated opposite H0(3eps,z) coefficients through epsilon2, before geometric and colour convolution.'))
            global_reads=[ref(run+f'/M_{k}_{j}.txt',f'M_{{{k},{j}}}','Loaded by the four-entry batch; adjacent geometry needs no M integral.') for k,j in [(0,1),(1,1),(2,1),(0,3)]]
            code.append(ref('scripts/chains_pi2.py','generate, gamma_jet, rational_coeffs, product','Checks exact original two-endpoint kernel; evaluates Gamma/H0 and all four boundary corners; regenerates U/inverseU and full external colour coefficients.'))
            budget={k:v for k,v in old.items() if 'budget' in k or k=='H0_jet_order'}
            budget['actual_conservative_batch_H0_jet_order']=budget.pop('H0_jet_order')
            budget['required_H0_jet_order_for_this_geometry']=None if oid.endswith('-01') else 2
            budget['actual_conservative_batch_U_inverseU_order']=budget.pop('Pi2_Gamma_ratio_budget_before_cancellation')
            budget['scope']='Conservative shared batch orders are recorded as implementation choices, not mathematical requirements of an adjacent H0-free geometry. Original complete shared metadata remains in archived_subintegrals.'
            source=read(own/'input.json')['source_bindings']['native_prepared']
            rel=source['path'].split('/runtime/',1)[1]
            if oid.startswith('2e-Pi2'):
                operators.append(ref(rel,'Full-D Pi2 A+B/U','Adopted original internal two-loop operator. Current code regenerates phi=4eps²(A+B/U) through the required finite budget and verifies its pole cancellation; it does not read an old operator expansion.',json_pointer=source['field']+'/operator'))
            else:operators.append({'object':'B1(eps)^2','role':'Two separate Pi1 insertions, producing AA/Af/ff weights; not the single two-loop 1PI Pi2 operator.'})
        code.append(ref('snapshots/workflow/checks/single-opposite/hlog_algebra.py','Exact HLOG parsing/shuffle arithmetic','Adopted expression arithmetic, not an integral evaluator.'))
        colors=[{'monomial':c.get('monomial',c.get('color_monomial')),'powers':list(c['coefficients']),'role':'Scalar/geometric I convolved with the complete source external operator.'} for c in result['color_channels']]
        data={'orbit_id':oid,'runtime_base':'../../../families/three-loop-insertions/runtime','paths_relative_to':'runtime_base; original_ordered_root is explicitly relative to this data/index.json file','original_ordered_root':'../input.json','archived_subintegrals':ref(f'integrals/three-loop/{oid}/subintegrals.json','Original source-specific intermediate kernel and order budget','Exact historical reduction metadata, not a substitute for current mathematical prose.'),'required_budget':budget,'mathematically_required_evaluated_functions':mathematical,'batch_global_function_read_closure':global_reads,'external_operator_inputs':operators,'actual_code':code,'color_assembly':colors,'adopted_result':ref(f'integrals/three-loop/{oid}/result.json','Complete curated scalar/geometric and coloured result','Accepted result; preserved separately from generation inputs.'),'comparison_targets':[ref(f'references/{oid}.json','Historical whole-I reference','Post-generation comparison only. Pi1 references contain scalar I, not a historical complete colour table.')] if (RUNTIME/f'references/{oid}.json').is_file() else [],'scope':'Original root, scalar/geometric quantity, evaluated functions, external operator, coloured result and comparison targets have distinct roles. Shared batch loads are not independent graph dependencies or new integrations.'}
        data['comparison_targets']=[ref(f'references/results/{oid}.json','Historical whole-I reference','Opened only after generation for exact comparison. Pi1 has historical scalar I only; its generated colours are a separate operator convolution.')]
        put(own/'data/index.json',data);count+=len(mathematical)
    put(HERE/'data-index-check.json',{'status':'bound','entries':15,'per_entry_evaluated_function_references':count,'scope':'Route-specific required functions are distinguished from globally loaded dictionaries; operator inputs have their own type.'});print({'entries':15,'evaluated_function_references':count})
if __name__=='__main__':main()
