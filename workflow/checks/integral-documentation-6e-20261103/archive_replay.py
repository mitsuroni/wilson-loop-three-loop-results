"""Retain actual compact replay evidence outside ignored transient work dirs."""
import shutil
from prepare import HERE,ROOT,read,put,sha

def copy(source,target):target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
def main():
    primary=HERE/'runs/primary-r001';isolated=HERE/'runs/isolation-r001'
    assert read(primary/'replay.json')['status']=='passed' and read(isolated/'replay.json')['status']=='passed'
    for rel in ('replay.json','generation-freeze.json','verification.json','independent-controls.json','generate.log','verify.log','controls.log','generate-command.json','verify-command.json','controls-command.json','generation/generated.json','generation/run.json','generation/primitive-replay.json'):
        copy(primary/rel,HERE/'evidence/primary'/rel)
    for rel in ('replay.json','generation-freeze.json','verification.json','generate.log','verify.log','generate-command.json','verify-command.json','generation/run.json','generation/primitive-replay.json'):
        copy(isolated/rel,HERE/'evidence/isolation'/rel)
    g=read(primary/'generation/run.json');p=read(primary/'generation/primitive-replay.json');v=read(primary/'verification.json');c=read(primary/'independent-controls.json');i=read(isolated/'replay.json')
    results=[]
    for e in read(HERE/'entries.json'):
        assert sha(ROOT/e['path']/'result.json')==e['result_sha256']
        results.append({'orbit_id':e['orbit_id'],'adopted_result_sha256':e['result_sha256'],'seven_coefficients_verified':True,'isolated_coefficient_payload_same':True})
    summary={'status':'passed','integrals':39,'coefficient_positions':v['coefficient_positions'],'function_back_substitutions':v['function_back_substitutions'],'new_Gamma_expansions':p['scalar_expansion_count'],'archived_primitive_events':len(p['replay_events']),'new_external_CAS_processes':g['new_cas_processes'],'old_final_I_read_by_generation':g['old_final_I_read'],'generation_seconds':g['elapsed_seconds'],'source_controls':len(c['source_checks']),'representative_polynomial_integrations':2*len(c['representative_polynomial_controls']),'original_finite_controls':len(c['finite_original_domain_controls']),'independent_control_negatives':len(c['negative_controls']),'verification_negatives':v['negative_controls'],'finite_alias_scope_control':c['finite_alias_scope_control'],'isolation':{'status':i['status'],'references_absent_before_generation':i['references_absent_during_generation'],'same_coefficient_payloads':i['same_coefficient_payloads']},'members':results,'primary':'evidence/primary','copied_closure':'evidence/isolation','scope':'Actual original-root recognition, fresh endpoint/Gamma expansion or adopted raw subintegral replay, and physical assembly, followed by exact comparisons and scoped controls. Same-environment copied generation excludes whole-I references until verification. Full raw-result JSON provenance may change after relocation. These checks do not grade prose, rerun all original CAS integrals, prove arbitrary complex branches or grant formal status.','formal':'not_run'}
    put(HERE/'replay.json',summary)
    print({k:summary[k] for k in ('status','integrals','coefficient_positions','function_back_substitutions','new_Gamma_expansions','archived_primitive_events')})
if __name__=='__main__':main()
