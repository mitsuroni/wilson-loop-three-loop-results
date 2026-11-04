"""Retain actual evidence outside transient output folders."""
import shutil
from prepare import HERE,ROOT,read,put,sha
def main():
    primary=HERE/'runs/primary-r001';iso=HERE/'runs/isolation-r001'
    assert read(primary/'replay.json')['status']=='passed' and read(iso/'replay.json')['status']=='passed'
    for rel in ('replay.json','generation-freeze.json','verification.json','representations.json','independent-controls.json','generate.log','verify.log','representations.log','controls.log','generate-command.json','verify-command.json','representations-command.json','controls-command.json','generation/generated.json','generation/run.json'):
        dst=HERE/'evidence/primary'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(primary/rel,dst)
    for rel in ('replay.json','generation-freeze.json','verification.json','representations.json','generate.log','verify.log','representations.log','generate-command.json','verify-command.json','representations-command.json','generation/run.json'):
        dst=HERE/'evidence/isolation'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(iso/rel,dst)
    g=read(primary/'generation/run.json');v=read(primary/'verification.json');r=read(primary/'representations.json');c=read(primary/'independent-controls.json');i=read(iso/'replay.json')
    members=[]
    for entry in read(HERE/'entries.json'):
        assert sha(ROOT/entry['path']/'result.json')==entry['result_sha256'];members.append({'orbit_id':entry['orbit_id'],'adopted_result_sha256':entry['result_sha256'],'scalar_and_colour_payload_identical_in_copy':True})
    summary={'status':'passed','entries':15,'scalar_or_geometric_positions':v['I_positions'],'generated_colour_positions':r['positions']-v['I_positions'],'historical_colour_positions_compared':v['colored_positions_compared'],'function_representation_positions':r['positions'],'new_external_CAS_processes':g['new_CAS_processes'],'old_final_tables_read_by_generation':g['old_final_tables_read'],'generation_seconds':g['seconds'],'source_orbits':c['source_orbits'],'paired_polynomial_integrals':c['paired_polynomial_integrals'],'Gamma_finite_controls':len(c['Gamma_finite_coefficient_controls']),'negative_controls':len(c['negative_controls']),'isolation':{'status':i['status'],'references_absent_during_generation':i['references_absent_during_generation'],'generated_json_same_bytes':i['generated_json_same_bytes']},'members':members,'scope':'Original source reductions and fresh scalar/geometric/operator assembly with adopted M/rho dictionaries and full-D Pi2 operator. 105 scalar and84 historically stored colour positions are compared;343 representations include238 generated colour positions. Same-environment isolation is not an independent integral calculation. Text continuity is reviewed separately; formal not_run.'}
    put(HERE/'replay.json',summary);print({k:summary[k] for k in ('status','entries','scalar_or_geometric_positions','generated_colour_positions','historical_colour_positions_compared','function_representation_positions')})
if __name__=='__main__':main()
