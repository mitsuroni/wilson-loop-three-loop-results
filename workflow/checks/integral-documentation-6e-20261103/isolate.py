"""Copy the operational closure, generate without final tables, then verify."""
import shutil,sys
from prepare import HERE,ROOT,RUNTIME,read,put,sha
from run_batch import run

def main():
    copy=HERE/'isolation-r001';assert not copy.exists(),'Use a new isolation revision'
    reader=copy/'integrals';runtime=reader/'families/three-loop-6e/runtime'
    for name in ('engine','primitives','scripts'):shutil.copytree(RUNTIME/name,runtime/name)
    for row in read(HERE/'entries.json'):
        p=ROOT/row['path']/'code/reproduce.py';dst=reader/'three-loop'/row['orbit_id']/'code/reproduce.py';dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dst)
    files={p.relative_to(reader).as_posix():sha(p) for p in reader.rglob('*') if p.is_file()}
    assert not (runtime/'references').exists() and not (runtime/'integrals').exists()
    put(HERE/'isolation-input.json',{'files':files,'reader_root':str(reader),'reference_tables_present_before_generation':False,'scope':'Real operational code/primitive copy. Reader manuscripts and accepted whole-I tables are not needed by generation.'})
    out=HERE/'runs/isolation-r001';out.mkdir(parents=True)
    first=run('generate',[reader/'three-loop/6e-01/code/reproduce.py','--output',out/'generation'],out)
    frozen=sha(out/'generation/generated.json');put(out/'generation-freeze.json',{'generated_sha256':frozen,'references_absent_during_generation':True})
    for name in ('inputs','references'):shutil.copytree(RUNTIME/name,runtime/name)
    second=run('verify',[runtime/'scripts/verify_batch.py','--results',out/'generation/generated.json','--out',out/'verification.json'],out)
    a=read(HERE/'runs/primary-r001/generation/generated.json');b=read(out/'generation/generated.json');members=[]
    for oid in a:
        for key in ('coefficients','complete_powers','lower_bound'):assert a[oid][key]==b[oid][key],(oid,key)
        members.append({'orbit_id':oid,'coefficient_payload_same':True,'complete_powers_same':True,'lower_bound_same':True})
    assert sha(out/'generation/generated.json')==frozen
    put(out/'replay.json',{'status':'passed','generation':first,'verification':second,'integrals':39,'coefficients':273,'members':members,'references_absent_during_generation':True,'same_coefficient_payloads':39,'scope':'Same installed Python environment, relocated real file closure. Full result JSON provenance paths may differ. No second analytic integration, OS sandbox or formal proof.'})
    print({'isolation':'passed','coefficient_payloads_same':39,'references_absent_before_generation':True})
if __name__=='__main__':main()
