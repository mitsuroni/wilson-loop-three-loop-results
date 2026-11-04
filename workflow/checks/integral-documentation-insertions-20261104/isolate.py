"""Copy real input/code closure, generate without final references, then check."""
import shutil,sys
from prepare import HERE,ROOT,RUNTIME,read,put,sha
from run_batch import run
def main():
    copy=HERE/'isolation-r001';assert not copy.exists();reader=copy/'integrals';runtime=reader/'families/three-loop-insertions/runtime'
    for name in ('scripts','snapshots'):shutil.copytree(RUNTIME/name,runtime/name)
    shutil.copyfile(RUNTIME/'inputs.json',runtime/'inputs.json')
    for entry in read(HERE/'entries.json'):
        src=ROOT/entry['path']/'code/reproduce.py';dst=reader/'three-loop'/entry['orbit_id']/'code/reproduce.py';dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
    assert not (runtime/'references').exists()
    files={p.relative_to(reader).as_posix():sha(p) for p in reader.rglob('*') if p.is_file()}
    put(HERE/'isolation-input.json',{'files':files,'references_present_before_generation':False,'scope':'Real relocated code/input/dictionary copy. Native A/B operator and source snapshots remain adopted inputs.'})
    out=HERE/'runs/isolation-r001';out.mkdir(parents=True)
    gen=run('generate',[reader/'three-loop/4e-Pi1-01/code/reproduce.py','--output',out/'generation'],out);frozen=sha(out/'generation/generated.json')
    put(out/'generation-freeze.json',{'generated_sha256':frozen,'references_absent_during_generation':True})
    for name in ('inputs','references'):shutil.copytree(RUNTIME/name,runtime/name)
    stages=[]
    for label,script,filename in [('verify','verify_batch.py','verification.json'),('representations','check_representations.py','representations.json')]:
        stages.append(run(label,[runtime/'scripts'/script,'--results',out/'generation/generated.json','--out',out/filename],out))
    primary=HERE/'runs/primary-r001/generation/generated.json';a=read(primary);b=read(out/'generation/generated.json');members=[]
    for oid in a:
        assert a[oid]['laurent']==b[oid]['laurent'] and a[oid]['color_channels']==b[oid]['color_channels'],oid
        members.append({'orbit_id':oid,'scalar_and_color_payload_same':True})
    assert sha(out/'generation/generated.json')==frozen
    put(out/'replay.json',{'status':'passed','generation':gen,'checking':stages,'entries':15,'members':members,'generated_json_same_bytes':sha(primary)==frozen,'references_absent_during_generation':True,'scope':'Same installed Python environment, real copied closure. Adopted source/function/operator inputs shared; no new CAS or formal proof.'})
    print({'isolated_entries':15,'same_payloads':15,'whole_generated_json_same_bytes':sha(primary)==frozen})
if __name__=='__main__':main()
