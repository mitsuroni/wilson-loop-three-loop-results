"""Preserve prior evidence and bind a new continuous-reading layer for39 roots."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,subprocess,shutil
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
SOURCE=ROOT/'workflow/checks/three-loop-6e-batch-20261103/package'
FAMILY=ROOT/'integrals/families/three-loop-6e';RUNTIME=FAMILY/'runtime'
PUBLIC=ROOT/'.wilson/publications/wilson-loop-three-loop-results'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(p,v):p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())
def git(*args,cwd=ROOT):return subprocess.check_output(['git','-C',str(cwd),'-c','safe.directory='+cwd.as_posix(),'-c','core.longpaths=true',*args],text=True,encoding='utf-8').strip()
def exact(src,dst):
    dst=Path(dst);dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists():assert dst.read_bytes()==src.read_bytes(),dst
    else:shutil.copyfile(src,dst)
def preflight():
    assert not (HERE/'preflight.json').exists()
    previous=ROOT/'workflow/checks/integral-documentation-low-loops-20261103'
    protected=dict(read(previous/'preflight.json')['protected'])
    navigation={'README.md','workflow/master-plan.md','workflow/project-map.md','workflow/data-flow.md','integrals/README.md'}
    for rel,row in read(previous/'manifest.json')['files'].items():
        if rel not in navigation:protected[rel]=row['sha256']
    protected[(previous/'manifest.json').relative_to(ROOT).as_posix()]=sha(previous/'manifest.json')
    for rel,h in protected.items():assert sha(ROOT/rel)==h,('Prior source drift',rel)
    public=read(previous/'preflight.json')['public_files']
    for rel,h in public.items():assert sha(PUBLIC/rel)==h,('Public drift',rel)
    put(HERE/'preflight.json',{'utc':datetime.now(timezone.utc).isoformat(),'owner':'codex-main','task':'WILSON-5.25 pending CLI creation','HEAD':git('rev-parse','HEAD'),'branch':git('branch','--show-current'),'staged':git('diff','--cached','--name-only').splitlines(),'tracked_dirty':git('status','--short','--untracked-files=no').splitlines(),'protected':protected,'public_files':public,'public_HEAD':git('rev-parse','HEAD',cwd=PUBLIC),'allowed_paths':['integrals/three-loop/6e-01..39','integrals/families/three-loop-6e',HERE.relative_to(ROOT).as_posix(),*sorted(navigation),'new task ledger entry'],'scope':'Curation and individual scientific reading audit; original science/results/low-loop reader/frozen packets/public protected. No Git/push or selector change.'})
    print({'protected':len(protected),'public_files':len(public)})
def stage():
    manifest=read(SOURCE/'manifest.json')
    for rel,row in manifest['files'].items():assert sha(SOURCE/rel)==row['sha256'],rel;exact(SOURCE/rel,RUNTIME/rel)
    exact(SOURCE/'manifest.json',RUNTIME/'manifest.json')
    entries=[]
    for row in read(SOURCE/'catalogue.json')['orbits']:
        oid=row['orbit_id'];old=SOURCE/'integrals/three-loop'/oid;own=ROOT/'integrals/three-loop'/oid;own.mkdir(parents=True,exist_ok=True)
        obj=read(old/'input.json');budget=read(old/'subintegrals.json');base='../../families/three-loop-6e/runtime/'
        def bind(ref):return {**ref,'path':base+ref['path']}
        relationship=obj.get('relationship',{}).copy()
        relationship['evidence']=[{**ref,'path':base+'evidence/historical/'+ref['path']} if (RUNTIME/'evidence/historical'/ref['path']).is_file() else {**{k:v for k,v in ref.items() if k!='path'},'historical_project_path':ref['path']} for ref in relationship.get('evidence',[])]
        put(own/'input.json',{'orbit_id':oid,'family':row['family'],'references_relative_to':'this input.json directory','source_bindings':{k:bind(v) for k,v in obj['source_bindings'].items()},'archived_input':{'path':base+f'integrals/three-loop/{oid}/input.json','sha256':sha(old/'input.json')},'original_ordered_J':budget['root_J'],'normalization':budget['normalization'],'relationship':relationship,'source_packet_manifest':{'path':base+'manifest.json','sha256':sha(SOURCE/'manifest.json')},'scope':'Initial source fields and exact ordered root are explicit; full archived metadata and historical provenance retain their original identity.'})
        exact(old/'result.json',own/'result.json')
        for sub in ('code','data','checks'):(own/sub).mkdir(exist_ok=True)
        (own/'code/reproduce.py').write_text('''"""Run the whole39-entry same-batch generation, then select this member."""
from pathlib import Path
import runpy,sys
runtime=Path(__file__).resolve().parents[3]/'families/three-loop-6e/runtime'
sys.path.insert(0,str(runtime/'scripts'))
runpy.run_path(str(runtime/'scripts/evaluate_batch.py'),run_name='__main__')
''',encoding='utf-8')
        put(own/'data/index.json',{'orbit_id':oid,'runtime_base':'../../../families/three-loop-6e/runtime','archived_subintegral_budget':{'path':f'integrals/three-loop/{oid}/subintegrals.json','sha256':sha(old/'subintegrals.json')},'scope':'Object-level primitive and code-role bindings are added after the actual new batch run; historical absolute paths are provenance, not live output destinations.'})
        put(own/'checks/index.json',{'orbit_id':oid,'historical_checks':{'path':f'../../../families/three-loop-6e/runtime/integrals/three-loop/{oid}/checks.json','sha256':sha(old/'checks.json')},'actual_replay':'../../../../workflow/checks/integral-documentation-6e-20261103/replay.json','scientific_reading_audit':'../../../../workflow/checks/integral-documentation-6e-20261103/review/README.md','scope':'Task, reading continuity, source calculation, replay and formal evidence remain separate.'})
        entries.append({'orbit_id':oid,'family':row['family'],'path':own.relative_to(ROOT).as_posix(),'result_sha256':sha(own/'result.json')})
    put(HERE/'entries.json',entries);put(HERE/'runtime-copy.json',{'files':len(manifest['files'])+1,'source_manifest_sha256':sha(SOURCE/'manifest.json'),'scope':'Full declared source packet copied without scientific patches; archived Markdown remains historical.'})
    print({'entries':len(entries),'runtime_files':len(manifest['files'])+1})
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('phase',choices=['preflight','stage']);a=p.parse_args();globals()[a.phase]()
