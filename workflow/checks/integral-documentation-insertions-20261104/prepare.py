"""Bind the15 current articles while preserving original scientific objects."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,shutil,subprocess
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
SOURCE=ROOT/'workflow/checks/three-loop-insertions-p5-20261103/package'
FAMILY=ROOT/'integrals/families/three-loop-insertions';RUNTIME=FAMILY/'runtime'
PUBLIC=ROOT/'.wilson/publications/wilson-loop-three-loop-results'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def put(p,v):p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())
def git(*args,cwd=ROOT):return subprocess.check_output(['git','-C',str(cwd),'-c','safe.directory='+cwd.as_posix(),'-c','core.longpaths=true',*args],text=True,encoding='utf-8').strip()
def exact(src,dst):
    dst=Path(dst);dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists():assert dst.read_bytes()==src.read_bytes(),dst
    else:shutil.copyfile(src,dst)
def ptr(value,path):
    for token in path.lstrip('/').split('/') if path else []:
        key=token.replace('~1','/').replace('~0','~');value=value[int(key)] if isinstance(value,list) else value[key]
    return value
def preflight():
    assert not (HERE/'preflight.json').exists()
    old=read(ROOT/'workflow/checks/integral-documentation-6e-20261103/preflight.json')['protected']
    paths={ROOT/rel for rel in old}
    for folder in ('integrals','workflow/checks/integral-documentation-low-loops-20261103','workflow/checks/integral-documentation-6e-20261103','workflow/checks/four-leg-3344-workflow-closed-20261104'):
        for p in (ROOT/folder).rglob('*'):
            if p.is_file() and 'runs' not in p.parts and '__pycache__' not in p.parts and not any(v.startswith('isolation-r') for v in p.parts):paths.add(p)
    for rel in read(ROOT/'workflow/checks/four-leg-3344-workflow-closed-20261104/preflight.json')['protected']:paths.add(ROOT/rel)
    for rel in read(SOURCE/'manifest.json')['files']:paths.add(SOURCE/rel)
    paths.add(SOURCE/'manifest.json')
    for rel in ('integrals/README.md','integrals/three-loop/README.md'):paths.discard(ROOT/rel)
    protected={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(paths) if p.is_file()}
    public={rel:sha(PUBLIC/rel) for rel in git('ls-files','-z',cwd=PUBLIC).split('\0') if rel}
    put(HERE/'preflight.json',{'utc':datetime.now(timezone.utc).isoformat(),'task':'WILSON-5.27','owner':'codex-main','HEAD':git('rev-parse','HEAD'),'branch':git('branch','--show-current'),'staged':git('diff','--cached','--name-only').splitlines(),'tracked_dirty':git('status','--short','--untracked-files=no').splitlines(),'public_HEAD':git('rev-parse','HEAD',cwd=PUBLIC),'protected':protected,'public_files':public,'allowed_paths':['15 new insertion entry directories','integrals/families/three-loop-insertions',HERE.relative_to(ROOT).as_posix(),'scoped root integral index/master/project-map/data-flow navigation','WILSON-5.27 task ledger'],'scope':'Starting bytes are preserved, not rebound as mathematical proof. Current four-leg work, earlier readers, original results and frozen/public packets are protected. No commit/push or scientific selector mutation.'})
    print({'protected':len(protected),'public_files':len(public)})
def stage():
    m=read(SOURCE/'manifest.json')
    for rel,row in m['files'].items():assert sha(SOURCE/rel)==row['sha256'],rel;exact(SOURCE/rel,RUNTIME/rel)
    exact(SOURCE/'manifest.json',RUNTIME/'manifest.json');entries=[]
    for row in read(SOURCE/'catalogue.json')['orbits']:
        oid=row['orbit_id'];old=SOURCE/'integrals/three-loop'/oid;own=ROOT/'integrals/three-loop'/oid;own.mkdir(parents=True,exist_ok=True)
        o=read(old/'input.json');result=read(old/'result.json');budget=read(old/'subintegrals.json');base='../../families/three-loop-insertions/runtime/'
        refs={k:{**o[k],'path':base+o[k]['path']} for k in ('graph','native_prepared','native_raw')}
        native=ptr(read(SOURCE/o['native_prepared']['path']),o['native_prepared']['field'])
        ordered=budget.get('raw_ordered_J',native.get('J_raw_s',native.get('J')))
        put(own/'input.json',{'orbit_id':oid,'family':row['family'],'references_relative_to':'this input.json directory','source_bindings':refs,'archived_input':{'path':base+f'integrals/three-loop/{oid}/input.json','sha256':sha(old/'input.json')},'original_ordered_J':ordered,'native_normalization':native.get('normalization'),'evaluated_quantity_normalization':result['normalization'],'native_operator':native.get('operator'),'quantity':result.get('quantity','normalized scalar I, external colour operator separate'),'relationship':o['relationship'],'scope':'Original ordered tensor kernel and source bindings are explicit. Scalar/geometric quantity and full outside colour operators retain distinct definitions; archived input retains full historical provenance.'})
        exact(old/'result.json',own/'result.json')
        for d in ('code','data','checks'):(own/d).mkdir(exist_ok=True)
        (own/'code/reproduce.py').write_text('''"""Generate the full15-entry batch; select this own member afterwards."""
from pathlib import Path
import runpy,sys
runtime=Path(__file__).resolve().parents[3]/'families/three-loop-insertions/runtime'
sys.path.insert(0,str(runtime/'scripts'))
runpy.run_path(str(runtime/'scripts/evaluate_batch.py'),run_name='__main__')
''',encoding='utf-8')
        put(own/'data/index.json',{'orbit_id':oid,'runtime_base':'../../../families/three-loop-insertions/runtime','archived_subintegrals':{'path':f'integrals/three-loop/{oid}/subintegrals.json','sha256':sha(old/'subintegrals.json')},'scope':'Exact kernel/order budget. Specific evaluated-object and operator/code roles are expanded after current execution.'})
        put(own/'checks/index.json',{'orbit_id':oid,'historical_checks':{'path':f'../../../families/three-loop-insertions/runtime/integrals/three-loop/{oid}/checks.json','sha256':sha(old/'checks.json')},'actual_replay':'../../../../workflow/checks/integral-documentation-insertions-20261104/replay.json','scientific_reading_audit':'../../../../workflow/checks/integral-documentation-insertions-20261104/review/README.md'})
        entries.append({'orbit_id':oid,'family':row['family'],'path':own.relative_to(ROOT).as_posix(),'result_sha256':sha(own/'result.json')})
    put(HERE/'entries.json',entries);put(HERE/'runtime-copy.json',{'files':len(m['files'])+1,'source_manifest_sha256':sha(SOURCE/'manifest.json'),'scope':'Declared original P5 packet copied byte-exactly; old Markdown retains historical status.'});print({'entries':len(entries),'runtime_files':len(m['files'])+1})
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('phase',choices=['preflight','stage']);a=p.parse_args();globals()[a.phase]()
