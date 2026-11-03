"""Preserve original science and prepare a current low-loop reading layer."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,subprocess,shutil
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];PUBLIC=ROOT/'.wilson/publications/wilson-loop-three-loop-results'
READING=ROOT/'integrals';RUNTIME=READING/'families/low-loop-runtime'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def put(p,v):p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())
def exact(src,dst):
    dst=Path(dst);dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists():assert dst.read_bytes()==src.read_bytes(),str(dst)
    else:dst.write_bytes(src.read_bytes())
def git(*args,cwd=ROOT):return subprocess.check_output(['git','-C',str(cwd),'-c','safe.directory='+cwd.as_posix(),'-c','core.longpaths=true',*args],text=True,encoding='utf-8').strip()
def preflight():
    assert not (HERE/'preflight.json').exists()
    p7=ROOT/'workflow/checks/three-loop-tgvlegs-p7-20261103';protected=dict(read(p7/'preflight.json')['protected'])
    for rel,row in read(p7/'package/manifest.json')['files'].items():protected[(p7/'package'/rel).relative_to(ROOT).as_posix()]=row['sha256']
    protected[(p7/'package/manifest.json').relative_to(ROOT).as_posix()]=sha(p7/'package/manifest.json')
    for folder in ('one-loop','two-loop'):
        m=read(PUBLIC/'integrals'/folder/'source-manifest.json')
        def visit(v):
            if isinstance(v,dict):
                path=v.get('original_path')
                if path and (ROOT/path).is_file():protected[path]=sha(ROOT/path)
                for x in v.values():visit(x)
            elif isinstance(v,list):
                for x in v:visit(x)
        visit(m)
    for p,h in protected.items():assert sha(ROOT/p)==h,('Original drift',p)
    public={p:sha(PUBLIC/p) for p in git('ls-files','-z',cwd=PUBLIC).split('\0') if p}
    put(HERE/'preflight.json',{'utc':datetime.now(timezone.utc).isoformat(),'owner':'codex-main','task':'WILSON-5.24','HEAD':git('rev-parse','HEAD'),'branch':git('branch','--show-current'),
        'staged':git('diff','--cached','--name-only').splitlines(),'tracked_dirty':git('status','--short','--untracked-files=no').splitlines(),'protected':protected,'public_files':public,'public_HEAD':git('rev-parse','HEAD',cwd=PUBLIC),
        'allowed_paths':['workflow/skills/wilson-integral-derivations','workflow/integral-documentation-standard.md','integrals',HERE.relative_to(ROOT).as_posix(),'scoped current navigation/task ledger','installed user skill identical to reviewed project source'],
        'policy':'New reader/standard/skill; preserve original science/results/proofs/frozen packets/public. No Git/push or coldCAS. Document review is separate from numerical/regression validation.'})
    print(json.dumps({'protected':len(protected),'public_files':len(public)}))
def stage():
    # Keep the original public relative context for unmodified replay programs.
    count=0
    for name in ('one-loop','two-loop'):
        src=PUBLIC/'integrals'/name;dst=RUNTIME/'integrals'/name
        for p in src.rglob('*'):
            if p.is_file():exact(p,dst/p.relative_to(src));count+=1
    upstream={}
    for oid in ('2e-01','2e-02'):
        obj=read(PUBLIC/f'integrals/one-loop/{oid}/input.json')
        def visit(v):
            if isinstance(v,dict):
                if 'path' in v and 'sha256' in v and isinstance(v['path'],str):
                    path=v['path']
                    if path.startswith('data/'):
                        src=PUBLIC/path;assert sha(src)==v['sha256'];exact(src,RUNTIME/path);upstream[path]=sha(src)
                for x in v.values():visit(x)
            elif isinstance(v,list):
                for x in v:visit(x)
        visit(obj)
    files={p.relative_to(RUNTIME).as_posix():sha(p) for p in RUNTIME.rglob('*') if p.is_file()}
    put(RUNTIME/'source-manifest.json',{'format':'wilson-lowloop-reader-runtime/1','source_public_commit':git('rev-parse','HEAD',cwd=PUBLIC),'files':files,'upstream_data':upstream,
        'scope':'Unmodified historical/public code+data context for actual replay. Archived MD are not current reader articles. No claim of new cold CAS integration.'})
    entries=[]
    for folder,ids in [('one-loop',['2e-01','2e-02']),('two-loop',[f'4e-{n:02}' for n in range(1,7)]+[f'3e-3g-{n:02}' for n in range(1,5)]+['2e-Pi1-01','2e-Pi1-02'])]:
        for oid in ids:
            old=PUBLIC/'integrals'/folder/oid;own=READING/folder/oid;own.mkdir(parents=True,exist_ok=True)
            exact(old/'result.json',own/'result.json')
            if (old/'input.json').exists():exact(old/'input.json',own/'input.json')
            else:
                cat=read(PUBLIC/'integrals/two-loop/catalog.json');rows=cat.get('integrals',cat.get('orbits',cat))
                bind=next(x for x in rows if x.get('orbit_id',x.get('id'))==oid) if isinstance(rows,list) else rows.get(oid,{})
                put(own/'input.json',{'orbit_id':oid,'source_public_commit':git('rev-parse','HEAD',cwd=PUBLIC),'accepted_result':{'path':f'../../families/low-loop-runtime/integrals/{folder}/{oid}/result.json','sha256':sha(old/'result.json')},
                    'family_catalog_entry':bind,'source_manifest':'../../families/low-loop-runtime/integrals/two-loop/source-manifest.json','scope':'Accepted current source and original ordered domain; retain source-manifest exact native/ordered refs.'})
            for sub in ('data','checks','code'):(own/sub).mkdir(exist_ok=True)
            put(own/'data/index.json',{'orbit_id':oid,'runtime_base':'../../../families/low-loop-runtime','family_source_manifest':f'integrals/{folder}/source-manifest.json',
                'accepted_result':{'path':f'integrals/{folder}/{oid}/result.json','sha256':sha(old/'result.json')},'scope':'Source/primitives and exact result kept distinct. Old whole-I data are comparison targets, not original integration code.'})
            put(own/'checks/index.json',{'orbit_id':oid,'historical_context':f'../../../families/low-loop-runtime/integrals/{folder}',
                'current_run_report':'../../../../workflow/checks/integral-documentation-low-loops-20261103/replay.json','scope':'Historical and new replay checks separate from current scientific-document review.'})
            code='''"""Run the shared same-loop calculation; read this integral's member."""
from pathlib import Path
import runpy,sys
root=Path(__file__).resolve().parents[3]
runtime=root/'families/low-loop-runtime/integrals/LOOP'
sys.path.insert(0,str(runtime))
runpy.run_path(str(runtime/'reproduce.py'),run_name='__main__')
'''.replace('LOOP',folder)
            (own/'code/reproduce.py').write_text(code,encoding='utf-8')
            entries.append({'orbit_id':oid,'loop_directory':folder,'path':own.relative_to(ROOT).as_posix(),'result_sha256':sha(own/'result.json')})
    put(HERE/'entries.json',entries);print(json.dumps({'entries':len(entries),'runtime_files':len(files),'extra_upstream_files':len(upstream),'copied_public_family_files':count}))
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('phase',choices=['preflight','stage']);a=p.parse_args();globals()[a.phase]()
