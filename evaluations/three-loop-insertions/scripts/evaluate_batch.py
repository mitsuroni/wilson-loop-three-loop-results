"""Evaluate P5 from frozen original kernels and admitted function dictionaries."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,importlib.metadata,json,os,sys,time
sys.dont_write_bytecode=True
PACK=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def save(p,d):Path(p).write_bytes((json.dumps(d,ensure_ascii=False,indent=2)+'\n').encode('utf-8'))
def guard(out):
    allowed=[PACK/'scripts',PACK/'snapshots',PACK/'inputs.json',Path(sys.prefix),Path(sys.base_prefix),out]
    reads=set();denied=[]
    def audit(event,args):
        if event=='open' and isinstance(args[0],(str,bytes,os.PathLike)):
            raw=os.fsdecode(args[0]);p=Path(raw).resolve()
            if raw.lower() in ('nul','/dev/null'):return
            mode=args[1];flags=args[2] if len(args)>2 and isinstance(args[2],int) else 0
            writing=isinstance(mode,str) and any(c in mode for c in 'wax+') or bool(flags&(os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC))
            permitted=p.is_relative_to(out) if writing else any(p==b or p.is_relative_to(b) for b in allowed)
            if '/results/' in p.as_posix().lower() or p.name=='operator-expansion-r001.json':permitted=False
            if not permitted:
                denied.append({'event':event,'path':str(p),'writing':writing});raise PermissionError('Outside P5 generation closure: '+str(p))
            if not writing and p.exists():reads.add(str(p))
        if event in ['subprocess.Popen','os.system','socket.connect','socket.bind']:
            denied.append({'event':event});raise PermissionError('External process/network forbidden')
    sys.addaudithook(audit);return reads,denied
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve()
    if out.exists():p.error('Choose a NEW output directory')
    if any(out.is_relative_to(PACK/n) for n in ['scripts','snapshots','references','integrals','evidence']):p.error('Output must be separate from inputs/code')
    source=read(PACK/'inputs.json')
    for rel,h in source['snapshot_files'].items():
        if sha(PACK/'snapshots'/rel)!=h:raise ValueError('Changed input/code '+rel)
    out.mkdir(parents=True)
    reads,denied=guard(out);start=time.perf_counter()
    import pi1,chains_pi2
    rows={}
    for module in [pi1,chains_pi2]:
        addition=module.generate()
        if set(rows)&set(addition):raise ValueError('Duplicate source result')
        rows.update(addition)
        save(out/'generated.json',rows)
        print(module.__name__,len(addition),'computed',flush=True)
    wanted={o['orbit_id'] for o in source['orbits']}
    if set(rows)!=wanted or len(rows)!=15:raise ValueError('Incomplete15 batch')
    for oid,r in rows.items():
        if set(r['laurent']['coefficients'])!={str(n) for n in range(-6,1)}:raise ValueError('Incomplete seven-position result '+oid)
    scripts={p.relative_to(PACK).as_posix():sha(p) for p in (PACK/'scripts').glob('*.py') if p.name not in ['independent_controls.py','verify_batch.py','verify_manifest.py']}
    report={'format':'wilson-p5-generation/1','status':'computed','orbit_count':15,'I_positions':105,
        'seconds':time.perf_counter()-start,'result_sha256':sha(out/'generated.json'),'created_utc':datetime.now(timezone.utc).isoformat(),
        'code':scripts,'snapshot_manifest_sha256':sha(PACK/'inputs.json'),'read_files':sorted(reads),'denied_events':denied,
        'new_CAS_processes':0,'old_final_tables_read':False,'python':sys.version,
        'versions':{n:importlib.metadata.version(n) for n in ['sympy','mpmath']},
        'scope':'New original-root reductions/Gamma/order assembly using existing evaluated M/HLOG dictionaries; Pi2 full-D literature operator adopted. No second independent CAS integration or formal proof.',
        'resource_scope':'Fixed finite15-input in-process computation; no newly enforced memory or wall bound. Python audit is not an OS sandbox.',
        'formal':'not_run'}
    save(out/'run.json',report)
    if denied:raise ValueError('Unexpected generation access denial')
    print(json.dumps({k:report[k] for k in ['status','orbit_count','I_positions','seconds','result_sha256']}))
if __name__=='__main__':main()
