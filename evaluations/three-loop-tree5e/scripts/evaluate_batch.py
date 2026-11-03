"""Reassemble the 44 P6 tree-three-gluon sources from admitted raw primitives."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, importlib.metadata, json, os, sys, time
sys.dont_write_bytecode = True
PACK = Path(__file__).resolve().parents[1]

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def save(p, v): Path(p).write_bytes((json.dumps(v, ensure_ascii=False, indent=2)+'\n').encode('utf-8'))

def guard(out):
    allowed = [PACK/'scripts', PACK/'snapshots', PACK/'inputs.json', Path(sys.prefix), Path(sys.base_prefix), out]
    reads, denied = set(), []
    def audit(event, args):
        if event == 'open' and isinstance(args[0], (str, bytes, os.PathLike)):
            raw = os.fsdecode(args[0]); p = Path(raw).resolve()
            if raw.lower() in ('nul', '/dev/null'): return
            mode = args[1]; flags = args[2] if len(args)>2 and isinstance(args[2], int) else 0
            writing = (isinstance(mode, str) and any(c in mode for c in 'wax+')) or bool(flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC))
            ok = p.is_relative_to(out) if writing else any(p == b or p.is_relative_to(b) for b in allowed)
            if '/results/' in p.as_posix().lower() or p.name in ('coefficients.json', 'construction.json', 'representation.json', 'representations.json'):
                ok = False
            if not ok:
                denied.append({'event':event, 'path':str(p), 'writing':writing})
                raise PermissionError('Outside P6 generation closure: '+str(p))
            if not writing and p.is_file(): reads.add(str(p))
        if event in ('subprocess.Popen','os.system','socket.connect','socket.bind'):
            denied.append({'event':event}); raise PermissionError('External process/network forbidden')
    sys.addaudithook(audit)
    return reads, denied

def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True); a = p.parse_args(); out = a.output.resolve()
    if out.exists(): p.error('Choose a NEW output directory')
    if any(out.is_relative_to(PACK/n) for n in ('scripts','snapshots','references','integrals','evidence')): p.error('Output must be separate from inputs/code')
    source = read(PACK/'inputs.json')
    for rel, h in source['snapshot_files'].items():
        if sha(PACK/'snapshots'/rel) != h: raise ValueError('Changed snapshot '+rel)
    out.mkdir(parents=True); reads, denied = guard(out); start = time.perf_counter()
    import t123, t4
    rows = {}
    for module in (t123, t4):
        addition = module.generate()
        if set(addition) & set(rows): raise ValueError('Duplicate result')
        rows.update(addition); save(out/'generated.json', rows)
        print(module.__name__, len(addition), 'computed', flush=True)
    if set(rows) != {x['orbit_id'] for x in source['orbits']} or len(rows) != 44: raise ValueError('Incomplete44 batch')
    for oid, r in rows.items():
        if set(r['laurent']['coefficients']) != {str(n) for n in range(-6,1)}: raise ValueError('Incomplete seven-position result '+oid)
    scripts = {p.relative_to(PACK).as_posix():sha(p) for p in (PACK/'scripts').glob('*.py') if p.name.startswith(('t123','t4','evaluate_'))}
    report = {'format':'wilson-P6-generation/1','status':'computed','orbit_count':44,'I_positions':308,
        'seconds':time.perf_counter()-start,'created_utc':datetime.now(timezone.utc).isoformat(),'result_sha256':sha(out/'generated.json'),
        'code':scripts,'snapshot_manifest_sha256':sha(PACK/'inputs.json'),'read_files':sorted(reads),'denied_events':denied,
        'new_CAS_processes':0,'old_final_tables_read':False,'python':sys.version,
        'versions':{n:importlib.metadata.version(n) for n in ('sympy','mpmath')},
        'scope':'New exact reassembly using existing raw HyperInt primitives and integrated global-IBP blocks, Gamma jets and explicit order budgets. Not 44 new independent integrations.',
        'resource_scope':'Finite in-process batch; no new wall/memory cap. Python access audit is not an OS sandbox. Archived snapshots contain historical formulas; this is not a blinded-human experiment.',
        'formal':'not_run'}
    save(out/'run.json',report)
    if denied: raise ValueError('Unexpected denied generation event')
    print(json.dumps({k:report[k] for k in ('status','orbit_count','I_positions','seconds','result_sha256')}))
if __name__ == '__main__': main()
