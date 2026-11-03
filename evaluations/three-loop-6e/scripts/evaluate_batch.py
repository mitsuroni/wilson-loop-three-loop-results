"""Rebuild all 39 original ordered 6e roots using byte-bound raw primitives.

No old final I is an input. No network, subprocess, installed Maple or library
file is used. The unchanged solver is given process-local adapters for archived
CAS primitives, historical path strings, direct Gamma expansion and output I/O.
"""
from __future__ import annotations
import argparse
from contextlib import ExitStack
from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import time
from unittest.mock import patch

sys.dont_write_bytecode = True
PACKAGE = Path(__file__).resolve().parents[1]
ENGINE = PACKAGE/'engine'

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def save(p, value): Path(p).write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n',encoding='utf-8')

def check_snapshot():
    snapshot=load(ENGINE/'snapshot.json')
    for row in snapshot['source_files']:
        p=(ENGINE/row['path']).resolve()
        if not p.is_relative_to(ENGINE.resolve()) or sha(p)!=row['sha256']:
            raise ValueError('Changed engine or ordered input: '+row['path'])
    profile=load(PACKAGE/'primitives/archive-profile.json')
    for row in profile['files']:
        p=(PACKAGE/'primitives'/row['path']).resolve()
        if not p.is_relative_to((PACKAGE/'primitives').resolve()) or sha(p)!=row['sha256']:
            raise ValueError('Changed archived primitive: '+row['path'])
    return profile

def install_guard(output):
    allowed_read=[ENGINE.resolve(),(PACKAGE/'primitives').resolve(),(PACKAGE/'scripts').resolve(),
                  output.resolve(),Path(sys.base_prefix).resolve(),Path(sys.prefix).resolve()]
    reads=set(); denied=[]
    def guard(event,args):
        if event=='open' and isinstance(args[0],(str,bytes,os.PathLike)):
            raw=os.fsdecode(args[0]); p=Path(raw).resolve()
            if raw.lower() in ('nul','/dev/null'):return
            mode=args[1]; flags=args[2] if len(args)>2 and isinstance(args[2],int) else 0
            writing=(isinstance(mode,str) and any(c in mode for c in 'wax+')) or bool(flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC))
            allowed=p.is_relative_to(output) if writing else any(p.is_relative_to(base) for base in allowed_read)
            if not allowed:
                denied.append({'event':event,'path':str(p),'writing':writing})
                raise PermissionError('Batch file access outside declared closure: '+str(p))
            if not writing and p.exists():reads.add(str(p))
        if event in ('subprocess.Popen','os.system','socket.connect','socket.bind'):
            denied.append({'event':event});raise PermissionError('No subprocess or network in raw replay')
    sys.addaudithook(guard)
    return reads,denied

def generate(output):
    profile=check_snapshot()
    sys.path.insert(0,str(ENGINE))
    reads,denied=install_guard(output)
    import sympy
    from wilson_calc import solve, backends, ordered_backend, single, ordered, solver
    from wilson_calc.runtime import CalculationObjects, engine_fingerprint
    from wilson_calc.assembly import expand_scalar
    from wilson_calc.algebra import ast_to_sympy
    from wilson_calc.exchange import EPS,Z1
    from primitive_replay import PrimitiveReplay,frozen_config,historical_maple_path

    class OutputStore(CalculationObjects):
        """Reuse immutable byte operations, with the CLI's output-only boundary.

        The original constructor requires output below its source tree. This
        explicit I/O adapter accepts the checked new output and opens no DB.
        """
        def __init__(self,path):
            self.path=Path(path).resolve()
            if not self.path.is_relative_to(output):raise ValueError('Output store escaped its output')
            (self.path/'objects').mkdir(parents=True)

    scalar_runs=[]
    def direct_scalar(payload,seconds,memory_mb):
        coefficients,lower=expand_scalar(ast_to_sympy(payload['J'],{'eps':EPS,'z1':Z1}),payload['normalization'],
                    payload['quantity'],payload['through'],EPS,Z1,payload['loop_order'])
        scalar_runs.append({'mode':'fresh_original_Gamma_expansion_in_process','seconds_argument':seconds,
                            'memory_mb_argument':memory_mb,'resource_enforcement':'none; finite fixed task'})
        return coefficients,lower
    def forbidden_process(*args,**kwargs):raise PermissionError('A new external process is forbidden')

    started=time.perf_counter(); rows={}; events_by_orbit={}
    with ExitStack() as stack:
        # Only filesystem/config/process providers are replaced. Mathematical
        # recognition, finite-order budgets, parsers and assembly are unchanged.
        for module in (backends,single,ordered):stack.enter_context(patch.object(module,'find_maple_config',frozen_config))
        for module in (backends,ordered_backend):
            stack.enter_context(patch.object(module,'_maple_path',historical_maple_path))
            stack.enter_context(patch.object(module,'run_bounded_process',forbidden_process))
        stack.enter_context(patch.object(solver,'scalar_worker',direct_scalar))
        replay=PrimitiveReplay(PACKAGE/'primitives',profile['approved_resource_modules'])
        for module,name,method in ((single,'compute_single_opposite_moments',replay.single),
                                  (ordered_backend,'compute_ordered_constant',replay.constant),
                                  (ordered,'compute_ordered_moments',replay.ordered)):
            stack.enter_context(patch.object(module,name,method))
        store=OutputStore(output/'objects')
        manifest=ENGINE/'benchmarks/2ne/prepared/manifest.json'
        for i in range(1,40):
            oid=f'6e-{i:02d}'; before=len(replay.events)
            result=solve({'format':'wilson-calculation-request/0.1',
                          'source':{'kind':'prepared','integral_id':oid,'manifest':str(manifest)},
                          'quantity':'I','through':0,'cache_policy':'recompute',
                          'limits':{'seconds':120,'cas_seconds':120,'max_order':8,'max_terms':100000,'memory_mb':2048}},store=store)
            result['provenance']['batch_execution_mode']='new_assembly_with_archived_raw_primitives_or_fresh_Gamma'
            rows[oid]=result; events_by_orbit[oid]=list(range(before,len(replay.events)))
            save(output/'generated.json',rows)
            print(oid,result['status'],result['reason']['code'],flush=True)
        report=replay.report()
        report.update(events_by_orbit=events_by_orbit,new_scalar_expansions=scalar_runs,
                      scalar_expansion_count=len(scalar_runs),elapsed_seconds=time.perf_counter()-started,
                      engine_fingerprint=engine_fingerprint(),
                      interpretation='Original raw CAS integrations reused, current request/script/parser/physical assembly reexecuted. Not independent CAS integration.')
        save(output/'primitive-replay.json',report)
    local={name:str(Path(module.__file__).resolve()) for name,module in sys.modules.items()
           if name.split('.')[0] in ('wilson_calc','preprocessing','interfaces','proof_store') and getattr(module,'__file__',None)}
    if any(not Path(p).is_relative_to(ENGINE.resolve()) for p in local.values()):raise ValueError('Project import escaped copied engine')
    result={'format':'wilson-6e-batch-generation/1','status_counts':dict(Counter(r['status'] for r in rows.values())),
            'orbit_count':len(rows),'elapsed_seconds':time.perf_counter()-started,
            'generated_sha256':sha(output/'generated.json'),'primitive_replay_sha256':sha(output/'primitive-replay.json'),
            'new_cas_processes':0,'old_final_I_read':False,'formal_state':'not_run',
            'python':sys.version,'versions':{p:importlib.metadata.version(p) for p in ('sympy','mpmath','jsonschema')},
            'project_modules':local,'read_files':sorted(reads),'denied_events':denied,
            'guard_scope':'Python open/subprocess/socket audit; no operating-system sandbox claim.',
            'source_checks':'All engine/input and archived primitive files hash checked before generation.',
            'resource_scope':'Finite existing 39-root batch; Gamma runs in process, no newly enforced memory/time bound.'}
    save(output/'run.json',result)
    if len(rows)!=39 or any(r['status']!='computed' for r in rows.values()):raise ValueError('Incomplete39 batch; inspect saved generated/run evidence')
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);args=p.parse_args()
    output=args.output.resolve()
    if output.exists():p.error('Choose a NEW output directory')
    if any(output.is_relative_to((PACKAGE/name).resolve()) for name in ('engine','primitives','scripts','inputs','references','integrals')):
        p.error('Output must be separate from source/input/reference/article directories')
    output.mkdir(parents=True)
    result=generate(output)
    print(json.dumps({k:result[k] for k in ('orbit_count','status_counts','elapsed_seconds','generated_sha256')}))

if __name__=='__main__':main()
