"""Execute supported39-entry replay and keep generation separate from checks."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import argparse,subprocess,sys,time
from prepare import HERE,ROOT,RUNTIME,sha,read,put

def run(label,args,logroot):
    started=time.perf_counter();log=logroot/(label+'.log')
    with log.open('wb') as f:r=subprocess.run([sys.executable,'-X','utf8','-B',*map(str,args)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    record={'stage':label,'arguments':[str(a) for a in args],'returncode':r.returncode,'elapsed_seconds':time.perf_counter()-started,'log_sha256':sha(log)}
    put(logroot/(label+'-command.json'),record)
    assert r.returncode==0,(label,log)
    return record

def main():
    p=argparse.ArgumentParser();p.add_argument('--output-root',required=True,type=Path);p.add_argument('--reader-root',default=ROOT/'integrals',type=Path);p.add_argument('--skip-controls',action='store_true');a=p.parse_args()
    reader=a.reader_root.resolve();runtime=reader/'families/three-loop-6e/runtime';out=a.output_root.resolve()
    assert not out.exists(),'Choose a new output directory';out.mkdir(parents=True)
    generation=run('generate',[reader/'three-loop/6e-01/code/reproduce.py','--output',out/'generation'],out)
    frozen={'utc':datetime.now(timezone.utc).isoformat(),'generated_sha256':sha(out/'generation/generated.json'),'run_sha256':sha(out/'generation/run.json'),'stage':'Actual generation finished; reference comparison is the following operation. No human-blindness or trusted timestamp claim.'}
    put(out/'generation-freeze.json',frozen)
    jobs=[('verify',[runtime/'scripts/verify_batch.py','--results',out/'generation/generated.json','--out',out/'verification.json'])]
    if not a.skip_controls:jobs.append(('controls',[runtime/'scripts/independent_controls.py','--input',runtime/'engine/benchmarks/2ne/prepared/6e.json','--results',out/'generation/generated.json','--out',out/'independent-controls.json']))
    with ThreadPoolExecutor(max_workers=2) as pool:checks=list(pool.map(lambda item:run(*item,out),jobs))
    assert sha(out/'generation/generated.json')==frozen['generated_sha256']
    verification=read(out/'verification.json');result={'status':'passed','integrals':verification['orbit_count'],'coefficients':verification['coefficient_positions'],'function_back_substitutions':verification['function_back_substitutions'],'generation':generation,'checks':checks,'controls_executed':not a.skip_controls,'scope':'Fresh original-root recognition, Gamma expansion or byte-bound historical primitive replay, physical assembly, and subsequent verification. No new external CAS or formal integral proof.'}
    put(out/'replay.json',result);print({k:result[k] for k in ('status','integrals','coefficients','function_back_substitutions','controls_executed')})
if __name__=='__main__':main()
