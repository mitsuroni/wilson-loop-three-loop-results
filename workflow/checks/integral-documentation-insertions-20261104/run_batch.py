"""Run the supported15-entry generation, then independent checking stages."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import argparse,subprocess,sys,time
from prepare import HERE,ROOT,RUNTIME,sha,read,put
def run(label,args,out):
    start=time.perf_counter();log=out/(label+'.log')
    with log.open('wb') as f:r=subprocess.run([sys.executable,'-X','utf8','-B',*map(str,args)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    record={'stage':label,'arguments':list(map(str,args)),'returncode':r.returncode,'seconds':time.perf_counter()-start,'log_sha256':sha(log)}
    put(out/(label+'-command.json'),record);assert r.returncode==0,(label,log);return record
def main():
    p=argparse.ArgumentParser();p.add_argument('--output-root',required=True,type=Path);a=p.parse_args();out=a.output_root.resolve();assert not out.exists();out.mkdir(parents=True)
    gen=run('generate',[ROOT/'integrals/three-loop/4e-Pi1-01/code/reproduce.py','--output',out/'generation'],out)
    frozen=sha(out/'generation/generated.json');put(out/'generation-freeze.json',{'utc':datetime.now(timezone.utc).isoformat(),'generated_sha256':frozen,'scope':'Generation completed before the following comparisons. No trusted timestamp or human-blindness claim.'})
    jobs=[('verify','verify_batch.py','verification.json'),('representations','check_representations.py','representations.json'),('controls','independent_controls.py','independent-controls.json')]
    with ThreadPoolExecutor(max_workers=3) as pool:checks=list(pool.map(lambda job:run(job[0],[RUNTIME/'scripts'/job[1],'--results',out/'generation/generated.json','--out',out/job[2]],out),jobs))
    assert sha(out/'generation/generated.json')==frozen
    put(out/'replay.json',{'status':'passed','generation':gen,'checks':checks,'scope':'Source-matched root reductions, adopted moment dictionaries and full-D Pi2 operator, fresh physical/colour assembly, then exact comparison, representations and scoped controls. No new CAS or formal proof.'})
    print({'status':'passed','entries':15,'scalar_or_geometric_positions':105,'scalar_and_color_representations':343})
if __name__=='__main__':main()
