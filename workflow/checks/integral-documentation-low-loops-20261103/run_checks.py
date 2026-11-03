"""Execute reader entrypoints and existing mathematical controls in fresh outputs."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse,subprocess,sys,time,json
from prepare import HERE,ROOT,RUNTIME,sha,read,put
def run(item):
    name,args,logroot=item;log=logroot/(name+'.log');started=time.perf_counter()
    with log.open('wb') as f:r=subprocess.run([sys.executable,'-X','utf8','-B',*map(str,args)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    record={'stage':name,'arguments':[str(a) for a in args],'returncode':r.returncode,'seconds':time.perf_counter()-started,'log_sha256':sha(log)}
    if r.returncode:raise RuntimeError('Reader check '+name+' failed; see '+str(log))
    return record
def main():
    p=argparse.ArgumentParser();p.add_argument('--output-root',required=True,type=Path);p.add_argument('--reader-root',type=Path,default=ROOT/'integrals');a=p.parse_args();out=a.output_root.resolve();reader=a.reader_root.resolve();runtime=reader/'families/low-loop-runtime'
    assert not out.exists(),'Choose a new output root';out.mkdir(parents=True)
    jobs=[('one-loop',[reader/'one-loop/2e-01/code/reproduce.py','--output-dir',out/'one-loop','--order','3','--dps','40']),
        ('two-loop',[reader/'two-loop/4e-01/code/reproduce.py','--output-dir',out/'two-loop']),
        ('two-loop-controls',[runtime/'integrals/two-loop/check_integrals.py','--output-dir',out/'two-loop-controls'])]
    jobs=[(name,args,out) for name,args in jobs]
    with ThreadPoolExecutor(max_workers=3) as pool:stages=list(pool.map(run,jobs))
    result={'status':'passed','stages':stages,'integrals':14,'scope':'Actual new reader entrypoints run existing source/mother/Gamma/HLOG replay and controls from byte-preserved runtime. This does not certify document continuity or reintegrate every historical CAS primitive.','formal':'not_run'}
    put(out/'replay.json',result);print(json.dumps({'status':'passed','stages':len(stages),'integrals':14}))
if __name__=='__main__':main()
