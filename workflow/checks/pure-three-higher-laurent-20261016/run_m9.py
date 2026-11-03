"""Compute every weight-nine M moment from its natural one-dimensional source."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib
import json
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
HYPER=Path('C:/Users/Lenovo/Documents/learning/work/hyperint')
MAPLE=Path('C:/Program Files/Maple 2025/bin.X86_64_WINDOWS/cmaple.exe')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    run=HERE/'candidate-runs'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    run.mkdir(parents=True,exist_ok=False)
    print(run.relative_to(ROOT),flush=True)
    lines=['restart:',f'_hyper_autoload_periods:=["{(HYPER/"periodLookups.m").as_posix()}"]:',
           f'read "{(HYPER/"HyperInt.mpl").as_posix()}":',
           '_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:',
           'FB:=proc(f) fibrationBasis(convert(f,HlogRegInf),[z]) end proc:',
           'Emit:=proc(tag,f) local fd;fd:=fopen(cat(tag,".txt"),WRITE):fprintf(fd,"%s = %a\\n",tag,f):fclose(fd):end proc:',
           'Wlog:=ln(x)+ln(1-x):Llog:=ln(1-a)-ln(a)+ln(x)-ln(1-x):']
    for j in range(1,9):
        k=8-j
        lines += [f'printf("START M_{k}_{j} time=%a\\n",time()):',
                  f'ans:=hyperInt(convert(Wlog^{k}*Llog^{j}/(a-x),HlogRegInf),[x=0..1]):',
                  'result:=FB(subs(a=1/(1+z),fibrationBasis(ans,[a]))):',
                  f'Emit("M_{k}_{j}",result):',f'printf("DONE M_{k}_{j} time=%a\\n",time()):gc():']
    lines += ['printf("M9_COMPLETE time=%a\\n",time()):','quit:']
    src=run/'compute.mpl';src.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    command=[str(MAPLE),'-q',str(src)]
    dependencies=[Path(__file__).resolve(),ROOT/'workflow/checks/4e04-high-order/compute.mpl',HYPER/'HyperInt.mpl',HYPER/'periodLookups.m']
    (run/'provenance.json').write_text(json.dumps({'command':command,'source_sha256':sha(src),
        'inputs_sha256':{str(p):sha(p) for p in dependencies},'scope':'All M_k,j with k+j=8, independently of new three-kernel targets.'},indent=2)+'\n',encoding='utf-8')
    with (run/'maple.log').open('wb') as log:
        try:code=subprocess.run(command,cwd=run,stdout=log,stderr=subprocess.STDOUT,timeout=1200).returncode
        except subprocess.TimeoutExpired:code='timeout'
    log=(run/'maple.log').read_text(encoding='utf-8',errors='replace')
    passed=code==0 and 'M9_COMPLETE' in log and 'Error,' not in log
    (run/'status.json').write_text(json.dumps({'passed':passed,'returncode':code,'formal':'not_run',
        'outputs_sha256':{p.name:sha(p) for p in run.glob('*.txt')}},indent=2)+'\n',encoding='utf-8')
    print(log[-2000:],flush=True)
    print(json.dumps({'passed':passed,'run':run.relative_to(ROOT).as_posix()}))
    if not passed:raise SystemExit(1)

if __name__=='__main__':main()
