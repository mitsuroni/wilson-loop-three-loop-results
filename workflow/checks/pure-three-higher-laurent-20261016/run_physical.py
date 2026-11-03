"""Integrate fresh physical projections of exact exported source coefficients."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import subprocess
import sys
from physical_sources import build,ROOT,HERE,OLD

HYPER=Path('C:/Users/Lenovo/Documents/learning/work/hyperint')
MAPLE=Path('C:/Program Files/Maple 2025/bin.X86_64_WINDOWS/cmaple.exe')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--order',type=int,default=3)
    args=ap.parse_args();assert 0<=args.order<=3
    run=HERE/'physical-runs'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    run.mkdir(parents=True,exist_ok=False)
    print(run.relative_to(ROOT),flush=True)
    paths=[Path(__file__).resolve(),HERE/'physical_sources.py',OLD/'w_source_check.py',OLD/'w_lambda_gpl_check.py',
           HERE/'runs/r001/W-parent-words.json',HERE/'runs/r001/Lambda-sector-words.json',
           HYPER/'HyperInt.mpl',HYPER/'periodLookups.m']
    inputs={str(p):sha(p) for p in paths}
    build(run,args.order+3)
    lines=['restart:',f'_hyper_autoload_periods:=["{(HYPER/"periodLookups.m").as_posix()}"]:',
           f'read "{(HYPER/"HyperInt.mpl").as_posix()}":',
           '_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:',
           'FB:=proc(f) fibrationBasis(convert(f,HlogRegInf),[z]) end proc:',
           'Inv:=proc(f) local g;g:=fibrationBasis(convert(subs(z=1/zz,f),HlogRegInf),[zz]);subs(zz=z,g) end proc:',
           'Emit:=proc(tag,f) local fd;fd:=fopen(cat(tag,".txt"),WRITE):fprintf(fd,"%s = %a\\n",tag,f):fclose(fd):end proc:',
           'read "source-inputs.mpl":',
           'rl:=ln(z)-2*ln(1+z):']
    for n in range(3,args.order+4):
        for kind in ('W','L'):
            lines += [f'printf("START {kind}_num{n} time=%a\\n",time()):',
                      f'{kind}n{n}:=FB({kind}n{n}):',f'Emit("{kind}_num{n}",{kind}n{n}):',
                      f'printf("DONE {kind}_num{n} length=%d time=%a\\n",length({kind}n{n}),time()):gc():']
    for q in range(args.order+1):
        for kind in ('W','L'):
            form='+'.join(f'(3*rl/2)^{q+3-n}/{q+3-n}!*{kind}n{n}' for n in range(3,q+4))
            lines += [f'{kind}{q}:=FB(-({form})):']
            if kind=='W':
                lines += [f'if FB(W{q}-Inv(W{q}))<>0 then error "W crossing failed" fi:',f'Emit("W{q}",W{q}):']
            else:lines += [f'Lsum{q}:=FB(L{q}+Inv(L{q})):Emit("Lsum{q}",Lsum{q}):']
        hform='+'.join(f'(rl/2)^{q+1-k}/{q+1-k}!*(Hlog(1,[0,(-1/z)${k}])+Hlog(1,[0,(-z)${k}])-2*Hlog(1,[0,1${k}]))'
                       for k in range(1,q+2))
        lines += [f'H{q}:=FB(-({hform})):Emit("H{q}",H{q}):']
    # Gamma/Euler and crossing prefactor through epsilon cubed, no linear term.
    lines += ['p0:=1:p1:=0:p2:=3*zeta[2]/2-9*ln(z)^2/8:p3:=zeta[3]:']
    for q in range(args.order+1):
        cube='+'.join(f'H{i}*H{j}*H{q-i-j}' for i in range(q+1) for j in range(q-i+1))
        lines += [f'S{q}:=FB({cube}):',f'A{q}:=FB(S{q}/2-W{q}-Lsum{q}):',f'B{q}:=FB(W{q}+Lsum{q}/2-S{q}/3):']
        ia='+'.join(f'p{r}*A{q-r}' for r in range(q+1))
        ib='+'.join(f'p{r}*B{q-r}' for r in range(q+1))
        lines += [f'I20_{q}:=FB(-({ia})/16):Emit("I20_{q}",I20_{q}):',
                  f'I21_{q}:=FB(-({ib})/32):Emit("I21_{q}",I21_{q}):',
                  f'Emit("A{q}",A{q}):Emit("B{q}",B{q}):',
                  f'printf("PHYSICAL q={q} time=%a\\n",time()):']
    lines += ['printf("PHYSICAL_SOURCE_COMPLETE time=%a\\n",time()):','quit:']
    source=run/'compute.mpl';source.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    command=[str(MAPLE),'-q',str(source)]
    (run/'provenance.json').write_text(json.dumps({'command':command,'inputs_sha256':inputs,
        'source_sha256':sha(source),'formal':'not_run'},indent=2)+'\n',encoding='utf-8')
    with (run/'maple.log').open('wb') as log:
        try:code=subprocess.run(command,cwd=run,stdout=log,stderr=subprocess.STDOUT,timeout=2400).returncode
        except subprocess.TimeoutExpired:code='timeout'
    log=(run/'maple.log').read_text(encoding='utf-8',errors='replace')
    passed=code==0 and 'PHYSICAL_SOURCE_COMPLETE' in log and 'Error,' not in log
    assert all(sha(Path(p))==expected for p,expected in inputs.items())
    status={'passed':passed,'returncode':code,'order':args.order,'formal':'not_run',
            'outputs_sha256':{p.name:sha(p) for p in run.glob('*.txt')},
            'inputs_preserved':True,'source':'Exact exported generic GPL words, complete physical projection; no old target coefficients used.'}
    (run/'status.json').write_text(json.dumps(status,indent=2)+'\n',encoding='utf-8')
    print(log[-3500:],flush=True)
    print(json.dumps({'passed':passed,'run':run.relative_to(ROOT).as_posix()}),flush=True)
    if not passed:raise SystemExit(1)

if __name__=='__main__':main()
