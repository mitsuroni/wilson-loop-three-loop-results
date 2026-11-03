"""Assemble native D-frame weighted kernels, then the unique physical I."""
from pathlib import Path
import sys,json
from check_math import ROOT,HERE,read,dump
ids=[int(v) for v in sys.argv[1:]] or [13,14,15,40]
runs={k:ROOT/read(HERE/('current-'+k+'.json'))['run'] for k in ['easy','hard','recoil','derivatives']}
runs['tw']=ROOT/read(HERE/'current-tw_vfirst.json')['run']
lines=['restart:','_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:', 'read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":','_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:',
 'FB:=proc(f) if type(f,list) then fibrationBasis(f,[z]) else fibrationBasis(convert(subs(zeta[4]=2/5*zeta[2]^2,f),HlogRegInf),[z]) fi end proc:',
 'Inv:=proc(f) local a;a:=fibrationBasis(convert(subs(z=1/zz,f),HlogRegInf),[zz]);subs(zz=z,a) end proc:',
 'Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\\n",n,f);fclose(fd) end proc:']
used=[]
def locate(kind,name):
 if kind=='tw':
  for tag in ['tw_words','tw_14','tw_15','tw_vfirst']:
   pointer=HERE/('current-'+tag+'.json')
   if pointer.exists():
    p=ROOT/read(pointer)['run']/(name+'.txt')
    if p.exists():return p
 return runs[kind]/(name+'.txt')
def expr(kind,name):
 p=locate(kind,name);used.append(p.relative_to(ROOT).as_posix());return '('+p.read_text(encoding='utf-8').split('=',1)[1].strip()+')'
available={i:(2 if i==40 else max(q for q in range(-2,3) if locate('tw',f'TW{i}q{q}').exists())) for i in ids}

def ser(kind,prefix,lo,hi):return '+'.join(f'{expr(kind,prefix+str(q))}*eps^({q})' for q in range(lo,hi+1))
for typ in ['F0','F1','G1']:lines.append(f'{typ}:={ser("easy",typ+"q",0,4)}:')
lines+=['br:=convert(series(exp(add((-1)^n*zeta[n]*(1+2^n-3^n)*eps^n/n,n=2..4)),eps=0,5),polynom):',
 'zp:=convert(series(exp(2*eps*ln(z)),eps=0,5),polynom):',
 'pref:=convert(series(exp(add(zeta[n]*(1+2^n)*eps^n/n,n=2..4))/cosh(3*eps*ln(z)/2),eps=0,5),polynom):']
for oid in ids:
 hy=ser('hard',f'H{oid}q',0,2)
 y0=ser('easy',f'Y0i{oid}q',0,3) if oid in (13,40) else '0'
 y1=ser('easy',f'Y1i{oid}q',0,3) if oid!=13 else '0'
 r=ser('recoil',f'R{oid}q',-1,2)
 cuspX={13:'zp*G1/(6*z*eps^2)',14:'zp*(F0-G1)/(6*z*eps^2)',15:'zp*(F0-G1)/(3*z*eps^2)',40:'0'}[oid]
 cuspW={13:'F1/(2*eps^2)',14:'(F0-F1)/(2*eps^2)',15:'0',40:'F0*(1-br)/(2*eps^2)'}[oid]
 hardW='zp/z*Inv('+ser('hard',f'H{100+oid}q',0,2)+')' if oid in (15,40) else '0'
 if oid==40:deriv='z*eps*('+ser('derivatives','D5q',-1,1)+')'
 else:
  dy1=ser('derivatives','D1q',-2,1);dy2=ser('derivatives','D2q',-2,1)
  if oid==13:deriv=f'eps*(1+z)*(({dy1})-({dy2}))'
  elif oid==14:deriv=f'-eps*(1+z)*(({dy1})-({dy2}))'
  else:deriv=f'-eps*(1+z)*({dy2})'
  deriv+='+(1-2*eps)*('+ser('tw',f'TW{oid}q',-2,available[oid])+')'
  cuspW='0';hardW='0'

 lines.append(f'Qraw:=expand(z*({cuspX})-z*({hy})+(1+z)*(({y1})-({y0}))/eps-({hardW})+({cuspW})+({deriv})+4*(1-2*eps)*({r})):')
 lines+=['Q:=add(FB(coeff(Qraw,eps,q))*eps^q,q=-2..2):',f'for q from -2 to {available[oid]} do Emit(cat("Q{oid}q",q),coeff(Q,eps,q)) od:',
 'evenQ:=expand((exp(-eps*ln(z)/2)*Q+exp(eps*ln(z)/2)*Inv(Q))/2):',
 'evenQ:=convert(series(evenQ,eps=0,3),polynom):',
 'out:=expand(pref*evenQ/(32*eps^2)):',
 f'for q from -6 to {available[oid]-2} do ans:=FB(coeff(out,eps,q)):Emit(cat("I{oid}q",q),ans):printf("I {oid} q%d len%d time%a\\n",q,length(ans),time()):if FB(Inv(ans)-ans)<>0 then error "I crossing failure" fi od:']
lines+=['printf("ASSEMBLY_COMPLETE\\n"):quit:']
name='assemble_'+('partial_' if min(available.values())<2 else '')+('_'.join(map(str,ids)))+'.mpl'
(HERE/name).write_text('\n'.join(lines)+'\n',encoding='utf-8')
dump(HERE/(name+'.inputs.json'),{'input_files':sorted(set(used)),'ids':ids,'highest_I_coefficient':{str(i):q-2 for i,q in available.items()}})
print(name)