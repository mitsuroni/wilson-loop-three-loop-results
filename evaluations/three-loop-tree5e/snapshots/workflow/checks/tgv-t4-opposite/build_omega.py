from pathlib import Path
import json
h=Path(__file__).resolve().parent;root=h.parents[2]
r=root/json.loads((h/'current-outer28.json').read_text())['run'];i28=(r/'I28zero.txt').read_text().split('=',1)[1].strip()
s='''restart:interface(warnlevel=0):
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f) if type(f,list) then fibrationBasis(f,[z]) else fibrationBasis(convert(f,HlogRegInf),[z]) fi end proc:
Inv:=proc(f) local g;g:=fibrationBasis(convert(subs(z=1/zz,f),HlogRegInf),[zz]);subs(zz=z,g) end proc:
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\\n",n,f);fclose(fd) end proc:
'''+f'Ifinite28:={i28}:\n'
s+='''targets:=select(h->op(1,h)=z and op(1,op(2,h))=-1/2,indets(Ifinite28,specfunc(anything,Hlog))):
if nops(targets)=0 then error "expected a new-letter sector" fi:
Tprimitive:=0:Phi5:=0:
for hh in targets do
 cc:=coeff(expand(Ifinite28),hh,1):
 if indets(cc,specfunc(anything,Hlog))<>{} then error "nonlinear Hlog coefficient" fi:
 tailword:=[op(2..nops(op(2,hh)),op(2,hh))]:
 Tprimitive:=Tprimitive+32*cc*hh:Phi5:=Phi5+32*cc*Hlog(z,tailword):
od:
Phi5:=FB(Phi5):Tprimitive:=FB(Tprimitive):
Omega6:=FB(Tprimitive+Inv(Tprimitive)):
residual28:=FB(Ifinite28-Omega6/32):
Emit("Phi5",Phi5):Emit("Tprimitive",Tprimitive):Emit("Omega6",Omega6):Emit("residual28",residual28):
if FB(Inv(Omega6)-Omega6)<>0 then error "Omega crossing failure" fi:
if FB(diff(Tprimitive,z)-Phi5/(z+1/2))<>0 then error "primitive derivative failed" fi:
save Phi5,Tprimitive,Omega6,residual28,"omega_basis.m":
printf("OMEGA_COMPLETE lengths %d %d %d\\n",length(Phi5),length(Omega6),length(residual28)):
quit:
'''
(h/'omega_basis.mpl').write_text(s,encoding='utf-8')