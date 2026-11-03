restart:interface(warnlevel=0):
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=1:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f,vs) if type(f,list) then fibrationBasis(f,vs) else fibrationBasis(convert(f,HlogRegInf),vs) fi end proc:
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\n",n,f);fclose(fd) end proc:
Tpol:=1+r+s:Xpol:=(1+t)*(1+u)*(1+v):
Fpol:=r*t*u*v*z+r*s*u*v+r*t*v*z+s*t*v*z+r*s*u+s*v*z+r*u+r:
bzero:=1/(Tpol*Xpol*Fpol^2):
RR:=(1+z)/2*(1+u)*(r*(1+v)+z*(r-s)*v*(1+t)):
Eu:=proc(f) r*diff(f,r)-u*diff(f,u) end proc:
ell:=-eps+(-1-eps)*r/Tpol+(1+2*eps)*u/(1+u)+(-2+2*eps)*Eu(Fpol)/Fpol:
Rg:=factor(normal((eps*RR-Eu(RR)-RR*ell)/eps)):
rm1:=normal(eval(eps*Rg,eps=0)):rz:=normal(eval(diff(eps*Rg,eps),eps=0)):
Lb:=-ln(r)-ln(s)-ln(Tpol)-2*ln(1+t)-2*ln(1+u)-2*ln(1+v)+2*ln(Fpol):
for q from -1 to 0 do
 if q=-1 then kern:=rm1*bzero else kern:=(rz+rm1*Lb)*bzero fi:
 printf("START partial R%d time%a\n",q,time()):
 ans:=hyperInt(kern,[r,u,s]):ans:=FB(ans,[t,v,z]):
 if q=-1 then Rm1:=ans else Rzero:=ans fi:
 Emit(cat("Rpartial",q),ans):save Rm1,Rzero,"partialR.m":
 printf("DONE partial R%d length%d time%a\n",q,length(ans),time()):gc():
od:
printf("START partial H0 time%a\n",time()):
Hdensity:=FB(hyperInt(1/(Tpol*Xpol*Fpol),[r,u,s]),[t,v,z]):
HH:=normal(Hdensity*(1+t)^2*(1+v)^2):
Emit("Hpartial",HH):save HH,"partialH.m":
printf("DONE partial H0 length%d time%a\n",length(HH),time()):
xx:=t/(1+t):ww:=v/(1+v):Jtv:=1/((1+t)^2*(1+v)^2):
Y:=table():
for face from 0 to 1 do
 if face=0 then dpol:=(1-b)*(1-xx)*(1-ww)+z*(b+(1-b)*xx)*ww
 else dpol:=(1-xx)*(1-(1-b)*ww)+z*xx*(1-b)*ww fi:
 for q from 0 to 1 do
  ans:=FB(hyperInt((2*ln(dpol)-ln(b)-ln(1-b))^q/(q!*dpol),[b=0..1]),[t,v,z]):
  Y[face,q]:=ans:Emit(cat("Y",face,"q",q),ans):
 od:
od:
Qm1:=FB((1+z)*Jtv*(Y[1,0]-Y[0,0])+4*Rm1,[t,v,z]):
Qzero:=FB(z/(1+v)^2*diff(HH,t)-1/(1+t)^2*diff(HH,v)+(1+z)*Jtv*(Y[1,1]-Y[0,1])+4*Rzero-8*Rm1,[t,v,z]):
Emit("Qpartial_m1",Qm1):Emit("Qpartial_0",Qzero):
save Qm1,Qzero,HH,Y,Rm1,Rzero,Rg,Fpol,"partial_kernel.m":
printf("PARTIAL_COMPLETE sizes%d %d time%a\n",length(Qm1),length(Qzero),time()):quit: