restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f) fibrationBasis(convert(f,HlogRegInf),[z]) end proc:
Inv:=proc(f) local g;g:=fibrationBasis(convert(subs(z=1/zz,f),HlogRegInf),[zz]);subs(zz=z,g) end proc:
Emit:=proc(tag,f)
 local fd;
 fd:=fopen(cat(tag,".txt"),WRITE):fprintf(fd,"%s = %a\n",tag,f):fclose(fd):
end proc:
KernelD:=a*(1-u)*(1-v)+(1-a)*u*v:
LL:=ln(KernelD)-ln(a)-ln(1-v):
LU:=ln(1-a)+ln(u)-ln(KernelD):
Kzero:=LL*LU/((v-a)*(u-a)):
Phi:=ln(KernelD)+(ln(1-v)+ln(u)-ln(a)-ln(1-a))/2:
for q from 0 to 4 do
 poly:=0:
 for ii from 0 to floor(q/2) do
  for jj from 0 to floor((q-2*ii)/2) do
   pp:=q-2*ii-2*jj:
   poly:=poly+Phi^pp*LL^(2*ii)*LU^(2*jj)/(pp!*2^(2*ii+2*jj)*(2*ii+1)!*(2*jj+1)!):
  od:
 od:
 printf("START Q%d time=%a\n",q,time()):
 ans:=hyperInt(Kzero*poly,[u=0..1,v=0..1]):
 f:=FB(subs(a=1/(1+z),fibrationBasis(ans,[a]))):
 Emit(sprintf("Q%d_raw",q),f):
 crossed:=Inv(f):
 Emit(sprintf("Q%d_plus",q),FB((f+crossed)/2)):
 Emit(sprintf("Q%d_minus",q),FB((f-crossed)/2)):
 printf("DONE Q%d length=%d time=%a\n",q,length(f),time()):gc():
od:
# Complete M moments of weight eight, computed rather than guessed.
for nn from 1 to 7 do
 kk:=7-nn:
 Wlog:=ln(x)+ln(1-x):Llog:=ln(1-a)-ln(a)+ln(x)-ln(1-x):
 printf("START M_%d_%d time=%a\n",kk,nn,time()):
 ans:=hyperInt(convert(Wlog^kk*Llog^nn/(a-x),HlogRegInf),[x=0..1]):
 f:=FB(subs(a=1/(1+z),fibrationBasis(ans,[a]))):
 Emit(sprintf("M_%d_%d",kk,nn),f):
 printf("DONE M_%d_%d length=%d time=%a\n",kk,nn,length(f),time()):gc():
od:
printf("FOURE04_HIGH_COMPLETE time=%a\n",time()):
quit:
