restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f) fibrationBasis(convert(f,HlogRegInf),[z]) end proc:
Inv:=proc(f) local g;g:=fibrationBasis(convert(subs(z=1/zz,f),HlogRegInf),[zz]);subs(zz=z,g) end proc:
base:=-(1+z)/((1+x)*(1+y)*(1+z*x*y)):
Ld:=ln(1+z*x*y)-ln(1+x)-ln(1+y)-ln(z)/2:
params:=table(["H0"=[0,0,0,0],"H1"=[0,0,1,0],"H2"=[0,0,2,0],
               "HW"=[0,0,1,1],"HD"=[1,0,1,0],"HX"=[1,0,0,1]]):
for kind in ["H0","H1","H2","HW","HD","HX"] do
 p:=params[kind]:
 Lq:=Ld+p[1]*(ln(x)-ln(1+x))-p[2]*ln(1+x)+p[3]*(ln(y)-ln(1+y))-p[4]*ln(1+y):
 for q from 0 to 4 do
  printf("START %s q=%d time=%a\n",kind,q,time()):
  ans:=hyperInt(convert(base*Lq^q/q!,HlogRegInf),[x,y]):
  Tq:=fibrationBasis(ans,[z]):
  Tplus:=FB((Tq+Inv(Tq))/2):Tminus:=FB((Tq-Inv(Tq))/2):
  for channel in ["plus","minus"] do
   if channel="plus" then val:=Tplus else val:=Tminus fi:
   fd:=fopen(sprintf("%s_%d_%s.txt",kind,q,channel),WRITE):
   fprintf(fd,"%s_%d_%s = %a\n",kind,q,channel,val):fclose(fd):
  od:
  save Tq,Tplus,Tminus,q,p,sprintf("%s_%d.m",kind,q):
  printf("RESULT %s q=%d length=%d time=%a\n",kind,q,length(Tq),time()):
 od:
od:
for wt from 2 to 6 do
 for nn from 1 to wt-1 do
  kk:=wt-nn-1:
  hh:=ln(x)+ln(1-x):LL:=ln(1-a)-ln(a)+ln(x)-ln(1-x):
  ans:=hyperInt(convert(hh^kk*LL^nn/(a-x),HlogRegInf),[x=0..1]):
  aa:=fibrationBasis(ans,[a]):Mn:=FB(subs(a=1/(1+z),aa)):
  fd:=fopen(sprintf("M_%d_%d.txt",kk,nn),WRITE):
  fprintf(fd,"M_%d_%d = %a\n",kk,nn,Mn):fclose(fd):
  printf("RESULT M %d %d\n",kk,nn):
 od:
od:
printf("SINGLE_OPPOSITE_MOMENTS_COMPLETE\n"):
quit:
