restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=1:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f) fibrationBasis(convert(f,HlogRegInf),[z]) end proc:
Inv:=proc(f) local g;g:=fibrationBasis(convert(subs(z=1/zz,f),HlogRegInf),[zz]);subs(zz=z,g) end proc:
# t=x/(1+x), s=y/(1+y). D=(1+z*x*y)/((1+x)*(1+y)).
base:=-(1+z)/((1+x)*(1+y)*(1+z*x*y)):
Ld:=ln(1+z*x*y)-ln(1+x)-ln(1+y)-ln(z)/2:
LAA:=Ld+ln(x)-ln(1+x)+ln(y)-ln(1+y):
Lcross:=Ld+ln(x)-ln(1+x)-ln(1+y):
for kind in ["AA","cross"] do
 for q from 0 to 5 do
  filename:=sprintf("gf51_%s%d.m",kind,q):
  if FileTools:-Exists(filename) then read filename:
  else
   printf("START two-end %s q=%d time=%a\n",kind,q,time()):
   Lq:=`if`(kind="AA",LAA,Lcross):
   ans:=hyperInt(convert(base*Lq^q/q!,HlogRegInf),[x,y]):
   Tq:=fibrationBasis(ans,[z]):
   Tplus:=FB((Tq+Inv(Tq))/2):Tminus:=FB((Tq-Inv(Tq))/2):
   save Tq,Tplus,Tminus,q,filename:
  fi:
  fd:=fopen(sprintf("gf51_%s%d_plus.txt",kind,q),WRITE):fprintf(fd,"T%s%dplus = %a\n",kind,q,Tplus):fclose(fd):
  fd:=fopen(sprintf("gf51_%s%d_minus.txt",kind,q),WRITE):fprintf(fd,"T%s%dminus = %a\n",kind,q,Tminus):fclose(fd):
  printf("RESULT two-end %s q=%d length=%d time=%a\n",kind,q,length(Tplus),time()):
 od:
od:
printf("GF51 COMPLETE\n"):quit:
