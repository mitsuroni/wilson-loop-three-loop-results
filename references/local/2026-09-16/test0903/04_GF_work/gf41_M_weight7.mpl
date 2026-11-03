restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f) fibrationBasis(convert(f,HlogRegInf),[z]) end proc:
MM:=table():
for wt from 2 to 7 do
 for nn from 1 to wt-1 do
  kk:=wt-nn-1:
  filename:=sprintf("gf41_M_%d_%d.m",kk,nn):
  if FileTools:-Exists(filename) then read filename:
  else
   h:=ln(x)+ln(1-x):L:=ln(1-a)-ln(a)+ln(x)-ln(1-x):
   ans:=hyperInt(convert(h^kk*L^nn/(a-x),HlogRegInf),[x=0..1]):
   aa:=fibrationBasis(ans,[a]):Mn:=FB(subs(a=1/(1+z),aa)):
   save Mn,filename:
  fi:
  MM[kk,nn]:=Mn:
  fd:=fopen(sprintf("gf41_M_%d_%d.txt",kk,nn),WRITE):fprintf(fd,"M_%d_%d = %a\n",kk,nn,Mn):fclose(fd):
  printf("RESULT M[%d,%d] weight=%d length=%d\n",kk,nn,wt,length(Mn)):
 od:
od:
save MM,"gf41_M_all.m":printf("GF41 COMPLETE\n"):quit:
