restart:
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_check_divergences := true:
_hyper_abort_on_divergence := true:
_hyper_verbosity := 0:
Digits := 35:
printf("BEGIN_T1\n"):
# F_2(e) = int_0^1 int_0^1 x/(1-x*u) exp(e*(-log(x)-log(1-x)+2*log(1-u))) dx du.
LF := -log(x)-log(1-x)+2*log(1-u):
Fcoeff := table():
for q from 0 to 5 do
  tic := time():
  Fcoeff[q] := hyperInt(expand(x/(1-x*u)*LF^q/q!),[u=0..1,x=0..1]):
  Fcoeff[q] := fibrationBasis(Fcoeff[q],[]):
  Fplain := simplify(convert(Fcoeff[q],Hlog)):
  printf("F %d = %a ; seconds %a\n",q,Fplain,time()-tic):
  fd := fopen(cat("F",q,".txt"),WRITE): fprintf(fd,"%a\n",Fplain): fclose(fd):
  save Fcoeff, "Fcoeff.m":
od:
# Z(e) = int x/(y*(1-x+r*x))*exp(e*LZ)*(1-exp(e*log(1-r*y))) dx dr dy.
LZ := log(r)-log(x)-log(1-x)+2*log(y):
Zcoeff := table():
for q from 1 to 4 do
  tic := time():
  Zkernel := -x/(y*(1-x+r*x))*add(log(1-r*y)^j*LZ^(q-j)/(j!*(q-j)!),j=1..q):
  Zcoeff[q] := hyperInt(expand(Zkernel),[y=0..1,x=0..1,r=0..1]):
  Zcoeff[q] := fibrationBasis(Zcoeff[q],[]):
  Zplain := simplify(convert(Zcoeff[q],Hlog)):
  printf("Z %d = %a ; seconds %a\n",q,Zplain,time()-tic):
  fd := fopen(cat("Z",q,".txt"),WRITE): fprintf(fd,"%a\n",Zplain): fclose(fd):
  save Zcoeff, "Zcoeff.m":
od:
printf("COMPLETE_T1\n"):
quit:
