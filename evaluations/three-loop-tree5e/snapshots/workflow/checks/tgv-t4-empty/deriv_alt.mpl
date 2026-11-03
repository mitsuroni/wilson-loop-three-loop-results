restart:interface(warnlevel=0):
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\n",n,f);fclose(fd) end proc:
read "D:/wilson-agent/workflow/checks/tgv-t4-empty/runs/20260918T014950651961Z-derivatives/regularized3.m":
Rp:=normal(Rg*eps^3):
for q from -3 to 1 do
 kern:=add(coeff(Rp,eps,j)*Lf^(q+3-j)/(q+3-j)!,j=0..q+3)*bzero:
 kern:=collect(expand(kern),indets(kern,specfunc(anything,ln)),normal):
 printf("START ALT D3 q%d length%d time%a\n",q,length(kern),time()):
 ans:=timelimit(300,hyperInt(kern,[t,r,u,s,v])):
 ans:=fibrationBasis(ans,[z]):Emit(cat("D3q",q),ans):
 printf("DONE ALT D3 q%d length%d time%a\n",q,length(ans),time()):gc():
od:
printf("ALT_D3_COMPLETE\n"):quit: