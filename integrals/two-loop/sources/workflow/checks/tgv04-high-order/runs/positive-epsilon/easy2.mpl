restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=1:
_hyper_check_divergences:=true:
_hyper_abort_on_divergence:=true:
tstart:=time():
AA:=a*(1-v)-x*(a-v):
la:=ln(AA)-ln(a):lxz:=ln(x)+ln(1-a)-ln(a):l0:=ln(1-x):lv:=ln(v)+ln(1-v):
Ecoeff:=proc(A,B) local j; add((-lv)^j/j!*2^(3-j)/(3-j+1)!*(A^(3-j+1)-B^(3-j+1)),j=0..3) end proc:
fnew:=a*(Ecoeff(la,lxz)/(1-v)-Ecoeff(l0,la)/v)/(a-x):
ans:=timelimit(600,hyperInt(convert(fnew,HlogRegInf),[v=0..1,x=0..1])):
Da:=fibrationBasis(ans,[a]):
D3:=fibrationBasis(convert(subs(a=1/(1+z),Da),HlogRegInf),[z]):
save D3,"D3.m":
fd:=fopen("D3.txt",WRITE):fprintf(fd,"D3 = %a\n",D3):fclose(fd):
printf("COMPLETE easy2 CPU=%a\n",time()-tstart):quit:
