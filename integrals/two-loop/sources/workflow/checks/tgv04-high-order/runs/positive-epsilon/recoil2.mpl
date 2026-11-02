restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=1:
_hyper_check_divergences:=true:
_hyper_abort_on_divergence:=true:
tstart:=time():
read "step09_optimized_input.m":
read "step10_optimized_integrands.m":
if numelems([indices(dnew)])>0 then error "recoil source divergent" fi:
if normal(diff(eps*Rreg,eps$2))<>0 then error "unexpected epsilon dependence" fi:
b0:=normal(eval(base,eps=0)): Rzero:=normal(eval(diff(eps*Rreg,eps),eps=0)):
Lb:=-ln(r)-ln(s)-ln(1+r+s)-2*ln(1+t)-2*ln(1+u)-2*ln(1+v)+2*ln(F):
fq:=normal(Rzero*b0)*Lb^2/2!+fm1*Lb^(2+1)/(2+1)!:
save fq,"O2-integrand.m":
ans:=timelimit(900,hyperInt(convert(fq,HlogRegInf),[r,u,s,t,v])):
O2:=fibrationBasis(ans,[z]):
save O2,"O2.m":
fd:=fopen("O2.txt",WRITE):fprintf(fd,"O2 = %a\n",O2):fclose(fd):
printf("COMPLETE recoil2 CPU=%a\n",time()-tstart):quit:
