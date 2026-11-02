restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=1:
_hyper_check_divergences:=true:
_hyper_abort_on_divergence:=true:
tstart:=time():
read "step04_hard_integrand.m":
if numelems([indices(dv)])>0 then error "hard source divergent" fi:
F:=r*(1+u)+r*s*u*(1+v)+z*s*v:
Lh:=-ln(r)-ln(s)-ln(1+r+s)-2*ln(1+u)-2*ln(1+v)+2*ln(F):
fh:=normal(eval(f,eps=0))*Lh^2/2!:
ans:=timelimit(700,hyperInt(convert(fh,HlogRegInf),[r,u,s,v])):
H2:=fibrationBasis(ans,[z]):
save H2,"H2.m":
fd:=fopen("H2.txt",WRITE):fprintf(fd,"H2 = %a\n",H2):fclose(fd):
printf("COMPLETE hard2 CPU=%a\n",time()-tstart):quit:
