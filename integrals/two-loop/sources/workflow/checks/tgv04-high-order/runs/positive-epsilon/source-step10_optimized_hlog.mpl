restart:
interface(warnlevel=0):
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity := 1:
_hyper_check_divergences := true:
_hyper_abort_on_divergence := true:
read "step09_optimized_input.m":
Eu := proc(f) r*diff(f,r)-u*diff(f,u) end proc:
ell := -eps+(-1-eps)*r/(1+r+s)+(1+2*eps)*u/(1+u)+(-2+2*eps)*Eu(F)/F:
Rreg := factor(normal((eps*RR-Eu(RR)-RR*ell)/eps)):
printf("regularized length=%d\n",length(Rreg));
dnew := findDivergences(Rreg*base,{z},[eps=0]):
printf("remaining=%a\n",eval(dnew));
if numelems([indices(dnew)])>0 then error "remaining divergences" fi:
b0 := normal(eval(base,eps=0)):
Rm1 := normal(eval(eps*Rreg,eps=0)):
Rzero := normal(eval(diff(eps*Rreg,eps),eps=0)):
Lb := -ln(r)-ln(s)-ln(1+r+s)-2*ln(1+t)-2*ln(1+u)-2*ln(1+v)+2*ln(F):
fm1 := normal(Rm1*b0):
fzero := normal(Rzero*b0)+fm1*Lb:
save Rreg,dnew,fm1,fzero,"step10_optimized_integrands.m":
printf("START optimized recoil pole length=%d\n",length(fm1));
outm1 := timelimit(240,hyperInt(fm1,[r,u,s,t,v])):
Om1 := fibrationBasis(outm1,[z]):
printf("RESULT Om1=%a\n",Om1);
save Om1,"step10_Om1.m":
printf("START optimized recoil finite length=%d\n",length(fzero));
out0 := timelimit(300,hyperInt(convert(fzero,HlogRegInf),[r,u,s,t,v])):
Ozero := fibrationBasis(out0,[z]):
printf("RESULT Ozero=%a\n",Ozero);
save Om1,Ozero,"step10_optimized_results.m":
printf("STEP10 COMPLETE\n");
quit:
