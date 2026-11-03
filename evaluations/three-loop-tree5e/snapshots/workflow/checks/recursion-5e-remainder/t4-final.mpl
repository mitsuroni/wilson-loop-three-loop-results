restart:
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:_hyper_verbosity:=0:
FB:=proc(f) if type(f,list) then fibrationBasis(f,[z]) else fibrationBasis(convert(f,HlogRegInf),[z]) fi end proc:
Emit:=proc(name,val) local fd;fd:=fopen(cat(name,".txt"),WRITE);fprintf(fd,"%s = %a\n",name,val);fclose(fd) end proc:
try
read "t4-final-inputs.mpl":
printf("BEGIN_NEW_T4_FINAL\n"):
r4C:=[1,0,(zeta[2]-ln(z)^2)/2,zeta[3]/3,9*zeta[2]^2/40+11*ln(z)^4/24-zeta[2]*ln(z)^2/4]:
for k from -6 to 0 do
 if k< -4 then r4low:=0 else r4low:=add(r4C[n+1]*r4P[k+2-n],n=0..k+4) fi:
 r4S:=FB(r4low-r4A[k]):Emit(cat("AdjacentEffective",k),r4S):Emit(cat("AdjacentLowerProduct",k),FB(r4low)):
od:
printf("COMPLETE_NEW_T4_FINAL\n"):
catch:
printf("FAILED_NEW_T4_FINAL %a\n",lasterror):
end try:
quit:
