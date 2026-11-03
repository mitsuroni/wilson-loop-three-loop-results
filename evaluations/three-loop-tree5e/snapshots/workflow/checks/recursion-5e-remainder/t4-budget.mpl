restart:
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:_hyper_verbosity:=0:
FB:=proc(f) if type(f,list) then fibrationBasis(f,[z]) else fibrationBasis(convert(f,HlogRegInf),[z]) fi end proc:
Inv:=proc(f) local g;g:=fibrationBasis(convert(subs(z=1/zz,f),HlogRegInf),[zz]);subs(zz=z,g) end proc:
Emit:=proc(name,val) local fd;fd:=fopen(cat(name,".txt"),WRITE);fprintf(fd,"%s = %a\n",name,val);fclose(fd) end proc:
try
read "t4-budget-inputs.mpl":
printf("BEGIN_NEW_T4_BUDGET\n"):
printf("SERIES_POLE_TEST %a\n",series(exp(eps)/eps^2,eps=0,3)):
r4Pc:=[1,0,5*zeta[2]/2-9*ln(z)^2/8,3*zeta[3],193*zeta[2]^2/40-45*zeta[2]*ln(z)^2/16+135*ln(z)^4/128]:
r4Correct:=table():r4Bug:=table():r4Even:=table():
for i in [13,14,15,40] do
 r4OldQ:=add(oldQ[i,q]*eps^q,q=-2..2):
 # Reproduce the actual historical call. Maple's order is relative to the
 # leading Laurent power; it differs between the four representatives.
 r4OldEven:=convert(series(expand((exp(-eps*ln(z)/2)*r4OldQ+exp(eps*ln(z)/2)*Inv(r4OldQ))/2),eps=0,3),polynom):
 r4OldOut:=expand(add(r4Pc[j+1]*eps^j,j=0..4)*r4OldEven/(32*eps^2)):
 printf("ACTUAL_OLD_SERIES_RANGE %d %d %d\n",i,ldegree(r4OldEven,eps),degree(r4OldEven,eps)):
 for k from -2 to 2 do
  r4Even[i,k]:=FB(add(((-ln(z)/2)^j*oldQ[i,k-j]+(ln(z)/2)^j*Inv(oldQ[i,k-j]))/(2*j!),j=0..k+2)):
 od:
 for k from -6 to 0 do
  if k< -4 then r4cv:=0
  else r4cv:=FB(add(r4Pc[j+1]*r4Even[i,k+2-j],j=0..k+4)/32):
  fi:
  r4bv:=FB(coeff(r4OldOut,eps,k)):
  r4Correct[i,k]:=r4cv:r4Bug[i,k]:=r4bv:
  Emit(cat("Correct",i,"q",k),r4cv):Emit(cat("OldBudget",i,"q",k),r4bv):
 od:
od:
for k from -6 to 0 do
 Emit(cat("CorrectEight",k),FB(2*add(r4Correct[i,k],i in [13,14,15,40]))):
 Emit(cat("OldBudgetEight",k),FB(2*add(r4Bug[i,k],i in [13,14,15,40]))):
od:
printf("COMPLETE_NEW_T4_BUDGET\n"):
catch:
printf("FAILED_NEW_T4_BUDGET %a\n",lasterror):
end try:
quit:
