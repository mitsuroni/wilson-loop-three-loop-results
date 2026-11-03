restart:interface(warnlevel=0):
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:_hyper_verbosity:=0:
try
read "t4-geometry.mpl":read "t4-paths.mpl":
printf("BEGIN_NEW_T4_RECOIL\n"):
r4out:=table():
for r4id in [1,2] do
 read cat(r4scanpath,"/regularized",r4id,".m"):
 r4T:=1+r+s:r4X:=(1+t)*(1+u)*(1+v):
 r4base:=(r*s)^(-eps)*r4T^(-1-eps)*r4X^(-1-2*eps)*F^(-2+2*eps):
 if r4id=1 then r4weight:=(1+u)^(-eps)*v^eps*(1+v)^(-eps);r4LW:=-ln(1+u)+ln(v)-ln(1+v)
 else r4weight:=(1+v)^(-eps);r4LW:=-ln(1+v) fi:
 r4poly:=normal(eps*r4fR/(r4base*r4weight)):
 if has(denom(r4poly),eps) then error "nonpolynomial regulator" fi:
 r4L:=-ln(r)-ln(s)-ln(r4T)-2*ln(1+t)-2*ln(1+u)-2*ln(1+v)+2*ln(F)+r4LW:
 r4zero:=1/(r4T*r4X*F^2):
 for r4q from -1 to 2 do
  r4kern:=r4zero*add(coeff(r4poly,eps,j)*r4L^(r4q+1-j)/(r4q+1-j)!,j=0..r4q+1):
  r4kern:=collect(expand(r4kern),indets(r4kern,specfunc(anything,ln)),normal):
  save r4kern,cat("inputR",r4id,"q",r4q,".m"):
  printf("RECOIL_START %d %d length=%d time=%a\n",r4id,r4q,length(r4kern),time()):
  r4out[r4id,r4q]:=timelimit(480,hyperInt(r4kern,[r,u,s,t,v])):
  save r4out,"recoil-raw.m":
  printf("RECOIL_DONE %d %d length=%d time=%a\n",r4id,r4q,length(r4out[r4id,r4q]),time()):gc():
 od:
od:
printf("COMPLETE_NEW_T4_RECOIL\n"):
catch:
printf("FAILED_NEW_T4_RECOIL %a\n",lasterror):
end try:
quit:
