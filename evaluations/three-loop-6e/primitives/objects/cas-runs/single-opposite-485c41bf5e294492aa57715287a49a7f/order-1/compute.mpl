restart:
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
_hyper_ignore_nonlinear_polynomials:=false:
try
  wbBase:=-(1+z)/((1+x)*(1+y)*(1+z*x*y)):
  wbLog:=ln(1+z*x*y)-ln(1+x)-ln(1+y)-ln(z)/2
    +(0)*(ln(x)-ln(1+x))-(2)*ln(1+x)
    +(0)*(ln(y)-ln(1+y))-(0)*ln(1+y):
  wbAnswer:=hyperInt(convert(wbBase*wbLog^(1)/(1)!,HlogRegInf),[x,y]):
  wbRaw:=fibrationBasis(wbAnswer,[z]):
  wbRef:=subs(zz=z,fibrationBasis(convert(subs(z=1/zz,wbRaw),HlogRegInf),[zz])):
  for wbChannel in ["raw","reflected"] do
    if wbChannel="raw" then wbValue:=wbRaw else wbValue:=wbRef fi:
    if type(wbValue,`+`) then wbCount:=nops(wbValue) else wbCount:=1 fi:
    if wbCount>100000 then error "HLOG word budget exceeded" fi:
    wbFile:=fopen(cat(wbChannel,".txt"),WRITE,TEXT):
    fprintf(wbFile,"%a\n",wbValue):fclose(wbFile):
  od:
  printf("WILSON_MOMENT_COMPLETE order=1\n"):
catch:
  printf("WILSON_MOMENT_FAILED %a\n",[lastexception]):
end try:
quit:
