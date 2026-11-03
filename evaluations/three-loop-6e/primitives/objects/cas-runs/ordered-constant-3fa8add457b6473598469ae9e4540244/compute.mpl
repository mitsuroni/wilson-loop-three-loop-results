restart:
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
_hyper_ignore_nonlinear_polynomials:=false:
try
  woInput:=(((-(((1)+(z)))))^(3))/(((((((1)-(wox1)))*(((1)-(wox2)))))+(((((z)*(wox1)))*(wox2)))))/(((((((1)-(wox3)))*(((1)-(wox4)))))+(((((z)*(wox3)))*(wox4)))))/(((((((1)-(wox5)))*(((1)-(wox6)))))+(((((z)*(wox5)))*(wox6))))):
  woAnswer:=hyperInt(woInput,[wox4=0..wox2,wox1=0..wox3,wox2=0..wox6,wox3=0..wox5,wox6=0..1,wox5=0..1]):
  woRaw:=fibrationBasis(woAnswer,[z]):
  woRef:=subs(wozz=z,fibrationBasis(convert(subs(z=1/wozz,woRaw),HlogRegInf),[wozz])):
  for woChannel in ["raw","reflected"] do
    if woChannel="raw" then woValue:=woRaw else woValue:=woRef fi:
    if type(woValue,`+`) then woCount:=nops(woValue) else woCount:=1 fi:
    if woCount>100000 then error "HLOG word budget exceeded" fi:
    woFile:=fopen(cat(woChannel,".txt"),WRITE,TEXT):
    fprintf(woFile,"%a\n",woValue):fclose(woFile):
  od:
  printf("WILSON_ORDERED_CONSTANT_COMPLETE\n"):
catch:
  printf("WILSON_ORDERED_CONSTANT_FAILED %a\n",[lastexception]):
end try:
quit:
