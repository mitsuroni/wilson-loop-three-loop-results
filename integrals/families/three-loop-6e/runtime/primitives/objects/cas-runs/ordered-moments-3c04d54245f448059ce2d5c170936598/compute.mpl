restart:
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
_hyper_ignore_nonlinear_polynomials:=false:
try
  woRegular:=(exp(woEps*ln(((((((1)-(wox1)))*(((1)-(wox2)))))+(((((z)*(wox1)))*(wox2))))))/(((((((1)-(wox1)))*(((1)-(wox2)))))+(((((z)*(wox1)))*(wox2))))))*(exp(woEps*ln(((((((1)-(wox3)))*(((1)-(wox4)))))+(((((z)*(wox3)))*(wox4))))))/(((((((1)-(wox3)))*(((1)-(wox4)))))+(((((z)*(wox3)))*(wox4))))))*((1-exp(woEps*ln(((1)-(wox1)))))*(1-0)):
  woTaylor:=convert(series(woRegular,woEps=0,3),polynom):
  for woQ from 0 to 2 do
    woInput:=expand(coeff(woTaylor,woEps,woQ)):
    if woInput=0 then woRaw:=0:woRef:=0 else
      woAnswer:=hyperInt(woInput,[wox2=0..wox4,wox1=0..wox3,wox4=0..1,wox3=0..1]):
      woRaw:=fibrationBasis(woAnswer,[z]):
      woRef:=subs(wozz=z,fibrationBasis(convert(subs(z=1/wozz,woRaw),HlogRegInf),[wozz])):
    fi:
    for woChannel in ["raw","reflected"] do
      if woChannel="raw" then woValue:=woRaw else woValue:=woRef fi:
      if type(woValue,`+`) then woCount:=nops(woValue) else woCount:=1 fi:
      if woCount>100000 then error "HLOG word budget exceeded" fi:
      woFile:=fopen(cat("q",woQ,".",woChannel,".txt"),WRITE,TEXT):
      fprintf(woFile,"%a\n",woValue):fclose(woFile):
    od:
    printf("WILSON_ORDERED_MOMENT_%d_COMPLETE\n",woQ):
    gc():
  od:
  printf("WILSON_ORDERED_MOMENTS_COMPLETE\n"):
catch:
  printf("WILSON_ORDERED_MOMENTS_FAILED %a\n",[lastexception]):
end try:
quit:
