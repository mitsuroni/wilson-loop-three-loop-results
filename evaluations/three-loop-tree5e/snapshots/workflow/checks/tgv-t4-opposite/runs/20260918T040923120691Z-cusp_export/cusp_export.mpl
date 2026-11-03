restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f) if type(f,list) then fibrationBasis(f,[z]) else fibrationBasis(convert(f,HlogRegInf),[z]) fi end proc:
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\n",n,f);fclose(fd) end proc:
Emit("CuspMnorm",FB(-(1+z)*(3*zeta[2]^2/(z+1)*Hlog(z,[-1])+zeta[2]/(z+1)*Hlog(z,[-1, 0, 0])+zeta[2]*zeta[3]/(z+1)-zeta[2]/(z+1)*Hlog(z,[0, 0, 0])-2*zeta[2]^2/(z+1)*Hlog(z,[0])))):
Emit("CuspZnorm",FB(-(1+z)*(4*zeta[2]/(z+1)*Hlog(z,[0, -1, 0, 0])+5*zeta[2]/(z+1)*Hlog(z,[-1, 0, 0, 0])-zeta[2]*zeta[3]/(z+1)*Hlog(z,[-1])+zeta[3]/(z+1)*Hlog(z,[-1, 0, 0])-1/10*(81*zeta[2]^3-10*zeta[3]^2)/(z+1)-4*zeta[2]/(z+1)*Hlog(z,[-1, -1, 0, 0])-12*zeta[2]^2/(z+1)*Hlog(z,[-1, -1])+9*zeta[2]^2/(z+1)*Hlog(z,[-1, 0])+12*zeta[2]^2/(z+1)*Hlog(z,[0, -1])-8*zeta[2]^2/(z+1)*Hlog(z,[0, 0])-zeta[3]/(z+1)*Hlog(z,[0, 0, 0])-zeta[2]*zeta[3]/(z+1)*Hlog(z,[0])-5*zeta[2]/(z+1)*Hlog(z,[0, 0, 0, 0])))):
printf("CUSP_EXPORT_COMPLETE\n"):quit:
