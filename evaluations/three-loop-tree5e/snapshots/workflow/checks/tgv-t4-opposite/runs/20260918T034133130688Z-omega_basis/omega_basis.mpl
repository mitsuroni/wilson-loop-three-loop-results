restart:interface(warnlevel=0):
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f) if type(f,list) then fibrationBasis(f,[z]) else fibrationBasis(convert(f,HlogRegInf),[z]) fi end proc:
Inv:=proc(f) local g;g:=fibrationBasis(convert(subs(z=1/zz,f),HlogRegInf),[zz]);subs(zz=z,g) end proc:
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\n",n,f);fclose(fd) end proc:
Ifinite28:=-31/32*ln(2)*zeta[5]+zeta[2]*zeta[1,-3]+1/2*zeta[1,-5]-547/840*zeta[2]^3-13/128*zeta[3]^2-7/4*ln(2)*zeta[2]*zeta[3]-11/16*zeta[2]*Hlog(z,[-1/2, -1, 0, 0])+1/4*zeta[3]*Hlog(z,[-2, -1, 0])-7/16*zeta[2]*Hlog(z,[-1/2, 0, -1, 0])-15/16*zeta[2]*Hlog(z,[-1/2, 0, 0, -1])-1/4*zeta[2]*Hlog(z,[-2, -1, -1, 0])+25/64*zeta[3]*Hlog(z,[0, 0, 0])-7/16*zeta[2]*Hlog(z,[0, 0, 0, 0])+25/32*zeta[2]^2*Hlog(z,[-1/2, 0])+11/160*zeta[2]^2*Hlog(z,[-2, -1])+1/2*zeta[2]*Hlog(z,[-2, -1, 0, 0])+3/16*zeta[2]*Hlog(z,[0, -1, -1, -1])+7/5*zeta[2]^2*Hlog(z,[0, -1])+3/16*zeta[2]*Hlog(z,[-1/2, 0, -1, -1])-1/16*zeta[3]*Hlog(z,[-1/2, -1, -1])-15/8*zeta[2]^2*Hlog(z,[-1, -1])-139/160*zeta[2]^2*Hlog(z,[-1, 0])-57/64*zeta[2]*Hlog(z,[0, 0, 0, -1])+17/16*zeta[2]*Hlog(z,[0, -1, 0, 0])+13/32*zeta[2]*Hlog(z,[-1, 0, -1, 0])-1/16*zeta[3]*Hlog(z,[-2, -1, -1])-3/16*zeta[3]*Hlog(z,[-1/2, -1, 0])+81/32*zeta[2]*Hlog(z,[-1, 0, 0, -1])-15/32*zeta[2]*Hlog(z,[0, -1, 0, -1])-3/16*zeta[2]*Hlog(z,[-1/2, -1, -1, -1])+11/16*zeta[2]*Hlog(z,[-1/2, 0, 0, 0])+1/16*zeta[3]*Hlog(z,[-1/2, 0, -1])+(-7/4*zeta[2]*zeta[3]-31/32*zeta[5])*Hlog(z,[-2])-51/160*zeta[2]^2*Hlog(z,[0, 0])-3/4*zeta[2]*Hlog(z,[-2, -1, 0, -1])-3/16*zeta[2]*Hlog(z,[-2, -1, -1, -1])+7/16*zeta[2]*Hlog(z,[-1/2, -1, -1, 0])-1/32*zeta[3]*Hlog(z,[0, -1, 0])-13/16*zeta[2]*Hlog(z,[-1, 0, 0, 0])-37/32*zeta[3]*Hlog(z,[-1, 0, 0])+15/32*zeta[2]*Hlog(z,[0, 0, -1, -1])+3/16*zeta[3]*Hlog(z,[-1/2, 0, 0])-25/32*zeta[2]^2*Hlog(z,[-1/2, -1])-9/64*zeta[2]*Hlog(z,[0, 0, -1, 0])+5/32*zeta[3]*Hlog(z,[0, 0, -1])-11/32*zeta[2]*Hlog(z,[0, -1, -1, 0])-21/16*zeta[2]*Hlog(z,[-1, -1, 0, 0])+1/16*zeta[3]*Hlog(z,[0, -1, -1])+15/16*zeta[2]*Hlog(z,[-1/2, -1, 0, -1])+(31/32*zeta[2]*zeta[3]+33/16*zeta[5])*Hlog(z,[-1])+(-1/2*zeta[2]*zeta[3]-7/32*zeta[5])*Hlog(z,[-1/2])+(7/64*zeta[2]*zeta[3]-21/32*zeta[5])*Hlog(z,[0])+5/16*Hlog(z,[-1/2, -1, 0, -1, 0, 0])-3/8*Hlog(z,[-1/2, 0, -1, 0, 0, 0])-5/32*Hlog(z,[0, -1, 0, -1, 0, 0])+1/16*Hlog(z,[-1, 0, 0, 0, 0, 0])+1/16*Hlog(z,[0, -1, -1, -1, 0, 0])-1/16*Hlog(z,[-2, -1, -1, -1, 0, 0])-5/16*Hlog(z,[-2, -1, -1, 0, 0, 0])+13/16*Hlog(z,[0, -1, 0, 0, 0, 0])+27/32*Hlog(z,[-1, 0, 0, -1, 0, 0])-1/4*Hlog(z,[-2, -1, 0, -1, 0, 0])-1/16*Hlog(z,[-1/2, -1, -1, -1, 0, 0])-5/16*Hlog(z,[-1/2, 0, 0, -1, 0, 0])-19/64*Hlog(z,[0, 0, 0, -1, 0, 0])-25/32*Hlog(z,[0, 0, 0, 0, 0, 0])+13/32*Hlog(z,[-1, 0, -1, 0, 0, 0])+3/8*Hlog(z,[-1/2, -1, -1, 0, 0, 0])+1/64*Hlog(z,[0, 0, -1, 0, 0, 0])-1/4*Hlog(z,[-1/2, -1, 0, 0, 0, 0])-11/8*Hlog(z,[-1, -1, 0, 0, 0, 0])-9/32*Hlog(z,[0, -1, -1, 0, 0, 0])+3/8*Hlog(z,[-2, -1, 0, 0, 0, 0])+1/4*Hlog(z,[-1/2, 0, 0, 0, 0, 0])+1/16*Hlog(z,[-1/2, 0, -1, -1, 0, 0])+5/32*Hlog(z,[0, 0, -1, -1, 0, 0]):
targets:=select(h->op(1,h)=z and op(1,op(2,h))=-1/2,indets(Ifinite28,specfunc(anything,Hlog))):
if nops(targets)=0 then error "expected a new-letter sector" fi:
Tprimitive:=0:Phi5:=0:
for hh in targets do
 cc:=coeff(expand(Ifinite28),hh,1):
 if indets(cc,specfunc(anything,Hlog))<>{} then error "nonlinear Hlog coefficient" fi:
 tailword:=[op(2..nops(op(2,hh)),op(2,hh))]:
 Tprimitive:=Tprimitive+32*cc*hh:Phi5:=Phi5+32*cc*Hlog(z,tailword):
od:
Phi5:=FB(Phi5):Tprimitive:=FB(Tprimitive):
Omega6:=FB(Tprimitive+Inv(Tprimitive)):
residual28:=FB(Ifinite28-Omega6/32):
Emit("Phi5",Phi5):Emit("Tprimitive",Tprimitive):Emit("Omega6",Omega6):Emit("residual28",residual28):
if FB(Inv(Omega6)-Omega6)<>0 then error "Omega crossing failure" fi:
if FB(diff(Tprimitive,z)-Phi5/(z+1/2))<>0 then error "primitive derivative failed" fi:
save Phi5,Tprimitive,Omega6,residual28,"omega_basis.m":
printf("OMEGA_COMPLETE lengths %d %d %d\n",length(Phi5),length(Omega6),length(residual28)):
quit:
