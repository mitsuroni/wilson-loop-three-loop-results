restart:interface(warnlevel=0):
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=1:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f,vs) if type(f,list) then fibrationBasis(f,vs) else fibrationBasis(convert(f,HlogRegInf),vs) fi end proc:
Inv:=proc(f) local g;g:=FB(subs(z=1/zz,f),[zz]);subs(zz=z,g) end proc:
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\n",n,f);fclose(fd) end proc:
read "D:/wilson-agent/workflow/checks/tgv-t4-opposite/runs/20260918T032241443101Z-partial/partial_kernel.m":
read "D:/wilson-agent/workflow/checks/tgv-t4-opposite/runs/20260918T032322455592Z-weights/weights.m":
Jtv:=1/((1+t)^2*(1+v)^2):
T4AtOrigin:=proc(f,var,vs) local hs,sub,h;
 hs:=select(h->op(1,h)=var,indets(f,specfunc(anything,Hlog))):sub:={}:
 for h in hs do if {op(op(2,h))} subset {0} then error "logarithmically singular zero limit" fi:sub:=sub union {h=0} od:
 FB(subs(var=0,subs(sub,f)),vs) end proc:
Hleft:=T4AtOrigin(HH,t,[v,z]):
Htop:=FB(subs({v=1/t,z=1/z},Hleft)/z,[t,z]):
Emit("Hleft",Hleft):Emit("Htop",Htop):
L0:=ln(1-p):L1:=ln(z)+ln(p):
Cb1:=FB(hyperInt(-2*L0*(L1-L0)/((1+z)*p-1),[p=0..1]),[z]):
Cb2:=FB(hyperInt((-L0*(L1^2-L0^2)-2*L0^2*(L1-L0))/((1+z)*p-1),[p=0..1]),[z]):
# F0 is rederived here from its finite two-dimensional definition.
Fa0:=FB(hyperInt(1/(1-b*a),[a=0..1,b=0..1]),[]):
Fa1:=FB(hyperInt((-ln(b)-ln(1-b)+2*ln(1-b*a))/(1-b*a),[a=0..1,b=0..1]),[]):
CuspM:=FB(Fa0*Cb1/2,[z]):CuspZ:=FB((Fa0*Cb2+Fa1*Cb1)/2,[z]):
Emit("CuspM",CuspM):Emit("CuspZ",CuspZ):
for id in [27,28] do
 printf("START global id%d pole time%a\n",id,time()):
 poleInteg:=((1+z)*Jtv*(Y[1,0]-Y[0,0])+4*Rm1)*W[id,0]:
 Vpole:=FB(hyperInt(poleInteg,[t,v]),[z]):
 if id=27 then Vpole:=FB(Vpole+CuspM,[z]) fi:
 Vm1:=FB(-(1+z)*Vpole,[z]):Emit(cat("V",id,"m1"),Vm1):
 printf("START global id%d finite time%a\n",id,time()):
 contact:=-z/(1+v)^2*HH*diff(W[id,0],t)+HH/(1+t)^2*diff(W[id,0],v):
 finiteInteg:=contact+(1+z)*Jtv*((Y[1,1]-Y[0,1])*W[id,0]+(Y[1,0]-Y[0,0])*W[id,1])+4*((Rzero-2*Rm1)*W[id,0]+Rm1*W[id,1]):
 bulk:=FB(hyperInt(finiteInteg,[t,v]),[z]):Emit(cat("bulk",id),bulk):
 Wleft:=T4AtOrigin(W[id,0],t,[v,z]):
 faceX:=FB(hyperInt(Hleft*Wleft/(1+v)^2,[v]),[z]):Emit(cat("faceX",id),faceX):
 if id=28 then
  Wtop:=FB(hyperInt(1/((1+p)*(1+q)*(1+z*p*q)),[q,p=t..infinity]),[t,z]):
  faceW:=FB(hyperInt(Htop*Wtop/(1+t)^2,[t]),[z]):Emit("faceW28",faceW):
  Vzero:=FB(-(1+z)*(bulk-z*faceX-faceW),[z]):
 else Vzero:=FB(-(1+z)*(bulk-z*faceX+CuspZ),[z]) fi:
 Emit(cat("V",id,"zero"),Vzero):
 Im1:=FB((Vm1+Inv(Vm1))/64,[z]):Izero:=FB((Vzero+Inv(Vzero))/64+3*ln(z)*(Inv(Vm1)-Vm1)/128,[z]):
 Emit(cat("I",id,"m1"),Im1):Emit(cat("I",id,"zero"),Izero):
 if FB(Inv(Im1)-Im1,[z])<>0 or FB(Inv(Izero)-Izero,[z])<>0 then error "I crossing failed" fi:
 if id=28 then if FB(Inv(Vm1)-Vm1,[z])<>0 or FB(Inv(Vzero)-Vzero+3*ln(z)*Vm1,[z])<>0 then error "28 self crossing failed" fi fi:
 printf("DONE global id%d sizes%d %d time%a\n",id,length(Im1),length(Izero),time()):gc():
od:
printf("GLOBAL_COMPLETE\n"):quit: