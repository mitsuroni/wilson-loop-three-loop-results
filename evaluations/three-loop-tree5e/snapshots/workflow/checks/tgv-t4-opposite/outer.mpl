restart:interface(warnlevel=0):
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=1:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f) if type(f,list) then fibrationBasis(f,[z]) else fibrationBasis(convert(f,HlogRegInf),[z]) fi end proc:
Inv:=proc(f) local g;g:=fibrationBasis(convert(subs(z=1/zz,f),HlogRegInf),[zz]);subs(zz=z,g) end proc:
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\n",n,f);fclose(fd) end proc:
read "D:/wilson-agent/workflow/checks/tgv-t4-opposite/runs/20260918T032241443101Z-partial/partial_kernel.m":
read "D:/wilson-agent/workflow/checks/tgv-t4-opposite/runs/20260918T032322455592Z-weights/weights.m":
for id in [27,28] do
 printf("START V%d pole time%a\n",id,time()):
 Vm1:=FB(hyperInt(-(1+z)*Qm1*W[id,0],[t,v])):
 Emit(cat("V",id,"m1"),Vm1):
 printf("DONE V%d pole length%d time%a\n",id,length(Vm1),time()):
 printf("START V%d finite time%a\n",id,time()):
 Vzero:=FB(hyperInt(-(1+z)*(Qzero*W[id,0]+Qm1*W[id,1]),[t,v])):
 Emit(cat("V",id,"zero"),Vzero):
 Im1:=FB((Vm1+Inv(Vm1))/64):
 Izero:=FB((Vzero+Inv(Vzero))/64+3*ln(z)*(Inv(Vm1)-Vm1)/128):
 Emit(cat("I",id,"m1"),Im1):Emit(cat("I",id,"zero"),Izero):
 if FB(Inv(Im1)-Im1)<>0 or FB(Inv(Izero)-Izero)<>0 then error "I crossing failed" fi:
 if id=28 then
  if FB(Inv(Vm1)-Vm1)<>0 or FB(Inv(Vzero)-Vzero+3*ln(z)*Vm1)<>0 then error "28 self crossing failed" fi:
 fi:
 printf("DONE I%d length%d time%a\n",id,length(Izero),time()):gc():
od:
printf("OUTER_COMPLETE\n"):quit:
