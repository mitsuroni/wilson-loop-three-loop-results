restart:interface(warnlevel=0):
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f,vs) if type(f,list) then fibrationBasis(f,vs) else fibrationBasis(convert(f,HlogRegInf),vs) fi end proc:
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\n",n,f);fclose(fd) end proc:
Bpure:=1/((1+p)*(1+q)*(1+z*p*q)):
Lpure:=ln(1+z*p*q)-ln(1+p)-ln(1+q):W:=table():
for id in [27,28] do for n from 0 to 1 do
 if id=27 then seqvars:=[p=t..infinity,q=v..infinity] else seqvars:=[p=t..infinity,q=0..v] fi:
 printf("START W%d q%d time%a\n",id,n,time()):
 W[id,n]:=FB(hyperInt(Bpure*Lpure^n/n!,seqvars),[t,v,z]):
 Emit(cat("W",id,"q",n),W[id,n]):save W,"weights.m":
 printf("DONE W%d q%d length%d time%a\n",id,n,length(W[id,n]),time()):
od od:
printf("WEIGHTS_COMPLETE\n"):quit: