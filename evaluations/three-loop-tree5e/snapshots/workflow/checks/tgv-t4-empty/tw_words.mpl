restart:interface(warnlevel=0):
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=1:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
read "D:/wilson-agent/workflow/checks/tgv-t4-empty/runs/20260918T021524768893Z-tw_vfirst/regularized_TW.m":
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\n",n,f);fclose(fd) end proc:
RP:=normal(Rg*eps^2):
PW:=table():RC:=table():
for av from 0 to 1 do for bv from 0 to 1 do
 lg:=convert(subs({aa=av,bb=bv},LL),HlogRegInf):PW[av,bv,0]:=[[1,[]]]:
 for k from 1 to 4 do PW[av,bv,k]:=scalarMul(shuffleSymbolic(PW[av,bv,k-1],lg),1/k) od:
 for j from 0 to 2 do RC[av,bv,j]:=normal(subs({aa=av,bb=bv},coeff(RP,eps,j))*bzero) od:
od od:
printf("WORDS_READY time%a\n",time()):
for id in [13,14,15] do
 if id=13 then choices:=[[1,1,0],[-1,1,1]]
 elif id=14 then choices:=[[1,0,0],[-1,1,0],[-1,0,1],[1,1,1]]
 else choices:=[[1,0,1],[-1,1,1]] fi:
 for q from -2 to 2 do
  if id=13 and q<1 then next fi:
  if id<>13 and q<2 then next fi:
  wl:=[]:
  for choice in choices do
   av:=choice[2]:bv:=choice[3]:
   for j from 0 to min(2,q+2) do
    wl:=[op(wl),op(scalarMul(PW[av,bv,q+2-j],choice[1]*RC[av,bv,j]))]:
   od:
  od:
  wl:=collectWords(wl):wl:=remove(vv->vv[1]=0,map(vv->[normal(vv[1]),vv[2]],wl)):
  printf("WORD_KERNEL %d q%d terms%d length%d time%a\n",id,q,nops(wl),length(wl),time()):save wl,cat("words",id,"q",q,".m"):
  ans:=hyperInt(wl,[v,t,u,r,s]):ans:=fibrationBasis(ans,[z]):Emit(cat("TW",id,"q",q),ans):
  printf("DONE WORD %d q%d length%d time%a\n",id,q,length(ans),time()):gc():
 od:
od:
printf("TW_WORDS_COMPLETE\n"):quit: