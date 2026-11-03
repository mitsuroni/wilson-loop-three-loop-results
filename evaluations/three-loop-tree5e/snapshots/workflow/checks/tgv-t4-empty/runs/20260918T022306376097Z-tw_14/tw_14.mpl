restart:interface(warnlevel=0):
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=1:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
read "D:/wilson-agent/workflow/checks/tgv-t4-empty/runs/20260918T014450574703Z-scan/scan.m":
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\n",n,f);fclose(fd) end proc:
FB:=proc(f) if type(f,list) then fibrationBasis(f,[z]) else fibrationBasis(convert(f,HlogRegInf),[z]) fi end proc:
Euler:=proc(f,sub) local var,ret;ret:=0;for var in sub do if type(var,name) then ret:=ret+var*diff(f,var) else ret:=ret-op(1,var)*diff(f,op(1,var)) fi od;ret end proc:
Nw:=(1+v)*(r*(z*t-1)*(1+u)+z*s*(1+t)):
# aa,bb are 0 or 1; both Euler degrees are changed by the w^eps weight.
LL:=-ln(r)-ln(s)-ln(1+r+s)-2*ln(1+t)-2*ln(1+u)-2*ln(1+v)+2*ln(F)-aa*ln(1+u)+bb*(ln(v)-ln(1+v)):
bzero:=normal(eval(base,eps=0)):
Rg:=Nw:
for pair in [[{r,v},(1+bb)*eps],[{v,1/t},(2+bb)*eps]] do
 ss:=pair[1]:sd:=pair[2]:dd:=numelems(select(type,ss,name))-numelems(remove(type,ss,name)):
 el:=normal(Euler(bzero,ss)/bzero+eps*Euler(LL,ss)):
 Rg:=factor(normal(((sd-dd)*Rg-Euler(Rg,ss)-Rg*el)/sd)):
od:
printf("REG_TW length%d time%a\n",length(Rg),time()):
save Rg,Nw,LL,bzero,base,"regularized_TW.m":
for av in [0,1] do for bv in [0,1] do
 ff:=base*(1+u)^(-av*eps)*v^(bv*eps)*(1+v)^(-bv*eps)*subs({aa=av,bb=bv},Rg):
 dv:=findDivergences(ff,{z},[eps=0]):printf("REMAIN %d %d %a\n",av,bv,eval(dv)):
 if numelems([indices(dv)])>0 then error "remaining TW divergences" fi:
od od:
Rp:=normal(Rg*eps^2):
for id in [14] do
 for q from -2 to 2 do
  mon:=bzero*add(coeff(Rp,eps,j)*LL^(q+2-j)/(q+2-j)!,j=0..q+2):
  if id=13 then kern:=subs({aa=1,bb=0},mon)-subs({aa=1,bb=1},mon)
  elif id=14 then kern:=subs({aa=0,bb=0},mon)-subs({aa=1,bb=0},mon)-subs({aa=0,bb=1},mon)+subs({aa=1,bb=1},mon)
  else kern:=subs({aa=0,bb=1},mon)-subs({aa=1,bb=1},mon) fi:
  kern:=collect(expand(kern),indets(kern,specfunc(anything,ln)),normal):
  printf("START TW%d q%d length%d time%a\n",id,q,length(kern),time()):save kern,cat("inputTW",id,"q",q,".m"):
  ans:=timelimit(360,hyperInt(kern,[v,t,u,r,s])):ans:=FB(ans):Emit(cat("TW",id,"q",q),ans):
  printf("DONE TW%d q%d length%d time%a\n",id,q,length(ans),time()):gc():
 od:
od:
printf("TW_COMPLETE\n"):quit: