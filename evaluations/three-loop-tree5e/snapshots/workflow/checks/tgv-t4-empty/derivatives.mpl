restart:interface(warnlevel=0):
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f) if type(f,list) then fibrationBasis(f,[z]) else fibrationBasis(convert(f,HlogRegInf),[z]) fi end proc:
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\n",n,f);fclose(fd) end proc:
Euler:=proc(f,sub) local var,ret;ret:=0;for var in sub do if type(var,name) then ret:=ret+var*diff(f,var) else ret:=ret-op(1,var)*diff(f,op(1,var)) fi od;ret end proc:
for id in [5,1,2,3,4] do
 read cat("D:/wilson-agent/workflow/checks/tgv-t4-empty/runs/20260918T014636406615Z-ibp_scan/deriv",id,".m"):
 bzero:=normal(eval(fw,eps=0)):
 Lf:=simplify(eval(diff(fw,eps)/fw,eps=0),ln):
 if id=5 then subslist:=[{t}]:sdds:=[eps]:nn:=1
 elif id<3 then subslist:=[{1/u},{r,1/u}]:sdds:=[eps,2*eps]:nn:=2
 else subslist:=[{v},{r,v},{v,1/t}]:sdds:=[eps,2*eps,3*eps]:nn:=3 fi:
 Rg:=1:
 for i to nn do
  ss:=subslist[i]:dd:=numelems(select(type,ss,name))-numelems(remove(type,ss,name)):
  elog:=normal(Euler(bzero,ss)/bzero+eps*Euler(Lf,ss)):
  Rg:=factor(normal(((sdds[i]-dd)*Rg-Euler(Rg,ss)-Rg*elog)/sdds[i])):
 od:
 printf("REG %d length%d time%a\n",id,length(Rg),time()):
 dnew:=findDivergences(Rg*fw,{z},[eps=0]):printf("remaining %a\n",eval(dnew)):
 if numelems([indices(dnew)])>0 then error "remaining divergences" fi:
 save Rg,bzero,Lf,fw,dnew,cat("regularized",id,".m"):
 Rpoly:=normal(Rg*eps^nn):
 for q from -nn to 1 do
  kern:=add(coeff(Rpoly,eps,j)*Lf^(q+nn-j)/(q+nn-j)!,j=0..q+nn)*bzero:
  kern:=collect(expand(kern),indets(kern,specfunc(anything,ln)),normal):
  save kern,cat("inputD",id,"q",q,".m"):
  printf("START D%d q%d length%d time%a\n",id,q,length(kern),time()):
  ans:=timelimit(420,hyperInt(kern,[r,u,s,t,v])):ans:=FB(ans):
  Emit(cat("D",id,"q",q),ans):
  printf("DONE D%d q%d length%d time%a\n",id,q,length(ans),time()):gc():
 od:
od:
printf("DERIV_COMPLETE\n"):quit: