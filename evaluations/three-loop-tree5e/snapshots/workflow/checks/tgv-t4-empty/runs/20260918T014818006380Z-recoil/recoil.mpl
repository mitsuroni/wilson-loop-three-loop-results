restart:interface(warnlevel=0):
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
read "D:/wilson-agent/workflow/checks/tgv-t4-empty/runs/20260918T014450574703Z-scan/scan.m":
FB:=proc(f) if type(f,list) then fibrationBasis(f,[z]) else fibrationBasis(convert(f,HlogRegInf),[z]) fi end proc:
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\n",n,f);fclose(fd) end proc:
Eu:=proc(f) r*diff(f,r)-u*diff(f,u) end proc:
RR:=(1+z)/2*(1+u)*(r*(1+v)+z*(r-s)*v*(1+t)):
ell:=-eps+(-1-eps)*r/(1+r+s)+(1+2*eps)*u/(1+u)+(-2+2*eps)*Eu(F)/F:
RegA:=factor(normal(((1+aa)*eps*RR-Eu(RR)-RR*(ell+aa*eps*u/(1+u)))/((1+aa)*eps))):
rm1:=normal(eval(eps*RegA,eps=0)):rr0:=normal(eval(diff(eps*RegA,eps),eps=0)):
bzero:=normal(eval(base,eps=0)):
Lb:=-ln(r)-ln(s)-ln(1+r+s)-2*ln(1+t)-2*ln(1+u)-2*ln(1+v)+2*ln(F):
LF:=Lb-aa*ln(1+u)+bb*(ln(v)-ln(1+v))+cc*(ln(t)-ln(1+t)):
Mon:=proc(q) if q=-1 then bzero*rm1 else bzero*(rm1*LF^(q+1)/(q+1)!+rr0*LF^q/q!) fi end proc:
Choice:=proc(q,id) local m; m:=Mon(q);
 if id=13 then subs({aa=1,bb=0,cc=0},m)-subs({aa=1,bb=1,cc=0},m)
 elif id=14 then subs({aa=0,bb=0,cc=0},m)-subs({aa=1,bb=0,cc=0},m)-subs({aa=0,bb=1,cc=0},m)+subs({aa=1,bb=1,cc=0},m)
 elif id=15 then subs({aa=0,bb=1,cc=0},m)-subs({aa=1,bb=1,cc=0},m)
 else subs({aa=0,bb=0,cc=0},m)-subs({aa=0,bb=0,cc=1},m) fi end proc:
for id in [40,13,14,15] do
 for q from -1 to 2 do
  kern:=collect(expand(Choice(q,id)),indets(Choice(q,id),specfunc(anything,ln)),normal):
  save kern,cat("input",id,"q",q,".m"):
  printf("START %d q%d length%d time%a\n",id,q,length(kern),time()):
  ans:=timelimit(420,hyperInt(kern,[r,u,s,t,v])):
  ans:=FB(ans):Emit(cat("R",id,"q",q),ans):
  printf("DONE %d q%d length%d time%a\n",id,q,length(ans),time()):gc():
 od:
od:
printf("RECOIL_COMPLETE\n"):quit: