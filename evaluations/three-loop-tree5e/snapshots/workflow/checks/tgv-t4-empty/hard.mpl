restart:interface(warnlevel=0):
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f) if type(f,list) then fibrationBasis(f,[z]) else fibrationBasis(convert(f,HlogRegInf),[z]) fi end proc:
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\n",n,f);fclose(fd) end proc:
F:=r*(1+u)+r*s*u*(1+v)+z*s*v:
bzero:=1/((1+r+s)*(1+u)*(1+v)*F):
Lb:=-ln(r)-ln(s)-ln(1+r+s)-2*ln(1+u)-2*ln(1+v)+2*ln(F):
lu:=-ln(1+u):lv:=ln(v)-ln(1+v):
for id in [40,13,14,15,115,140] do
 for q from 0 to 2 do
  if id=40 then pol:=Lb^q
  elif id=13 then pol:=(Lb+lu)^q-(Lb+lu+lv)^q
  elif id=14 then pol:=Lb^q-(Lb+lu)^q-(Lb+lv)^q+(Lb+lu+lv)^q
  elif id=15 then pol:=(Lb+lv)^q-(Lb+lu+lv)^q
  elif id=115 then pol:=Lb^q-(Lb+ln(u)-ln(1+u))^q
  else pol:=Lb^q-(Lb-ln(1+v))^q fi:
  kern:=bzero*expand(pol)/q!:
  printf("START H%d q%d length%d time%a\n",id,q,length(kern),time()):
  ans:=timelimit(420,hyperInt(kern,[r,u,s,v])):ans:=FB(ans):
  Emit(cat("H",id,"q",q),ans):printf("DONE H%d q%d length%d time%a\n",id,q,length(ans),time()):gc():
 od:
od:
printf("HARD_COMPLETE\n"):quit: