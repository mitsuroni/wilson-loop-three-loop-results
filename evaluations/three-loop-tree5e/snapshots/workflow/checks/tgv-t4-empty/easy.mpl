restart:interface(warnlevel=0):
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f,vs) if type(f,list) then fibrationBasis(f,vs) else fibrationBasis(convert(f,HlogRegInf),vs) fi end proc:
Emit:=proc(n,f) local fd;fd:=fopen(cat(n,".txt"),WRITE);fprintf(fd,"%s = %a\n",n,f);fclose(fd) end proc:
# The two cusp integrals after the exact alpha-scale and x/w integral.
for typ in ["F","G"] do
 for aa from 0 to 1 do
  if typ="G" and aa=0 then next fi:
  if typ="F" then LL:=-ln(b)-ln(1-b)+aa*ln(u)+2*ln(1-b*u)
  else LL:=-ln(b)-ln(1-b)+aa*ln(1-u)+2*ln(1-b*u) fi:
  for q from 0 to 4 do
   ans:=FB(hyperInt(LL^q/(q!*(1-b*u)),[u=0..1,b=0..1]),[]):
   Emit(cat(typ,aa,"q",q),ans):printf("CUSP %s %d q%d time%a\n",typ,aa,q,time()):
  od:
 od:
od:
# y=0 and y=1: alpha scale is 1/eps, remaining integrals are finite.
for face from 0 to 1 do
 if face=0 then Q:=(1-b)*(1-x)*(1-w)+z*(b+(1-b)*x)*w
 else Q:=(1-x)*(1-(1-b)*w)+z*x*(1-b)*w fi:
 Lb:=-ln(b)-ln(1-b)+2*ln(Q):
 for id in [13,14,15,40] do
  if (face=0 and id in {14,15}) or (face=1 and id=13) then next fi:
  for q from 0 to 3 do
   if id=40 then pol:=(Lb^q-(Lb+ln(x))^q)/q!
   elif id in {13,14} then pol:=(Lb^q-(Lb+ln(w))^q)/q!
   else pol:=(Lb+ln(w))^q/q! fi:
   printf("START Y%d id%d q%d time%a\n",face,id,q,time()):
   ans:=timelimit(300,hyperInt(pol/Q,[x=0..1,w=0..1,b=0..1])):
   ans:=FB(ans,[z]):Emit(cat("Y",face,"i",id,"q",q),ans):
   printf("DONE Y%d id%d q%d length%d time%a\n",face,id,q,length(ans),time()):gc():
  od:
 od:
od:
printf("EASY_COMPLETE\n"):quit: