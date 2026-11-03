restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f) if type(f,list) then fibrationBasis(f,[z]) else fibrationBasis(convert(f,HlogRegInf),[z]) fi end proc:
Inv:=proc(f) local g;g:=fibrationBasis(convert(subs(z=1/zz,f),HlogRegInf),[zz]);subs(zz=z,g) end proc:
Emit:=proc(name,val) local fd;fd:=fopen(cat(name,".txt"),WRITE);fprintf(fd,"%s = %a\n",name,val);fclose(fd) end proc:
printf("BEGIN_T1_OPPOSITE\n"):
Lb:=-ln(b)-ln(1-b)+2*ln(1-b+b*w):
Hc:=table():
for q from 0 to 3 do
  Hc[q]:=fibrationBasis(hyperInt((1-w)/(1-b+b*w)*Lb^q/q!,[b=0..1]),[w]):
  Emit(cat("H",q),Hc[q]):
  save Hc,"Hcoeff.m":
  printf("H %d length=%d time=%a\n",q,length(Hc[q]),time()):
od:
Tc:=table():Ac:=table():Kc:=table():Tc[0]:=0:Kc[0]:=0:
Bcoef:=[1,2,4-zeta[2],8-2*zeta[2]-2*zeta[3]]:
pole:=1/(1+z):
for q from 1 to 4 do
  kernelT:=-add((2*ln(w))^j/j!*((ln(z)+ln(w))^(q+1-j)-ln(1-w)^(q+1-j))/(q+1-j)!,j=1..q)/(w-pole):
  Tc[q]:=FB(hyperInt(kernelT,[w=0..1])):
  Emit(cat("T",q),Tc[q]):
  printf("T %d length=%d time=%a\n",q,length(Tc[q]),time()):
  kernelA:=add(((ln(z)+ln(w))^(j+1)-ln(1-w)^(j+1))/(j+1)!*Hc[q-1-j],j=0..q-1)/(w-pole):
  Ac[q-1]:=FB(hyperInt(kernelA,[w=0..1])):
  Emit(cat("A",q-1),Ac[q-1]):
  Kc[q]:=FB(2*Ac[q-1]-add(Bcoef[j+1]*Tc[q-j],j=0..q-1)):
  Emit(cat("K",q),Kc[q]):
  save Hc,Tc,Ac,Kc,"moments.m":
  printf("K %d length=%d time=%a\n",q,length(Kc[q]),time()):
od:
if Kc[1]<>0 then error "K linear coefficient must cancel" fi:
for q from 1 to 4 do
  ET:=FB(add((ln(z)/2)^j/j!*Tc[q-j]+(-ln(z)/2)^j/j!*Inv(Tc[q-j]),j=0..q-1)/2):
  EK:=FB(add((ln(z)/2)^j/j!*Kc[q-j]+(-ln(z)/2)^j/j!*Inv(Kc[q-j]),j=0..q-1)/2):
  Emit(cat("ET",q),ET):Emit(cat("EK",q),EK):
  printf("EVEN %d time=%a\n",q,time()):
od:
printf("T1_OPPOSITE_COMPLETE\n"):
quit:
