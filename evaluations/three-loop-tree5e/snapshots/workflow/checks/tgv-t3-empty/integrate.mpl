restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f,vars) if type(f,list) then fibrationBasis(f,vars) else fibrationBasis(convert(f,HlogRegInf),vars) fi end proc:
Inv:=proc(f) local g;g:=FB(subs(z=1/zz,f),[zz]);subs(zz=z,g) end proc:
Emit:=proc(name,val) local fd;fd:=fopen(cat(name,".txt"),WRITE);fprintf(fd,"%s = %a\n",name,val);fclose(fd) end proc:
printf("BEGIN_T3_EMPTY\n"):
Ac:=table():Rc:=table():Lc:=table():Hc:=table():
Bcoef:=[1,2,4-zeta[2],8-2*zeta[2]-2*zeta[3],16-4*zeta[2]-4*zeta[3]-zeta[2]^2/10]:
for q from 0 to 3 do
 Ac[q]:=FB(hyperInt((ln(u)-2*ln(v)-ln(1-v))^q/(v*q!),[v=u..1]),[u]):
 Rc[q]:=FB(hyperInt(((1-2*u*v)/(1-u*v)*(ln(1-v)-2*ln(1-u*v)-ln(v))^q-(-ln(v))^q)/(v*q!),[v=0..1]),[u]):
 Emit(cat("A",q),Ac[q]):Emit(cat("R",q),Rc[q]):
 printf("INNER %d time=%a\n",q,time()):
od:
for q from 1 to 4 do
 Lc[9,q]:=FB(add(Bcoef[j+1]*ln(u)^(q-j)/(q-j)!,j=0..q-1)+Ac[q-1],[u]):
 Lc[10,q]:=FB((-ln(u)-ln(1-u))^q/q!-Ac[q-1]+Rc[q-1],[u]):
 Emit(cat("L9q",q),Lc[9,q]):Emit(cat("L10q",q),Lc[10,q]):
od:
if Lc[9,1]<>0 or Lc[10,1]<>0 then error "linear cancellation failed" fi:
for id in [9,10] do
 Hc[id,0]:=0:Hc[id,1]:=0:
 for q from 2 to 4 do
  kern:=add((ln(u)^(j+1)-(ln(z)+ln(1-u))^(j+1))*2^j/(j+1)!*Lc[id,q-j],j=0..q-2)/(u-z/(1+z)):
  Hc[id,q]:=FB(hyperInt(kern,[u=0..1]),[z]):
  Emit(cat("H",id,"q",q),Hc[id,q]):save Ac,Rc,Lc,Hc,"moments.m":
  printf("H %d %d length=%d time=%a\n",id,q,length(Hc[id,q]),time()):
 od:
od:
Fpol:=z*(1-t)*(1-u)+t*u:
for q from 2 to 4 do
 kern:=(1+z)/Fpol*add(add(-ln(t)^j/j!*((-ln(u)-ln(1-u))^k/k!-Bcoef[k+1])*(2*ln(Fpol))^(q-j-k)/(q-j-k)!,k=1..q-j),j=1..q-1):
 inner26:=FB(hyperInt(kern,[t=0..1]),[u,z]):
 Emit(cat("inner2626q",q),inner26):
 Hc[26,q]:=FB(hyperInt(inner26,[u=0..1]),[z]):
 Emit(cat("H26q",q),Hc[26,q]):save Ac,Rc,Lc,Hc,"moments.m":
 printf("H 26 %d length=%d time=%a\n",q,length(Hc[26,q]),time()):
od:
for id in [9,10,26] do
 for q from 2 to 4 do
  EC:=FB(add((-3*ln(z)/2)^j/j!*Hc[id,q-j]+(3*ln(z)/2)^j/j!*Inv(Hc[id,q-j]),j=0..q-2)/2,[z]):
  Emit(cat("E",id,"q",q),EC):
  printf("EVEN %d %d time=%a\n",id,q,time()):
 od:
od:
printf("T3_EMPTY_COMPLETE\n"):
quit:
