restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f,vars) if type(f,list) then fibrationBasis(f,vars) else fibrationBasis(convert(f,HlogRegInf),vars) fi end proc:
Inv:=proc(f) local g;g:=FB(subs(z=1/zz,f),[zz]);subs(zz=z,g) end proc:
Emit:=proc(name,val) local fd;fd:=fopen(cat(name,".txt"),WRITE);fprintf(fd,"%s = %a\n",name,val);fclose(fd) end proc:
printf("BEGIN_T3_OPPOSITE\n"):
pole:=z/(1+z):
K0:=v->(ln(v)-ln(z)-ln(1-v))/(v-pole):
K1:=v->(ln(v)^2-(ln(z)+ln(1-v))^2)/(2*(v-pole)):
L1:=-ln(u)-ln(1-u):L2:=-ln(q)-ln(1-q):
V1:=2*ln(q)-ln(u)-ln(q-u):
V2:=2*ln(1-u)-ln(q-u)-ln(1-q):
dL1:=L1-V1:dH1:=L2-V2:
dL2:=(L1^2-V1^2)/2:dH2:=(L2^2-V2^2)/2:
kern1:=K0(u)*K0(q)*(dL1+dH1):
kern2:=K0(u)*K0(q)*(dL2+dH2)+(2*K1(u)*K0(q)+K0(u)*K1(q))*dL1+(2*K1(q)*K0(u)+K0(q)*K1(u))*dH1:
KC:=table():
for n from 1 to 2 do
 k:=eval(cat(kern,n)):
 firstu:=FB(hyperInt(k,[u=0..q]),[q,z]):
 Emit(cat("mid_inner",n),firstu):
 KC[n]:=FB(hyperInt(firstu,[q=0..1]),[z]):
 Emit(cat("Kmid",n),KC[n]):save KC,"mid.m":
 printf("MID %d length=%d time=%a\n",n,length(KC[n]),time()):
od:
Emit("E23q1",KC[1]):
Emit("E23q2",FB(KC[2]-3*ln(z)*KC[1]/2,[z])):
FTu:=z*(1-t)*(1-u)+t*u:FTq:=z*(1-t)*(1-q)+t*q:
P0:=(ln(q)-ln(FTq))/(q-pole):
P1:=(ln(q)^2-ln(FTq)^2)/(2*(q-pole)):
W1:=V1-2:W2:=V1^2/2-4+zeta[2]:
kern1:=(1+z)/FTu*P0*W1:
kern2:=(1+z)/FTu*(P0*(W2+2*ln(FTu)*W1)+P1*W1):
HC:=table():
for n from 1 to 2 do
 k:=eval(cat(kern,n)):
 firstt:=FB(hyperInt(k,[t=0..1]),[u,q,z]):
 Emit(cat("low_t",n),firstt):
 secondu:=FB(hyperInt(firstt,[u=0..q]),[q,z]):
 Emit(cat("low_u",n),secondu):
 HC[n]:=FB(hyperInt(secondu,[q=0..1]),[z]):
 Emit(cat("H22q",n),HC[n]):save HC,"low.m":
 printf("LOW %d length=%d time=%a\n",n,length(HC[n]),time()):
od:
EC1:=FB((HC[1]+Inv(HC[1]))/2,[z]):
EC2:=FB((HC[2]+Inv(HC[2]))/2+3*ln(z)*(Inv(HC[1])-HC[1])/4,[z]):
Emit("E22q1",EC1):Emit("E22q2",EC2):
printf("T3_OPPOSITE_COMPLETE\n"):
quit:
