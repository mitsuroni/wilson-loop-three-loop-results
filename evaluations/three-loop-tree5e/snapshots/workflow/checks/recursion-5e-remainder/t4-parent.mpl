restart:interface(warnlevel=0):
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:_hyper_verbosity:=0:
FB:=proc(f,vs) if type(f,list) then fibrationBasis(f,vs) else fibrationBasis(convert(f,HlogRegInf),vs) fi end proc:
Inv:=proc(f) local g;g:=FB(subs(z=1/zz,f),[zz]);subs(zz=z,g) end proc:
Emit:=proc(name,val) local fd;fd:=fopen(cat(name,".txt"),WRITE);fprintf(fd,"%s = %a\n",name,val);fclose(fd) end proc:
Scale:=proc(wl,c) scalarMul(wl,c) end proc:
Join:=proc() local wl,i;wl:=[];for i to nargs do wl:=[op(wl),op(args[i])] od;collectWords(wl) end proc:
try
read "t4-geometry.mpl":read "t4-paths.mpl":read "t4-adj-paths.mpl":
printf("BEGIN_NEW_T4_PARENT\n"):
r4T:=1+r+s:r4X:=(1+t)*(1+u)*(1+v):
r4base:=(r*s)^(-eps)*r4T^(-1-eps)*r4X^(-1-2*eps)*F^(-2+2*eps):
r4pReg:=dimregPartial(NR*r4base,{r,1/u},eps):
r4pDiv:=findDivergences(r4pReg,{z},[eps=0]):if numelems([indices(r4pDiv)])>0 then error "parent R divergence remains" fi:
r4pPoly:=factor(normal(eps*r4pReg/r4base)):
r4pL:=-ln(r)-ln(s)-ln(r4T)-2*ln(1+t)-2*ln(1+u)-2*ln(1+v)+2*ln(F):
r4pZero:=1/(r4T*r4X*F^2):r4pR:=table():r4pH:=table():r4pY:=table():
r4pPW:=table():r4pPW[0]:=[[1,[]]]:r4pLW:=convert(r4pL,HlogRegInf):
for n from 1 to 3 do r4pPW[n]:=scalarMul(shuffleSymbolic(r4pPW[n-1],r4pLW),1/n) od:
for k from -1 to 2 do
 r4wl:=[]:for j from 0 to min(1,k+1) do r4wl:=[op(r4wl),op(scalarMul(r4pPW[k+1-j],normal(r4pZero*coeff(r4pPoly,eps,j))))] od:
 r4wl:=collectWords(r4wl):printf("PARENT_R_START %d words=%d time=%a\n",k,nops(r4wl),time()):
 r4pR[k]:=timelimit(720,hyperInt(r4wl,[v,t,u,r,s])):save r4pR,"parent-recoil.m":
 printf("PARENT_R_DONE %d time=%a\n",k,time()):gc():
od:
r4FH:=r*(1+u)+r*s*u*(1+v)+z*s*v:
r4pL:=-ln(r)-ln(s)-ln(1+r+s)-2*ln(1+u)-2*ln(1+v)+2*ln(r4FH):
for k from 0 to 2 do
 r4pH[k]:=timelimit(480,hyperInt(expand(r4pL^k/(k!*(1+r+s)*(1+u)*(1+v)*r4FH)),[r,u,s,v])):
 printf("PARENT_H_DONE %d time=%a\n",k,time()):
od:
r4Q0:=(1-b)*(1-x)*(1-w)+z*(b+(1-b)*x)*w:
r4Q1:=(1-x)*(1-(1-b)*w)+z*x*(1-b)*w:
for f from 0 to 1 do
 if f=0 then r4Q:=r4Q0 else r4Q:=r4Q1 fi:
 r4pL:=-ln(b)-ln(1-b)+2*ln(r4Q):
 for k from 0 to 3 do
  r4pY[f,k]:=timelimit(480,hyperInt(expand(r4pL^k/(k!*r4Q)),[x=0..1,b=0..1,w=0..1])):
  printf("PARENT_Y_DONE %d %d time=%a\n",f,k,time()):
 od:
od:
r4pPacket:=table():r4pProjected:=table():
for k from -2 to 2 do
 r4cus:=add((2*ln(z))^n/n!*FB(r4Fmom[0,k+2-n],[]),n=0..k+2)/2+FB(r4Fmom[0,k+2],[])/2:
 r4row:=convert(r4cus,HlogRegInf):
 if k>=0 then r4row:=Join(r4row,Scale(r4pH[k],-z));r4h:=FB(r4pH[k],[z]);r4hm:=Inv(r4h);
  r4wface:=add((2*ln(z))^n/n!*Inv(FB(r4pH[k-n],[z])),n=0..k)/z:
  r4row:=Join(r4row,convert(-r4wface,HlogRegInf)) fi:
 if k>=-1 then r4row:=Join(r4row,Scale(Join(r4pY[1,k+1],Scale(r4pY[0,k+1],-1)),1+z),Scale(r4pR[k],4)) fi:
 if k>=0 then r4row:=Join(r4row,Scale(r4pR[k-1],-8)) fi:
 r4pPacket[k]:=FB(r4row,[z]):Emit(cat("ParentPacket",k),r4pPacket[k]):
 r4pProjected[k]:=(r4pPacket[k]+Inv(r4pPacket[k]))/2:
 if k>=-1 then r4pProjected[k]:=r4pProjected[k]+ln(z)*(Inv(r4pPacket[k-1])-r4pPacket[k-1])/2 fi:
 if k>=1 then r4pProjected[k]:=r4pProjected[k]-ln(z)^3*(Inv(r4pPacket[k-3])-r4pPacket[k-3])/6 fi:
od:
r4g2:=[1,0,2*zeta[2],8*zeta[3]/3,18*zeta[2]^2/5]:
for k from -2 to 2 do
 r4pVal:=FB(add(r4g2[j+1]*r4pProjected[k-j],j=0..k+2)/8,[z]):
 Emit(cat("ParentI",k),r4pVal):
od:
save r4pPacket,r4pProjected,"parent-assembled.m":
printf("COMPLETE_NEW_T4_PARENT\n"):
catch:
printf("FAILED_NEW_T4_PARENT %a\n",lasterror):
end try:
quit:
