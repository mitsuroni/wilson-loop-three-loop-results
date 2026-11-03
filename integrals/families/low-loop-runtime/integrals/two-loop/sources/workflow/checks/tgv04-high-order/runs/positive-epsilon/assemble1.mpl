restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
FB:=proc(f) fibrationBasis(convert(f,HlogRegInf),[z]) end proc:
Inv:=proc(f) local w;w:=fibrationBasis(convert(subs(z=1/zz,f),HlogRegInf),[zz]);subs(zz=z,w) end proc:
read "step03_easy_results.m":read "step05_H0.m":read "step10_optimized_results.m":
read "step11_hlog_results.m":
read "D2.m":read "H1.m":read "O1.m":
ell:=Hlog(z,[0]):
Bm2:=zeta[2]:
Bm1:=FB(zeta[2]*ell+zeta[3]-(1+z)*D0z):
B0:=FB(zeta[2]*ell^2+zeta[3]*ell+4/5*zeta[2]^2-(1+z)*D1z-z*H0-Inv(H0)/z):
B1:=FB(2/3*zeta[2]*ell^3+zeta[3]*ell^2+4/5*zeta[2]^2*ell+3*zeta[5]-zeta[2]*zeta[3]-(1+z)*D2-z*H1-(Inv(H1)+2*ell*Inv(H0))/z):
J1:=FB((B1+2*B0+4*Bm1+8*Bm2)/4+O1):
normSeries:=convert(series((1-2*eps)*exp(2*zeta[2]*eps^2+8/3*zeta[3]*eps^3+4*(2/5*zeta[2]^2)*eps^4)/(1+exp(2*eps*ell)),eps=0,5),polynom):
I1:=FB(J1/2+coeff(normSeries,eps,1)*J0+coeff(normSeries,eps,2)*Jm1+coeff(normSeries,eps,3)*Jm2):
cm1:=FB(Inv(J1)-J1+2*ell*J0-2*ell^2*Jm1+4/3*ell^3*Jm2):
ci1:=FB(Inv(I1)-I1):
if cm1<>0 or ci1<>0 then error "first-order crossing failed" fi:
save J1,I1,B1,normSeries,"assembled1.m":
fd:=fopen("I1.txt",WRITE):fprintf(fd,"I1 = %a\n",I1):fclose(fd):
printf("ASSEMBLY1 COMPLETE\n"):
quit:
