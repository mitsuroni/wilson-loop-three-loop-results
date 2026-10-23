restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:
FB:=proc(f) fibrationBasis(convert(subs(Pi^6=216*zeta[2]^3,Pi^4=36*zeta[2]^2,Pi^2=6*zeta[2],f),HlogRegInf),[z]) end proc:
read "D:/wilson-agent/workflow/checks/w3-ca-nf/source_inputs.mpl":
read "D:/wilson-agent/workflow/checks/w3-ca-nf/physical.m":
CompactShape:=-M13/36-N12/24+11*M21/72-11*M01^2/432+289*M11/216
 +(2051/1296-zeta[2]/8-zeta[3]/2)*M01
 +412765/104976-56071*zeta[2]/11664-76*zeta[3]/81
 +2239*zeta[2]^2/2160+7*zeta[2]*zeta[3]/18+2*zeta[5]:
CompactFinite:=11*L^4/108+(289/324-zeta[2]/9)*L^3
 +(2051/648-43*zeta[2]/18-11*rho^2/72)*L^2
 +(11*M11/18+(zeta[2]/12-289/432)*rho^2+31313/5832
 -2807*zeta[2]/324+17*zeta[2]^2/10-113*zeta[3]/54)*L+CompactShape:
Rules:={M01=Mom[0,1],M11=Mom[1,1],M21=Mom[2,1],M13=Mom[1,3],N12=Nom[1,2],rho=ln(z)}:
if FB(subs(Rules,CompactShape)-Shape)<>0 then error "compact shape differs" fi:
if FB(subs(Rules,CompactFinite)-Finite)<>0 then error "compact finite differs" fi:
fd:=fopen("compact-coefficients.txt",WRITE,TEXT):
fprintf(fd,"shape\t0\t%a\n",CompactShape):fprintf(fd,"finite\t0\t%a\n",CompactFinite):
for q from -6 to -1 do fprintf(fd,"pole\t%d\t%a\n",q,Post[q]) od:
fprintf(fd,"k3\t0\t%a\n",k3):fprintf(fd,"h3\t0\t%a\n",h3):
fclose(fd):printf("COMPACT_MN_SHAPE_AND_FINITE_EXACT_CHECKS_PASS=2\n"):
quit:
