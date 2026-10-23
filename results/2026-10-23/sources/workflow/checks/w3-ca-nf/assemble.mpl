restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
interface(warnlevel=0):
FB:=proc(f) fibrationBasis(convert(subs(Pi^6=216*zeta[2]^3,Pi^4=36*zeta[2]^2,Pi^2=6*zeta[2],f),HlogRegInf),[z]) end proc:
read "D:/wilson-agent/workflow/checks/w3-ca-nf/source_inputs.mpl":
partnames:=["Pi1","chain","Pi2","legs","Gamma"]:
Raw:=table():Bare:=table():Post:=table():Checks:=[]:
fd:=fopen("assembly-coefficients.txt",WRITE,TEXT):
for q from -6 to 0 do
 for p in partnames do
  parts[p,q]:=FB(subs(r=1/(1+z),parts[p,q])):
  fprintf(fd,"part_%s\t%d\t%a\n",p,q,parts[p,q]):
 od:
 Raw[q]:=FB(add(parts[p,q],p in partnames)):
 fprintf(fd,"raw\t%d\t%a\n",q,Raw[q]):
 printf("RAW_CHANNEL eps%d length%d\n",q,length(Raw[q])):
od:
Scale:=proc(n,j) expand(n^j*((L+ln(z)/2)^j+(L-ln(z)/2)^j)/j!) end proc:
for q from -6 to 0 do
 Bare[q]:=FB(add(Raw[j]*Scale(3,q-j),j=-6..q)):
 Post[q]:=FB(Bare[q]+subs(rho=ln(z),CT[q])):
 fprintf(fd,"bare\t%d\t%a\n",q,Bare[q]):
 fprintf(fd,"post\t%d\t%a\n",q,Post[q]):
 printf("POST_COUNTERTERM eps%d length%d\n",q,length(Post[q])):
od:
Check:=proc(name,expr) local v;global Checks,fd;
 v:=FB(expr):Checks:=[op(Checks),[name,evalb(v=0),v]]:
 fprintf(fd,"check_%s\t0\t%a\n",name,v):
 printf("CHECK %s %a residual_length%d\n",name,evalb(v=0),length(v)):
 return evalb(v=0)
end proc:
Check("superleading_minus6",Post[-6]):
Check("superleading_minus5",Post[-5]):
k3:=FB(-3/2*diff(Post[-1],L)):
h3:=FB(-3*subs(L=0,Post[-1])):
k3ref:=-209/216+5*zeta[2]/9-7*zeta[3]/6:
fprintf(fd,"k3\t0\t%a\n",k3):fprintf(fd,"h3\t0\t%a\n",h3):
Check("cusp_geometry_independence",diff(k3,z)):
Check("cusp_scale_independence",diff(k3,L)):
Check("noncusp_geometry_independence",diff(h3,z)):
Check("cusp_GHKM_eq134",k3-k3ref):
b0A:=11/12:b0f:=-1/3:b1Af:=-5/12:
CrossB0Squared:=2*b0A*b0f:
CrossB0k2:=b0A*k2f+b0f*k2A:
CrossB0h2:=b0A*h2f+b0f*h2A:
Expected:=table():
Expected[-4]:=-11/9*CrossB0Squared:
Expected[-3]:=-2/3*CrossB0Squared*L+5/9*CrossB0k2+8/9*b1Af:
Expected[-2]:=2/3*(CrossB0k2+b1Af)*L+CrossB0h2/3-2*k3/9:
Expected[-1]:=-(2*k3*L+h3)/3:
for q from -4 to -1 do Check(cat("pole_template_",q),Post[q]-Expected[q]) od:
Finite:=Post[0]:Shape:=subs(L=0,Finite):
LowerRG:=subs(rho=ln(z),2*b0A*F2f+2*b0f*F2A+b1Af*F1):
Df:=FB(diff(Finite,L)-LowerRG):
Check("finite_RG",Df+2*k3*L+h3):
Check("finite_cusp_extraction",-diff(Df,L)/2-k3):
Check("finite_noncusp_extraction",-subs(L=0,Df)-h3):
for q from -6 to 0 do Check(cat("crossing_",q),Post[q]-subs(z=1/z,Post[q])) od:
fprintf(fd,"finite\t0\t%a\n",Finite):fprintf(fd,"shape\t0\t%a\n",Shape):
fprintf(fd,"lower_RG\t0\t%a\n",FB(LowerRG)):
for controlLabel in ["minus6","minus5","minus4"] do
 q:=`if`(controlLabel="minus6",-6,`if`(controlLabel="minus5",-5,-4)):
 fprintf(fd,"omit_Gamma\t%d\t%a\n",q,FB(Post[q]-add(parts["Gamma",j]*Scale(3,q-j),j=-6..q))):
od:
fclose(fd):
save Raw,Bare,Post,k3,h3,Finite,Shape,Checks,parts,"assembly.m":
printf("CHANNEL_ASSEMBLY_COMPLETED checks%d passed%d\n",nops(Checks),nops(select(x->x[2],Checks))):
quit:
