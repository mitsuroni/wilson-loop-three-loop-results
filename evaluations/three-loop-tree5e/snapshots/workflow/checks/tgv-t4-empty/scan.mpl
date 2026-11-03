restart:
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
T:=1+r+s:X:=(1+t)*(1+u)*(1+v):
SS:=a1*a2*y*(1-x)+a1*a3*((1-x)*(1-w)+z*x*w)+z*a2*a3*w*(1-y):
PP:=(a1*a2*(-(1+z)*(1-x)+z*y)+a1*a3*(1+3*z-(1+z)*x+z*(1+z)*w)+z*a2*a3*(1-y-(1+z)*w))/4:
mapping:={a1=r/T,a2=s/T,a3=1/T,x=t/(1+t),y=u/(1+u),w=v/(1+v)}:
F:=factor(normal(subs(mapping,SS)*T^2*X)):
N:=factor(normal(subs(mapping,PP)*T^2*X)):
base:=(r*s)^(-eps)*T^(-1-eps)*X^(-1-2*eps)*F^(-2+2*eps):
f:=N*base:
printf("F=%a\nN=%a\n",F,N):
dv:=findDivergences(f,{z},[eps=0]):
printf("RAW divergences=%a\n",eval(dv)):
for sub in [indices(dv,'nolist')] do
 printf("SUB %a SDD %a\n",sub,dv[sub]):
od:
save F,N,base,f,dv,"scan.m":
printf("SCAN_COMPLETE time=%a\n",time()):quit: