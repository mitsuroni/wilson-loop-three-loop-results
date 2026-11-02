restart:
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity := 0:
# Positive-domain simplex chart alpha1=r/(1+r+s), alpha2=s/(1+r+s), alpha3=1/(1+r+s).
# Hard face: x2=u/(1+u), x3=v/(1+v).
F := r*(1+u)+r*s*u*(1+v)+z*s*v:
polys := {r,s,u,v,1+r+s,1+u,1+v,F}:
L := table():
L[{}] := [polys,{op(combinat[choose](polys,2))}]:
printf("START hard-face cgReduction\n");
cgReduction(L,{z}):
ord := suggestIntegrationOrder(L):
printf("hard-face order=%a\n",ord);
printf("hard-face complete=%a\n",assigned(L[{r,s,u,v}]));
save L,ord,"step04_hard_reduction.m":
f := (r*s)^(-eps)*(1+r+s)^(-1-eps)*(1+u)^(-1-2*eps)*(1+v)^(-1-2*eps)*F^(-1+2*eps):
dv := findDivergences(f,{z},[eps=0]):
printf("hard-face divergences=%a\n",eval(dv));
save f,dv,"step04_hard_integrand.m":
quit:
