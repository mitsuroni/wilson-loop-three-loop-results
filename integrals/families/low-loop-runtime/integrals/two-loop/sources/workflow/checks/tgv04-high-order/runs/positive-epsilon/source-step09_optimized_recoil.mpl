restart:
interface(warnlevel=0):
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity := 0:
read "step06_recoil_scan.m":
SS := a1*a2*x2*(1-x1)+a1*a3*((1-x1)*(1-x3)+z*x1*x3)+z*a2*a3*x3*(1-x2):
PP := (a1*a2*(-(1+z)*(1-x1)+z*x2)+a1*a3*(1+3*z-(1+z)*x1+z*(1+z)*x3)+z*a2*a3*(1-x2-(1+z)*x3))/4:
PRopt := (1+z)/2*a3*(a1*(1-x1+z*x3)-z*a2*x3):
cert := expand(PP-(-z*diff(SS,x1)-(1+z)*diff(SS,x2)+diff(SS,x3))/4-PRopt):
printf("optimized P certificate residual=%a\n",cert);
if cert<>0 then error "invalid optimized numerator" fi:
RR := (1+z)/2*(1+u)*(r*(1+v)+z*(r-s)*v*(1+t)):
fo := RR*base:
dopt := findDivergences(fo,{z},[eps=0]):
printf("optimized divergences=%a\n",eval(dopt));
save fo,RR,base,F,dopt,"step09_optimized_input.m":
quit:
