restart:
_hyper_autoload_periods := ["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_check_divergences := true:
_hyper_abort_on_divergence := true:
_hyper_verbosity := 0:
for q from 0 to 5 do
    kernel := expand(x/(1-x*u)*(-log(x)-log(1-x)+2*log(1-u))^q/q!):
    alternate := simplify(convert(fibrationBasis(hyperInt(kernel,[x=0..1,u=0..1]),[]),Hlog)):
    fd := fopen(cat("alternateF",q,".txt"),WRITE):fprintf(fd,"%a\n",alternate):fclose(fd):
    printf("ALTERNATE_F %d %a\n",q,alternate):
    barekernel := expand(x/(1-x*u)*(-log(x)-log(1-x)+log(1-u))^q/q!):
    bare := simplify(convert(fibrationBasis(hyperInt(barekernel,[u=0..1,x=0..1]),[]),Hlog)):
    fd := fopen(cat("bareF",q,".txt"),WRITE):fprintf(fd,"%a\n",bare):fclose(fd):
od:
printf("CONTROL_COMPLETE\n"):
quit:
