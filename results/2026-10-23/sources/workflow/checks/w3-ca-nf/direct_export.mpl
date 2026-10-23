restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:
FB:=proc(f) fibrationBasis(convert(subs(Pi^6=216*zeta[2]^3,Pi^4=36*zeta[2]^2,Pi^2=6*zeta[2],f),HlogRegInf),[z]) end proc:
read "D:/wilson-agent/workflow/checks/w3-ca-nf/source_inputs.mpl":

fd:=fopen("direct-export-checks.txt",WRITE,TEXT):
direct:=-1/8*Hlog(z,[-1, -1, 0, 0])-1/48*(12*z*zeta[2]-3*z*zeta[3]+13*zeta[2]-3*zeta[3])/(1+z)*Hlog(z,[0])+1/16*Hlog(z,[-1, 0, 0, 0])-3/8*zeta[2]*Hlog(z,[-1, -1])+3/16*zeta[2]*Hlog(z,[-1, 0])+3/16*zeta[2]*Hlog(z,[0, -1])+(-1/24*zeta[2]-53/288)*Hlog(z,[0, 0])-1/144*(14+11*z)/(1+z)*Hlog(z,[0, 0, 0])-1/24*Hlog(z,[0, 0, 0, 0])+1/16*Hlog(z,[0, -1, 0, 0])+25/144*Hlog(z,[-1, 0, 0])+1/16*zeta[3]-1/40*zeta[2]^2-445/864*zeta[2]+1/72+(25/48*zeta[2]-1/8*zeta[3])*Hlog(z,[-1]):
v:=FB(diagrams["3e-Gamma1-03",-1]-direct):
fprintf(fd,"3e-Gamma1-03 export_diff\t%a\n",v):
dd:=FB(direct-subs(z=1/z,direct)):
fprintf(fd,"3e-Gamma1-03 direct_crossing\t%a\n",dd):
direct:=-1/144*(8*z+5)/(1+z)*Hlog(z,[0, 0, 0])-1/4*zeta[2]*Hlog(z,[0, 0])+3/16*Hlog(z,[0, -1, 0, 0])+13/144*Hlog(z,[-1, 0, 0])-3/8*Hlog(z,[-1, -1, 0, 0])-1/8*Hlog(z,[0, 0, 0, 0])-9/8*zeta[2]*Hlog(z,[-1, -1])+(13/48*zeta[2]-3/8*zeta[3])*Hlog(z,[-1])+9/16*zeta[2]*Hlog(z,[-1, 0])-1/48*(7*z*zeta[2]-9*z*zeta[3]+6*zeta[2]-9*zeta[3])/(1+z)*Hlog(z,[0])-11/15*zeta[2]^2-13/162*zeta[2]+83/432*zeta[3]+9/16*zeta[2]*Hlog(z,[0, -1])+3/16*Hlog(z,[-1, 0, 0, 0]):
v:=FB(diagrams["3e-Gamma1-04",-1]-direct):
fprintf(fd,"3e-Gamma1-04 export_diff\t%a\n",v):
dd:=FB(direct-subs(z=1/z,direct)):
fprintf(fd,"3e-Gamma1-04 direct_crossing\t%a\n",dd):
direct:=coeff(expand(convert(series(-8*(1-eps)/(3-2*eps),eps=0,5),polynom)*((0)*eps^(-4)+(-1/192*zeta[2]+1/384)*eps^(-3)+(1/96*zeta[3]+1/64*Hlog(z,[-1, 0, 0])+3/64*zeta[2]*Hlog(z,[-1])-3/128*zeta[2]*Hlog(z,[0])-1/128*Hlog(z,[0, 0, 0])+3/256*Hlog(z,[0, 0])+1/384+19/768*zeta[2])*eps^(-2)+(1/512*(11*z*zeta[2]+28*z*zeta[3]+13*zeta[2]+28*zeta[3])/(z+1)*Hlog(z,[0])+7/128*Hlog(z,[0, -1, 0, 0])+21/128*zeta[2]*Hlog(z,[-1, 0])+(-5/64*zeta[2]+3/256)*Hlog(z,[0, 0])-1/64*Hlog(z,[-1, 0, 0])+(-7/64*zeta[3]-3/64*zeta[2])*Hlog(z,[-1])+1/24*zeta[2]+3/128-89/384*zeta[2]^2-5/192*zeta[3]+1/512*(3*z+5)/(z+1)*Hlog(z,[0, 0, 0])-7/64*Hlog(z,[-1, -1, 0, 0])-21/64*zeta[2]*Hlog(z,[-1, -1])-3/64*Hlog(z,[0, 0, 0, 0])+7/128*Hlog(z,[-1, 0, 0, 0])+21/128*zeta[2]*Hlog(z,[0, -1]))*eps^(-1)+(17/384*zeta[2]+17/384+5/128*Hlog(z,[-1, 0, 0, 0])-45/256*Hlog(z,[0, 0, 0, 0])+5/128*Hlog(z,[0, -1, 0, 0])-1/64*Hlog(z,[-1, 0, 0, 0, 0])+7/64*Hlog(z,[0, 0, -1, 0, 0])+5/128*Hlog(z,[0, 0, 0, 0, 0])+5/128*Hlog(z,[0, -1, 0, 0, 0])-5/64*Hlog(z,[-1, -1, 0, 0])+27/64*Hlog(z,[-1, -1, -1, 0, 0])-27/128*Hlog(z,[-1, 0, -1, 0, 0])-27/128*Hlog(z,[0, -1, -1, 0, 0])-27/128*Hlog(z,[-1, -1, 0, 0, 0])-81/128*zeta[2]*Hlog(z,[-1, 0, -1])+1/4*zeta[2]*Hlog(z,[0, -1, 0])+81/64*zeta[2]*Hlog(z,[-1, -1, -1])-81/128*zeta[2]*Hlog(z,[0, -1, -1])+21/64*zeta[2]*Hlog(z,[0, 0, -1])-81/128*zeta[2]*Hlog(z,[-1, -1, 0])+41/48*zeta[2]*zeta[3]-997/3840*zeta[2]^2+17/384*zeta[3]+5/12*zeta[5]-1/2560*(788*z*zeta[2]^2+65*z*zeta[2]-100*z*zeta[3]+788*zeta[2]^2+55*zeta[2]-100*zeta[3])/(z+1)*Hlog(z,[0])-1/512*(36*z*zeta[2]+5*z+36*zeta[2]+3)/(z+1)*Hlog(z,[0, 0, 0])+(-27/128*zeta[3]+15/128*zeta[2])*Hlog(z,[-1, 0])+(-59/512*zeta[2]+11/64*zeta[3]+9/256)*Hlog(z,[0, 0])+(25/128*zeta[2]+1/64)*Hlog(z,[-1, 0, 0])+(-5/64*zeta[3]+309/640*zeta[2]^2+3/64*zeta[2])*Hlog(z,[-1])+(27/64*zeta[3]-15/64*zeta[2])*Hlog(z,[-1, -1])+(-27/128*zeta[3]+15/128*zeta[2])*Hlog(z,[0, -1]))*eps^(0))),eps,-1):
v:=FB(diagrams["3e-3g+Pi1-11",-1]-direct):
fprintf(fd,"3e-3g+Pi1-11 export_diff\t%a\n",v):
dd:=FB(direct-subs(z=1/z,direct)):
fprintf(fd,"3e-3g+Pi1-11 direct_crossing\t%a\n",dd):
fclose(fd):printf("DIRECT_EXPORT_CHECK_COMPLETE\n"):quit:
