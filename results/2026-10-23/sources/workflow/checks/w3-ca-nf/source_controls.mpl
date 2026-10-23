restart:
_hyper_autoload_periods:=["C:/Users/Lenovo/Documents/learning/work/hyperint/periodLookups.m"]:
read "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl":
_hyper_verbosity:=0:
FB:=proc(f) fibrationBasis(convert(subs(Pi^6=216*zeta[2]^3,Pi^4=36*zeta[2]^2,Pi^2=6*zeta[2],f),HlogRegInf),[z]) end proc:
read "D:/wilson-agent/workflow/checks/w3-ca-nf/source_inputs.mpl":
fd:=fopen("source-crossing.txt",WRITE,TEXT):
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["4e-Pi1-01",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"4e-Pi1-01\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 4e-Pi1-01 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["4e-Pi1-02",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"4e-Pi1-02\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 4e-Pi1-02 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["4e-Pi1-03",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"4e-Pi1-03\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 4e-Pi1-03 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["4e-Pi1-04",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"4e-Pi1-04\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 4e-Pi1-04 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["4e-Pi1-05",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"4e-Pi1-05\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 4e-Pi1-05 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["4e-Pi1-06",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"4e-Pi1-06\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 4e-Pi1-06 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["4e-Pi1-07",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"4e-Pi1-07\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 4e-Pi1-07 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["4e-Pi1-08",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"4e-Pi1-08\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 4e-Pi1-08 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["4e-Pi1-09",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"4e-Pi1-09\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 4e-Pi1-09 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["4e-Pi1-10",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"4e-Pi1-10\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 4e-Pi1-10 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["4e-Pi1-11",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"4e-Pi1-11\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 4e-Pi1-11 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["2e-Pi1+Pi1-01",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"2e-Pi1+Pi1-01\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 2e-Pi1+Pi1-01 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["2e-Pi1+Pi1-02",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"2e-Pi1+Pi1-02\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 2e-Pi1+Pi1-02 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["2e-Pi2-01",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"2e-Pi2-01\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 2e-Pi2-01 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["2e-Pi2-02",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"2e-Pi2-02\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 2e-Pi2-02 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-3g+Pi1-01",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-3g+Pi1-01\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-3g+Pi1-01 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-3g+Pi1-02",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-3g+Pi1-02\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-3g+Pi1-02 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-3g+Pi1-03",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-3g+Pi1-03\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-3g+Pi1-03 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-3g+Pi1-04",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-3g+Pi1-04\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-3g+Pi1-04 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-3g+Pi1-05",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-3g+Pi1-05\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-3g+Pi1-05 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-3g+Pi1-06",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-3g+Pi1-06\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-3g+Pi1-06 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-3g+Pi1-07",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-3g+Pi1-07\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-3g+Pi1-07 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-3g+Pi1-08",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-3g+Pi1-08\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-3g+Pi1-08 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-3g+Pi1-09",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-3g+Pi1-09\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-3g+Pi1-09 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-3g+Pi1-10",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-3g+Pi1-10\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-3g+Pi1-10 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-3g+Pi1-11",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-3g+Pi1-11\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-3g+Pi1-11 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-3g+Pi1-12",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-3g+Pi1-12\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-3g+Pi1-12 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-Gamma1-01",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-Gamma1-01\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-Gamma1-01 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-Gamma1-02",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-Gamma1-02\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-Gamma1-02 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-Gamma1-03",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-Gamma1-03\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-Gamma1-03 eps%d length%d\n",q,length(dd)) fi:
 od:
for q from -6 to 0 do
 v:=FB(subs(r=1/(1+z),diagrams["3e-Gamma1-04",q])):
 dd:=FB(v-subs(z=1/z,v)):
 fprintf(fd,"3e-Gamma1-04\t%d\t%a\n",q,dd):
 if dd<>0 then printf("SOURCE_CROSSING_FAILURE 3e-Gamma1-04 eps%d length%d\n",q,length(dd)) fi:
 od:
fclose(fd):printf("SOURCE_CROSSING_AUDIT_COMPLETE\n"):quit:
