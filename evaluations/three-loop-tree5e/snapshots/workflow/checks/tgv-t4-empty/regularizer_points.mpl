restart:
read "D:/wilson-agent/workflow/checks/tgv-t4-empty/runs/20260918T021524768893Z-tw_vfirst/regularized_TW.m":
fd:=fopen("regularizer-points.txt",WRITE):
for k from 1 to 3 do for bv from 0 to 1 do
 pointValue:=subs({r=k,s=k+1,t=k+2,u=k+3,v=k+4,z=2/5,eps=13/20,aa=0,bb=bv},Rg):
 fprintf(fd,"%d %d %a\n",k,bv,pointValue):
od od:fclose(fd):quit: