restart:
read "D:/wilson-agent/workflow/checks/tgv-t4-empty/runs/20260918T021524768893Z-tw_vfirst/regularized_TW.m":
fd:=fopen("regularized_TW.txt",WRITE):fprintf(fd,"%a\n",Rg):fclose(fd):quit: