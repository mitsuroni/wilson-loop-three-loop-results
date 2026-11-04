"""Refresh management metadata after final task integration; not mathematical proof."""
import subprocess,sys,json
from prepare import HERE,ROOT,put
def main():
    for name in ('build','check'):
        args=[sys.executable,'-X','utf8','-B','-m','workflow.project_index',name,'--json'];r=subprocess.run(args,cwd=ROOT,text=True,encoding='utf-8',capture_output=True)
        record={'arguments':args,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
        if r.returncode==0:record['parsed']=json.loads(r.stdout)
        put(HERE/('project-index-'+name+'.json'),record);assert r.returncode==0,(name,r.stderr,r.stdout)
    print({'project_index_ok':record['parsed']['ok'],'counts':record['parsed']['counts'],'scope':'Management consistency only; historical warnings are not rebound.'})
if __name__=='__main__':main()
