"""Refresh the required project-management index; no mathematical proof claim."""
import subprocess,sys
from prepare import HERE,ROOT,read,put

def main():
    for name in ('build','check'):
        command=[sys.executable,'-X','utf8','-B','-m','workflow.project_index',name,'--json']
        r=subprocess.run(command,cwd=ROOT,text=True,encoding='utf-8',capture_output=True)
        record={'arguments':command,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
        if r.returncode==0:record['parsed']=__import__('json').loads(r.stdout)
        put(HERE/('project-index-'+name+'.json'),record);assert r.returncode==0,(name,r.stderr,r.stdout)
    result=read(HERE/'project-index-check.json')['parsed'];print({'project_index_ok':result['ok'],'counts':result['counts'],'scope':'Management consistency only. Historical warnings are preserved; no mathematical result is promoted.'})
if __name__=='__main__':main()
