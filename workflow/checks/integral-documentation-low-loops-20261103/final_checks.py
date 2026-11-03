"""Record actual skill and management checks, separate from mathematical audit."""
from pathlib import Path
import subprocess,sys
from prepare import HERE,ROOT,read,put

def run(label,args,json_output=False):
    r=subprocess.run([sys.executable,'-X','utf8','-B',*map(str,args)],cwd=ROOT,text=True,encoding='utf-8',capture_output=True)
    record={'arguments':[str(a) for a in args],'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
    if json_output and r.returncode==0:record['parsed']=__import__('json').loads(r.stdout)
    put(HERE/(label+'.json'),record);assert r.returncode==0,(label,r.stderr,r.stdout)
    return record

def main():
    validator=Path('C:/Users/Lenovo/.codex/skills/.system/skill-creator/scripts/quick_validate.py')
    run('skill-validation-project',[validator,ROOT/'workflow/skills/wilson-integral-derivations'])
    run('skill-validation-installed',[validator,Path('C:/Users/Lenovo/.codex/skills/wilson-integral-derivations')])
    run('project-index-build',['-m','workflow.project_index','build','--json'],True)
    result=run('project-index-check',['-m','workflow.project_index','check','--json'],True)['parsed']
    print({'skills_validated':2,'project_index_ok':result['ok'],'counts':result['counts'],'scope':'Management/schema checks only; existing historical warnings retain their original status.'})
if __name__=='__main__':main()
