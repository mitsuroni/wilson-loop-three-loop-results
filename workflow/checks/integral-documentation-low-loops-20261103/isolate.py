"""Replay the new reader entrypoints from a real, self-contained file copy."""
from pathlib import Path
import shutil,subprocess,sys
from prepare import HERE,ROOT,READING,sha,put

def main():
    copy=HERE/'isolation-r001';assert not copy.exists(),'Choose a new isolation revision'
    target=copy/'integrals'
    shutil.copytree(READING,target,ignore=shutil.ignore_patterns('__pycache__'))
    files={p.relative_to(target).as_posix():sha(p) for p in target.rglob('*') if p.is_file()}
    put(HERE/'isolation-input.json',{'reader_root':str(target),'files':files,'scope':'Real copied reader and historical runtime. Mathematical code/data closure can run without the original publication checkout. Article text may receive later editorial fixes; this manifest binds the actual copied version.'})
    out=HERE/'runs/isolation-r001'
    command=[sys.executable,'-X','utf8','-B',str(HERE/'run_checks.py'),'--reader-root',str(target),'--output-root',str(out)]
    subprocess.run(command,cwd=ROOT,check=True)
    print('Isolated replay completed: '+str(out))
if __name__=='__main__':main()
