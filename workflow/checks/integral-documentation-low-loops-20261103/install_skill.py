"""Install only the reviewed project skill into the user's discoverable skills."""
from pathlib import Path
import shutil,json
from prepare import HERE,ROOT,sha,put
def main():
    source=ROOT/'workflow/skills/wilson-integral-derivations'
    destination=Path('C:/Users/Lenovo/.codex/skills/wilson-integral-derivations')
    assert destination.parent.resolve()==Path('C:/Users/Lenovo/.codex/skills').resolve()
    files={p.relative_to(source).as_posix():sha(p) for p in source.rglob('*') if p.is_file()}
    if destination.exists():
        assert {p.relative_to(destination).as_posix():sha(p) for p in destination.rglob('*') if p.is_file()}==files,'Existing different skill; do not overwrite silently'
    else:shutil.copytree(source,destination)
    for p,h in files.items():assert sha(destination/p)==h
    report={'status':'installed_and_byte_verified','name':'wilson-integral-derivations','source':source.as_posix(),'destination':destination.as_posix(),'files':files,
        'scope':'Only this newly requested skill; no unrelated configuration, automatic commit, upload, or selector mutation.'}
    put(HERE/'skill-install.json',report);print(json.dumps({'status':report['status'],'name':report['name'],'files':len(files)}))
if __name__=='__main__':main()
