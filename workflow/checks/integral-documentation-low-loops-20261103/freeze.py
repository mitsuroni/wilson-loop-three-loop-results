"""Inventory the final new reader/skill/evidence; never edit historical manifests."""
from prepare import HERE,ROOT,READING,read,put,sha

def main():
    files=set(p for p in READING.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    files.update(p for p in (ROOT/'workflow/skills/wilson-integral-derivations').rglob('*') if p.is_file())
    files.update(p for p in HERE.rglob('*') if p.is_file() and p!=HERE/'manifest.json' and 'runs' not in p.parts and 'isolation-r001' not in p.parts and '__pycache__' not in p.parts)
    for rel in ('README.md','workflow/master-plan.md','workflow/project-map.md','workflow/data-flow.md','workflow/integral-documentation-standard.md','workflow/backlog/tasks/wilson-5.24 - Establish-integral-derivation-standard-and-audit-the-low-loop-pilot.md'):
        files.add(ROOT/rel)
    manifest={p.relative_to(ROOT).as_posix():{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(files)}
    put(HERE/'manifest.json',{'format':'wilson-integral-reader-trial/1','scope':'Exact new reading/skill/evidence and scoped integration files; temporary runs/copies excluded. Historical originals are independently protected by preflight/artifact-check. This manifest is not mathematical evidence.','files':manifest})
    print({'manifest_files':len(manifest)})
if __name__=='__main__':main()
