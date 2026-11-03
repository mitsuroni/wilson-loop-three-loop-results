"""Inventory the new current reader/metadata/evidence after final integration."""
from prepare import HERE,ROOT,FAMILY,read,put,sha

def main():
    files={p for p in FAMILY.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    for entry in read(HERE/'entries.json'):files.update(p for p in (ROOT/entry['path']).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    files.update(p for p in HERE.rglob('*') if p.is_file() and p!=HERE/'manifest.json' and 'runs' not in p.parts and not any(v.startswith('isolation-') for v in p.relative_to(HERE).parts[:-1]) and '__pycache__' not in p.parts)
    for rel in ('integrals/README.md','integrals/three-loop/README.md','workflow/master-plan.md','workflow/project-map.md','workflow/data-flow.md','workflow/backlog/tasks/wilson-5.25 - Complete-continuous-derivations-and-audit-all39-three-loop6e-integrals.md'):files.add(ROOT/rel)
    inventory={p.relative_to(ROOT).as_posix():{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(files)}
    put(HERE/'manifest.json',{'format':'wilson-6e-continuous-reader/1','entries':39,'scope':'Current reader, exact copied runtime, object indices and actual evidence. Temporary runs/copies excluded. This manifest is a binding inventory, not mathematical proof.','files':inventory})
    print({'manifest_files':len(inventory)})
if __name__=='__main__':main()
